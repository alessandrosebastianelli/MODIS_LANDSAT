# ============================================
# Landsat Download - Full Bands + LST + SST
# 15-day composites at 30m resolution
# ============================================

import ee
from datetime import datetime, timedelta
import time
from pathlib import Path
import requests
import urllib.parse
import yaml
import logging
import rasterio
from rasterio.merge import merge

# ============================================
# LOGGING SETUP
# ============================================

def setup_logging(config):
    log_level = getattr(logging, config['logging']['level'])
    log_file = config['logging']['log_file']
    
    logger = logging.getLogger('Landsat')
    logger.setLevel(log_level)
    logger.handlers = []
    
    fh = logging.FileHandler(log_file, mode='a')
    fh.setLevel(log_level)
    fh.setFormatter(logging.Formatter(
        '%(asctime)s | %(levelname)-8s | %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    ))
    logger.addHandler(fh)
    
    ch = logging.StreamHandler()
    ch.setLevel(log_level)
    
    if config['logging']['console_colors']:
        class ColoredFormatter(logging.Formatter):
            COLORS = {'DEBUG': '\033[36m', 'INFO': '\033[32m', 'WARNING': '\033[33m',
                     'ERROR': '\033[31m', 'CRITICAL': '\033[35m'}
            RESET = '\033[0m'
            
            def format(self, record):
                color = self.COLORS.get(record.levelname, self.RESET)
                record.levelname = f"{color}{record.levelname:<8}{self.RESET}"
                return super().format(record)
        
        ch.setFormatter(ColoredFormatter('%(asctime)s | %(levelname)s | %(message)s',
                                        datefmt='%Y-%m-%d %H:%M:%S'))
    else:
        ch.setFormatter(logging.Formatter('%(asctime)s | %(levelname)-8s | %(message)s',
                                         datefmt='%Y-%m-%d %H:%M:%S'))
    
    logger.addHandler(ch)
    return logger

# ============================================
# INITIALIZE
# ============================================

with open('config_refactored.yaml', 'r') as f:
    config = yaml.safe_load(f)

logger = setup_logging(config)

try:
    ee.Initialize(project=config['ee_project'])
    logger.info("Earth Engine initialized")
except:
    ee.Authenticate()
    ee.Initialize(project=config['ee_project'])
    logger.info("Earth Engine authenticated and initialized")

# ============================================
# PARAMETERS
# ============================================

geometry = ee.Geometry.Rectangle([
    config['geometry']['min_lon'],
    config['geometry']['min_lat'],
    config['geometry']['max_lon'],
    config['geometry']['max_lat']
])

start_date = config['temporal']['start_date']
end_date = config['temporal']['end_date']
interval_days = config['temporal']['interval_days']
TARGET_RESOLUTION = config['spatial']['target_resolution']
TILE_SIZE = config['spatial']['tile_size']

base_dir = Path(config['output']['base_dir'])
landsat_dir = base_dir / config['output']['subdirs']['landsat']
tiles_dir = landsat_dir / 'tiles_temp'

landsat_dir.mkdir(parents=True, exist_ok=True)
tiles_dir.mkdir(parents=True, exist_ok=True)

logger.info("="*70)
logger.info("Landsat Download - Full Bands + LST + SST")
logger.info("="*70)
logger.info(f"Date range: {start_date} to {end_date}")
logger.info(f"Interval: {interval_days} days")
logger.info(f"Resolution: {TARGET_RESOLUTION}m")

# ============================================
# TILING
# ============================================

def create_tiles(bounds, tile_size):
    coords = bounds.coordinates().getInfo()[0]
    min_lon = min(c[0] for c in coords)
    max_lon = max(c[0] for c in coords)
    min_lat = min(c[1] for c in coords)
    max_lat = max(c[1] for c in coords)
    
    tiles = []
    lat = min_lat
    tile_id = 0
    while lat < max_lat:
        lon = min_lon
        while lon < max_lon:
            tile_bounds = ee.Geometry.Rectangle([
                lon, lat,
                min(lon + tile_size, max_lon),
                min(lat + tile_size, max_lat)
            ])
            tiles.append({'id': tile_id, 'bounds': tile_bounds})
            tile_id += 1
            lon += tile_size
        lat += tile_size
    
    return tiles

tiles = create_tiles(geometry, TILE_SIZE)
logger.info(f"Created {len(tiles)} tiles")

# ============================================
# LST/SST CALCULATION
# ============================================

def calculate_lst_sst(image):
    """Calculate LST and SST from Landsat thermal band"""
    # NDVI
    ndvi = image.normalizedDifference(['SR_B5', 'SR_B4']).rename('NDVI')
    
    # Fractional Vegetation
    fv = ndvi.subtract(0.2).divide(0.5).pow(2).clamp(0, 1).rename('FV')
    
    # Emissivity
    emissivity = fv.multiply(0.004).add(0.986).rename('emissivity')
    
    # Brightness Temperature
    BT = image.select('ST_B10').multiply(0.00341802).add(149.0)
    
    # LST calculation
    lambda_um = 10.895
    rho = 14388.0
    par = BT.multiply(lambda_um).divide(rho)
    eps_safe = emissivity.where(emissivity.lte(0), 0.001).where(emissivity.gt(1), 1)
    lst = BT.divide(par.multiply(eps_safe.log()).add(1)).rename('LST_K')
    
    # SST (same as LST, will be masked to water later)
    sst = lst.rename('SST_K')
    
    # Cloud mask
    qa = image.select('QA_PIXEL')
    cloud_bit = 1 << 3
    shadow_bit = 1 << 4
    snow_bit = 1 << 5
    cloud_mask = (qa.bitwiseAnd(cloud_bit).eq(0)
                  .And(qa.bitwiseAnd(shadow_bit).eq(0))
                  .And(qa.bitwiseAnd(snow_bit).eq(0))).rename('cloud_mask').uint8()
    
    # LST with cloud mask
    lst_masked = lst.updateMask(cloud_mask).rename('LST_K_masked')
    
    return (image
            .addBands(lst)
            .addBands(lst_masked)
            .addBands(sst)
            .addBands(cloud_mask)
            .addBands(ndvi)
            .addBands(emissivity))

# ============================================
# LOAD LANDSAT COLLECTIONS
# ============================================

def load_landsat_sensor(sensor_config):
    """Load Landsat collection for one sensor"""
    collection = (ee.ImageCollection(sensor_config['collection'])
                  .filterBounds(geometry)
                  .filterDate(sensor_config['start_date'], sensor_config['end_date'])
                  .filter(ee.Filter.lt('CLOUD_COVER', config['landsat']['max_cloud_cover'])))
    
    # Path/Row filtering
    if config['landsat']['path_row_filter']['enabled']:
        paths = config['landsat']['path_row_filter']['paths']
        rows = config['landsat']['path_row_filter']['rows']
        
        if paths:
            collection = collection.filter(ee.Filter.inList('WRS_PATH', paths))
        if rows:
            collection = collection.filter(ee.Filter.inList('WRS_ROW', rows))
    
    return collection.map(calculate_lst_sst)

logger.info("Loading Landsat collections...")

all_landsat = ee.ImageCollection([])
for sensor in config['landsat']['sensors']:
    coll = load_landsat_sensor(sensor)
    count = coll.size().getInfo()
    logger.info(f"{sensor['name']}: {count} images")
    all_landsat = all_landsat.merge(coll)

total_images = all_landsat.size().getInfo()
logger.info(f"Total Landsat images: {total_images}")

# ============================================
# CREATE COMPOSITES
# ============================================

def create_composites(collection, start, end, interval):
    """Create temporal composites"""
    composites = []
    current = datetime.fromisoformat(start)
    end_dt = datetime.fromisoformat(end)
    
    comp_num = 0
    while current < end_dt:
        start_str = current.strftime('%Y-%m-%d')
        end_str = (current + timedelta(days=interval)).strftime('%Y-%m-%d')
        
        period_coll = collection.filterDate(start_str, end_str)
        
        # Composite method
        if config['compositing']['method'] == 'median':
            composite = period_coll.median()
        elif config['compositing']['method'] == 'mean':
            composite = period_coll.mean()
        else:  # mosaic
            if config['compositing']['quality_mosaic']:
                composite = period_coll.sort('CLOUD_COVER').mosaic()
            else:
                composite = period_coll.mosaic()
        
        composite = (composite
                     .set('system:time_start', ee.Date(start_str).millis())
                     .set('date_string', start_str)
                     .set('composite_number', comp_num))
        
        composites.append(composite)
        current += timedelta(days=interval)
        comp_num += 1
    
    return ee.ImageCollection(composites)

landsat_composites = create_composites(all_landsat, start_date, end_date, interval_days)
num_composites = landsat_composites.size().getInfo()
logger.info(f"Created {num_composites} composites")

# ============================================
# ADD COORDINATES
# ============================================

def add_coordinates(image):
    """Add lat/lon/time"""
    latlon = ee.Image.pixelLonLat().reproject(crs='EPSG:4326', scale=TARGET_RESOLUTION)
    
    time_img = ee.Image.constant(
        ee.Number(image.get('system:time_start')).divide(86400000)
    ).rename('time').float()
    
    return (image
            .addBands(latlon.select('longitude').rename('lon'))
            .addBands(latlon.select('latitude').rename('lat'))
            .addBands(time_img))

# ============================================
# DOWNLOAD FUNCTIONS
# ============================================

def download_tile(image, tile_bounds, filename):
    try:
        img_prepared = (image
                        .clip(tile_bounds)
                        .reproject(crs='EPSG:4326', scale=TARGET_RESOLUTION))
        
        url = img_prepared.getDownloadURL({
            'scale': TARGET_RESOLUTION,
            'crs': 'EPSG:4326',
            'region': tile_bounds,
            'format': 'GEO_TIFF',
            'filePerBand': False
        })
        
        url_safe = urllib.parse.quote(url, safe=':/?&=')
        r = requests.get(url_safe, timeout=config['download']['timeout_seconds'])
        
        if r.status_code == 200:
            with open(filename, 'wb') as f:
                f.write(r.content)
            if filename.exists() and filename.stat().st_size > 1000:
                return True
        return False
    except Exception as e:
        logger.error(f"Tile error: {str(e)[:80]}")
        return False

def mosaic_tiles(tile_files, output_file):
    src_files = [rasterio.open(f) for f in tile_files if f.exists()]
    
    if len(src_files) == 0:
        return False
    
    mosaic, out_trans = merge(src_files)
    
    out_meta = src_files[0].meta.copy()
    out_meta.update({
        "driver": "GTiff",
        "height": mosaic.shape[1],
        "width": mosaic.shape[2],
        "transform": out_trans,
        "compress": "deflate",
        "predictor": 2,
        "tiled": True,
        "blockxsize": 256,
        "blockysize": 256
    })
    
    with rasterio.open(output_file, "w", **out_meta) as dest:
        dest.write(mosaic)
    
    for src in src_files:
        src.close()
    
    return output_file.exists()

# ============================================
# BATCH DOWNLOAD
# ============================================

composites_list = landsat_composites.toList(num_composites)

logger.info("\n" + "="*70)
logger.info(f"Downloading {num_composites} Landsat composites")
logger.info("="*70)

stats = {'downloaded': 0, 'skipped': 0, 'errors': 0}

for i in range(num_composites):
    composite = ee.Image(composites_list.get(i))
    date_str = composite.get('date_string').getInfo()
    
    final_file = landsat_dir / f"Landsat_{date_str}_{i:04d}.tif"
    
    if final_file.exists():
        logger.info(f"[{i+1:3d}/{num_composites}] SKIPPED: {final_file.name}")
        stats['skipped'] += 1
        continue
    
    logger.info(f"[{i+1:3d}/{num_composites}] DOWNLOADING: {final_file.name}")
    
    composite_final = add_coordinates(composite)
    
    # Download tiles
    tile_files = []
    tiles_success = 0
    
    for tile in tiles:
        tile_file = tiles_dir / f"img{i:04d}_t{tile['id']:02d}.tif"
        
        if download_tile(composite_final, tile['bounds'], tile_file):
            tile_files.append(tile_file)
            tiles_success += 1
        
        time.sleep(config['download']['delay_between_tiles'])
    
    logger.info(f"[{i+1:3d}/{num_composites}] Tiles: {tiles_success}/{len(tiles)}")
    
    # Check minimum tiles
    min_required = len(tiles) * config['download']['min_tiles_success_ratio']
    if tiles_success < min_required:
        logger.error(f"[{i+1:3d}/{num_composites}] INSUFFICIENT tiles")
        stats['errors'] += 1
        for f in tile_files:
            if f.exists():
                f.unlink()
        continue
    
    # Mosaic
    if mosaic_tiles(tile_files, final_file):
        logger.info(f"[{i+1:3d}/{num_composites}] SUCCESS")
        stats['downloaded'] += 1
        
        # Cleanup tiles
        for f in tile_files:
            if f.exists():
                f.unlink()
    else:
        logger.error(f"[{i+1:3d}/{num_composites}] MOSAIC FAILED")
        stats['errors'] += 1
    
    time.sleep(config['download']['delay_between_images'])

# Cleanup
try:
    tiles_dir.rmdir()
except:
    pass

# ============================================
# SUMMARY
# ============================================

logger.info("\n" + "="*70)
logger.info("LANDSAT DOWNLOAD SUMMARY")
logger.info("="*70)
logger.info(f"Downloaded: {stats['downloaded']}")
logger.info(f"Skipped: {stats['skipped']}")
logger.info(f"Errors: {stats['errors']}")
logger.info("="*70)

# Band structure info
logger.info("\nBand structure in each file:")
logger.info("  Original Landsat bands:")
logger.info("    SR_B1-SR_B7: Surface Reflectance")
logger.info("    ST_B10: Thermal band")
logger.info("    QA_PIXEL: Quality flags")
logger.info("  Derived products:")
logger.info("    LST_K: Land Surface Temperature")
logger.info("    LST_K_masked: LST with cloud mask")
logger.info("    SST_K: Sea Surface Temperature")
logger.info("    cloud_mask: Binary mask")
logger.info("    NDVI, emissivity")
logger.info("  Coordinates:")
logger.info("    lon, lat, time")
logger.info("="*70)
