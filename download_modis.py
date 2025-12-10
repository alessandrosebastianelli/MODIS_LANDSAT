# ============================================
# MODIS Download - Full Bands + LST + SST
# 15-day composites, resampled 1km -> 30m
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
    
    logger = logging.getLogger('MODIS')
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
    #ee.Authenticate()
    #ee.Initialize(project=config['ee_project'])
    KEY_FILE = config['service_account']['key_file']
    SERVICE_ACCOUNT = config['service_account']['email']
    credentials = ee.ServiceAccountCredentials(SERVICE_ACCOUNT, KEY_FILE)
    ee.Initialize(credentials)
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
modis_dir = base_dir / config['output']['subdirs']['modis']
tiles_dir = modis_dir / 'tiles_temp'

modis_dir.mkdir(parents=True, exist_ok=True)
tiles_dir.mkdir(parents=True, exist_ok=True)

logger.info("="*70)
logger.info("MODIS Download - Full Bands + LST + SST")
logger.info("="*70)
logger.info(f"Date range: {start_date} to {end_date}")
logger.info(f"Interval: {interval_days} days")
logger.info(f"Resolution: 1km -> {TARGET_RESOLUTION}m (resampled)")

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
# LAND COVER AND EMISSIVITY
# ============================================

def get_land_cover_for_date(date_str):
    """Get MODIS land cover for a specific date (yearly data)"""
    year = datetime.fromisoformat(date_str).year
    
    lc_collection = ee.ImageCollection('MODIS/061/MCD12Q1')
    lc_image = lc_collection.filterDate(f'{year}-01-01', f'{year}-12-31').first()
    
    if lc_image:
        return lc_image.select('LC_Type1').rename('land_cover')
    else:
        return lc_collection.sort('system:time_start', False).first().select('LC_Type1').rename('land_cover')

def calculate_emissivity_from_lc(land_cover_image):
    """Calculate emissivity from MODIS land cover"""
    emis_lookup = config['dynamic_layers']['emissivity']['lookup']
    
    emis = ee.Image(0.0)
    for lc_class, emis_value in emis_lookup.items():
        mask = land_cover_image.eq(int(lc_class))
        emis = emis.where(mask, emis_value)
    
    return emis.rename('emissivity_dynamic')

# ============================================
# PROCESS MODIS IMAGE
# ============================================

def process_modis(image):
    """Process MODIS: scale values, add LST/SST"""
    # Scale LST (Kelvin * 50 -> Kelvin)
    lst_day = image.select('LST_Day_1km').multiply(0.02).rename('MODIS_LST_Day_K')
    lst_night = image.select('LST_Night_1km').multiply(0.02).rename('MODIS_LST_Night_K')
    
    # SST (same as LST for MODIS, will be masked to water later)
    sst_day = lst_day.rename('MODIS_SST_Day_K')
    sst_night = lst_night.rename('MODIS_SST_Night_K')
    
    # QC bands
    qc_day = image.select('QC_Day').rename('MODIS_QC_Day')
    qc_night = image.select('QC_Night').rename('MODIS_QC_Night')
    
    # Emissivity bands (scale: value * 0.002 + 0.49)
    emis_31 = image.select('Emis_31').multiply(0.002).add(0.49).rename('MODIS_Emis_31')
    emis_32 = image.select('Emis_32').multiply(0.002).add(0.49).rename('MODIS_Emis_32')
    
    # View angles (scale: value * 1.0, degrees)
    view_angle_day = image.select('View_angle_Day').rename('MODIS_View_angle_Day')
    view_angle_night = image.select('View_angle_Night').rename('MODIS_View_angle_Night')
    
    # View times (scale: value * 0.1, hours)
    view_time_day = image.select('View_time_Day').multiply(0.1).rename('MODIS_View_time_Day')
    view_time_night = image.select('View_time_Night').multiply(0.1).rename('MODIS_View_time_Night')
    
    # Combine all
    result = (lst_day
              .addBands(lst_night)
              .addBands(sst_day)
              .addBands(sst_night)
              .addBands(qc_day)
              .addBands(qc_night)
              .addBands(emis_31)
              .addBands(emis_32)
              .addBands(view_angle_day)
              .addBands(view_angle_night)
              .addBands(view_time_day)
              .addBands(view_time_night))
    
    return result.copyProperties(image, ['system:time_start'])

# ============================================
# LOAD MODIS COLLECTIONS
# ============================================

logger.info("Loading MODIS collections...")

# Terra (morning overpass ~10:30 AM)
terra = (ee.ImageCollection(config['modis']['collections']['terra'])
         .filterBounds(geometry)
         .filterDate(start_date, end_date)
         .map(process_modis))

# Aqua (afternoon overpass ~1:30 PM)
aqua = (ee.ImageCollection(config['modis']['collections']['aqua'])
        .filterBounds(geometry)
        .filterDate(start_date, end_date)
        .map(process_modis))

# Merge
modis_all = terra.merge(aqua)

terra_count = terra.size().getInfo()
aqua_count = aqua.size().getInfo()
total_count = modis_all.size().getInfo()

logger.info(f"Terra (MOD11A1): {terra_count} images")
logger.info(f"Aqua (MYD11A1): {aqua_count} images")
logger.info(f"Total MODIS: {total_count} images")

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

modis_composites = create_composites(modis_all, start_date, end_date, interval_days)
num_composites = modis_composites.size().getInfo()
logger.info(f"Created {num_composites} MODIS composites")

# ============================================
# ADD COORDINATES AND DYNAMIC LAYERS
# ============================================

def add_coordinates_and_dynamic_layers(image):
    """Add lat/lon/time + land cover + emissivity"""
    latlon = ee.Image.pixelLonLat().reproject(crs='EPSG:4326', scale=TARGET_RESOLUTION)
    
    time_img = ee.Image.constant(
        ee.Number(image.get('system:time_start')).divide(86400000)
    ).rename('time').float()
    
    # Get date for this composite
    date_str = ee.String(image.get('date_string')).getInfo()
    
    # Add land cover for this date
    land_cover = get_land_cover_for_date(date_str)
    land_cover = land_cover.reproject(crs='EPSG:4326', scale=TARGET_RESOLUTION)
    
    # Calculate emissivity from land cover
    emissivity_dyn = calculate_emissivity_from_lc(land_cover)
    
    return (image
            .addBands(latlon.select('longitude').rename('lon'))
            .addBands(latlon.select('latitude').rename('lat'))
            .addBands(time_img)
            .addBands(land_cover)
            .addBands(emissivity_dyn))

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

composites_list = modis_composites.toList(num_composites)

logger.info("\n" + "="*70)
logger.info(f"Downloading {num_composites} MODIS composites")
logger.info("="*70)

stats = {'downloaded': 0, 'skipped': 0, 'errors': 0}

for i in range(num_composites):
    composite = ee.Image(composites_list.get(i))
    date_str = composite.get('date_string').getInfo()
    
    final_file = modis_dir / f"MODIS_{date_str}_{i:04d}.tif"
    
    if final_file.exists():
        logger.info(f"[{i+1:3d}/{num_composites}] SKIPPED: {final_file.name}")
        stats['skipped'] += 1
        continue
    
    logger.info(f"[{i+1:3d}/{num_composites}] DOWNLOADING: {final_file.name}")
    
    composite_final = add_coordinates_and_dynamic_layers(composite)
    
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
logger.info("MODIS DOWNLOAD SUMMARY")
logger.info("="*70)
logger.info(f"Downloaded: {stats['downloaded']}")
logger.info(f"Skipped: {stats['skipped']}")
logger.info(f"Errors: {stats['errors']}")
logger.info("="*70)

logger.info("\nBand structure in each file:")
logger.info("  MODIS bands (resampled 1km -> 30m):")
logger.info("    MODIS_LST_Day_K, MODIS_LST_Night_K")
logger.info("    MODIS_SST_Day_K, MODIS_SST_Night_K")
logger.info("    MODIS_QC_Day, MODIS_QC_Night")
logger.info("    MODIS_Emis_31, MODIS_Emis_32")
logger.info("    MODIS_View_angle_Day, MODIS_View_angle_Night")
logger.info("    MODIS_View_time_Day, MODIS_View_time_Night")
logger.info("  Dynamic layers (15-day frequency):")
logger.info("    land_cover: MODIS land cover")
logger.info("    emissivity_dynamic: Emissivity from land cover")
logger.info("  Coordinates:")
logger.info("    lon, lat, time")
logger.info(f"  Total: ~17 bands per file")
logger.info("="*70)