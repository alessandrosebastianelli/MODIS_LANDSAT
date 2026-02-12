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
from tqdm import tqdm
import argparse
import sys
import subprocess

# ============================================
# COMMAND LINE ARGUMENTS
# ============================================

parser = argparse.ArgumentParser(description='Download MODIS data')
parser.add_argument('--year', type=int, default=None,
                    help='Download data for specific year only')
parser.add_argument('--start-year', type=int, default=None,
                    help='Start year for filtering')
parser.add_argument('--end-year', type=int, default=None,
                    help='End year for filtering')
args = parser.parse_args()

# ============================================
# LOGGING SETUP
# ============================================

def setup_logging(config, year_suffix=None):
    """Setup clean, consistent logging"""
    log_level = getattr(logging, config['logging']['level'])
    log_file = config['logging']['log_file']
    
    # Add year suffix if provided
    if year_suffix:
        log_file = log_file.replace('.log', f'_{year_suffix}.log')
    
    # Move to logs/ directory
    log_file = Path('logs') / Path(log_file).name
    log_file.parent.mkdir(parents=True, exist_ok=True)
    
    logger = logging.getLogger('MODIS')
    logger.setLevel(log_level)
    logger.handlers = []
    logger.propagate = False
    
    # File handler
    fh = logging.FileHandler(log_file, mode='a')
    fh.setLevel(log_level)
    fh.setFormatter(logging.Formatter(
        '%(asctime)s | %(name)-15s | %(levelname)-8s | %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    ))
    logger.addHandler(fh)
    
    # Console handler
    ch = logging.StreamHandler()
    ch.setLevel(log_level)
    
    if config['logging'].get('console_colors', False):
        class ColoredFormatter(logging.Formatter):
            COLORS = {
                'DEBUG': '\033[36m',
                'INFO': '\033[32m',
                'WARNING': '\033[33m',
                'ERROR': '\033[31m',
                'CRITICAL': '\033[35m'
            }
            RESET = '\033[0m'
            
            def format(self, record):
                color = self.COLORS.get(record.levelname, self.RESET)
                levelname_colored = f"{color}{record.levelname[0]}{self.RESET}"
                return f"{levelname_colored} | {record.getMessage()}"
        
        ch.setFormatter(ColoredFormatter())
    else:
        ch.setFormatter(logging.Formatter('%(levelname)s | %(message)s'))
    
    logger.addHandler(ch)
    return logger

# ============================================
# INITIALIZE
# ============================================

with open('config_refactored.yaml', 'r') as f:
    config = yaml.safe_load(f)

logger = setup_logging(config, year_suffix=args.year if args.year else None)

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

# Override dates if year filtering requested
if args.year:
    start_date = f'{args.year}-01-01'
    end_date = f'{args.year}-12-31'
    logger.info(f"Filtering to year: {args.year}")
elif args.start_year or args.end_year:
    if args.start_year:
        start_date = f'{args.start_year}-01-01'
    if args.end_year:
        end_date = f'{args.end_year}-12-31'
    logger.info(f"Filtering to range: {start_date} to {end_date}")

logger.info(f"Date range: {start_date} to {end_date}")

base_dir = Path(config['output']['base_dir'])
modis_dir = base_dir / config['output']['subdirs']['modis']
tiles_dir = modis_dir / f'tiles_temp_{args.year}' if args.year else modis_dir / 'tiles_temp'

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
    view_angle_day = image.select('Day_view_angle').rename('MODIS_View_angle_Day')
    view_angle_night = image.select('Night_view_angle').rename('MODIS_View_angle_Night')
    
    # View times (scale: value * 0.1, hours)
    view_time_day = image.select('Day_view_time').multiply(0.1).rename('MODIS_View_time_Day')
    view_time_night = image.select('Night_view_time').multiply(0.1).rename('MODIS_View_time_Night')
    
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
    """Mosaic tiles into single file using GDAL VRT"""
    valid_tiles = [f for f in tile_files if f.exists() and f.stat().st_size > 1000]
    
    if len(valid_tiles) == 0:
        logger.error("No valid tiles to mosaic")
        return False
    
    try:
        vrt_file = output_file.parent / f"{output_file.stem}_temp.vrt"
        
        vrt_cmd = ['gdalbuildvrt', str(vrt_file)] + [str(f) for f in valid_tiles]
        result = subprocess.run(vrt_cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            logger.error(f"VRT creation failed: {result.stderr}")
            return False
        
        tif_cmd = [
            'gdal_translate',
            '-co', 'COMPRESS=DEFLATE',
            '-co', 'PREDICTOR=2',
            '-co', 'TILED=YES',
            '-co', 'BLOCKXSIZE=256',
            '-co', 'BLOCKYSIZE=256',
            str(vrt_file),
            str(output_file)
        ]
        
        result = subprocess.run(tif_cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            logger.error(f"GeoTIFF creation failed: {result.stderr}")
            return False
        
        if vrt_file.exists():
            vrt_file.unlink()
        
        return output_file.exists() and output_file.stat().st_size > 1000
        
    except Exception as e:
        logger.error(f"Mosaic error: {e}")
        return False

# ============================================
# BATCH DOWNLOAD WITH CLEAN LOGGING
# ============================================

composites_list = modis_composites.toList(num_composites)

logger.info("="*70)
logger.info(f"Starting MODIS download: {num_composites} composites")
logger.info(f"Estimated tiles per composite: {len(tiles)}")
logger.info(f"Total downloads: {num_composites * len(tiles)} tiles")
logger.info(f"Estimated time: {num_composites * len(tiles) * config['download']['delay_between_tiles'] / 60:.1f} min (delays only)")
logger.info("="*70)

stats = {'downloaded': 0, 'skipped': 0, 'errors': 0, 'tiles_downloaded': 0}

for i in range(num_composites):
    composite = ee.Image(composites_list.get(i))
    date_str = composite.get('date_string').getInfo()
    
    final_file = modis_dir / f"MODIS_{date_str}_{i:04d}.tif"
    
    progress_pct = (i + 1) / num_composites * 100
    logger.info(f"[{i+1:3d}/{num_composites}] ({progress_pct:5.1f}%) Composite: {date_str}")
    
    if final_file.exists():
        logger.info("Status: SKIPPED (file exists)")
        stats['skipped'] += 1
        continue
    
    # Add coordinates and dynamic layers
    composite_final = add_coordinates_and_dynamic_layers(composite)
    
    # Download all tiles for this composite
    tile_files = []
    tiles_success = 0
    
    logger.info(f"Downloading {len(tiles)} tiles")
    
    for tile_idx, tile in enumerate(tqdm(tiles, desc="Tiles", unit="tile", leave=False, ncols=80)):
        tile_file = tiles_dir / f"{date_str.replace('-', '')}_c{i:03d}_t{tile['id']:02d}.tif"
        
        logger.debug(f"  Tile {tile_idx+1}/{len(tiles)}: bounds={tile['bounds']}")
        
        if download_tile(composite_final, tile['bounds'], tile_file):
            tile_files.append(tile_file)
            tiles_success += 1
            stats['tiles_downloaded'] += 1
            file_size = tile_file.stat().st_size / 1024  # KB
            logger.debug(f"    SUCCESS: {file_size:.1f} KB")
        else:
            logger.warning(f"    FAILED: Tile {tile_idx} download error")
        
        time.sleep(config['download']['delay_between_tiles'])
    
    logger.info(f"Tiles completed: {tiles_success}/{len(tiles)} successful ({tiles_success/len(tiles)*100:.1f}%)")
    
    # Check minimum tiles requirement
    min_required = len(tiles) * config['download']['min_tiles_success_ratio']
    if tiles_success < min_required:
        logger.error(f"Status: FAILED - Insufficient tiles ({tiles_success} < {min_required:.0f} required)")
        stats['errors'] += 1
        for tile_file in tile_files:
            if tile_file.exists():
                tile_file.unlink()
        continue
    
    # Mosaic tiles
    logger.info(f"Mosaicking {tiles_success} tiles")
    if mosaic_tiles(tile_files, final_file):
        file_size_mb = final_file.stat().st_size / (1024**2)
        logger.info(f"Status: SUCCESS - {file_size_mb:.1f} MB")
        stats['downloaded'] += 1
        
        # Cleanup tile files
        for tile_file in tile_files:
            if tile_file.exists():
                tile_file.unlink()
    else:
        logger.error("Status: FAILED - Mosaic error")
        stats['errors'] += 1
    
    time.sleep(config['download']['delay_between_images'])

# Cleanup tiles directory
try:
    tiles_dir.rmdir()
except:
    logger.warning("Could not remove tiles directory (may contain files)")

# ============================================
# SUMMARY
# ============================================

logger.info("\n" + "="*70)
logger.info("MODIS DOWNLOAD SUMMARY")
logger.info("="*70)
logger.info(f"Downloaded: {stats['downloaded']}")
logger.info(f"Skipped: {stats['skipped']}")
logger.info(f"Errors: {stats['errors']}")
logger.info(f"Total tiles downloaded: {stats['tiles_downloaded']}")
logger.info("="*70)

logger.info("\nBand structure in each file:")
logger.info("  MODIS bands (resampled 1km -> 30m):")
logger.info("    MODIS_LST_Day_K, MODIS_LST_Night_K")
logger.info("    MODIS_SST_Day_K, MODIS_SST_Night_K")
logger.info("    MODIS_QC_Day, MODIS_QC_Night")
logger.info("    MODIS_Emis_31, MODIS_Emis_32")
logger.info("    MODIS_View_angle_Day, MODIS_View_angle_Night")
logger.info("    MODIS_View_time_Day, MODIS_View_time_Night")
logger.info(f"  Dynamic layers ({interval_days}-day frequency):")
logger.info("    land_cover: MODIS land cover")
logger.info("    emissivity_dynamic: Emissivity from land cover")
logger.info("  Coordinates:")
logger.info("    lon, lat, time")
logger.info(f"  Total: ~17 bands per file")
logger.info("="*70)

# Cleanup temporary tiles directory
import shutil
if tiles_dir.exists():
    try:
        shutil.rmtree(tiles_dir)
        logger.info(f"Cleaned up temporary tiles: {tiles_dir}")
    except Exception as e:
        logger.warning(f"Could not remove tiles directory: {e}")