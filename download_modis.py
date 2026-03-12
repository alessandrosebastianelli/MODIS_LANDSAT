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
parser.add_argument('--start-year', type=int, default=None)
parser.add_argument('--end-year', type=int, default=None)
args = parser.parse_args()

# ============================================
# LOGGING SETUP
# ============================================

def setup_logging(config, year_suffix=None):
    log_level = getattr(logging, config['logging']['level'])
    log_file = config['logging']['log_file']
    if year_suffix:
        log_file = log_file.replace('.log', f'_{year_suffix}.log')
    log_file = Path('logs') / Path(log_file).name
    log_file.parent.mkdir(parents=True, exist_ok=True)

    logger = logging.getLogger('MODIS')
    logger.setLevel(log_level)
    logger.handlers = []
    logger.propagate = False

    fh = logging.FileHandler(log_file, mode='a')
    fh.setLevel(log_level)
    fh.setFormatter(logging.Formatter(
        '%(asctime)s | %(name)-15s | %(levelname)-8s | %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    ))
    logger.addHandler(fh)

    ch = logging.StreamHandler()
    ch.setLevel(log_level)
    if config['logging'].get('console_colors', False):
        class ColoredFormatter(logging.Formatter):
            COLORS = {'DEBUG': '\033[36m', 'INFO': '\033[32m',
                      'WARNING': '\033[33m', 'ERROR': '\033[31m', 'CRITICAL': '\033[35m'}
            RESET = '\033[0m'
            def format(self, record):
                color = self.COLORS.get(record.levelname, self.RESET)
                return f"{color}{record.levelname[0]}{self.RESET} | {record.getMessage()}"
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
except Exception:
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

if args.year:
    start_date = f'{args.year}-01-01'
    end_date = f'{args.year}-12-31'
    logger.info(f"Filtering to year: {args.year}")
elif args.start_year or args.end_year:
    if args.start_year:
        start_date = f'{args.start_year}-01-01'
    if args.end_year:
        end_date = f'{args.end_year}-12-31'

base_dir = Path(config['output']['base_dir'])
modis_dir = base_dir / config['output']['subdirs']['modis']
tiles_dir = modis_dir / (f'tiles_temp_{args.year}' if args.year else 'tiles_temp')

modis_dir.mkdir(parents=True, exist_ok=True)
tiles_dir.mkdir(parents=True, exist_ok=True)

logger.info("="*70)
logger.info("MODIS Download - Full Bands + LST + SST")
logger.info("="*70)
logger.info(f"Date range: {start_date} to {end_date}")
logger.info(f"Collections: {config['modis']['collections']['terra']} + {config['modis']['collections']['aqua']}")

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
logger.info(f"Tiles per composite: {len(tiles)}")

# ============================================
# LAND COVER AND EMISSIVITY
# ============================================

def get_land_cover_for_date(date_str):
    year = datetime.fromisoformat(date_str).year
    if year == 2021:
        return ee.Image('ESA/WorldCover/v200/2021').select('Map').rename('land_cover')
    elif year == 2020:
        return ee.Image('ESA/WorldCover/v100/2020').select('Map').rename('land_cover')
    else:
        lc_col = ee.ImageCollection('MODIS/061/MCD12Q1')
        lc_img = lc_col.filterDate(f'{year}-01-01', f'{year}-12-31').first()
        lc_img = ee.Algorithms.If(
            lc_img,
            lc_img,
            lc_col.sort('system:time_start', False).first()
        )
        return ee.Image(lc_img).select('LC_Type1').rename('land_cover')

def calculate_emissivity_from_lc(land_cover_image):
    emis_lookup = config['dynamic_layers']['emissivity']['lookup']
    emis = ee.Image(0.986)
    for lc_class, emis_value in emis_lookup.items():
        mask = land_cover_image.eq(int(lc_class))
        emis = emis.where(mask, emis_value)
    return emis.rename('emissivity_dynamic')

# ============================================
# PROCESS MODIS IMAGE
# ============================================
#
# Scaling factors for MOD11A1 / MYD11A1 (Collection 6.1):
#   LST_Day_1km, LST_Night_1km : DN * 0.02         -> Kelvin
#   Emis_31, Emis_32           : DN * 0.002 + 0.49 -> emissivity [0,1]
#   Day_view_time,Night_view_time: DN * 0.1         -> hours (local solar time)
#   Day_view_angle,Night_view_angle: DN - 65        -> degrees (signed)
#   QC_Day, QC_Night           : bitmask, no scaling
#
# NOTE: SST is NOT simply a copy of LST.
# MOD11A1 does not provide a separate SST product; it provides LST for all
# land AND water pixels (the algorithm uses the split-window method with emissivity
# from bands 31/32 which is valid over water too).
# We therefore compute:
#   - LST_Day/Night   : raw scaled value (valid everywhere including water)
#   - SST_Day/Night   : same value, but the distinction is semantic — in the
#                       netcdf step the WATERMASK will be used to select which
#                       pixels are physically SST vs LST.
#   - LST_MODIS       : mean of Day + Night (NaN if either is missing)
#   - SST_MODIS       : same as LST_MODIS (will be masked to water in netcdf)
# ============================================

def process_modis(image):
    """Scale MODIS bands to physical units and derive composite products."""

    # --- Temperature (Kelvin) ---
    lst_day   = image.select('LST_Day_1km').multiply(0.02).rename('MODIS_LST_Day_K')
    lst_night = image.select('LST_Night_1km').multiply(0.02).rename('MODIS_LST_Night_K')

    # Daily mean LST: average of day and night where both are available.
    # unmask(0) replaces masked pixels with 0; we track validity with the
    # individual masks and divide only where both observations exist.
    day_valid   = lst_day.mask()
    night_valid = lst_night.mask()
    both_valid  = day_valid.And(night_valid)

    lst_mean = (lst_day.unmask(0).add(lst_night.unmask(0))
                       .divide(2)
                       .updateMask(both_valid)
                       .rename('LST_MODIS'))

    # SST: semantically same data — water masking applied in make_netcdf step
    sst_day   = lst_day.rename('MODIS_SST_Day_K')
    sst_night = lst_night.rename('MODIS_SST_Night_K')
    sst_mean  = lst_mean.rename('SST_MODIS')

    # --- Quality flags (no scaling) ---
    qc_day   = image.select('QC_Day').rename('MODIS_QC_Day')
    qc_night = image.select('QC_Night').rename('MODIS_QC_Night')

    # --- Emissivity (dimensionless) ---
    emis_31 = image.select('Emis_31').multiply(0.002).add(0.49).rename('MODIS_Emis_31')
    emis_32 = image.select('Emis_32').multiply(0.002).add(0.49).rename('MODIS_Emis_32')

    # --- View geometry ---
    # View angles: stored as uint8 with offset=-65, range [-65, 65] degrees
    view_angle_day   = image.select('Day_view_angle').subtract(65).rename('MODIS_View_angle_Day')
    view_angle_night = image.select('Night_view_angle').subtract(65).rename('MODIS_View_angle_Night')

    # View times: scale=0.1, range [0, 24] hours local solar time
    view_time_day   = image.select('Day_view_time').multiply(0.1).rename('MODIS_View_time_Day')
    view_time_night = image.select('Night_view_time').multiply(0.1).rename('MODIS_View_time_Night')

    return (lst_day
            .addBands(lst_night)
            .addBands(lst_mean)
            .addBands(sst_day)
            .addBands(sst_night)
            .addBands(sst_mean)
            .addBands(qc_day)
            .addBands(qc_night)
            .addBands(emis_31)
            .addBands(emis_32)
            .addBands(view_angle_day)
            .addBands(view_angle_night)
            .addBands(view_time_day)
            .addBands(view_time_night)
            .copyProperties(image, ['system:time_start']))

# ============================================
# LOAD MODIS COLLECTIONS
# ============================================

logger.info("Loading MODIS collections...")

terra = (ee.ImageCollection(config['modis']['collections']['terra'])
         .filterBounds(geometry)
         .filterDate(start_date, end_date)
         .map(process_modis))

aqua = (ee.ImageCollection(config['modis']['collections']['aqua'])
        .filterBounds(geometry)
        .filterDate(start_date, end_date)
        .map(process_modis))

modis_all = terra.merge(aqua)

terra_count = terra.size().getInfo()
aqua_count  = aqua.size().getInfo()
logger.info(f"Terra (MOD11A1): {terra_count} images")
logger.info(f"Aqua  (MYD11A1): {aqua_count} images")
logger.info(f"Total merged   : {modis_all.size().getInfo()} images")

# ============================================
# CREATE COMPOSITES
# ============================================

def create_composites(collection, start, end, interval):
    composites = []
    current = datetime.fromisoformat(start)
    end_dt = datetime.fromisoformat(end)
    comp_num = 0
    while current < end_dt:
        start_str = current.strftime('%Y-%m-%d')
        end_str = (current + timedelta(days=interval)).strftime('%Y-%m-%d')
        period_coll = collection.filterDate(start_str, end_str)
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
logger.info(f"Created {num_composites} composites")

# ============================================
# ADD COORDINATES AND DYNAMIC LAYERS
# ============================================

def add_coordinates_and_dynamic_layers(image, date_str):
    latlon = ee.Image.pixelLonLat().reproject(crs='EPSG:4326', scale=TARGET_RESOLUTION)
    time_img = ee.Image.constant(
        ee.Number(image.get('system:time_start')).divide(86400000)
    ).rename('time_days').float()

    try:
        land_cover = get_land_cover_for_date(date_str)
        land_cover = land_cover.reproject(crs='EPSG:4326', scale=TARGET_RESOLUTION)
    except Exception as e:
        logger.warning(f"Land cover unavailable for {date_str}: {e}")
        land_cover = ee.Image.constant(0).rename('land_cover')

    try:
        emissivity_dyn = calculate_emissivity_from_lc(land_cover)
    except Exception as e:
        logger.warning(f"Emissivity failed for {date_str}: {e}")
        emissivity_dyn = ee.Image.constant(0.95).rename('emissivity_dynamic')

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
    max_retries = config['download']['max_retries']
    timeout = config['download']['timeout_seconds']
    for attempt in range(max_retries):
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
            r = requests.get(url_safe, timeout=timeout)
            if r.status_code == 200:
                with open(filename, 'wb') as f:
                    f.write(r.content)
                if filename.exists() and filename.stat().st_size > 1000:
                    return True
            else:
                logger.error(f"HTTP {r.status_code}: {r.text[:200]}")
            return False
        except (requests.Timeout, requests.ConnectionError) as e:
            logger.warning(f"Network error attempt {attempt+1}/{max_retries}: {type(e).__name__}")
            if attempt < max_retries - 1:
                time.sleep(5 * (attempt + 1))
                continue
            return False
        except Exception as e:
            logger.error(f"Tile error: {str(e)[:200]}")
            return False
    return False

def mosaic_tiles(tile_files, output_file):
    valid_tiles = [f for f in tile_files if f.exists() and f.stat().st_size > 1000]
    if not valid_tiles:
        logger.error("No valid tiles to mosaic")
        return False
    try:
        vrt_file = output_file.parent / f"{output_file.stem}_temp.vrt"
        result = subprocess.run(
            ['gdalbuildvrt', str(vrt_file)] + [str(f) for f in valid_tiles],
            capture_output=True, text=True
        )
        if result.returncode != 0:
            logger.error(f"VRT failed: {result.stderr}")
            return False
        result = subprocess.run([
            'gdal_translate',
            '-co', 'COMPRESS=DEFLATE', '-co', 'PREDICTOR=2',
            '-co', 'TILED=YES', '-co', 'BLOCKXSIZE=256', '-co', 'BLOCKYSIZE=256',
            str(vrt_file), str(output_file)
        ], capture_output=True, text=True)
        if result.returncode != 0:
            logger.error(f"gdal_translate failed: {result.stderr}")
            return False
        if vrt_file.exists():
            vrt_file.unlink()
        return output_file.exists() and output_file.stat().st_size > 1000
    except Exception as e:
        logger.error(f"Mosaic error: {e}")
        return False

# ============================================
# BATCH DOWNLOAD
# ============================================

composites_list = modis_composites.toList(num_composites)

logger.info("="*70)
logger.info(f"Starting MODIS download: {num_composites} composites x {len(tiles)} tiles")
logger.info("="*70)

stats = {'downloaded': 0, 'skipped': 0, 'errors': 0, 'tiles_downloaded': 0}

for i in range(num_composites):
    composite = ee.Image(composites_list.get(i))
    date_str = composite.get('date_string').getInfo()

    final_file = modis_dir / f"MODIS_{date_str}_{i:04d}.tif"
    progress_pct = (i + 1) / num_composites * 100
    logger.info(f"[{i+1:3d}/{num_composites}] ({progress_pct:5.1f}%) {date_str}")

    if final_file.exists():
        logger.info("  SKIPPED (exists)")
        stats['skipped'] += 1
        continue

    composite_final = add_coordinates_and_dynamic_layers(composite, date_str)

    try:
        band_names = composite_final.bandNames().getInfo()
        logger.info(f"  Bands ({len(band_names)}): {', '.join(band_names)}")
    except Exception:
        pass

    tile_files = []
    tiles_success = 0

    for tile in tqdm(tiles, desc="  Tiles", unit="tile", leave=False, ncols=80):
        tile_file = tiles_dir / f"{date_str.replace('-','')}_c{i:03d}_t{tile['id']:02d}.tif"
        if download_tile(composite_final, tile['bounds'], tile_file):
            tile_files.append(tile_file)
            tiles_success += 1
            stats['tiles_downloaded'] += 1
        else:
            logger.warning(f"  Tile {tile['id']} failed")
        time.sleep(config['download']['delay_between_tiles'])

    logger.info(f"  Tiles: {tiles_success}/{len(tiles)} ({tiles_success/len(tiles)*100:.1f}%)")

    min_required = len(tiles) * config['download']['min_tiles_success_ratio']
    if tiles_success < min_required:
        logger.error(f"  FAILED: {tiles_success} < {min_required:.0f} required")
        stats['errors'] += 1
        for tf in tile_files:
            if tf.exists():
                tf.unlink()
        continue

    if mosaic_tiles(tile_files, final_file):
        size_mb = final_file.stat().st_size / (1024**2)
        logger.info(f"  SUCCESS: {size_mb:.1f} MB")
        stats['downloaded'] += 1
        for tf in tile_files:
            if tf.exists():
                tf.unlink()
    else:
        logger.error("  FAILED: mosaic error")
        stats['errors'] += 1

    time.sleep(config['download']['delay_between_images'])

# Cleanup
import shutil
if tiles_dir.exists():
    try:
        shutil.rmtree(tiles_dir)
        logger.info(f"Cleaned up: {tiles_dir}")
    except Exception as e:
        logger.warning(f"Could not remove tiles dir: {e}")

# ============================================
# SUMMARY
# ============================================

logger.info("\n" + "="*70)
logger.info("MODIS DOWNLOAD SUMMARY")
logger.info("="*70)
logger.info(f"Downloaded : {stats['downloaded']}")
logger.info(f"Skipped    : {stats['skipped']}")
logger.info(f"Errors     : {stats['errors']}")
logger.info(f"Tiles total: {stats['tiles_downloaded']}")
logger.info("="*70)
logger.info("Band layout per file:")
logger.info("  MODIS_LST_Day_K, MODIS_LST_Night_K  [Kelvin]")
logger.info("  LST_MODIS  [day+night mean, Kelvin]")
logger.info("  MODIS_SST_Day_K, MODIS_SST_Night_K  [Kelvin, water pixels only in netcdf]")
logger.info("  SST_MODIS  [day+night mean, Kelvin, water pixels only in netcdf]")
logger.info("  MODIS_QC_Day, MODIS_QC_Night  [bitmask]")
logger.info("  MODIS_Emis_31, MODIS_Emis_32  [dimensionless]")
logger.info("  MODIS_View_angle_Day/Night  [degrees, -65 to +65]")
logger.info("  MODIS_View_time_Day/Night   [hours, local solar time]")
logger.info("  lon, lat, time_days, land_cover, emissivity_dynamic")
logger.info(f"  Total: 19 bands per file")
logger.info("="*70)