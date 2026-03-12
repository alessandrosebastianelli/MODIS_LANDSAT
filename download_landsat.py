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
from tqdm import tqdm
import argparse
import sys
import subprocess

# ============================================
# COMMAND LINE ARGUMENTS
# ============================================

parser = argparse.ArgumentParser(description='Download Landsat data')
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
    
    if year_suffix:
        log_file = log_file.replace('.log', f'_{year_suffix}.log')
    
    log_file = Path('logs') / Path(log_file).name
    log_file.parent.mkdir(parents=True, exist_ok=True)
    
    logger = logging.getLogger('Landsat')
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
    logger.info(f"Filtering to range: {start_date} to {end_date}")

logger.info(f"Date range: {start_date} to {end_date}")

base_dir = Path(config['output']['base_dir'])
landsat_dir = base_dir / config['output']['subdirs']['landsat']
tiles_dir = landsat_dir / f'tiles_temp_{args.year}' if args.year else landsat_dir / 'tiles_temp'

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
# LAND COVER AND EMISSIVITY
# ============================================

def get_land_cover_for_date(date_str):
    """Get land cover for a specific date - ESA WorldCover (2020-2021) or MODIS fallback"""
    year = datetime.fromisoformat(date_str).year
    
    if year == 2021:
        lc_image = ee.Image('ESA/WorldCover/v200/2021').select('Map').rename('land_cover')
        return lc_image
    elif year == 2020:
        lc_image = ee.Image('ESA/WorldCover/v100/2020').select('Map').rename('land_cover')
        return lc_image
    else:
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
# ADDITIONAL SATELLITE DATA FOR T2m ML
# ============================================

def get_aerosol_data_for_date(date_str):
    """Get aerosol optical depth from MODIS for atmospheric correction"""
    date = datetime.fromisoformat(date_str)
    
    aod_collection = ee.ImageCollection('MODIS/061/MCD19A2_GRANULES')
    
    start_date = date.strftime('%Y-%m-%d')
    end_date = (date + timedelta(days=15)).strftime('%Y-%m-%d')
    
    aod_data = aod_collection.filterDate(start_date, end_date).mean()
    aod = aod_data.select(['Optical_Depth_047', 'Optical_Depth_055']).rename(['AOD_047', 'AOD_055'])
    
    return aod

def get_precipitation_for_date(date_str):
    """Get precipitation data from GPM for soil moisture proxy"""
    date = datetime.fromisoformat(date_str)
    
    gpm_collection = ee.ImageCollection('NASA/GPM_L3/IMERG_V07')
    
    start_date = date.strftime('%Y-%m-%d')
    end_date = (date + timedelta(days=15)).strftime('%Y-%m-%d')
    
    precip_data = gpm_collection.filterDate(start_date, end_date).sum()
    precip = precip_data.select('precipitationCal').rename('GPM_precipitation')
    
    return precip

# ============================================
# LANDSAT BAND HARMONIZATION
# ============================================

def harmonize_landsat(image):
    """Harmonize Landsat 5/7 bands to match Landsat 8/9 naming and cast to Float32"""
    satellite = ee.String(image.get('SPACECRAFT_ID'))
    
    def rename_l57_bands(img):
        sr_bands = img.select(
            ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7'],
            ['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7']
        ).toFloat()
        
        sr_b1 = img.select('SR_B1').multiply(0).rename('SR_B1').toFloat()
        thermal = img.select('ST_B6').rename('ST_B10').toFloat()
        qa = img.select('QA_PIXEL')
        
        return sr_b1.addBands(sr_bands).addBands(thermal).addBands(qa).copyProperties(img, img.propertyNames())
    
    def keep_l89_bands(img):
        sr_thermal = img.select(['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7', 'ST_B10']).toFloat()
        qa = img.select('QA_PIXEL')
        return sr_thermal.addBands(qa).copyProperties(img, img.propertyNames())
    
    is_l57 = satellite.compareTo(ee.String('LANDSAT_8')).lt(0)
    
    return ee.Algorithms.If(
        is_l57,
        rename_l57_bands(image),
        keep_l89_bands(image)
    )

# ============================================
# LST/SST CALCULATION
# ============================================

def calculate_lst_sst(image):
    """Calculate LST, SST, and spectral indices from Landsat"""
    
    # Spectral Indices
    ndvi = image.normalizedDifference(['SR_B5', 'SR_B4']).rename('NDVI')
    ndwi = image.normalizedDifference(['SR_B3', 'SR_B5']).rename('NDWI')
    ndbi = image.normalizedDifference(['SR_B6', 'SR_B5']).rename('NDBI')
    mndwi = image.normalizedDifference(['SR_B3', 'SR_B6']).rename('MNDWI')
    
    nir   = image.select('SR_B5')
    red   = image.select('SR_B4')
    blue  = image.select('SR_B2')
    swir1 = image.select('SR_B6')
    swir2 = image.select('SR_B7')
    
    savi = nir.subtract(red).divide(nir.add(red).add(0.5)).multiply(1.5).rename('SAVI')
    
    evi = nir.subtract(red).divide(
        nir.add(red.multiply(6)).subtract(blue.multiply(7.5)).add(1.0)
    ).multiply(2.5).rename('EVI')
    
    bsi = swir1.add(red).subtract(nir.add(blue)).divide(
        swir1.add(red).add(nir.add(blue))
    ).rename('BSI')
    
    ui = swir2.subtract(nir).divide(swir2.add(nir)).rename('UI')
    
    albedo = (blue.multiply(0.356)
              .add(red.multiply(0.130))
              .add(nir.multiply(0.373))
              .add(swir1.multiply(0.085))
              .add(swir2.multiply(0.072))
              .subtract(0.0018)).rename('Albedo')
    
    # Fractional Vegetation (FIX: removed duplicate definition)
    fv = ndvi.subtract(0.2).divide(0.5).pow(2).clamp(0, 1).rename('FV')
    
    # Emissivity
    emissivity = fv.multiply(0.004).add(0.986).rename('emissivity')
    
    # Brightness Temperature (DN -> Kelvin, Collection 2 scale factors)
    BT = image.select('ST_B10').multiply(0.00341802).add(149.0)
    
    # LST (Jimenez-Munoz & Sobrino 2003)
    # LST = BT / (1 + (lambda * BT / rho) * ln(emissivity))
    lambda_um = 10.895  # effective wavelength B10 [um]
    rho = 14388.0       # h*c/k_B [um*K]
    eps_safe = emissivity.where(emissivity.lte(0), 0.001).where(emissivity.gt(1), 1)
    # FIX: corrected formula - was BT/(par*ln(eps)+1), should be BT/(1+(lambda*BT/rho)*ln(eps))
    lst = BT.divide(
        BT.multiply(lambda_um / rho).multiply(eps_safe.log()).add(1)
    ).rename('LST_K')
    
    # SST (same as LST, masked to water later)
    sst = lst.rename('SST_K')
    
    # Cloud mask
    qa = image.select('QA_PIXEL')
    cloud_bit  = 1 << 3
    shadow_bit = 1 << 4
    snow_bit   = 1 << 5
    cloud_mask = (qa.bitwiseAnd(cloud_bit).eq(0)
                  .And(qa.bitwiseAnd(shadow_bit).eq(0))
                  .And(qa.bitwiseAnd(snow_bit).eq(0))).rename('cloud_mask').uint8()
    
    lst_masked = lst.updateMask(cloud_mask).rename('LST_K_masked')
    
    return (image
            .addBands(lst)
            .addBands(lst_masked)
            .addBands(sst)
            .addBands(cloud_mask)
            .addBands(ndvi)
            .addBands(ndwi)
            .addBands(ndbi)
            .addBands(mndwi)
            .addBands(savi)
            .addBands(evi)
            .addBands(bsi)
            .addBands(ui)
            .addBands(albedo)
            .addBands(fv)
            .addBands(emissivity))

# ============================================
# LOAD LANDSAT COLLECTIONS
# ============================================

def load_landsat_sensor(sensor_config, date_start, date_end):
    """Load Landsat collection for one sensor with date filtering"""
    sensor_start = max(datetime.fromisoformat(sensor_config['start_date']),
                      datetime.fromisoformat(date_start))
    sensor_end = min(datetime.fromisoformat(sensor_config['end_date']),
                    datetime.fromisoformat(date_end))
    
    if sensor_start > sensor_end:
        return ee.ImageCollection([])
    
    collection_all = (ee.ImageCollection(sensor_config['collection'])
                      .filterBounds(geometry)
                      .filterDate(sensor_start.isoformat(), sensor_end.isoformat()))
    
    if config['landsat']['path_row_filter']['enabled']:
        paths = config['landsat']['path_row_filter']['paths']
        rows = config['landsat']['path_row_filter']['rows']
        
        if paths:
            collection_all = collection_all.filter(ee.Filter.inList('WRS_PATH', paths))
        if rows:
            collection_all = collection_all.filter(ee.Filter.inList('WRS_ROW', rows))
    
    collection_clean = collection_all.filter(ee.Filter.lt('CLOUD_COVER', config['landsat']['max_cloud_cover']))
    
    count_clean = collection_clean.size().getInfo()
    if count_clean < 20:
        relaxed_threshold = min(config['landsat']['max_cloud_cover'] * 2, 95)
        collection = collection_all.filter(ee.Filter.lt('CLOUD_COVER', relaxed_threshold))
    else:
        collection = collection_clean
    
    return collection.map(lambda img: calculate_lst_sst(ee.Image(harmonize_landsat(img))))

logger.info("Loading Landsat collections...")

all_landsat = ee.ImageCollection([])
for sensor in config['landsat']['sensors']:
    coll = load_landsat_sensor(sensor, start_date, end_date)
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
        
        if config['compositing']['method'] == 'median':
            composite = period_coll.median()
        elif config['compositing']['method'] == 'mean':
            composite = period_coll.mean()
        else:
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
# ADD COORDINATES AND DYNAMIC LAYERS
# ============================================

def add_coordinates_and_dynamic_layers(image, date_str):
    """Add lat/lon/time + land cover + emissivity + pure satellite data (no reanalysis)"""
    latlon = ee.Image.pixelLonLat().reproject(crs='EPSG:4326', scale=TARGET_RESOLUTION)
    
    time_img = ee.Image.constant(
        ee.Number(image.get('system:time_start')).divide(86400000)
    ).rename('time').float()
    
    try:
        land_cover = get_land_cover_for_date(date_str)
        land_cover = land_cover.reproject(crs='EPSG:4326', scale=TARGET_RESOLUTION)
    except Exception as e:
        logger.warning(f"Land cover unavailable for {date_str}: {e}")
        land_cover = ee.Image.constant(0).rename('land_cover')
    
    try:
        emissivity_dyn = calculate_emissivity_from_lc(land_cover)
    except Exception as e:
        logger.warning(f"Emissivity calculation failed for {date_str}: {e}")
        emissivity_dyn = ee.Image.constant(0.95).rename('emissivity_dynamic')
    
    try:
        aod_data = get_aerosol_data_for_date(date_str)
        aod_data = aod_data.reproject(crs='EPSG:4326', scale=TARGET_RESOLUTION)
    except Exception as e:
        logger.debug(f"Aerosol data unavailable for {date_str}")
        aod_data = ee.Image.constant([0, 0]).rename(['AOD_047', 'AOD_055'])
    
    return (image
            .addBands(latlon.select('longitude').rename('lon'))
            .addBands(latlon.select('latitude').rename('lat'))
            .addBands(time_img)
            .addBands(land_cover)
            .addBands(emissivity_dyn)
            .addBands(aod_data))

# ============================================
# DOWNLOAD FUNCTIONS
# ============================================

def download_tile(image, tile_bounds, filename):
    max_retries = 3
    timeout = 600
    
    for attempt in range(max_retries):
        try:
            logger.debug(f"    Attempt {attempt+1}/{max_retries}: Preparing image...")
            img_prepared = (image
                            .clip(tile_bounds)
                            .reproject(crs='EPSG:4326', scale=TARGET_RESOLUTION))
            
            logger.debug(f"    Getting download URL...")
            url = img_prepared.getDownloadURL({
                'scale': TARGET_RESOLUTION,
                'crs': 'EPSG:4326',
                'region': tile_bounds,
                'format': 'GEO_TIFF',
                'filePerBand': False
            })
            
            url_safe = urllib.parse.quote(url, safe=':/?&=')
            logger.debug(f"    Downloading tile data...")
            r = requests.get(url_safe, timeout=timeout)
            
            logger.debug(f"    HTTP Status: {r.status_code}")
            if r.status_code == 200:
                with open(filename, 'wb') as f:
                    f.write(r.content)
                if filename.exists() and filename.stat().st_size > 1000:
                    logger.debug(f"    File size: {filename.stat().st_size} bytes")
                    return True
                else:
                    logger.warning(f"    File too small: {filename.stat().st_size if filename.exists() else 0} bytes")
            else:
                logger.error(f"    HTTP {r.status_code}: {r.text[:200]}")
            return False
            
        except (requests.Timeout, requests.ConnectionError) as e:
            logger.warning(f"    Network error (attempt {attempt+1}/{max_retries}): {type(e).__name__}")
            if attempt < max_retries - 1:
                wait_time = 5 * (attempt + 1)
                logger.debug(f"    Retrying in {wait_time}s...")
                time.sleep(wait_time)
                continue
            logger.error(f"Tile timeout after {max_retries} attempts: {str(e)}")
            return False
        except Exception as e:
            import traceback
            error_msg = f"Tile error: {type(e).__name__}: {str(e)}"
            logger.error(error_msg[:300])
            logger.debug(f"Full traceback:\n{traceback.format_exc()}")
            return False
    
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

composites_list = landsat_composites.toList(num_composites)

# Auxiliary bands (not Landsat data) - used to detect empty composites
AUX_BANDS = {'lon', 'lat', 'time', 'land_cover', 'emissivity_dynamic', 'AOD_047', 'AOD_055'}

logger.info("="*70)
logger.info(f"Starting Landsat download: {num_composites} composites")
logger.info(f"Estimated tiles per composite: {len(tiles)}")
logger.info(f"Total downloads: {num_composites * len(tiles)} tiles")
logger.info(f"Estimated time: {num_composites * len(tiles) * config['download']['delay_between_tiles'] / 60:.1f} min (delays only)")
logger.info("="*70)

stats = {'downloaded': 0, 'skipped': 0, 'errors': 0, 'tiles_downloaded': 0}

for i in range(num_composites):
    composite = ee.Image(composites_list.get(i))
    date_str = composite.get('date_string').getInfo()
    
    final_file = landsat_dir / f"Landsat_{date_str}_{i:04d}.tif"
    
    progress_pct = (i + 1) / num_composites * 100
    logger.info(f"[{i+1:3d}/{num_composites}] ({progress_pct:5.1f}%) Composite: {date_str}")
    
    if final_file.exists():
        logger.info("Status: SKIPPED (file exists)")
        stats['skipped'] += 1
        continue
    
    composite_final = add_coordinates_and_dynamic_layers(composite, date_str)
    
    try:
        band_names = composite_final.bandNames().getInfo()
        logger.info(f"Composite has {len(band_names)} bands: {', '.join(band_names[:10])}{'...' if len(band_names) > 10 else ''}")
        
        # Skip composites with no Landsat data (empty period)
        landsat_bands = [b for b in band_names if b not in AUX_BANDS]
        if not landsat_bands:
            logger.warning(f"Status: SKIPPED (no Landsat data in period {date_str})")
            stats['skipped'] += 1
            continue
    except Exception as e:
        logger.warning(f"Cannot retrieve band names: {e}")
    
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
            file_size = tile_file.stat().st_size / 1024
            logger.debug(f"    SUCCESS: {file_size:.1f} KB")
        else:
            logger.warning(f"    FAILED: Tile {tile_idx} download error")
        
        time.sleep(config['download']['delay_between_tiles'])
    
    logger.info(f"Tiles completed: {tiles_success}/{len(tiles)} successful ({tiles_success/len(tiles)*100:.1f}%)")
    
    min_required = len(tiles) * config['download']['min_tiles_success_ratio']
    if tiles_success < min_required:
        logger.error(f"Status: FAILED - Insufficient tiles ({tiles_success} < {min_required:.0f} required)")
        stats['errors'] += 1
        for tile_file in tile_files:
            if tile_file.exists():
                tile_file.unlink()
        continue
    
    logger.info(f"Mosaicking {tiles_success} tiles")
    if mosaic_tiles(tile_files, final_file):
        file_size_mb = final_file.stat().st_size / (1024**2)
        logger.info(f"Status: SUCCESS - {file_size_mb:.1f} MB")
        stats['downloaded'] += 1
        for tile_file in tile_files:
            if tile_file.exists():
                tile_file.unlink()
    else:
        logger.error("Status: FAILED - Mosaic error")
        stats['errors'] += 1
    
    time.sleep(config['download']['delay_between_images'])

try:
    tiles_dir.rmdir()
except:
    logger.warning("Could not remove tiles directory (may contain files)")

# ============================================
# SUMMARY
# ============================================

logger.info("\n" + "="*70)
logger.info("LANDSAT DOWNLOAD SUMMARY")
logger.info("="*70)
logger.info(f"Downloaded: {stats['downloaded']}")
logger.info(f"Skipped: {stats['skipped']}")
logger.info(f"Errors: {stats['errors']}")
logger.info(f"Total tiles downloaded: {stats['tiles_downloaded']}")
logger.info("="*70)

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
logger.info("    NDVI, NDWI, NDBI, MNDWI, SAVI, EVI, BSI, UI, Albedo, FV, emissivity")
logger.info(f"  Dynamic layers ({interval_days}-day frequency):")
logger.info("    land_cover, emissivity_dynamic")
logger.info("  Coordinates: lon, lat, time")
logger.info("  Aerosol: AOD_047, AOD_055")
logger.info(f"  Total: ~26 bands per file")
logger.info("="*70)

logger.info("="*70)
logger.info("Download Summary")
logger.info("="*70)
logger.info(f"Composites processed: {stats['downloaded'] + stats['skipped'] + stats['errors']}")
logger.info(f"  Downloaded: {stats['downloaded']}")
logger.info(f"  Skipped (existing): {stats['skipped']}")
logger.info(f"  Errors: {stats['errors']}")
logger.info(f"Total tiles downloaded: {stats['tiles_downloaded']}")
if stats['downloaded'] > 0:
    logger.info(f"Average tiles per composite: {stats['tiles_downloaded'] / stats['downloaded']:.1f}")
logger.info("="*70)

import shutil
if tiles_dir.exists():
    try:
        shutil.rmtree(tiles_dir)
        logger.info(f"Cleaned up temporary tiles: {tiles_dir}")
    except Exception as e:
        logger.warning(f"Could not remove tiles directory: {e}")