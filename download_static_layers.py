# ============================================
# Download Static Layers
# DEM + Land Cover + Water Mask + Emissivity
# ============================================

import ee
from pathlib import Path
import requests
import urllib.parse
import yaml
import logging
import rasterio
from rasterio.merge import merge
import numpy as np

# ============================================
# LOGGING SETUP
# ============================================

def setup_logging(config):
    log_level = getattr(logging, config['logging']['level'])
    log_file = config['logging']['log_file']
    
    logger = logging.getLogger('Static_Layers')
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

TARGET_RESOLUTION = config['spatial']['target_resolution']
TILE_SIZE = config['spatial']['tile_size']

base_dir = Path(config['output']['base_dir'])
static_dir = base_dir / config['output']['subdirs']['static']
tiles_dir = static_dir / 'tiles_temp'

static_dir.mkdir(parents=True, exist_ok=True)
tiles_dir.mkdir(parents=True, exist_ok=True)

logger.info("="*70)
logger.info("Static Layers Download - Rome Area")
logger.info("="*70)
logger.info(f"Resolution: {TARGET_RESOLUTION}m")
logger.info(f"Output: {static_dir}")

# ============================================
# TILING FUNCTIONS
# ============================================

def create_tiles(bounds, tile_size):
    """Create tile grid"""
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

def download_tile(image, tile_bounds, filename):
    """Download single tile"""
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
        r = requests.get(url_safe, timeout=300)
        
        if r.status_code == 200:
            with open(filename, 'wb') as f:
                f.write(r.content)
            if filename.exists() and filename.stat().st_size > 1000:
                return True
        return False
    except Exception as e:
        logger.error(f"Tile download error: {str(e)[:80]}")
        return False

def mosaic_tiles(tile_files, output_file):
    """Mosaic tiles into single file"""
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
# 1. DOWNLOAD DEM
# ============================================

if config['static_layers']['dem']['enabled']:
    logger.info("\n" + "="*70)
    logger.info("1. DIGITAL ELEVATION MODEL (SRTM 30m)")
    logger.info("="*70)
    
    dem_file = static_dir / 'DEM.tif'
    
    if dem_file.exists():
        logger.info("DEM already exists, skipping")
    else:
        dem = ee.Image(config['static_layers']['dem']['collection'])
        dem = dem.select(config['static_layers']['dem']['bands'])
        
        # Download in tiles
        tile_files = []
        for tile in tiles:
            tile_file = tiles_dir / f"dem_t{tile['id']:02d}.tif"
            if download_tile(dem, tile['bounds'], tile_file):
                tile_files.append(tile_file)
        
        if mosaic_tiles(tile_files, dem_file):
            logger.info(f"✓ DEM saved: {dem_file.name}")
            # Cleanup
            for f in tile_files:
                if f.exists():
                    f.unlink()
        else:
            logger.error("✗ DEM mosaic failed")

# ============================================
# 2. DOWNLOAD LAND COVER
# ============================================

if config['static_layers']['land_cover']['enabled']:
    logger.info("\n" + "="*70)
    logger.info("2. LAND COVER (ESA WorldCover)")
    logger.info("="*70)
    
    lc_file = static_dir / 'LandCover.tif'
    
    if lc_file.exists():
        logger.info("Land Cover already exists, skipping")
    else:
        lc = ee.ImageCollection(config['static_layers']['land_cover']['collection'])
        lc = lc.first().select(config['static_layers']['land_cover']['bands'])
        
        tile_files = []
        for tile in tiles:
            tile_file = tiles_dir / f"lc_t{tile['id']:02d}.tif"
            if download_tile(lc, tile['bounds'], tile_file):
                tile_files.append(tile_file)
        
        if mosaic_tiles(tile_files, lc_file):
            logger.info(f"✓ Land Cover saved: {lc_file.name}")
            for f in tile_files:
                if f.exists():
                    f.unlink()
        else:
            logger.error("✗ Land Cover mosaic failed")

# ============================================
# 3. CREATE WATER MASK
# ============================================

if config['static_layers']['water_mask']['enabled']:
    logger.info("\n" + "="*70)
    logger.info("3. WATER MASK (Combined)")
    logger.info("="*70)
    
    water_file = static_dir / 'WaterMask.tif'
    
    if water_file.exists():
        logger.info("Water Mask already exists, skipping")
    else:
        # Source 1: JRC Global Surface Water
        gsw = ee.Image('JRC/GSW1_4/GlobalSurfaceWater')
        water_occurrence = gsw.select('occurrence').gt(50)  # >50% water occurrence
        
        # Source 2: ESA WorldCover water classes
        lc = ee.ImageCollection('ESA/WorldCover/v200').first()
        water_classes = lc.select('Map')
        water_from_lc = (water_classes.eq(80)  # Permanent water
                        .Or(water_classes.eq(90))  # Wetlands
                        .Or(water_classes.eq(95)))  # Mangroves
        
        # Combined water mask (binary: 1=water, 0=land)
        water_mask = water_occurrence.Or(water_from_lc).rename('water_mask')
        
        tile_files = []
        for tile in tiles:
            tile_file = tiles_dir / f"water_t{tile['id']:02d}.tif"
            if download_tile(water_mask, tile['bounds'], tile_file):
                tile_files.append(tile_file)
        
        if mosaic_tiles(tile_files, water_file):
            logger.info(f"✓ Water Mask saved: {water_file.name}")
            for f in tile_files:
                if f.exists():
                    f.unlink()
        else:
            logger.error("✗ Water Mask mosaic failed")

# ============================================
# 4. CREATE EMISSIVITY MAP
# ============================================

if config['static_layers']['emissivity']['enabled']:
    logger.info("\n" + "="*70)
    logger.info("4. EMISSIVITY (from Land Cover)")
    logger.info("="*70)
    
    emis_file = static_dir / 'Emissivity.tif'
    
    if emis_file.exists():
        logger.info("Emissivity already exists, skipping")
    else:
        lc = ee.ImageCollection('ESA/WorldCover/v200').first().select('Map')
        
        # Create emissivity map from lookup table
        emis_lookup = config['static_layers']['emissivity']['lookup']
        emis = ee.Image(0.0)
        
        for lc_class, emis_value in emis_lookup.items():
            mask = lc.eq(int(lc_class))
            emis = emis.where(mask, emis_value)
        
        emis = emis.rename('emissivity')
        
        tile_files = []
        for tile in tiles:
            tile_file = tiles_dir / f"emis_t{tile['id']:02d}.tif"
            if download_tile(emis, tile['bounds'], tile_file):
                tile_files.append(tile_file)
        
        if mosaic_tiles(tile_files, emis_file):
            logger.info(f"✓ Emissivity saved: {emis_file.name}")
            for f in tile_files:
                if f.exists():
                    f.unlink()
        else:
            logger.error("✗ Emissivity mosaic failed")

# Cleanup tiles directory
try:
    tiles_dir.rmdir()
except:
    pass

# ============================================
# SUMMARY
# ============================================

logger.info("\n" + "="*70)
logger.info("STATIC LAYERS DOWNLOAD COMPLETE")
logger.info("="*70)

files_created = list(static_dir.glob('*.tif'))
logger.info(f"Files created: {len(files_created)}")
for f in files_created:
    size_mb = f.stat().st_size / (1024**2)
    logger.info(f"  • {f.name}: {size_mb:.2f} MB")

logger.info("="*70)
