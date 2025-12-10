# ============================================
# Create Combined NetCDF
# Landsat + MODIS + Static Layers + LST/SST Hybrid
# ============================================

import re
from pathlib import Path
import pandas as pd
import xarray as xr
import rioxarray
import dask
import warnings
import yaml
import logging
import numpy as np

dask.config.set(scheduler='threads')
warnings.filterwarnings('ignore', category=RuntimeWarning)

# ============================================
# LOGGING SETUP
# ============================================

def setup_logging(config):
    log_level = getattr(logging, config['logging']['level'])
    log_file = config['logging']['log_file']
    
    logger = logging.getLogger('NetCDF')
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
# LOAD CONFIG
# ============================================

with open('config_refactored.yaml', 'r') as f:
    config = yaml.safe_load(f)

logger = setup_logging(config)

base_dir = Path(config['output']['base_dir'])
landsat_dir = base_dir / config['output']['subdirs']['landsat']
modis_dir = base_dir / config['output']['subdirs']['modis']
static_dir = base_dir / config['output']['subdirs']['static']
out_nc = config['output']['netcdf_filename']

logger.info("="*70)
logger.info("Creating Combined NetCDF")
logger.info("="*70)
logger.info(f"Landsat directory: {landsat_dir}")
logger.info(f"MODIS directory: {modis_dir}")
logger.info(f"Static directory: {static_dir}")
logger.info(f"Output: {out_nc}")

# ============================================
# DATE EXTRACTION
# ============================================

date_re = re.compile(r"(\d{4}-\d{2}-\d{2})")

def extract_date(stem: str):
    m = date_re.search(stem)
    if m:
        return pd.to_datetime(m.group(1), errors="coerce")
    return None

# ============================================
# LOAD STATIC LAYERS
# ============================================

logger.info("\n" + "="*70)
logger.info("Loading static layers...")
logger.info("="*70)

static_layers = {}

# DEM
dem_file = static_dir / 'DEM.tif'
if dem_file.exists():
    logger.info(f"Loading {dem_file.name}")
    static_layers['dem'] = rioxarray.open_rasterio(dem_file).squeeze()
    logger.info(f"  ✓ DEM loaded")

# Land Cover
lc_file = static_dir / 'LandCover.tif'
if lc_file.exists():
    logger.info(f"Loading {lc_file.name}")
    static_layers['land_cover'] = rioxarray.open_rasterio(lc_file).squeeze()
    logger.info(f"  ✓ Land Cover loaded")

# Water Mask
water_file = static_dir / 'WaterMask.tif'
if water_file.exists():
    logger.info(f"Loading {water_file.name}")
    static_layers['water_mask'] = rioxarray.open_rasterio(water_file).squeeze()
    logger.info(f"  ✓ Water Mask loaded")

# Emissivity
emis_file = static_dir / 'Emissivity.tif'
if emis_file.exists():
    logger.info(f"Loading {emis_file.name}")
    static_layers['emissivity'] = rioxarray.open_rasterio(emis_file).squeeze()
    logger.info(f"  ✓ Emissivity loaded")

logger.info(f"Total static layers: {len(static_layers)}")

# ============================================
# LOAD LANDSAT TIME SERIES
# ============================================

logger.info("\n" + "="*70)
logger.info("Loading Landsat time series...")
logger.info("="*70)

landsat_files = sorted(landsat_dir.glob("*.tif"))
logger.info(f"Found {len(landsat_files)} Landsat files")

landsat_datasets = []

for tif in landsat_files:
    logger.debug(f"Processing: {tif.name}")
    date = extract_date(tif.stem)
    if date is None or pd.isna(date):
        logger.warning(f"Cannot parse date from {tif.name}")
        continue
    
    ds = rioxarray.open_rasterio(tif)
    
    # Drop spatial_ref
    for dim in list(ds.dims):
        if dim not in ("band", "y", "x"):
            ds = ds.squeeze(dim, drop=True)
    
    # Add time dimension
    ds = ds.squeeze(drop=True).expand_dims(time=[pd.Timestamp(date)])
    
    # Keep essential coords
    coords_to_keep = ["time", "x", "y", "band"]
    for c in list(ds.coords):
        if c not in coords_to_keep:
            ds = ds.drop_vars(c)
    
    # Convert to dataset with band names
    if "band" in ds.dims:
        # Landsat expected bands (from download_landsat.py)
        # Original: SR_B1-SR_B7, ST_B10, QA_PIXEL
        # Derived: LST_K, LST_K_masked, SST_K, cloud_mask, NDVI, emissivity
        # Coords: lon, lat, time
        # Total: 9 + 6 + 3 = 18 bands
        
        actual_bands = ds.sizes.get("band", 0)
        logger.debug(f"  Bands: {actual_bands}")
        
        if actual_bands >= 18:
            band_names = [
                'SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7',
                'ST_B10', 'QA_PIXEL',
                'LST_K', 'LST_K_masked', 'SST_K', 'cloud_mask', 'NDVI', 'emissivity',
                'lon', 'lat', 'time'
            ][:actual_bands]
            
            ds = ds.assign_coords(band=band_names)
            ds = ds.to_dataset(dim="band")
            
            # Add prefix to differentiate
            ds = ds.rename({var: f'landsat_{var}' for var in ds.data_vars})
            
            landsat_datasets.append(ds)
        else:
            logger.warning(f"  Unexpected band count: {actual_bands}, skipping")

logger.info(f"Loaded {len(landsat_datasets)} Landsat datasets")

# ============================================
# LOAD MODIS TIME SERIES
# ============================================

logger.info("\n" + "="*70)
logger.info("Loading MODIS time series...")
logger.info("="*70)

modis_files = sorted(modis_dir.glob("*.tif"))
logger.info(f"Found {len(modis_files)} MODIS files")

modis_datasets = []

for tif in modis_files:
    logger.debug(f"Processing: {tif.name}")
    date = extract_date(tif.stem)
    if date is None or pd.isna(date):
        logger.warning(f"Cannot parse date from {tif.name}")
        continue
    
    ds = rioxarray.open_rasterio(tif)
    
    for dim in list(ds.dims):
        if dim not in ("band", "y", "x"):
            ds = ds.squeeze(dim, drop=True)
    
    ds = ds.squeeze(drop=True).expand_dims(time=[pd.Timestamp(date)])
    
    coords_to_keep = ["time", "x", "y", "band"]
    for c in list(ds.coords):
        if c not in coords_to_keep:
            ds = ds.drop_vars(c)
    
    if "band" in ds.dims:
        # MODIS expected bands (from download_modis.py)
        # LST/SST: Day, Night (4)
        # QC: Day, Night (2)
        # Emis: 31, 32 (2)
        # View: angle Day/Night, time Day/Night (4)
        # Coords: lon, lat, time (3)
        # Total: 15 bands
        
        actual_bands = ds.sizes.get("band", 0)
        logger.debug(f"  Bands: {actual_bands}")
        
        if actual_bands >= 15:
            band_names = [
                'MODIS_LST_Day_K', 'MODIS_LST_Night_K',
                'MODIS_SST_Day_K', 'MODIS_SST_Night_K',
                'MODIS_QC_Day', 'MODIS_QC_Night',
                'MODIS_Emis_31', 'MODIS_Emis_32',
                'MODIS_View_angle_Day', 'MODIS_View_angle_Night',
                'MODIS_View_time_Day', 'MODIS_View_time_Night',
                'lon', 'lat', 'time'
            ][:actual_bands]
            
            ds = ds.assign_coords(band=band_names)
            ds = ds.to_dataset(dim="band")
            
            modis_datasets.append(ds)
        else:
            logger.warning(f"  Unexpected band count: {actual_bands}, skipping")

logger.info(f"Loaded {len(modis_datasets)} MODIS datasets")

# ============================================
# COMBINE TIME SERIES
# ============================================

logger.info("\n" + "="*70)
logger.info("Combining time series...")
logger.info("="*70)

if len(landsat_datasets) == 0 and len(modis_datasets) == 0:
    logger.error("No time series data to process")
    exit(1)

# Concatenate Landsat
if len(landsat_datasets) > 0:
    landsat_datasets = sorted(landsat_datasets, key=lambda d: pd.Timestamp(d.time.values[0]))
    landsat_combined = xr.concat(
        landsat_datasets,
        dim="time",
        data_vars="all",
        coords="minimal",
        compat='override',
        join='override'
    )
    logger.info(f"Landsat: {len(landsat_combined.time)} timesteps")
else:
    landsat_combined = None

# Concatenate MODIS
if len(modis_datasets) > 0:
    modis_datasets = sorted(modis_datasets, key=lambda d: pd.Timestamp(d.time.values[0]))
    modis_combined = xr.concat(
        modis_datasets,
        dim="time",
        data_vars="all",
        coords="minimal",
        compat='override',
        join='override'
    )
    logger.info(f"MODIS: {len(modis_combined.time)} timesteps")
else:
    modis_combined = None

# Merge Landsat + MODIS
if landsat_combined is not None and modis_combined is not None:
    combined = xr.merge([landsat_combined, modis_combined], compat='override')
elif landsat_combined is not None:
    combined = landsat_combined
else:
    combined = modis_combined

logger.info(f"Combined time series: {len(combined.time)} timesteps")

# ============================================
# ADD STATIC LAYERS
# ============================================

logger.info("\n" + "="*70)
logger.info("Adding static layers...")
logger.info("="*70)

for name, data in static_layers.items():
    logger.info(f"Adding {name}")
    # Broadcast static layer to all timesteps
    expanded = data.expand_dims(time=combined.time)
    combined[name] = expanded

# ============================================
# CREATE HYBRID LST/SST LAYER
# ============================================

logger.info("\n" + "="*70)
logger.info("Creating hybrid LST/SST layer...")
logger.info("="*70)

if 'water_mask' in combined and 'landsat_LST_K' in combined and 'landsat_SST_K' in combined:
    water_mask = combined['water_mask']
    
    # Hybrid layer: LST over land (water_mask=0), SST over water (water_mask=1)
    hybrid_temp = xr.where(
        water_mask == 1,
        combined['landsat_SST_K'],  # Water pixels: use SST
        combined['landsat_LST_K']   # Land pixels: use LST
    )
    hybrid_temp = hybrid_temp.rename('LST_SST_hybrid')
    combined['LST_SST_hybrid'] = hybrid_temp
    
    logger.info("✓ Created LST_SST_hybrid layer")
    logger.info("  Land pixels: LST from Landsat")
    logger.info("  Water pixels: SST from Landsat")

# ============================================
# TIME FORMATTING
# ============================================

combined["time"] = pd.to_datetime(combined.time.values)
combined["time"].encoding.update({
    "units": "days since 1970-01-01 00:00:00",
    "calendar": "gregorian"
})

# ============================================
# METADATA
# ============================================

logger.info("\n" + "="*70)
logger.info("Adding metadata...")
logger.info("="*70)

# Global attributes
combined.attrs.update({
    "title": "Multi-Sensor LST/SST Dataset - Rome, Italy",
    "institution": "CMCC Foundation",
    "source": "Landsat 8/9, MODIS Terra/Aqua, SRTM, ESA WorldCover",
    "history": f"Created on {pd.Timestamp.now().isoformat()}",
    "comment": "15-day composites at 30m resolution with static layers",
    
    # Temporal
    "time_coverage_start": str(combined.time.min().values),
    "time_coverage_end": str(combined.time.max().values),
    "time_coverage_resolution": f"P{config['temporal']['interval_days']}D",
    
    # Spatial
    "geospatial_lat_min": float(config['geometry']['min_lat']),
    "geospatial_lat_max": float(config['geometry']['max_lat']),
    "geospatial_lon_min": float(config['geometry']['min_lon']),
    "geospatial_lon_max": float(config['geometry']['max_lon']),
    "spatial_resolution": f"{config['spatial']['target_resolution']} meters",
    
    # Conventions
    "Conventions": "CF-1.8",
    "creation_date": pd.Timestamp.now().isoformat(),
    "creator_name": "Alessandro Sebastianelli",
    "creator_institution": "CMCC Foundation",
    
    # Data description
    "landsat_sensors": "Landsat 8 (2020-2021), Landsat 9 (2021-2025)",
    "modis_sensors": "Terra (MOD11A1), Aqua (MYD11A1)",
    "static_layers": "DEM (SRTM), Land Cover (ESA WorldCover), Water Mask, Emissivity",
    "hybrid_layer": "LST_SST_hybrid combines LST over land and SST over water"
})

# Variable attributes
var_attrs = {
    'LST_SST_hybrid': {
        'long_name': 'Hybrid Land/Sea Surface Temperature',
        'units': 'K',
        'description': 'LST over land pixels, SST over water pixels',
        'source': 'Landsat 8/9 with water mask applied'
    },
    'water_mask': {
        'long_name': 'Water Mask',
        'description': 'Combined water mask (seas, oceans, rivers, lakes)',
        'flag_values': [0, 1],
        'flag_meanings': 'land water',
        'source': 'JRC Global Surface Water + ESA WorldCover'
    },
    'dem': {
        'long_name': 'Digital Elevation Model',
        'units': 'm',
        'source': 'SRTM 30m'
    },
    'land_cover': {
        'long_name': 'Land Cover Classification',
        'source': 'ESA WorldCover 2021'
    },
    'emissivity': {
        'long_name': 'Land Surface Emissivity',
        'description': 'Derived from land cover classification'
    }
}

for var, attrs in var_attrs.items():
    if var in combined:
        combined[var].attrs.update(attrs)

# ============================================
# SAVE NETCDF
# ============================================

logger.info("\n" + "="*70)
logger.info(f"Saving NetCDF: {out_nc}")
logger.info("="*70)

# Encoding with compression
encoding = {}
for var in combined.data_vars:
    chunk_y = min(256, combined.sizes.get('y', 256))
    chunk_x = min(256, combined.sizes.get('x', 256))
    
    encoding[var] = {
        "zlib": True,
        "complevel": 4,
        "chunksizes": (10, chunk_y, chunk_x)
    }

logger.info("Writing NetCDF (this may take several minutes)...")
with dask.config.set(scheduler='threads', num_workers=2):
    combined.to_netcdf(
        out_nc,
        encoding=encoding,
        mode='w',
        compute=True
    )

file_size_mb = Path(out_nc).stat().st_size / (1024**2)

logger.info("\n" + "="*70)
logger.info("NETCDF CREATION COMPLETE")
logger.info("="*70)
logger.info(f"File: {out_nc}")
logger.info(f"Size: {file_size_mb:.2f} MB")
logger.info(f"Timesteps: {len(combined.time)}")
logger.info(f"Variables: {len(combined.data_vars)}")
logger.info("="*70)

# ============================================
# SUMMARY
# ============================================

logger.info("\nDataset contents:")
logger.info("\n1. Landsat variables (time series):")
for var in sorted([v for v in combined.data_vars if 'landsat_' in v]):
    logger.info(f"  • {var}")

logger.info("\n2. MODIS variables (time series):")
for var in sorted([v for v in combined.data_vars if 'MODIS_' in v]):
    logger.info(f"  • {var}")

logger.info("\n3. Static layers:")
for var in sorted([v for v in combined.data_vars if v in static_layers]):
    logger.info(f"  • {var}")

logger.info("\n4. Hybrid layer:")
if 'LST_SST_hybrid' in combined:
    logger.info("  • LST_SST_hybrid (LST over land, SST over water)")

logger.info("\n" + "="*70)
logger.info("✓ All processing complete!")
logger.info("="*70)
