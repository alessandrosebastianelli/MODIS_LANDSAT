#!/usr/bin/env python3
# ============================================
# Create Yearly NetCDF with Proper Band Names
# ============================================

import re
from pathlib import Path
import pandas as pd
import xarray as xr
import rioxarray
import rasterio
import dask
import warnings
import yaml
import logging
import numpy as np
from collections import defaultdict
import argparse
import sys


# Optimize Dask for MINIMAL memory usage
dask.config.set({
    'scheduler': 'threads',
    'num_workers': 2,  # Few workers to minimize memory
    'threads_per_worker': 1,
    'array.chunk-size': '64MB',  # Small chunks to stay under memory limits
    'distributed.worker.memory.target': 0.70,  # Conservative memory usage
    'distributed.worker.memory.spill': 0.80,
    'distributed.worker.memory.pause': 0.85,
    'distributed.worker.memory.terminate': 0.90
})

warnings.filterwarnings('ignore', category=RuntimeWarning)

# ============================================
# ARGUMENTS
# ============================================

parser = argparse.ArgumentParser(description='Create NetCDF for specific year(s)')
parser.add_argument('--year', type=int, default=None, help='Process specific year only')
parser.add_argument('--start-year', type=int, default=None, help='Start year')
parser.add_argument('--end-year', type=int, default=None, help='End year')
args = parser.parse_args()

# ============================================
# LOGGING
# ============================================

def setup_logging(config, year_suffix=None):
    log_level = getattr(logging, config['logging']['level'])
    log_file = config['logging']['log_file']
    
    # Add year suffix if provided
    if year_suffix:
        log_file = log_file.replace('.log', f'_{year_suffix}.log')
    
    # Move to logs/ directory
    log_file = Path('logs') / Path(log_file).name
    log_file.parent.mkdir(parents=True, exist_ok=True)
    
    logger = logging.getLogger('NetCDF_Yearly')
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
    
    class ColoredFormatter(logging.Formatter):
        COLORS = {'DEBUG': '\033[36m', 'INFO': '\033[32m', 'WARNING': '\033[33m',
                  'ERROR': '\033[31m', 'CRITICAL': '\033[35m'}
        RESET = '\033[0m'
        
        def format(self, record):
            color = self.COLORS.get(record.levelname, self.RESET)
            levelname_colored = f"{color}{record.levelname[0]}{self.RESET}"
            return f"{levelname_colored} | {record.getMessage()}"
    
    ch.setFormatter(ColoredFormatter())
    logger.addHandler(ch)
    return logger

# ============================================
# LOAD CONFIG
# ============================================

with open('config_refactored.yaml', 'r') as f:
    config = yaml.safe_load(f)

logger = setup_logging(config, year_suffix=args.year if args.year else None)

base_dir = Path(config['output']['base_dir'])
landsat_dir = base_dir / config['output']['subdirs']['landsat']
modis_dir = base_dir / config['output']['subdirs']['modis']
static_dir = base_dir / config['output']['subdirs']['static']
out_dir = base_dir / 'netcdf_yearly'
out_dir.mkdir(parents=True, exist_ok=True)

logger.info("="*70)
logger.info("Creating Yearly NetCDF Files with Proper Band Names")
logger.info("="*70)
logger.info(f"Output directory: {out_dir}")

# ============================================
# BAND NAME DEFINITIONS
# ============================================

# Expected Landsat bands in order
LANDSAT_BANDS = [
    'SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7',  # Reflectance
    'ST_B10',  # Thermal
    'QA_PIXEL',  # Quality
    'NDVI', 'NDWI', 'NDBI', 'MNDWI', 'SAVI', 'EVI', 'BSI', 'UI', 'Albedo', 'FV',  # Indices (10)
    'LST_K', 'LST_K_masked', 'SST_K', 'LSSTSST',  # Temperature (4)
    'emissivity', 'emissivity_dynamic',  # Emissivity (2)
    'cloud_mask',  # Cloud (1)
    'lon', 'lat', 'time_days',  # Coordinates (3)
    'land_cover',  # Land cover (1)
    'AOD_047', 'AOD_055'  # Aerosol (2)
]

# Expected MODIS bands
MODIS_BANDS = [
    'LST_Day_K', 'LST_Night_K', 'SST_Day_K', 'SST_Night_K',
    'LST_MODIS', 'SST_MODIS',
    'QC_Day', 'QC_Night',
    'Emis_31', 'Emis_32',
    'View_angle_Day', 'View_angle_Night',
    'View_time_Day', 'View_time_Night'
]

# ============================================
# HELPER FUNCTIONS
# ============================================

date_re = re.compile(r"(\d{4})-(\d{2})-(\d{2})")

def extract_date_and_year(stem: str):
    m = date_re.search(stem)
    if m:
        year = int(m.group(1))
        date = pd.to_datetime(f"{m.group(1)}-{m.group(2)}-{m.group(3)}", errors="coerce")
        return date, year
    return None, None

def read_tiff_with_band_names(tif_file, expected_bands, prefix):
    """Read GeoTIFF with minimal memory footprint using lazy Dask arrays"""
    with rasterio.open(tif_file) as src:
        transform = src.transform
        crs = src.crs
        n_bands = src.count
        shape = (src.height, src.width)
    
    # Create coordinates
    y_coords = np.arange(shape[0]) * transform[4] + transform[5]
    x_coords = np.arange(shape[1]) * transform[0] + transform[2]
    
    # Create dataset with Dask arrays for MINIMAL memory usage
    import dask.array as da
    dataset_dict = {}
    
    # Small chunks to minimize memory: ~64MB per chunk
    # For typical 5000x5000 @ 30m = 25M pixels × 4 bytes = 100MB per band
    # Chunk into 256x256 pieces = ~256KB per chunk
    chunk_size = 256
    
    for i in range(min(n_bands, len(expected_bands))):
        band_name = f"{prefix}_{expected_bands[i]}"
        
        # Use delayed loading - data not read until needed
        def read_band(band_idx, filename):
            with rasterio.open(filename) as src:
                return src.read(band_idx + 1)
        
        # Create lazy Dask array that reads from file only when computed
        lazy_array = da.from_delayed(
            dask.delayed(read_band)(i, tif_file),
            shape=shape,
            dtype=np.float32
        )
        
        # Rechunk to small pieces
        dask_array = lazy_array.rechunk((chunk_size, chunk_size))
        dataset_dict[band_name] = (['y', 'x'], dask_array)
    
    ds = xr.Dataset(
        dataset_dict,
        coords={'y': y_coords, 'x': x_coords}
    )
    
    return ds

# ============================================
# SCAN FILES
# ============================================

logger.info("Scanning files...")

landsat_by_year = defaultdict(list)
for f in sorted(landsat_dir.glob("*.tif")):
    date, year = extract_date_and_year(f.stem)
    if year:
        landsat_by_year[year].append((f, date))

modis_by_year = defaultdict(list)
for f in sorted(modis_dir.glob("*.tif")):
    date, year = extract_date_and_year(f.stem)
    if year:
        modis_by_year[year].append((f, date))

all_years = sorted(set(landsat_by_year.keys()) | set(modis_by_year.keys()))

# Filter years
if args.year:
    all_years = [args.year] if args.year in all_years else []
    if not all_years:
        logger.error(f"No data for year {args.year}")
        sys.exit(1)
elif args.start_year or args.end_year:
    start_y = args.start_year if args.start_year else min(all_years)
    end_y = args.end_year if args.end_year else max(all_years)
    all_years = [y for y in all_years if start_y <= y <= end_y]

logger.info(f"Years to process: {all_years}")

# Load static layers
logger.info("Loading static layers...")
dem_file = static_dir / 'DEM.tif'
water_file = static_dir / 'WaterMask.tif'

static_data = {}
if dem_file.exists():
    try:
        dem_raw = rioxarray.open_rasterio(dem_file)
        # Remove band dimension if present
        if 'band' in dem_raw.dims:
            dem_raw = dem_raw.isel(band=0)
        dem_raw = dem_raw.squeeze()
        static_data['DEM'] = dem_raw
        logger.info(f"Loaded DEM: dims={list(dem_raw.dims)}, shape={dem_raw.shape}")
    except Exception as e:
        logger.error(f"Failed to load DEM: {e}")

if water_file.exists():
    try:
        water_raw = rioxarray.open_rasterio(water_file)
        if 'band' in water_raw.dims:
            water_raw = water_raw.isel(band=0)
        water_raw = water_raw.squeeze()
        static_data['WATERMASK'] = water_raw
        logger.info(f"Loaded WATERMASK: dims={list(water_raw.dims)}, shape={water_raw.shape}")
    except Exception as e:
        logger.error(f"Failed to load WATERMASK: {e}")

if not static_data:
    logger.warning(f"No static files found in {static_dir}")

# ============================================
# PROCESS EACH YEAR
# ============================================

for year in all_years:
    logger.info("="*70)
    logger.info(f"Processing year: {year}")
    logger.info("="*70)
    
    out_file = out_dir / f"Rome_LST_SST_{year}_30m_15day.nc"
    
    if out_file.exists():
        logger.info(f"File exists, skipping: {out_file.name}")
        continue
    
    # Process Landsat
    landsat_datasets = []
    if year in landsat_by_year:
        logger.info(f"Processing {len(landsat_by_year[year])} Landsat files")
        
        for tif_file, date in landsat_by_year[year]:
            try:
                ds = read_tiff_with_band_names(tif_file, LANDSAT_BANDS, 'LANDSAT')
                ds = ds.expand_dims(time=[pd.Timestamp(date)])
                landsat_datasets.append(ds)
            except Exception as e:
                logger.warning(f"Error loading {tif_file.name}: {e}")
    
    # Process MODIS
    modis_datasets = []
    if year in modis_by_year:
        logger.info(f"Processing {len(modis_by_year[year])} MODIS files")
        
        for tif_file, date in modis_by_year[year]:
            try:
                ds = read_tiff_with_band_names(tif_file, MODIS_BANDS, 'MODIS')
                ds = ds.expand_dims(time=[pd.Timestamp(date)])
                modis_datasets.append(ds)
            except Exception as e:
                logger.warning(f"Error loading {tif_file.name}: {e}")
    
    if not landsat_datasets and not modis_datasets:
        logger.warning(f"No data for {year}")
        continue
    
    # Combine datasets
    combined_vars = {}
    
    if landsat_datasets:
        landsat_datasets = sorted(landsat_datasets, key=lambda d: d.time.values[0])
        landsat_combined = xr.concat(landsat_datasets, dim="time")
        for var in landsat_combined.data_vars:
            combined_vars[var] = landsat_combined[var]
        logger.info(f"Landsat: {len(landsat_combined.time)} timesteps")
        time_coord = landsat_combined.time
        x_coord = landsat_combined.x
        y_coord = landsat_combined.y
    
    if modis_datasets:
        modis_datasets = sorted(modis_datasets, key=lambda d: d.time.values[0])
        modis_combined = xr.concat(modis_datasets, dim="time")
        for var in modis_combined.data_vars:
            combined_vars[var] = modis_combined[var]
        logger.info(f"MODIS: {len(modis_combined.time)} timesteps")
        if not landsat_datasets:
            time_coord = modis_combined.time
            x_coord = modis_combined.x
            y_coord = modis_combined.y
    
    # Create dataset
    ds_year = xr.Dataset(combined_vars, coords={'time': time_coord, 'x': x_coord, 'y': y_coord})
    
    # Add static layers (broadcast and interpolate to data grid)
    if 'DEM' in static_data:
        try:
            logger.info(f"DEM coords: x={static_data['DEM'].x.shape}, y={static_data['DEM'].y.shape}")
            logger.info(f"Data coords: x={x_coord.shape}, y={y_coord.shape}")
            logger.info(f"DEM x range: [{static_data['DEM'].x.min().values:.4f}, {static_data['DEM'].x.max().values:.4f}]")
            logger.info(f"Data x range: [{x_coord.min().values:.4f}, {x_coord.max().values:.4f}]")
            
            dem_interp = static_data['DEM'].interp(x=x_coord, y=y_coord, method='linear', kwargs={'fill_value': 'extrapolate'})
            ds_year['DEM'] = dem_interp.expand_dims(time=time_coord)
            logger.info(f"Added DEM: shape={dem_interp.shape}")
        except Exception as e:
            logger.error(f"Failed to add DEM: {e}")
            import traceback
            logger.error(traceback.format_exc())
    else:
        logger.warning("DEM not in static_data")
    
    if 'WATERMASK' in static_data:
        try:
            water_interp = static_data['WATERMASK'].interp(x=x_coord, y=y_coord, method='nearest', kwargs={'fill_value': 'extrapolate'})
            ds_year['WATERMASK'] = water_interp.expand_dims(time=time_coord)
            logger.info(f"Added WATERMASK: shape={water_interp.shape}")
        except Exception as e:
            logger.error(f"Failed to add WATERMASK: {e}")
            import traceback
            logger.error(traceback.format_exc())
    else:
        logger.warning("WATERMASK not in static_data")
    
    # Format time
    ds_year["time"] = pd.to_datetime(ds_year.time.values)
    ds_year["time"].encoding.update({
        "units": "days since 1970-01-01 00:00:00",
        "calendar": "gregorian"
    })
    
    # Coordinate attributes
    ds_year["x"].attrs.update({
        "long_name": "longitude",
        "standard_name": "longitude",
        "units": "degrees_east"
    })
    
    ds_year["y"].attrs.update({
        "long_name": "latitude",
        "standard_name": "latitude",
        "units": "degrees_north"
    })
    
    # Add comprehensive metadata with config values
    ds_year.attrs.update({
        "title": f"Multi-Sensor LST/SST Dataset - Rome, Italy - {year}",
        "institution": "CMCC Foundation",
        "creator_name": "Alessandro Sebastianelli",
        "year": year,
        "temporal_resolution": f"{config['temporal']['interval_days']} days",
        "spatial_resolution": f"{config['spatial']['target_resolution']} meters",
        "spatial_resolution_original": "30m (Landsat), 1km (MODIS)",
        "resampled_resolution": f"{config['spatial']['target_resolution']}m",
        "time_coverage_start": str(ds_year.time.min().values),
        "time_coverage_end": str(ds_year.time.max().values),
        "area_of_interest": f"Rome: {config['geometry']}",
        "Conventions": "CF-1.8, ACDD-1.3",
        "history": f"Created {pd.Timestamp.now().isoformat()}",
        "source": "Landsat 5/7/8, MODIS Terra/Aqua",
        "processing_level": "Level-3",
        "tile_size_degrees": config['spatial']['tile_size']
    })
    
    # Save with MINIMAL memory footprint
    logger.info(f"Saving: {out_file.name}")
    logger.info(f"Variables: {len(ds_year.data_vars)}")
    logger.info(f"Dimensions: time={len(ds_year.time)}, y={ds_year.sizes['y']}, x={ds_year.sizes['x']}")
    
    encoding = {}
    for var in ds_year.data_vars:
        encoding[var] = {
            "zlib": True,
            "complevel": 4,  # Good compression/speed balance
            "shuffle": True,  # Improves compression
            # Small chunks to minimize memory during write
            "chunksizes": (1,  # One timestep at a time
                          min(256, ds_year.sizes.get('y', 256)),
                          min(256, ds_year.sizes.get('x', 256)))
        }
    
    # Write with MINIMAL memory footprint (single worker, small chunks)
    logger.info("Writing NetCDF with minimal memory usage (may take 5-10 minutes)...")
    with dask.config.set(scheduler='threads', num_workers=1):
        ds_year.to_netcdf(out_file, encoding=encoding, mode='w', compute=True)
    
    file_size_mb = out_file.stat().st_size / (1024**2)
    logger.info(f"✓ Created: {out_file.name} ({file_size_mb:.2f} MB)")


logger.info("="*70)
logger.info("COMPLETE")
logger.info("="*70)