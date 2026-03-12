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

dask.config.set({
    'scheduler': 'threads',
    'num_workers': 2,
    'threads_per_worker': 1,
    'array.chunk-size': '64MB',
})

warnings.filterwarnings('ignore', category=RuntimeWarning)

# ============================================
# ARGUMENTS
# ============================================

parser = argparse.ArgumentParser(description='Create NetCDF for specific year(s)')
parser.add_argument('--year', type=int, default=None)
parser.add_argument('--start-year', type=int, default=None)
parser.add_argument('--end-year', type=int, default=None)
args = parser.parse_args()

# ============================================
# LOGGING
# ============================================

def setup_logging(config, year_suffix=None):
    log_level = getattr(logging, config['logging']['level'])
    log_file = config['logging']['log_file']
    if year_suffix:
        log_file = log_file.replace('.log', f'_{year_suffix}.log')
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
        datefmt='%Y-%m-%d %H:%M:%S'))
    logger.addHandler(fh)

    ch = logging.StreamHandler()
    ch.setLevel(log_level)

    class ColoredFormatter(logging.Formatter):
        COLORS = {'DEBUG': '\033[36m', 'INFO': '\033[32m',
                  'WARNING': '\033[33m', 'ERROR': '\033[31m', 'CRITICAL': '\033[35m'}
        RESET = '\033[0m'
        def format(self, record):
            color = self.COLORS.get(record.levelname, self.RESET)
            return f"{color}{record.levelname[0]}{self.RESET} | {record.getMessage()}"

    ch.setFormatter(ColoredFormatter())
    logger.addHandler(ch)
    return logger

# ============================================
# LOAD CONFIG
# ============================================

with open('config_refactored.yaml', 'r') as f:
    config = yaml.safe_load(f)

logger = setup_logging(config, year_suffix=args.year if args.year else None)

base_dir   = Path(config['output']['base_dir'])
landsat_dir = base_dir / config['output']['subdirs']['landsat']
modis_dir   = base_dir / config['output']['subdirs']['modis']
static_dir  = base_dir / config['output']['subdirs']['static']
out_dir     = base_dir / config['output'].get('netcdf_yearly_subdir', 'netcdf_yearly')
out_dir.mkdir(parents=True, exist_ok=True)

logger.info("="*70)
logger.info("Creating Yearly NetCDF Files")
logger.info("="*70)

# ============================================
# BAND NAME DEFINITIONS
# These must match EXACTLY the band order in download_landsat.py
# and download_modis.py
# ============================================

# Landsat bands in the order they appear in the output TIF
# (matches the addBands() sequence in download_landsat.py)
LANDSAT_BANDS = [
    # Original SR bands (after harmonization to L8 naming)
    'SR_B1',        # Coastal/aerosol (zero for L5/L7)
    'SR_B2',        # Blue
    'SR_B3',        # Green
    'SR_B4',        # Red
    'SR_B5',        # NIR
    'SR_B6',        # SWIR1
    'SR_B7',        # SWIR2
    'ST_B10',       # Thermal brightness temperature (K)
    'QA_PIXEL',     # Quality flags bitmask
    # Spectral indices
    'NDVI',
    'NDWI',
    'NDBI',
    'MNDWI',
    'SAVI',
    'EVI',
    'BSI',
    'UI',
    'Albedo',
    'FV',           # Fractional vegetation cover
    # Emissivity
    'emissivity',   # FVC-based emissivity (used for LST)
    # Cloud
    'cloud_mask',
    # Temperature products
    'LST_K',        # Land Surface Temperature (Kelvin)
    'LST_K_masked', # LST with cloud mask applied
    'SST_K',        # Sea Surface Temperature (= LST; water masking in netcdf step)
    'LSSTSST',      # Merged land-sea temperature field
    # Coordinates and auxiliary
    'lon',
    'lat',
    'time_days',    # Days since 1970-01-01
    'land_cover',
    'emissivity_dynamic',  # Emissivity from land cover lookup
    'AOD_047',
    'AOD_055',
]

# MODIS bands in output order (matches download_modis.py addBands sequence)
MODIS_BANDS = [
    'MODIS_LST_Day_K',
    'MODIS_LST_Night_K',
    'LST_MODIS',           # Day+Night mean
    'MODIS_SST_Day_K',
    'MODIS_SST_Night_K',
    'SST_MODIS',           # Day+Night mean (water pixels in netcdf step)
    'MODIS_QC_Day',
    'MODIS_QC_Night',
    'MODIS_Emis_31',
    'MODIS_Emis_32',
    'MODIS_View_angle_Day',
    'MODIS_View_angle_Night',
    'MODIS_View_time_Day',
    'MODIS_View_time_Night',
    # Auxiliary
    'lon',
    'lat',
    'time_days',
    'land_cover',
    'emissivity_dynamic',
]

# ============================================
# HELPER FUNCTIONS
# ============================================

date_re = re.compile(r"(\d{4})-(\d{2})-(\d{2})")

def extract_date_and_year(stem):
    m = date_re.search(stem)
    if m:
        year = int(m.group(1))
        date = pd.to_datetime(f"{m.group(1)}-{m.group(2)}-{m.group(3)}", errors='coerce')
        return date, year
    return None, None

def read_tiff_with_band_names(tif_file, expected_bands, prefix):
    """
    Read a GeoTIFF, assigning band names from expected_bands.
    Bands whose name already starts with the prefix (e.g. MODIS_LST_Day_K)
    are used as-is; others are prefixed with 'LANDSAT_' or 'MODIS_'.
    Auxiliary bands (lon, lat, time_days, land_cover, emissivity_dynamic)
    are NOT prefixed.
    """
    AUX_BANDS = {'lon', 'lat', 'time_days', 'land_cover', 'emissivity_dynamic'}

    with rasterio.open(tif_file) as src:
        transform = src.transform
        n_bands   = src.count
        shape     = (src.height, src.width)

    y_coords = np.arange(shape[0]) * transform[4] + transform[5]
    x_coords = np.arange(shape[1]) * transform[0] + transform[2]

    import dask.array as da
    dataset_dict = {}
    chunk_size = 256

    for i in range(min(n_bands, len(expected_bands))):
        raw_name = expected_bands[i]

        if raw_name in AUX_BANDS:
            band_name = raw_name
        elif raw_name.startswith('MODIS_') or raw_name.startswith('LST_MODIS') or raw_name.startswith('SST_MODIS'):
            band_name = raw_name
        else:
            band_name = f"{prefix}_{raw_name}"

        def read_band(band_idx, filename=tif_file):
            with rasterio.open(filename) as src:
                return src.read(band_idx + 1).astype(np.float32)

        lazy_array = da.from_delayed(
            dask.delayed(read_band)(i),
            shape=shape,
            dtype=np.float32
        ).rechunk((chunk_size, chunk_size))

        dataset_dict[band_name] = (['y', 'x'], lazy_array)

    return xr.Dataset(dataset_dict, coords={'y': y_coords, 'x': x_coords})

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

if args.year:
    all_years = [args.year] if args.year in all_years else []
    if not all_years:
        logger.error(f"No data for year {args.year}")
        sys.exit(1)
elif args.start_year or args.end_year:
    sy = args.start_year or min(all_years)
    ey = args.end_year or max(all_years)
    all_years = [y for y in all_years if sy <= y <= ey]

logger.info(f"Years to process: {all_years}")

# ============================================
# LOAD STATIC LAYERS
# ============================================

logger.info("Loading static layers...")
dem_file   = static_dir / 'DEM.tif'
water_file = static_dir / 'WaterMask.tif'

static_data = {}
for name, path in [('DEM', dem_file), ('WATERMASK', water_file)]:
    if path.exists():
        try:
            arr = rioxarray.open_rasterio(path)
            if 'band' in arr.dims:
                arr = arr.isel(band=0)
            arr = arr.squeeze()
            static_data[name] = arr
            logger.info(f"  {name}: shape={arr.shape}, "
                        f"x=[{float(arr.x.min()):.4f}, {float(arr.x.max()):.4f}]")
        except Exception as e:
            logger.error(f"  Failed to load {name}: {e}")
    else:
        logger.warning(f"  {name} not found: {path}")

# ============================================
# PROCESS EACH YEAR
# ============================================

for year in all_years:
    logger.info("="*70)
    logger.info(f"Processing year: {year}")
    logger.info("="*70)

    out_file = out_dir / f"Rome_LST_SST_{year}_30m_15day.nc"
    if out_file.exists():
        logger.info(f"Exists, skipping: {out_file.name}")
        continue

    # --- Landsat ---
    landsat_datasets = []
    if year in landsat_by_year:
        logger.info(f"Loading {len(landsat_by_year[year])} Landsat files...")
        for tif_file, date in landsat_by_year[year]:
            try:
                ds = read_tiff_with_band_names(tif_file, LANDSAT_BANDS, 'LANDSAT')
                ds = ds.expand_dims(time=[pd.Timestamp(date)])
                landsat_datasets.append(ds)
            except Exception as e:
                logger.warning(f"  Skip {tif_file.name}: {e}")

    # --- MODIS ---
    modis_datasets = []
    if year in modis_by_year:
        logger.info(f"Loading {len(modis_by_year[year])} MODIS files...")
        for tif_file, date in modis_by_year[year]:
            try:
                ds = read_tiff_with_band_names(tif_file, MODIS_BANDS, 'MODIS')
                ds = ds.expand_dims(time=[pd.Timestamp(date)])
                modis_datasets.append(ds)
            except Exception as e:
                logger.warning(f"  Skip {tif_file.name}: {e}")

    if not landsat_datasets and not modis_datasets:
        logger.warning(f"No data for {year}, skipping")
        continue

    # --- Combine ---
    combined_vars = {}
    time_coord = x_coord = y_coord = None

    if landsat_datasets:
        landsat_datasets.sort(key=lambda d: d.time.values[0])
        lc = xr.concat(landsat_datasets, dim="time")
        for var in lc.data_vars:
            combined_vars[var] = lc[var]
        time_coord, x_coord, y_coord = lc.time, lc.x, lc.y
        logger.info(f"Landsat: {len(lc.time)} timesteps, {len(lc.data_vars)} vars")

    if modis_datasets:
        modis_datasets.sort(key=lambda d: d.time.values[0])
        mc = xr.concat(modis_datasets, dim="time")
        for var in mc.data_vars:
            combined_vars[var] = mc[var]
        if time_coord is None:
            time_coord, x_coord, y_coord = mc.time, mc.x, mc.y
        logger.info(f"MODIS: {len(mc.time)} timesteps, {len(mc.data_vars)} vars")

    ds_year = xr.Dataset(combined_vars,
                         coords={'time': time_coord, 'x': x_coord, 'y': y_coord})

    # --- Static layers ---
    for layer_name, interp_method in [('DEM', 'linear'), ('WATERMASK', 'nearest')]:
        if layer_name in static_data:
            try:
                interpolated = static_data[layer_name].interp(
                    x=x_coord, y=y_coord,
                    method=interp_method,
                    kwargs={'fill_value': 'extrapolate'}
                )
                ds_year[layer_name] = interpolated.expand_dims(time=time_coord)
                logger.info(f"Added {layer_name}: shape={interpolated.shape}")
            except Exception as e:
                logger.error(f"Failed to add {layer_name}: {e}")

    # --- Apply water mask to SST variables ---
    # SST_K, MODIS_SST_Day_K, MODIS_SST_Night_K, SST_MODIS are physically only
    # meaningful over water pixels. Apply the water mask here.
    if 'WATERMASK' in ds_year:
        water_bool = ds_year['WATERMASK'] > 0  # 1=water, 0=land
        for sst_var in ['LANDSAT_SST_K', 'LANDSAT_LSSTSST',
                        'MODIS_SST_Day_K', 'MODIS_SST_Night_K', 'SST_MODIS']:
            if sst_var in ds_year:
                ds_year[sst_var] = ds_year[sst_var].where(water_bool)
                logger.info(f"Applied water mask to {sst_var}")
    else:
        logger.warning("WATERMASK not available — SST variables not masked to water pixels")

    # --- Time encoding ---
    ds_year["time"] = pd.to_datetime(ds_year.time.values)
    ds_year["time"].encoding.update({
        "units": "days since 1970-01-01 00:00:00",
        "calendar": "gregorian"
    })

    ds_year["x"].attrs.update({"long_name": "longitude",
                                "standard_name": "longitude",
                                "units": "degrees_east"})
    ds_year["y"].attrs.update({"long_name": "latitude",
                                "standard_name": "latitude",
                                "units": "degrees_north"})

    # --- Variable attributes ---
    var_attrs = {
        # Landsat SR
        'LANDSAT_SR_B1':  {'long_name': 'Landsat SR Band 1 (Coastal/Aerosol; zero for L5/L7)', 'units': 'reflectance'},
        'LANDSAT_SR_B2':  {'long_name': 'Landsat SR Band 2 (Blue)', 'units': 'reflectance'},
        'LANDSAT_SR_B3':  {'long_name': 'Landsat SR Band 3 (Green)', 'units': 'reflectance'},
        'LANDSAT_SR_B4':  {'long_name': 'Landsat SR Band 4 (Red)', 'units': 'reflectance'},
        'LANDSAT_SR_B5':  {'long_name': 'Landsat SR Band 5 (NIR)', 'units': 'reflectance'},
        'LANDSAT_SR_B6':  {'long_name': 'Landsat SR Band 6 (SWIR1)', 'units': 'reflectance'},
        'LANDSAT_SR_B7':  {'long_name': 'Landsat SR Band 7 (SWIR2)', 'units': 'reflectance'},
        'LANDSAT_ST_B10': {'long_name': 'Landsat Brightness Temperature (thermal)', 'units': 'K'},
        'LANDSAT_QA_PIXEL': {'long_name': 'Landsat QA_PIXEL bitmask'},
        # Indices
        'LANDSAT_NDVI':  {'long_name': 'Normalized Difference Vegetation Index', 'valid_range': [-1, 1]},
        'LANDSAT_NDWI':  {'long_name': 'Normalized Difference Water Index (McFeeters 1996)', 'valid_range': [-1, 1]},
        'LANDSAT_NDBI':  {'long_name': 'Normalized Difference Built-up Index (Zha 2003)', 'valid_range': [-1, 1]},
        'LANDSAT_MNDWI': {'long_name': 'Modified NDWI (Xu 2006)', 'valid_range': [-1, 1]},
        'LANDSAT_SAVI':  {'long_name': 'Soil Adjusted Vegetation Index (L=0.5)', 'valid_range': [-1.5, 1.5]},
        'LANDSAT_EVI':   {'long_name': 'Enhanced Vegetation Index'},
        'LANDSAT_BSI':   {'long_name': 'Bare Soil Index'},
        'LANDSAT_UI':    {'long_name': 'Urban Index'},
        'LANDSAT_Albedo': {'long_name': 'Surface Albedo (Liang 2001)', 'valid_range': [0, 1]},
        'LANDSAT_FV':    {'long_name': 'Fractional Vegetation Cover', 'valid_range': [0, 1]},
        'LANDSAT_emissivity': {'long_name': 'Surface emissivity (FVC method, Sobrino 2004)', 'valid_range': [0.9, 1.0]},
        'LANDSAT_cloud_mask': {'long_name': 'Cloud/shadow/snow mask (1=clear)', 'flag_values': [0, 1]},
        # Temperature
        'LANDSAT_LST_K':        {'long_name': 'Land Surface Temperature (Planck inversion)', 'units': 'K',
                                   'standard_name': 'surface_temperature',
                                   'comment': 'Formula: BT / (1 + (lambda*BT/rho)*ln(emissivity)); lambda=10.895µm(L8) or 11.45µm(L5/L7)'},
        'LANDSAT_LST_K_masked': {'long_name': 'LST with cloud mask applied', 'units': 'K'},
        'LANDSAT_SST_K':        {'long_name': 'Sea Surface Temperature (LST over water pixels)', 'units': 'K',
                                   'comment': 'Masked to water pixels using WATERMASK'},
        'LANDSAT_LSSTSST':      {'long_name': 'Merged land+sea surface temperature', 'units': 'K',
                                   'comment': 'LST for land, SST for water pixels'},
        # Auxiliary
        'lon':  {'long_name': 'pixel longitude', 'units': 'degrees_east'},
        'lat':  {'long_name': 'pixel latitude', 'units': 'degrees_north'},
        'time_days': {'long_name': 'days since 1970-01-01', 'units': 'days'},
        'land_cover': {'long_name': 'Land cover class (ESA WorldCover or MODIS IGBP)'},
        'emissivity_dynamic': {'long_name': 'Emissivity from land cover lookup table'},
        'AOD_047': {'long_name': 'Aerosol Optical Depth at 470nm (MODIS MCD19A2)', 'units': '1'},
        'AOD_055': {'long_name': 'Aerosol Optical Depth at 550nm (MODIS MCD19A2)', 'units': '1'},
        # MODIS
        'MODIS_LST_Day_K':   {'long_name': 'MODIS Daytime LST (MOD/MYD11A1 scaled)', 'units': 'K'},
        'MODIS_LST_Night_K': {'long_name': 'MODIS Nighttime LST', 'units': 'K'},
        'LST_MODIS':         {'long_name': 'MODIS LST day+night mean', 'units': 'K'},
        'MODIS_SST_Day_K':   {'long_name': 'MODIS Daytime SST (LST masked to water)', 'units': 'K'},
        'MODIS_SST_Night_K': {'long_name': 'MODIS Nighttime SST (LST masked to water)', 'units': 'K'},
        'SST_MODIS':         {'long_name': 'MODIS SST day+night mean (water pixels)', 'units': 'K'},
        'MODIS_QC_Day':   {'long_name': 'MODIS Daytime QC flags'},
        'MODIS_QC_Night': {'long_name': 'MODIS Nighttime QC flags'},
        'MODIS_Emis_31':  {'long_name': 'MODIS Band 31 emissivity (10.78-11.28µm)', 'units': '1'},
        'MODIS_Emis_32':  {'long_name': 'MODIS Band 32 emissivity (11.77-12.27µm)', 'units': '1'},
        'MODIS_View_angle_Day':   {'long_name': 'MODIS Daytime view zenith angle', 'units': 'degrees'},
        'MODIS_View_angle_Night': {'long_name': 'MODIS Nighttime view zenith angle', 'units': 'degrees'},
        'MODIS_View_time_Day':    {'long_name': 'MODIS Daytime observation time (local solar)', 'units': 'hours'},
        'MODIS_View_time_Night':  {'long_name': 'MODIS Nighttime observation time (local solar)', 'units': 'hours'},
        # Static
        'DEM':       {'long_name': 'SRTM Digital Elevation Model', 'units': 'm', 'source': 'USGS/SRTMGL1_003'},
        'WATERMASK': {'long_name': 'Water mask (1=water, 0=land)',
                      'source': 'JRC GSW1_4 (occurrence>50%) + ESA WorldCover classes 80/90/95',
                      'flag_values': [0, 1]},
    }

    for var, attrs in var_attrs.items():
        if var in ds_year:
            ds_year[var].attrs.update(attrs)

    # --- Global metadata ---
    ds_year.attrs.update({
        "title": f"Multi-Sensor LST/SST Dataset - Rome Metropolitan Area - {year}",
        "institution": "CMCC Foundation",
        "creator_name": "Alessandro Sebastianelli",
        "creator_email": "asebastianelli@ieee.org",
        "year": year,
        "temporal_resolution": f"{config['temporal']['interval_days']} days",
        "spatial_resolution": f"{config['spatial']['target_resolution']} meters",
        "time_coverage_start": str(ds_year.time.min().values),
        "time_coverage_end": str(ds_year.time.max().values),
        "geospatial_lon_min": config['geometry']['min_lon'],
        "geospatial_lon_max": config['geometry']['max_lon'],
        "geospatial_lat_min": config['geometry']['min_lat'],
        "geospatial_lat_max": config['geometry']['max_lat'],
        "area_of_interest": "Rome metropolitan area and Tyrrhenian coast, Italy",
        "Conventions": "CF-1.8, ACDD-1.3",
        "history": f"Created {pd.Timestamp.now().isoformat()}",
        "source_landsat": "Landsat 5/7/8 Collection 2 Level-2 (USGS)",
        "source_modis": f"{config['modis']['collections']['terra']} + {config['modis']['collections']['aqua']}",
        "source_dem": "USGS SRTMGL1_003",
        "source_watermask": "JRC/GSW1_4/GlobalSurfaceWater + ESA/WorldCover/v200",
        "source_aerosol": "MODIS/061/MCD19A2_GRANULES",
        "processing_platform": "Google Earth Engine",
        "tile_size_degrees": config['spatial']['tile_size'],
        "compositing_method": config['compositing']['method'],
        "lst_formula": "LST = BT / (1 + (lambda*BT/rho) * ln(emissivity)); rho=14388 µm·K",
        "lst_wavelength_l5l7": "11.45 µm (Band 6)",
        "lst_wavelength_l8": "10.895 µm (Band 10)",
        "emissivity_method": "FVC method (Sobrino et al. 2004): eps = FV*0.004 + 0.986",
        "sr_scaling": "DN * 2.75e-5 + (-0.2) = surface reflectance",
        "bt_scaling": "DN * 0.00341802 + 149.0 = brightness temperature [K]",
        "sst_note": "SST variables are LST values masked to water pixels (WATERMASK=1)",
        "landsat_harmonization": "L5/L7 bands renamed to L8-style; SR_B1 zero-filled for L5/L7",
        "variables_count": len(ds_year.data_vars),
    })

    # --- Save ---
    logger.info(f"Saving: {out_file.name}")
    logger.info(f"  Variables: {len(ds_year.data_vars)}")
    logger.info(f"  Dimensions: time={len(ds_year.time)}, y={ds_year.sizes['y']}, x={ds_year.sizes['x']}")

    encoding = {}
    for var in ds_year.data_vars:
        encoding[var] = {
            "zlib": True,
            "complevel": 4,
            "shuffle": True,
            "chunksizes": (
                1,
                min(256, ds_year.sizes.get('y', 256)),
                min(256, ds_year.sizes.get('x', 256))
            )
        }

    logger.info("Writing NetCDF (single-threaded, minimal memory)...")
    with dask.config.set(scheduler='threads', num_workers=1):
        ds_year.to_netcdf(out_file, encoding=encoding, mode='w', compute=True)

    size_mb = out_file.stat().st_size / (1024**2)
    logger.info(f"✓ Created: {out_file.name} ({size_mb:.2f} MB)")

logger.info("="*70)
logger.info("COMPLETE")
logger.info("="*70)