# ============================================
# Create Combined NetCDF with Specific Structure
# All fields broadcast to same temporal dimension
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
    """Setup clean, consistent logging"""
    log_level = getattr(logging, config['logging']['level'])
    log_file = config['logging']['log_file']
    
    logger = logging.getLogger('NetCDF')
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
logger.info("Creating Combined NetCDF - Standardized Structure")
logger.info("="*70)

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

logger.info("\nLoading static layers...")

dem_file = static_dir / 'DEM.tif'
water_mask_file = static_dir / 'WaterMask.tif'

dem_data = None
water_mask_data = None

if dem_file.exists():
    logger.info(f"  ✓ Loading {dem_file.name}")
    dem_data = rioxarray.open_rasterio(dem_file).squeeze()

if water_mask_file.exists():
    logger.info(f"  ✓ Loading {water_mask_file.name}")
    water_mask_data = rioxarray.open_rasterio(water_mask_file).squeeze()

# ============================================
# LOAD LANDSAT
# ============================================

logger.info("\nLoading Landsat data...")

landsat_files = sorted(landsat_dir.glob("*.tif"))
logger.info(f"Found {len(landsat_files)} Landsat files")

landsat_datasets = []

for tif in landsat_files:
    date = extract_date(tif.stem)
    if date is None or pd.isna(date):
        logger.warning(f"Cannot parse date from {tif.name}")
        continue
    
    ds = rioxarray.open_rasterio(tif)
    
    # Clean dimensions
    for dim in list(ds.dims):
        if dim not in ("band", "y", "x"):
            ds = ds.squeeze(dim, drop=True)
    
    ds = ds.squeeze(drop=True).expand_dims(time=[pd.Timestamp(date)])
    
    # Keep essential coords
    for c in list(ds.coords):
        if c not in ["time", "x", "y", "band"]:
            ds = ds.drop_vars(c)
    
    if "band" in ds.dims:
        actual_bands = ds.sizes.get("band", 0)
        
        # Expected ~20 bands
        if actual_bands >= 18:
            band_names = [
                'SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7',
                'ST_B10', 'QA_PIXEL',
                'LST_K', 'LST_K_masked', 'SST_K', 'cloud_mask', 'NDVI', 'emissivity',
                'lon', 'lat', 'time',
                'land_cover', 'emissivity_dynamic'
            ][:actual_bands]
            
            ds = ds.assign_coords(band=band_names)
            ds = ds.to_dataset(dim="band")
            landsat_datasets.append(ds)

logger.info(f"Loaded {len(landsat_datasets)} valid Landsat datasets")

# ============================================
# LOAD MODIS
# ============================================

logger.info("\nLoading MODIS data...")

modis_files = sorted(modis_dir.glob("*.tif"))
logger.info(f"Found {len(modis_files)} MODIS files")

modis_datasets = []

for tif in modis_files:
    date = extract_date(tif.stem)
    if date is None or pd.isna(date):
        logger.warning(f"Cannot parse date from {tif.name}")
        continue
    
    ds = rioxarray.open_rasterio(tif)
    
    for dim in list(ds.dims):
        if dim not in ("band", "y", "x"):
            ds = ds.squeeze(dim, drop=True)
    
    ds = ds.squeeze(drop=True).expand_dims(time=[pd.Timestamp(date)])
    
    for c in list(ds.coords):
        if c not in ["time", "x", "y", "band"]:
            ds = ds.drop_vars(c)
    
    if "band" in ds.dims:
        actual_bands = ds.sizes.get("band", 0)
        
        # Expected ~17 bands
        if actual_bands >= 15:
            band_names = [
                'LST_Day_K', 'LST_Night_K',
                'SST_Day_K', 'SST_Night_K',
                'QC_Day', 'QC_Night',
                'Emis_31', 'Emis_32',
                'View_angle_Day', 'View_angle_Night',
                'View_time_Day', 'View_time_Night',
                'lon', 'lat', 'time',
                'land_cover', 'emissivity_dynamic'
            ][:actual_bands]
            
            ds = ds.assign_coords(band=band_names)
            ds = ds.to_dataset(dim="band")
            modis_datasets.append(ds)

logger.info(f"Loaded {len(modis_datasets)} valid MODIS datasets")

# ============================================
# CONCATENATE TIME SERIES
# ============================================

logger.info("\nConcatenating time series...")

if len(landsat_datasets) > 0:
    landsat_datasets = sorted(landsat_datasets, key=lambda d: pd.Timestamp(d.time.values[0]))
    landsat_combined = xr.concat(landsat_datasets, dim="time", data_vars="all",
                                  coords="minimal", compat='override', join='override')
    logger.info(f"Landsat: {len(landsat_combined.time)} timesteps")
else:
    landsat_combined = None
    logger.warning("No Landsat data")

if len(modis_datasets) > 0:
    modis_datasets = sorted(modis_datasets, key=lambda d: pd.Timestamp(d.time.values[0]))
    modis_combined = xr.concat(modis_datasets, dim="time", data_vars="all",
                                coords="minimal", compat='override', join='override')
    logger.info(f"MODIS: {len(modis_combined.time)} timesteps")
else:
    modis_combined = None
    logger.warning("No MODIS data")

# ============================================
# CREATE FINAL DATASET WITH REQUIRED STRUCTURE
# ============================================

logger.info("\nCreating final dataset with required structure...")

# Determine common time coordinate (union of all times)
all_times = []
if landsat_combined is not None:
    all_times.extend(landsat_combined.time.values)
if modis_combined is not None:
    all_times.extend(modis_combined.time.values)

all_times = sorted(list(set(all_times)))
time_coord = pd.to_datetime(all_times)

logger.info(f"Total unique timesteps: {len(time_coord)}")

# Create empty dataset
final = xr.Dataset(coords={'time': time_coord})

# Get spatial coordinates from first available source
if landsat_combined is not None:
    final.coords['y'] = landsat_combined.y
    final.coords['x'] = landsat_combined.x
elif modis_combined is not None:
    final.coords['y'] = modis_combined.y
    final.coords['x'] = modis_combined.x

# ============================================
# 1. LANDSAT - ALL BANDS
# ============================================

logger.info("\n1. Adding LANDSAT bands...")

if landsat_combined is not None:
    for var in ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7',
                'ST_B10', 'QA_PIXEL']:
        if var in landsat_combined:
            final[f'LANDSAT_{var}'] = landsat_combined[var].reindex(time=time_coord)
            logger.info(f"  ✓ LANDSAT_{var}")

# ============================================
# 2. MODIS - ALL BANDS
# ============================================

logger.info("\n2. Adding MODIS bands...")

if modis_combined is not None:
    for var in ['QC_Day', 'QC_Night', 'Emis_31', 'Emis_32',
                'View_angle_Day', 'View_angle_Night',
                'View_time_Day', 'View_time_Night']:
        if var in modis_combined:
            final[f'MODIS_{var}'] = modis_combined[var].reindex(time=time_coord)
            logger.info(f"  ✓ MODIS_{var}")

# ============================================
# 3. LST LANDSAT
# ============================================

logger.info("\n3. Adding LST_LANDSAT...")

if landsat_combined is not None and 'LST_K' in landsat_combined:
    final['LST_LANDSAT'] = landsat_combined['LST_K'].reindex(time=time_coord)
    logger.info("  ✓ LST_LANDSAT")

# ============================================
# 4. LST MODIS
# ============================================

logger.info("\n4. Adding LST_MODIS...")

if modis_combined is not None:
    # Average of day and night
    if 'LST_Day_K' in modis_combined and 'LST_Night_K' in modis_combined:
        lst_modis = (modis_combined['LST_Day_K'] + modis_combined['LST_Night_K']) / 2.0
        final['LST_MODIS'] = lst_modis.reindex(time=time_coord)
        logger.info("  ✓ LST_MODIS (average of day/night)")

# ============================================
# 5. SST LANDSAT
# ============================================

logger.info("\n5. Adding SST_LANDSAT...")

if landsat_combined is not None and 'SST_K' in landsat_combined:
    final['SST_LANDSAT'] = landsat_combined['SST_K'].reindex(time=time_coord)
    logger.info("  ✓ SST_LANDSAT")

# ============================================
# 6. SST MODIS
# ============================================

logger.info("\n6. Adding SST_MODIS...")

if modis_combined is not None:
    # Average of day and night
    if 'SST_Day_K' in modis_combined and 'SST_Night_K' in modis_combined:
        sst_modis = (modis_combined['SST_Day_K'] + modis_combined['SST_Night_K']) / 2.0
        final['SST_MODIS'] = sst_modis.reindex(time=time_coord)
        logger.info("  ✓ SST_MODIS (average of day/night)")

# ============================================
# 7. WATERMASK - FOR ALL DATES
# ============================================

logger.info("\n7. Adding WATERMASK...")

if water_mask_data is not None:
    # Broadcast to all timesteps
    water_broadcast = water_mask_data.expand_dims(time=time_coord)
    final['WATERMASK'] = water_broadcast
    logger.info("  ✓ WATERMASK (broadcast to all dates)")

# ============================================
# 8. LSSTSST - USING WATERMASK
# ============================================

logger.info("\n8. Creating LSSTSST (hybrid LST/SST using watermask)...")

if 'WATERMASK' in final and 'LST_LANDSAT' in final and 'SST_LANDSAT' in final:
    # Where watermask==1 (water): use SST, else: use LST
    lsstsst = xr.where(
        final['WATERMASK'] == 1,
        final['SST_LANDSAT'],  # Water: SST
        final['LST_LANDSAT']   # Land: LST
    )
    final['LSSTSST'] = lsstsst
    logger.info("  ✓ LSSTSST (LST on land, SST on water)")

# ============================================
# 9. DEM - FOR ALL DATES
# ============================================

logger.info("\n9. Adding DEM...")

if dem_data is not None:
    # Broadcast to all timesteps
    dem_broadcast = dem_data.expand_dims(time=time_coord)
    final['DEM'] = dem_broadcast
    logger.info("  ✓ DEM (broadcast to all dates)")

# ============================================
# 10. EMISSIVITY - FOR ALL DATES
# ============================================

logger.info("\n10. Adding EMISSIVITY...")

if landsat_combined is not None and 'emissivity_dynamic' in landsat_combined:
    final['EMISSIVITY'] = landsat_combined['emissivity_dynamic'].reindex(time=time_coord)
    logger.info("  ✓ EMISSIVITY (from Landsat, dynamic)")
elif modis_combined is not None and 'emissivity_dynamic' in modis_combined:
    final['EMISSIVITY'] = modis_combined['emissivity_dynamic'].reindex(time=time_coord)
    logger.info("  ✓ EMISSIVITY (from MODIS, dynamic)")

# ============================================
# 11. LANDCOVER - FOR ALL DATES
# ============================================

logger.info("\n11. Adding LANDCOVER...")

if landsat_combined is not None and 'land_cover' in landsat_combined:
    final['LANDCOVER'] = landsat_combined['land_cover'].reindex(time=time_coord)
    logger.info("  ✓ LANDCOVER (from Landsat, dynamic)")
elif modis_combined is not None and 'land_cover' in modis_combined:
    final['LANDCOVER'] = modis_combined['land_cover'].reindex(time=time_coord)
    logger.info("  ✓ LANDCOVER (from MODIS, dynamic)")

# ============================================
# TIME FORMATTING
# ============================================

final["time"].encoding.update({
    "units": "days since 1970-01-01 00:00:00",
    "calendar": "gregorian"
})

# ============================================
# METADATA
# ============================================

logger.info("\nAdding comprehensive metadata...")

final.attrs.update({
    # Dataset identification
    "title": "Multi-Sensor Land and Sea Surface Temperature Dataset - Rome Metropolitan Area",
    "summary": "High-resolution (30m) time series of land and sea surface temperatures combining Landsat 8/9 and MODIS Terra/Aqua observations with ancillary geospatial data",
    "keywords": "Land Surface Temperature, Sea Surface Temperature, LST, SST, Landsat, MODIS, Remote Sensing, Earth Observation, Rome, Italy, Mediterranean, Thermal Remote Sensing",
    "id": f"LST_SST_Rome_{pd.Timestamp.now().strftime('%Y%m%d')}",
    
    # Temporal coverage
    "time_coverage_start": str(final.time.min().values),
    "time_coverage_end": str(final.time.max().values),
    "time_coverage_resolution": f"P{config['temporal']['interval_days']}D",
    "time_coverage_duration": f"P{(pd.Timestamp(final.time.max().values) - pd.Timestamp(final.time.min().values)).days}D",
    "temporal_resolution": "15 days",
    "temporal_description": "15-day median composites to reduce cloud contamination and provide consistent temporal sampling",
    
    # Spatial coverage
    "geospatial_lat_min": float(config['geometry']['min_lat']),
    "geospatial_lat_max": float(config['geometry']['max_lat']),
    "geospatial_lon_min": float(config['geometry']['min_lon']),
    "geospatial_lon_max": float(config['geometry']['max_lon']),
    "geospatial_lat_units": "degrees_north",
    "geospatial_lon_units": "degrees_east",
    "geospatial_vertical_min": 0.0,
    "geospatial_vertical_max": 0.0,
    "geospatial_vertical_units": "meters",
    "geospatial_vertical_positive": "up",
    "spatial_resolution": "30 meters",
    "spatial_description": "Rome metropolitan area including urban center, suburban areas, agricultural lands, and coastal Mediterranean waters (Tyrrhenian Sea)",
    
    # Data sources
    "source": "Satellite remote sensing: Landsat 8/9 Collection 2 Level-2, MODIS Terra/Aqua MOD11A1/MYD11A1 Collection 6.1, SRTM DEM v3, MODIS MCD12Q1 Land Cover, JRC Global Surface Water",
    "source_landsat": "LANDSAT/LC08/C02/T1_L2 (2020-2021), LANDSAT/LC09/C02/T1_L2 (2021-2025)",
    "source_modis": "MODIS/006/MOD11A1 (Terra), MODIS/006/MYD11A1 (Aqua)",
    "source_dem": "USGS/SRTMGL1_003 (Shuttle Radar Topography Mission, 30m resolution)",
    "source_land_cover": "MODIS/061/MCD12Q1 (MODIS Land Cover Type Yearly, IGBP classification)",
    "source_water_mask": "JRC/GSW1_4/GlobalSurfaceWater (JRC Global Surface Water) + ESA/WorldCover/v200",
    "platform": "Landsat-8, Landsat-9, Terra, Aqua",
    "sensor": "OLI/TIRS (Landsat-8), OLI-2/TIRS-2 (Landsat-9), MODIS (Terra/Aqua)",
    
    # Processing information
    "processing_level": "Level-3 derived product",
    "processing_description": """Multi-step processing pipeline:
1. Data acquisition from Google Earth Engine
2. Cloud masking using QA_PIXEL band (Landsat) and QC flags (MODIS)
3. Temporal compositing: 15-day median to reduce cloud contamination
4. LST calculation from thermal bands using emissivity correction
5. SST derived from thermal bands for water bodies
6. MODIS data resampled from 1km to 30m using bilinear interpolation
7. Integration of ancillary datasets (DEM, land cover, water mask)
8. Creation of hybrid LSSTSST product using water mask classification""",
    
    "processing_software": "Python 3.x with Google Earth Engine Python API, xarray, rasterio, rioxarray, dask",
    "processing_institution": "CMCC Foundation (Centro Euro-Mediterraneo sui Cambiamenti Climatici)",
    "processing_date": pd.Timestamp.now().isoformat(),
    
    # LST/SST methodology
    "lst_algorithm": "Split-window algorithm with NDVI-based emissivity correction",
    "lst_description": "Land Surface Temperature derived from Landsat thermal band (10.6-11.19 µm) using brightness temperature conversion with emissivity correction based on NDVI and fractional vegetation cover",
    "sst_description": "Sea Surface Temperature derived using same thermal algorithm as LST, applied to water bodies identified by water mask",
    "emissivity_method": "Calculated from NDVI-derived fractional vegetation cover for Landsat; from MODIS land cover classification lookup table for MODIS data",
    "cloud_masking_method": "Landsat: QA_PIXEL band filtering (cloud, cloud shadow, snow); MODIS: QC flags filtering",
    "compositing_method": "15-day median composite to minimize cloud contamination and atmospheric effects",
    
    # Quality information
    "quality_control": "Cloud masking, data validation, corrupted file removal, spatial consistency checks",
    "data_quality_description": "Products include both cloud-masked and non-masked versions for quality assessment",
    "uncertainties": "LST/SST accuracy affected by atmospheric conditions, viewing geometry, surface emissivity variations, and cloud contamination in composites",
    
    # Variables description
    "variables_landsat": "LANDSAT_SR_B1-B7 (surface reflectance), LANDSAT_ST_B10 (thermal), LANDSAT_QA_PIXEL (quality)",
    "variables_modis": "MODIS_QC_Day/Night (quality), MODIS_Emis_31/32 (emissivity), MODIS_View_angle_Day/Night, MODIS_View_time_Day/Night",
    "variables_temperature": "LST_LANDSAT, LST_MODIS, SST_LANDSAT, SST_MODIS, LSSTSST (hybrid)",
    "variables_ancillary": "DEM (elevation), WATERMASK (land/water classification), LANDCOVER (MODIS MCD12Q1), EMISSIVITY (surface emissivity)",
    
    # Hybrid product description
    "lsstsst_description": "Hybrid temperature product combining LST over land surfaces and SST over water bodies, seamlessly integrated using water mask classification",
    "lsstsst_methodology": "Conditional selection: LSSTSST = SST_LANDSAT where WATERMASK=1 (water), LST_LANDSAT where WATERMASK=0 (land)",
    
    # Conventions and standards
    "Conventions": "CF-1.8, ACDD-1.3",
    "standard_name_vocabulary": "CF Standard Name Table v79",
    "metadata_conventions": "Unidata Dataset Discovery v1.0, CF-1.8, ACDD-1.3",
    
    # Creator information
    "creator_name": "Alessandro Sebastianelli",
    "creator_email": "alessandro.sebastianelli@cmcc.it",
    "creator_url": "https://alessandrosebastianelli.github.io/",
    "creator_type": "person",
    "creator_institution": "CMCC Foundation - Centro Euro-Mediterraneo sui Cambiamenti Climatici",
    "institution": "CMCC Foundation",
    "institution_url": "https://www.cmcc.it/",
    
    # Project information
    "project": "Multi-Sensor LST/SST Analysis",
    "program": "Earth Observation and Remote Sensing",
    
    # Publication information
    "date_created": pd.Timestamp.now().isoformat(),
    "date_modified": pd.Timestamp.now().isoformat(),
    "date_issued": pd.Timestamp.now().isoformat(),
    "history": f"""Created {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}
Data downloaded from Google Earth Engine
Processing: Cloud masking, temporal compositing, LST/SST calculation, resampling, integration
Final product: NetCDF-4 with CF-1.8 conventions""",
    
    # License and usage
    "license": "Data are provided for scientific research purposes. Users should acknowledge the data sources (Landsat, MODIS, SRTM) and cite this dataset appropriately.",
    "acknowledgment": "Landsat data courtesy of the U.S. Geological Survey. MODIS data courtesy of NASA. SRTM data courtesy of NASA/JPL. Processed using Google Earth Engine.",
    
    # Data access
    "distribution_statement": "Data are available for research and educational purposes",
    
    # Technical specifications
    "cdm_data_type": "Grid",
    "data_structure": "Multidimensional array (time, y, x) with uniform spatial grid and regular temporal sampling",
    "coordinate_reference_system": "EPSG:4326 (WGS84)",
    "compression": "NetCDF-4 with zlib compression (level 4)",
    
    # References
    "references": """Landsat Collection 2: https://www.usgs.gov/landsat-missions/landsat-collection-2
MODIS Land Surface Temperature: https://lpdaac.usgs.gov/products/mod11a1v061/
SRTM: https://www2.jpl.nasa.gov/srtm/
MODIS Land Cover: https://lpdaac.usgs.gov/products/mcd12q1v061/
JRC Global Surface Water: https://global-surface-water.appspot.com/
Google Earth Engine: https://earthengine.google.com/""",
    
    # Additional comments
    "comment": """This dataset provides comprehensive multi-sensor surface temperature observations for the Rome metropolitan area. 
The hybrid LSSTSST product offers seamless temperature coverage across land and sea surfaces, useful for urban heat island studies, 
coastal dynamics, and land-ocean interactions. Dynamic land cover and emissivity products capture seasonal variations. 
All data are provided at consistent 30m spatial resolution with 15-day temporal sampling."""
})

# Variable-specific metadata
var_metadata = {
    # Landsat bands
    "LANDSAT_SR_B1": {
        "long_name": "Landsat Surface Reflectance Band 1 (Coastal/Aerosol)",
        "standard_name": "surface_bidirectional_reflectance",
        "units": "1",
        "wavelength": "0.43-0.45 µm",
        "valid_range": [0.0, 1.0]
    },
    "LANDSAT_SR_B2": {
        "long_name": "Landsat Surface Reflectance Band 2 (Blue)",
        "standard_name": "surface_bidirectional_reflectance",
        "units": "1",
        "wavelength": "0.45-0.51 µm",
        "valid_range": [0.0, 1.0]
    },
    "LANDSAT_SR_B3": {
        "long_name": "Landsat Surface Reflectance Band 3 (Green)",
        "standard_name": "surface_bidirectional_reflectance",
        "units": "1",
        "wavelength": "0.53-0.59 µm",
        "valid_range": [0.0, 1.0]
    },
    "LANDSAT_SR_B4": {
        "long_name": "Landsat Surface Reflectance Band 4 (Red)",
        "standard_name": "surface_bidirectional_reflectance",
        "units": "1",
        "wavelength": "0.64-0.67 µm",
        "valid_range": [0.0, 1.0]
    },
    "LANDSAT_SR_B5": {
        "long_name": "Landsat Surface Reflectance Band 5 (NIR)",
        "standard_name": "surface_bidirectional_reflectance",
        "units": "1",
        "wavelength": "0.85-0.88 µm",
        "valid_range": [0.0, 1.0]
    },
    "LANDSAT_SR_B6": {
        "long_name": "Landsat Surface Reflectance Band 6 (SWIR1)",
        "standard_name": "surface_bidirectional_reflectance",
        "units": "1",
        "wavelength": "1.57-1.65 µm",
        "valid_range": [0.0, 1.0]
    },
    "LANDSAT_SR_B7": {
        "long_name": "Landsat Surface Reflectance Band 7 (SWIR2)",
        "standard_name": "surface_bidirectional_reflectance",
        "units": "1",
        "wavelength": "2.11-2.29 µm",
        "valid_range": [0.0, 1.0]
    },
    "LANDSAT_ST_B10": {
        "long_name": "Landsat Surface Temperature Band 10 (Thermal Infrared)",
        "units": "K",
        "wavelength": "10.6-11.19 µm",
        "valid_range": [200.0, 400.0]
    },
    "LANDSAT_QA_PIXEL": {
        "long_name": "Landsat Quality Assessment Pixel Flags",
        "description": "Bit-packed quality flags including cloud, cloud shadow, snow/ice detection"
    },
    
    # Temperature products
    "LST_LANDSAT": {
        "long_name": "Land Surface Temperature from Landsat",
        "standard_name": "surface_temperature",
        "units": "K",
        "valid_range": [250.0, 350.0],
        "description": "Land surface temperature derived from Landsat thermal band with emissivity correction"
    },
    "LST_MODIS": {
        "long_name": "Land Surface Temperature from MODIS",
        "standard_name": "surface_temperature",
        "units": "K",
        "valid_range": [250.0, 350.0],
        "description": "Land surface temperature from MODIS (average of day and night observations)"
    },
    "SST_LANDSAT": {
        "long_name": "Sea Surface Temperature from Landsat",
        "standard_name": "sea_surface_temperature",
        "units": "K",
        "valid_range": [270.0, 320.0],
        "description": "Sea surface temperature derived from Landsat thermal band"
    },
    "SST_MODIS": {
        "long_name": "Sea Surface Temperature from MODIS",
        "standard_name": "sea_surface_temperature",
        "units": "K",
        "valid_range": [270.0, 320.0],
        "description": "Sea surface temperature from MODIS (average of day and night observations)"
    },
    "LSSTSST": {
        "long_name": "Hybrid Land and Sea Surface Temperature",
        "standard_name": "surface_temperature",
        "units": "K",
        "valid_range": [250.0, 350.0],
        "description": "Seamless temperature product: LST over land surfaces, SST over water bodies, integrated using water mask classification",
        "methodology": "Conditional: LSSTSST = SST_LANDSAT where WATERMASK=1, LST_LANDSAT where WATERMASK=0"
    },
    
    # Ancillary data
    "WATERMASK": {
        "long_name": "Water Mask",
        "flag_values": [0, 1],
        "flag_meanings": "land water",
        "description": "Binary land/water classification: 0=land, 1=water (seas, oceans, rivers, lakes, wetlands)",
        "source": "JRC Global Surface Water + ESA WorldCover"
    },
    "DEM": {
        "long_name": "Digital Elevation Model",
        "standard_name": "height_above_reference_ellipsoid",
        "units": "m",
        "description": "Elevation above WGS84 ellipsoid from SRTM 30m",
        "source": "NASA Shuttle Radar Topography Mission (SRTM) v3"
    },
    "EMISSIVITY": {
        "long_name": "Land Surface Emissivity",
        "units": "1",
        "valid_range": [0.90, 1.0],
        "description": "Surface emissivity in thermal infrared, derived from MODIS land cover classification",
        "note": "Dynamic product - varies with land cover changes over time"
    },
    "LANDCOVER": {
        "long_name": "Land Cover Classification",
        "standard_name": "land_cover_lccs",
        "description": "MODIS Land Cover Type (IGBP classification scheme)",
        "source": "MODIS MCD12Q1 Collection 6.1",
        "classification": "1=Evergreen Needleleaf, 2=Evergreen Broadleaf, 3=Deciduous Needleleaf, 4=Deciduous Broadleaf, 5=Mixed Forests, 6=Closed Shrublands, 7=Open Shrublands, 8=Woody Savannas, 9=Savannas, 10=Grasslands, 11=Permanent Wetlands, 12=Croplands, 13=Urban and Built-up, 14=Cropland/Natural Vegetation Mosaics, 15=Snow and Ice, 16=Barren, 17=Water Bodies",
        "note": "Dynamic product - captures seasonal and inter-annual land cover changes"
    }
}

# Apply variable metadata
for var, attrs in var_metadata.items():
    if var in final:
        final[var].attrs.update(attrs)

# MODIS variables metadata
modis_vars = {
    "MODIS_QC_Day": "MODIS Daytime Quality Control Flags",
    "MODIS_QC_Night": "MODIS Nighttime Quality Control Flags",
    "MODIS_Emis_31": "MODIS Band 31 Emissivity (10.78-11.28 µm)",
    "MODIS_Emis_32": "MODIS Band 32 Emissivity (11.77-12.27 µm)",
    "MODIS_View_angle_Day": "MODIS Daytime View Zenith Angle (degrees)",
    "MODIS_View_angle_Night": "MODIS Nighttime View Zenith Angle (degrees)",
    "MODIS_View_time_Day": "MODIS Daytime Observation Time (hours)",
    "MODIS_View_time_Night": "MODIS Nighttime Observation Time (hours)"
}

for var, long_name in modis_vars.items():
    if var in final:
        final[var].attrs['long_name'] = long_name
        final[var].attrs['source'] = "MODIS MOD11A1/MYD11A1 Collection 6.1"

logger.info("  ✓ Comprehensive metadata added")

# ============================================
# SAVE NETCDF
# ============================================

logger.info(f"\nSaving NetCDF: {out_nc}")

encoding = {}
for var in final.data_vars:
    chunk_y = min(256, final.sizes.get('y', 256))
    chunk_x = min(256, final.sizes.get('x', 256))
    
    encoding[var] = {
        "zlib": True,
        "complevel": 4,
        "chunksizes": (10, chunk_y, chunk_x)
    }

logger.info("Writing NetCDF...")
with dask.config.set(scheduler='threads', num_workers=2):
    final.to_netcdf(out_nc, encoding=encoding, mode='w', compute=True)

file_size_mb = Path(out_nc).stat().st_size / (1024**2)

# ============================================
# SUMMARY
# ============================================

logger.info("\n" + "="*70)
logger.info("NETCDF CREATION COMPLETE")
logger.info("="*70)
logger.info(f"File: {out_nc}")
logger.info(f"Size: {file_size_mb:.2f} MB")
logger.info(f"Timesteps: {len(final.time)}")
logger.info(f"Variables: {len(final.data_vars)}")
logger.info("\nFinal NetCDF Structure:")
logger.info("  1. LANDSAT_* - All Landsat bands (SR_B1-B7, ST_B10, QA_PIXEL)")
logger.info("  2. MODIS_* - All MODIS bands (QC, Emis, View angles/times)")
logger.info("  3. LST_LANDSAT - Landsat LST")
logger.info("  4. LST_MODIS - MODIS LST (average day/night)")
logger.info("  5. SST_LANDSAT - Landsat SST")
logger.info("  6. SST_MODIS - MODIS SST (average day/night)")
logger.info("  7. LSSTSST - Hybrid (LST on land, SST on water)")
logger.info("  8. WATERMASK - Water mask (all dates)")
logger.info("  9. DEM - Elevation (all dates)")
logger.info(" 10. EMISSIVITY - Surface emissivity (all dates)")
logger.info(" 11. LANDCOVER - Land cover classification (all dates)")
logger.info("="*70)
logger.info("✓ All processing complete!")
logger.info("="*70)