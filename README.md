# Multi-Sensor LST/SST Pipeline - Rome Area

Complete refactored pipeline for downloading and processing Land Surface Temperature (LST) and Sea Surface Temperature (SST) data from multiple sensors for the Rome metropolitan area (including coastal waters).

## Overview

This pipeline downloads and combines:

1. **Landsat 8/9** - 15-day composites at 30m with ALL bands + derived LST/SST
2. **MODIS Terra/Aqua** - 15-day composites resampled from 1km to 30m with ALL bands
3. **Static Layers** (time-independent):
   - DEM (SRTM 30m)
   - Land Cover (ESA WorldCover 2021)
   - Water Mask (combined from JRC Global Surface Water + ESA)
   - Emissivity (derived from land cover)
4. **Hybrid LST/SST Layer** - Uses LST over land pixels, SST over water pixels

**NO daily data** - only 15-day composites for consistency.

## Quick Start

```bash
# 1. Install dependencies
pip install -r requirements.txt

# 2. Authenticate with Earth Engine
earthengine authenticate

# 3. Edit configuration
nano config_refactored.yaml

# 4. Run complete pipeline
python run_pipeline.py
```

This will:
- Download static layers (DEM, land cover, water mask, emissivity)
- Download Landsat 15-day composites with all bands
- Download MODIS 15-day composites with all bands
- Validate files (with 2 re-download cycles)
- Create combined NetCDF with hybrid LST/SST layer

## Output Structure

```
LST_Rome_30m_15day/
├── static/
│   ├── DEM.tif
│   ├── LandCover.tif
│   ├── WaterMask.tif
│   └── Emissivity.tif
│
├── landsat/
│   ├── Landsat_2020-01-01_0000.tif  (18 bands each)
│   ├── Landsat_2020-01-16_0001.tif
│   └── ...
│
├── modis/
│   ├── MODIS_2020-01-01_0000.tif    (15 bands each)
│   ├── MODIS_2020-01-16_0001.tif
│   └── ...
│
└── Rome_MultiSensor_LST_SST_2020_2025.nc  (Combined NetCDF)
```

## Data Products

### Landsat Files (18 bands)

**Original Bands:**
- `SR_B1-B7`: Surface Reflectance (Coastal, Blue, Green, Red, NIR, SWIR1, SWIR2)
- `ST_B10`: Thermal Infrared
- `QA_PIXEL`: Quality flags

**Derived Products:**
- `LST_K`: Land Surface Temperature (K)
- `LST_K_masked`: LST with cloud mask applied
- `SST_K`: Sea Surface Temperature (K)
- `cloud_mask`: Binary mask (1=clear, 0=cloud/shadow/snow)
- `NDVI`: Normalized Difference Vegetation Index
- `emissivity`: Surface emissivity

**Coordinates:**
- `lon`, `lat`, `time`

### MODIS Files (15 bands)

**All MODIS bands resampled from 1km to 30m:**
- `MODIS_LST_Day_K`, `MODIS_LST_Night_K`
- `MODIS_SST_Day_K`, `MODIS_SST_Night_K`
- `MODIS_QC_Day`, `MODIS_QC_Night`
- `MODIS_Emis_31`, `MODIS_Emis_32` (Band 31/32 emissivity)
- `MODIS_View_angle_Day`, `MODIS_View_angle_Night`
- `MODIS_View_time_Day`, `MODIS_View_time_Night`
- `lon`, `lat`, `time`

### Static Layers

- **DEM**: SRTM 30m elevation (meters)
- **Land Cover**: ESA WorldCover 2021 classification
  - 10: Tree cover, 20: Shrubland, 30: Grassland, 40: Cropland
  - 50: Built-up, 60: Bare/sparse vegetation, 70: Snow and ice
  - 80: Permanent water, 90: Wetland, 95: Mangroves, 100: Moss/lichen
- **Water Mask**: Binary (0=land, 1=water)
  - Combines JRC Global Surface Water + ESA water classes
  - Includes: seas, oceans, rivers, lakes, wetlands
- **Emissivity**: Surface emissivity derived from land cover (0.93-0.99)

### NetCDF Combined Dataset

```python
import xarray as xr
ds = xr.open_dataset('Rome_MultiSensor_LST_SST_2020_2025.nc')

# Hybrid layer (KEY FEATURE)
ds['LST_SST_hybrid']  # LST over land, SST over water

# Landsat variables (with 'landsat_' prefix)
ds['landsat_SR_B1']  # Blue band
ds['landsat_LST_K']  # Land Surface Temperature
ds['landsat_SST_K']  # Sea Surface Temperature
# ... etc for all 18 Landsat bands

# MODIS variables (with 'MODIS_' prefix)
ds['MODIS_LST_Day_K']
ds['MODIS_SST_Night_K']
# ... etc for all 15 MODIS bands

# Static layers (broadcast to all timesteps)
ds['dem']
ds['land_cover']
ds['water_mask']
ds['emissivity']
```

## Key Features

### 1. Hybrid LST/SST Layer

The pipeline creates a seamless temperature layer:
- **Land pixels**: Uses Landsat LST
- **Water pixels**: Uses Landsat SST  
- **Separation**: Based on water mask (JRC + ESA)

This gives you:
- Rome city temperature (LST)
- Tiber River temperature (SST)
- Tyrrhenian Sea temperature (SST)
- All in one continuous layer!

### 2. Comprehensive Band Coverage

Unlike simplified pipelines:
- **All** Landsat bands preserved (SR_B1-7, ST_B10, QA_PIXEL)
- **All** MODIS bands included (LST, QC, emissivity, view geometry)
- Full spectral information for advanced analysis

### 3. Multi-Source Water Mask

Water mask combines:
- JRC Global Surface Water (>50% occurrence)
- ESA WorldCover water classes (permanent water, wetlands, mangroves)
- Captures: Mediterranean Sea, Tiber River, lakes, wetlands

## Scripts

- **config_refactored.yaml** - Main configuration file
- **download_static_layers.py** - Download DEM, land cover, water mask, emissivity
- **download_landsat.py** - Download Landsat 15-day composites
- **download_modis.py** - Download MODIS 15-day composites
- **validate_tiffs.py** - Validate and clean corrupted files
- **make_netcdf.py** - Create combined NetCDF with hybrid LST/SST
- **run_pipeline.py** - Master orchestrator
- **requirements.txt** - Python dependencies

## Requirements

- Python 3.8+
- Google Earth Engine account
- ~50GB disk space
- 8GB RAM minimum (16GB recommended)

## Contact

Alessandro Sebastianelli - CMCC Foundation
