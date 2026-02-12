# Multi-Sensor LST/SST Download Pipeline

## Overview

Automated pipeline for downloading and processing multi-sensor satellite data (Landsat + MODIS) over Rome metropolitan area (2007-2021) and converting to analysis-ready NetCDF format.

**Key Features:**
- Multi-sensor integration (Landsat 5/7/8, MODIS Terra/Aqua)
- 30m spatial resolution with 15-day temporal composites
- 32+ Landsat features + 19+ MODIS features for ML applications
- Yearly NetCDF output for manageable file sizes
- Comprehensive CF-1.8 compliant metadata
- Pure satellite data (no reanalysis) for operational deployment
- Parallel processing support via LSF job scheduler

**Current Status:** Production-ready with optimized tile-based download strategy

---

## Table of Contents

1. [Requirements](#requirements)
2. [Installation](#installation)
3. [Quick Start](#quick-start)
4. [Pipeline Components](#pipeline-components)
5. [Output Structure](#output-structure)
6. [Configuration](#configuration)
7. [Data Products](#data-products)
8. [Usage Examples](#usage-examples)
9. [Metadata Structure](#metadata-structure)
10. [Troubleshooting](#troubleshooting)
11. [Performance Notes](#performance-notes)

---

## Requirements

### Software Dependencies
- Python 3.8+
- Google Earth Engine account with project setup
- 100-150 GB free disk space
- 16-32 GB RAM (recommended)

### Python Packages
See `requirements.txt`:
```
earthengine-api>=0.1.350
numpy>=1.21.0
pandas>=1.3.0
xarray>=2023.1.0
rioxarray>=0.13.0
pyyaml>=6.0
requests>=2.27.0
rasterio>=1.3.0
netCDF4>=1.6.0
dask>=2023.1.0
tqdm>=4.65.0
```

---

## Installation

### 1. Clone Repository
```bash
git clone <repository-url>
cd LST_download_pipeline
```

### 2. Install Dependencies
```bash
pip install -r requirements.txt
```

### 3. Configure Google Earth Engine
```bash
# Authenticate (first time only)
earthengine authenticate

# Set project (update with your project ID)
# Edit config_refactored.yaml:
# ee_project: 'your-project-id'
```

### 4. Verify Installation
```bash
python -c "import ee; ee.Initialize(project='your-project-id'); print('GEE initialized successfully')"
```

---

## Quick Start

### Parallel Processing (Recommended for HPC)
```bash
# Submit parallel jobs for all years (2007-2021)
chmod +x submit_parallel.sh
./submit_parallel.sh

# Monitor progress
bjobs
tail -f logs/multisensor_pipeline_2021.log
```

### Single Year Processing
```bash
# Process specific year
python run_pipeline_yearly.py --year 2021

# Skip static layers (after first run)
python run_pipeline_yearly.py --year 2020 --skip-static --skip-validation
```

### Step-by-Step Manual Execution
```bash
# 1. Download static layers (run once)
python download_static_layers.py

# 2. Download Landsat data for specific year
python download_landsat.py --year 2021

# 3. Download MODIS data for specific year
python download_modis.py --year 2021

# 4. Validate files (optional but recommended)
python validate_tiffs.py --year 2021 --delete

# 5. Create NetCDF for specific year
python make_netcdf_yearly_v2.py --year 2021

# 6. Visualize results
python plot_yearly.py --year 2021 --mode overview
```

**What This Does:**
1. Downloads static layers (DEM, water mask) - **once**
2. Downloads Landsat composites (25 per year) - **~7 min/year** with parallel
3. Downloads MODIS composites (25 per year) - **~7 min/year** with parallel
4. Validates and removes corrupted files
5. Creates yearly NetCDF files (~1-2 GB each)
6. Optional: Creates visualization plots

---

## Pipeline Components

### Core Scripts

#### 1. `config_refactored.yaml`
Central configuration file containing:
- Temporal range (2007-2021, 15-day composites)
- Spatial extent (Rome area: 11.6-13.2°E, 41.4-42.4°N)
- Sensor configurations (Landsat 5/7/8, MODIS Terra/Aqua)
- Download parameters (tile_size: 0.2°, delays: 0.1s)
- Processing parameters (median compositing, cloud masking)

**Key Settings:**
```yaml
spatial:
  tile_size: 0.2  # degrees (~22km) - optimized for GEE limits
download:
  delay_between_tiles: 0.1  # seconds - 10x faster than default
  timeout_seconds: 300
logging:
  level: 'DEBUG'  # Detailed monitoring enabled
```

#### 2. `download_static_layers.py`
Downloads time-invariant data:
- SRTM DEM (30m elevation)
- Combined water mask (JRC + GLCF + ESA)

**Runtime:** 10-20 minutes  
**Output:** `static/DEM_30m.tif`, `static/WaterMask_30m.tif`

#### 3. `download_landsat.py`
Downloads Landsat data with full processing chain:
- **Sensors:** Landsat 5 (2007-2011), 7 (2007-2021), 8 (2013-2021)
- **Processing:** Band harmonization, cloud masking, LST/SST calculation
- **Features:** 32 bands including spectral indices, temperature products
- **Strategy:** Tile-based download (~170 tiles/composite with tile_size=0.2°)
- **Year filter:** `--year YYYY` for parallel processing

**Runtime:** ~7-10 min/year (with optimized delays)  
**Output:** ~25 GeoTIFF files/year (15-day composites)

**Key Features:**
- Adaptive cloud cover filtering (relaxes threshold if <20 images available)
- Automatic retry on timeout
- Skip existing files
- Detailed tile-level logging

#### 4. `download_modis.py`
Downloads MODIS LST products:
- **Collections:** MOD11A1 (Terra), MYD11A1 (Aqua)
- **Bands:** Day/Night LST, emissivity, QC, view geometry
- **Features:** 19 bands total
- **Resampling:** 1km → 30m for consistency with Landsat

**Runtime:** ~7-10 min/year  
**Output:** ~25 GeoTIFF files/year

#### 5. `validate_tiffs.py`
Validates GeoTIFF integrity:
- File size checks (>1KB)
- Band count validation
- Readability tests
- Optional year filtering: `--year 2021`

```bash
# List corrupted files for specific year
python validate_tiffs.py --year 2021

# Remove corrupted files
python validate_tiffs.py --year 2021 --delete
```

#### 6. `make_netcdf_yearly_v2.py`
Creates yearly NetCDF files with proper band naming:
- **Naming:** Descriptive band names (e.g., `LANDSAT_NDVI`, `MODIS_LST_Day_K`)
- **Interpolation:** DEM/WATERMASK interpolated to match data grid
- **Metadata:** Complete CF-1.8/ACDD-1.3 attributes
- **Compression:** Level 4 for optimal size/speed

**Runtime:** 5-10 minutes/year  
**Output:** `netcdf_yearly/Rome_LST_SST_YYYY_30m_15day.nc` (~1-2 GB)

**Key Improvements:**
- Band names preserved from GeoTIFF descriptions
- Static layers included in all timesteps
- Config parameters in metadata

#### 7. `run_pipeline_yearly.py`
Year-specific pipeline orchestrator:
- Executes download → validate → netcdf for one year
- Options: `--year`, `--skip-static`, `--skip-validation`
- Separate logs per year
- Error handling and recovery

**Usage:**
```bash
python run_pipeline_yearly.py --year 2021
python run_pipeline_yearly.py --year 2020 --skip-static
```

#### 8. `submit_parallel.sh`
LSF job submission script for HPC:
- Submits 15 parallel jobs (one per year)
- First job downloads static layers
- Subsequent jobs skip static download
- Isolated logs: `logs/out/pipeline_YYYY.out`

**Customization:**
```bash
# Edit year range in script
for year in {2007..2021}; do
  # Submit job
done
```

#### 9. `plot_yearly.py`
Visualization toolkit:
- **Modes:** map (single variable), overview (4-panel), timeseries
- **Variables:** Any NetCDF variable (e.g., `LANDSAT_LST_K`, `LANDSAT_NDVI`)
- **Output:** PDF files

```bash
python plot_yearly.py --year 2021 --mode map --variable LANDSAT_LST_K
python plot_yearly.py --year 2021 --mode overview
python plot_yearly.py --year 2021 --mode timeseries --variable LANDSAT_NDVI
```

#### 10. `cleanup_tiles.py`
Removes temporary tile directories:
```bash
python cleanup_tiles.py  # Removes all tiles_temp_* directories
```

---

## Download Strategy

### Tile-Based Approach
The pipeline uses a tile-based download strategy to work within GEE memory limits:

**Configuration:**
- `tile_size: 0.2°` → ~170 tiles per composite
- Each tile: ~22km × 22km at 30m resolution
- Total pixels per tile: ~730 × 730 = 533,000 pixels
- Request size per tile: ~20-30 MB (within 50 MB GEE limit)

**Why Not Larger Tiles?**
- `tile_size: 0.5°` → 45 tiles but **179 MB/tile** → exceeds GEE 50MB limit
- `tile_size: 1.0°` → 2-4 tiles but **>200 MB/tile** → fails

**Why Not Band-by-Band?**
Median/mean compositing loses original band names. The composite only preserves:
- Coordinates (lon, lat, time)
- Dynamic layers (land_cover, emissivity, AOD)

All Landsat bands (SR_B1-B7, ST_B10, indices, LST) are lost after compositing.

**Current Solution:**
Download complete composite (all 32 bands) in small tiles, then mosaic locally with GDAL.

---

## Output Structure

```
LST_Rome_30m_15day/
├── landsat/                    # Landsat GeoTIFF files
│   ├── Landsat_2007-01-01_0000.tif
│   ├── Landsat_2007-01-16_0001.tif
│   └── ... (~365 files)
│
├── modis/                      # MODIS GeoTIFF files
│   ├── MODIS_2007-01-01_0000.tif
│   └── ... (~365 files)
│
├── static/                     # Static layers
│   ├── DEM_30m.tif
│   └── WaterMask_30m.tif
│
├── netcdf_yearly/              # Yearly NetCDF files
│   ├── Rome_LST_SST_2007_30m_15day.nc (~1-2 GB)
│   ├── Rome_LST_SST_2008_30m_15day.nc
│   ├── ...
│   └── Rome_LST_SST_2021_30m_15day.nc
│
└── Rome_MultiSensor_LST_SST_2007_2021.nc  # Combined file (~15-20 GB)
```

---

## Configuration

### Key Parameters in `config_refactored.yaml`

#### Temporal Configuration
```yaml
temporal:
  start_date: '2007-01-01'
  end_date: '2021-12-31'
  interval_days: 15  # 15-day composites → 25 per year
```

#### Spatial Configuration
```yaml
geometry:
  min_lon: 11.6  # Western boundary (Tyrrhenian Sea)
  min_lat: 41.4  # Southern boundary
  max_lon: 13.2  # Eastern boundary
  max_lat: 42.4  # Northern boundary

spatial:
  target_resolution: 30  # meters
  tile_size: 0.2  # degrees - optimized for GEE memory limits
```

**Area Coverage:** ~176 km (E-W) × ~111 km (N-S) = 19,536 km²

#### Landsat Sensors
```yaml
landsat:
  sensors:
    - collection: 'LANDSAT/LT05/C02/T1_L2'
      start_date: '2007-01-01'
      end_date: '2011-12-31'
      name: 'L5'
    
    - collection: 'LANDSAT/LE07/C02/T1_L2'
      start_date: '2007-01-01'
      end_date: '2021-12-31'
      name: 'L7'
    
    - collection: 'LANDSAT/LC08/C02/T1_L2'
      start_date: '2013-04-01'
      end_date: '2021-12-31'
      name: 'L8'
  
  max_cloud_cover: 80  # Relaxed to 160 if <20 images available
  
  path_row_filter:
    enabled: true
    paths: [191]  # WRS-2 path covering Rome
    rows: []
```

#### MODIS Configuration
```yaml
modis:
  collections:
    terra: 'MODIS/061/MOD11A1'  # Collection 6.1 (updated from 006)
    aqua: 'MODIS/061/MYD11A1'
```

#### Download Parameters (Optimized)
```yaml
download:
  timeout_seconds: 600  # Increased for larger tiles
  delay_between_tiles: 0.1  # Reduced from 1s → 10x faster
  delay_between_images: 0.1  # Reduced from 2s
  min_tiles_success_ratio: 0.10  # At least 10% tiles required
  max_retries: 3
```

**Performance Impact:**
- Old: 1s delay × 170 tiles × 25 composites = 71 min/year
- New: 0.1s delay × 170 tiles × 25 composites = 7 min/year
- **~10x speedup** (delays only, actual download time additional)

#### Dynamic Layers
```yaml
dynamic_layers:
  land_cover:
    enabled: false  # Disabled to reduce memory usage
    # When enabled: ESA WorldCover (2020-2021) + MODIS fallback (2007-2019)
  
  emissivity:
    enabled: true  # Calculated from MODIS land cover
```

#### Logging (Enhanced)
```yaml
logging:
  level: 'DEBUG'  # Detailed tile-level monitoring
  log_file: 'multisensor_pipeline.log'
  console_colors: true
```

**Log Output Includes:**
- Estimated download time and tile count
- Per-tile success/failure with file sizes
- Band count and names in each composite
- Download statistics (tiles/composite, success rate)
- GEE error messages with full context

---

## Data Products

### Landsat Variables (32 total)

**Variable naming:** `LANDSAT_<band>` (e.g., `LANDSAT_SR_B1`, `LANDSAT_NDVI`)

#### Spectral Bands (9)
- `LANDSAT_SR_B1`: Coastal/Aerosol (0.43-0.45 µm)
- `LANDSAT_SR_B2`: Blue (0.45-0.51 µm)
- `LANDSAT_SR_B3`: Green (0.53-0.59 µm)
- `LANDSAT_SR_B4`: Red (0.64-0.67 µm)
- `LANDSAT_SR_B5`: NIR (0.85-0.88 µm)
- `LANDSAT_SR_B6`: SWIR1 (1.57-1.65 µm)
- `LANDSAT_SR_B7`: SWIR2 (2.11-2.29 µm)
- `LANDSAT_ST_B10`: Thermal Infrared (10.6-11.19 µm)
- `LANDSAT_QA_PIXEL`: Quality assessment flags

#### Spectral Indices (10)
- `LANDSAT_NDVI`: Normalized Difference Vegetation Index
- `LANDSAT_NDWI`: Normalized Difference Water Index
- `LANDSAT_NDBI`: Normalized Difference Built-up Index
- `LANDSAT_MNDWI`: Modified NDWI (water detection)
- `LANDSAT_SAVI`: Soil Adjusted Vegetation Index
- `LANDSAT_EVI`: Enhanced Vegetation Index
- `LANDSAT_BSI`: Bare Soil Index
- `LANDSAT_UI`: Urban Index
- `LANDSAT_Albedo`: Surface albedo
- `LANDSAT_FV`: Fractional Vegetation cover

#### Temperature Products (4)
- `LANDSAT_LST_K`: Land Surface Temperature (Kelvin)
- `LANDSAT_LST_K_masked`: LST with cloud mask applied
- `LANDSAT_SST_K`: Sea Surface Temperature (Kelvin)
- `LANDSAT_LSSTSST`: Hybrid land-sea temperature

#### Additional Products (9)
- `LANDSAT_emissivity`: Surface emissivity (calculated)
- `LANDSAT_emissivity_dynamic`: From land cover lookup
- `LANDSAT_cloud_mask`: Binary cloud/shadow/snow mask
- `LANDSAT_lon`: Longitude (degrees)
- `LANDSAT_lat`: Latitude (degrees)
- `LANDSAT_time_days`: Days since composite start
- `LANDSAT_land_cover`: MODIS land cover (when enabled)
- `LANDSAT_AOD_047`: Aerosol Optical Depth at 470 nm
- `LANDSAT_AOD_055`: Aerosol Optical Depth at 550 nm
- `MODIS_TPW`: Total Precipitable Water

### MODIS Variables (14+)

#### Temperature (4)
- `MODIS_LST_Day_K`: Daytime LST
- `MODIS_LST_Night_K`: Nighttime LST
- `MODIS_SST_Day_K`: Daytime SST
- `MODIS_SST_Night_K`: Nighttime SST

#### Derived Products (2)
- `LST_MODIS`: Mean LST (day+night average)
- `SST_MODIS`: Mean SST (day+night average)

#### Quality and Geometry (8)
- `MODIS_QC_Day`: Daytime quality flags
- `MODIS_QC_Night`: Nighttime quality flags
- `MODIS_Emis_31`: Band 31 emissivity
- `MODIS_Emis_32`: Band 32 emissivity
- `MODIS_View_angle_Day`: Daytime view angle
- `MODIS_View_angle_Night`: Nighttime view angle
- `MODIS_View_time_Day`: Daytime observation time
- `MODIS_View_time_Night`: Nighttime observation time

### Static Layers (2)

- `DEM`: Elevation (SRTM 30m)
- `WATERMASK`: Land/water classification (binary)

### Dynamic Layers (2)

- `LANDCOVER`: ESA WorldCover 10m (2020-2021) + MODIS fallback (2007-2019)
- `EMISSIVITY`: Surface emissivity from land cover

---

## Usage Examples

### Loading Yearly Data

```python
import xarray as xr
import matplotlib.pyplot as plt

# Load single year
ds = xr.open_dataset('netcdf_yearly/Rome_LST_SST_2015_30m_15day.nc')

# Access variables
lst = ds['LST_LANDSAT']
ndvi = ds['NDVI']
dem = ds['DEM']

# Temporal subset
summer = ds.sel(time=slice('2015-06-01', '2015-08-31'))

# Spatial subset (latitude/longitude)
urban_core = ds.sel(x=slice(12.4, 12.6), y=slice(41.8, 42.0))

# Plot mean LST
ds['LST_LANDSAT'].mean(dim='time').plot()
plt.title('Mean LST 2015')
plt.show()
```

### Loading Combined Data

```python
import xarray as xr

# Load entire time series (requires sufficient RAM)
ds_all = xr.open_dataset('Rome_MultiSensor_LST_SST_2007_2021.nc')

# Calculate long-term trend
lst_trend = ds_all['LST_LANDSAT'].mean(dim=['x', 'y'])
lst_trend.plot()

# Multi-year statistics
lst_mean = ds_all['LST_LANDSAT'].mean(dim='time')
lst_std = ds_all['LST_LANDSAT'].std(dim='time')
```

### Lazy Loading with Dask

```python
import xarray as xr

# Open with dask chunks (memory efficient)
ds = xr.open_dataset(
    'Rome_MultiSensor_LST_SST_2007_2021.nc',
    chunks={'time': 24, 'y': 256, 'x': 256}
)

# Operations are lazy
mean_lst = ds['LST_LANDSAT'].mean(dim=['x', 'y'])

# Compute only when needed
result = mean_lst.compute()
```

### Loading Multiple Years

```python
import xarray as xr
from pathlib import Path

# Load all yearly files
nc_files = sorted(Path('netcdf_yearly').glob('Rome_LST_SST_*.nc'))
datasets = [xr.open_dataset(f) for f in nc_files]

# Concatenate
ds_combined = xr.concat(datasets, dim='time')

# Time series analysis
lst_timeseries = ds_combined['LST_LANDSAT'].mean(dim=['x', 'y'])
```

### Accessing Metadata

```python
import xarray as xr

ds = xr.open_dataset('netcdf_yearly/Rome_LST_SST_2015_30m_15day.nc')

# Global metadata
print(ds.attrs['title'])
print(ds.attrs['processing_step_5'])  # LST calculation method
print(ds.attrs['source_landsat_8'])

# Variable metadata
print(ds['LST_LANDSAT'].attrs['long_name'])
print(ds['LST_LANDSAT'].attrs['units'])

# Spectral index formulas
print(ds.attrs['index_ndvi'])
print(ds.attrs['index_albedo'])
```

---

## Metadata Structure

Each NetCDF file contains comprehensive metadata following CF-1.8 and ACDD-1.3 conventions.

### Global Attributes

#### Identification
- `title`: Dataset title with year
- `institution`: CMCC Foundation
- `creator_name`, `creator_email`, `creator_institution`
- `project`: Project description

#### Temporal Coverage
- `time_coverage_start`, `time_coverage_end`
- `temporal_resolution`: "15 days"
- `time_coverage_resolution`: "P15D" (ISO 8601)

#### Spatial Coverage
- `geospatial_lat_min/max`, `geospatial_lon_min/max`
- `spatial_resolution`: "30 meters"
- `area_of_interest`: "Rome metropolitan area and Tyrrhenian coast, Italy"

#### Data Sources
- `source_landsat_5/7/8`: Landsat collection details
- `source_modis_terra/aqua`: MODIS collection details
- `source_precipitation`: GPM IMERG
- `source_aerosol`: MODIS AOD
- `source_solar_radiation`: CERES
- `source_water_vapor`: MODIS TPW
- `source_dem`: SRTM
- `source_watermask`: JRC + ESA
- `source_landcover`: ESA WorldCover + MODIS

#### Processing Steps (10 documented)
- `processing_step_1`: Data acquisition from GEE
- `processing_step_2`: Band harmonization
- `processing_step_3`: Cloud masking
- `processing_step_4`: Spectral indices calculation
- `processing_step_5`: LST calculation
- `processing_step_6`: Temporal compositing
- `processing_step_7`: Resampling
- `processing_step_8`: Multi-sensor integration
- `processing_step_9`: Tiling and mosaicking
- `processing_step_10`: NetCDF creation

#### Methodology
- `methodology_lst`: LST algorithm formula
- `methodology_emissivity`: Emissivity calculation
- `methodology_compositing`: Temporal aggregation method
- `index_*`: Formula for each spectral index

#### References
- Direct URLs to all data sources
- Algorithm references
- Processing documentation

### Variable Attributes

Each variable includes:
- `long_name`: Descriptive name
- `standard_name`: CF standard name (when applicable)
- `units`: Physical units
- `valid_range`: Expected value range
- `description`: Detailed explanation
- `methodology`: Calculation method (for derived products)
- `source`: Original data source

---

## Troubleshooting

### Common Issues

#### 1. Earth Engine Authentication Error
```bash
# Re-authenticate
earthengine authenticate

# Check project
python -c "import ee; ee.Initialize(project='your-project-id')"
```

#### 2. Out of Memory
For combined NetCDF:
```python
# Use dask chunking
ds = xr.open_dataset(file, chunks={'time': 10, 'y': 256, 'x': 256})
```

For yearly processing:
```bash
# Process years individually
python make_netcdf_yearly.py  # Creates one file at a time
```

#### 3. Corrupted GeoTIFF Files
```bash
# Identify and remove
python validate_tiffs.py --delete

# Re-download missing composites
python download_landsat.py  # Skips existing files
```

#### 4. Download Failures
- Check internet connection
- Verify GEE quota limits
- Increase timeout in config
- Reduce tile size if memory issues

#### 5. Missing Years in NetCDF
```bash
# Check if GeoTIFFs exist for that year
ls LST_Rome_30m_15day/landsat/ | grep "2015"

# Re-run NetCDF creation
python make_netcdf_yearly.py  # Skips existing files
```

---

## Performance Notes

### Disk Space Requirements
- GeoTIFF files (raw): 50-100 GB
- Yearly NetCDF: 10-30 GB (compressed)
- Combined NetCDF: 15-20 GB
- Total: 75-150 GB

### Runtime Estimates
- Static layers: 10-20 minutes
- Landsat download: 24-48 hours
- MODIS download: 24-48 hours
- Validation: 5-10 minutes
- Yearly NetCDF: 30-60 minutes
- Combined NetCDF: 60-120 minutes
- **Total pipeline: 2-4 days**

### Memory Requirements
- Download scripts: 8-16 GB RAM
- NetCDF creation: 16-32 GB RAM
- Data loading (with dask): 4-8 GB RAM
- Combined NetCDF (full load): 20-30 GB RAM

### Optimization Tips

1. **Parallel Downloads**: Not implemented but possible
2. **Chunk Size**: Adjust in config for memory vs speed
3. **Compression**: Level 4 provides good balance
4. **Dask Workers**: 2-4 workers recommended
5. **Tile Size**: Larger tiles = fewer requests but more memory

### Network Considerations
- Download speed depends on GEE server load
- Retry logic implemented for transient failures
- Rate limiting prevents quota issues
- Progress bars show real-time status

---

## Citation

If using this dataset, please cite:

```
Sebastianelli, A. (2025). Multi-Sensor Land and Sea Surface Temperature Dataset 
for Rome Metropolitan Area (2007-2021). CMCC Foundation. 
Processed using Google Earth Engine.
```

And acknowledge the data sources:
- Landsat: U.S. Geological Survey
- MODIS: NASA
- GPM: NASA/JAXA
- SRTM: NASA/JPL
- ESA WorldCover: ESA

---

## License

Individual data sources have their own licenses:
- Landsat: U.S. public domain
- MODIS: NASA open data policy
- GPM: NASA open data policy
- SRTM: Public domain
- ESA WorldCover: CC-BY 4.0

Processed dataset and code: Check with CMCC Foundation.

---

## Contact

**Alessandro Sebastianelli**  
Email: asebastianelli@ieee.org
Institution: CMCC Foundation  
Website: https://alessandrosebastianelli.github.io/

For issues, questions, or contributions, please contact the author.

---

## Acknowledgments

This work uses data from:
- U.S. Geological Survey (Landsat)
- NASA (MODIS, GPM, CERES, SRTM)
- European Space Agency (WorldCover)
- JRC (Global Surface Water)

Processing performed using Google Earth Engine platform.