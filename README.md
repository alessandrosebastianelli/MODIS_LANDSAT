# Multi-Sensor LST/SST Pipeline - Rome Metropolitan Area

**Author:** Alessandro Sebastianelli  
**Institution:** CMCC Foundation (Centro Euro-Mediterraneo sui Cambiamenti Climatici)  
**Website:** https://alessandrosebastianelli.github.io/  
**Contact:** alessandro.sebastianelli@cmcc.it

---

## Overview

This pipeline provides a comprehensive system for downloading, processing, and integrating multi-sensor Land Surface Temperature (LST) and Sea Surface Temperature (SST) data for the Rome metropolitan area, including coastal Mediterranean waters.

### Key Features

- ✅ **Multi-sensor integration**: Landsat 8/9 + MODIS Terra/Aqua
- ✅ **High resolution**: 30m spatial resolution, 15-day temporal frequency
- ✅ **Complete spectral coverage**: All Landsat and MODIS bands preserved
- ✅ **Hybrid LST/SST product**: Seamless land-water temperature integration
- ✅ **Dynamic ancillary data**: Land cover and emissivity at 15-day frequency
- ✅ **CF-compliant**: NetCDF-4 following CF-1.8 and ACDD-1.3 standards
- ✅ **Comprehensive metadata**: Full data provenance and processing documentation
- ✅ **Automated workflow**: One-command execution with validation cycles

---

## Final NetCDF Structure

The pipeline produces a single NetCDF file with the following variables, **all aligned to the same temporal dimension**:

### 1. LANDSAT - All Spectral Bands (time, y, x)
- `LANDSAT_SR_B1` - Coastal/Aerosol (0.43-0.45 µm)
- `LANDSAT_SR_B2` - Blue (0.45-0.51 µm)
- `LANDSAT_SR_B3` - Green (0.53-0.59 µm)
- `LANDSAT_SR_B4` - Red (0.64-0.67 µm)
- `LANDSAT_SR_B5` - NIR (0.85-0.88 µm)
- `LANDSAT_SR_B6` - SWIR1 (1.57-1.65 µm)
- `LANDSAT_SR_B7` - SWIR2 (2.11-2.29 µm)
- `LANDSAT_ST_B10` - Thermal Infrared (10.6-11.19 µm)
- `LANDSAT_QA_PIXEL` - Quality assessment flags

### 2. MODIS - All Relevant Bands (time, y, x)
- `MODIS_QC_Day` / `MODIS_QC_Night` - Quality control flags
- `MODIS_Emis_31` / `MODIS_Emis_32` - Band 31/32 emissivity
- `MODIS_View_angle_Day` / `MODIS_View_angle_Night` - View zenith angles
- `MODIS_View_time_Day` / `MODIS_View_time_Night` - Observation times

### 3. Temperature Products (time, y, x)
- `LST_LANDSAT` - Land Surface Temperature from Landsat (K)
- `LST_MODIS` - Land Surface Temperature from MODIS (K) [average day/night]
- `SST_LANDSAT` - Sea Surface Temperature from Landsat (K)
- `SST_MODIS` - Sea Surface Temperature from MODIS (K) [average day/night]
- `LSSTSST` - **Hybrid product**: LST over land + SST over water (K)

### 4. Ancillary Data (time, y, x)
- `WATERMASK` - Binary land/water mask (0=land, 1=water) [broadcast to all dates]
- `DEM` - Digital Elevation Model (m) [broadcast to all dates]
- `EMISSIVITY` - Surface emissivity (0.93-0.99) [dynamic, 15-day frequency]
- `LANDCOVER` - MODIS Land Cover classification [dynamic, 15-day frequency]

---

## Data Sources

### Satellite Data
- **Landsat 8/9**: Surface Reflectance + Thermal (Collection 2 Level-2)
  - `LANDSAT/LC08/C02/T1_L2` (2020-2021)
  - `LANDSAT/LC09/C02/T1_L2` (2021-2025)
  - Native resolution: 30m
  - Platforms: Landsat-8, Landsat-9
  - Sensors: OLI/TIRS, OLI-2/TIRS-2

- **MODIS Terra/Aqua**: Land Surface Temperature (Collection 6.1)
  - `MODIS/006/MOD11A1` (Terra, ~10:30 AM overpass)
  - `MODIS/006/MYD11A1` (Aqua, ~1:30 PM overpass)
  - Native resolution: 1km (resampled to 30m)

### Ancillary Data
- **DEM**: SRTM v3 (`USGS/SRTMGL1_003`, 30m)
- **Land Cover**: MODIS MCD12Q1 (`MODIS/061/MCD12Q1`, IGBP classification)
- **Water Mask**: JRC Global Surface Water + ESA WorldCover
  - `JRC/GSW1_4/GlobalSurfaceWater`
  - `ESA/WorldCover/v200`

---

## Processing Workflow

### 1. Data Acquisition
- Access via Google Earth Engine API
- Service account authentication (no interactive login)
- Spatial filtering: Rome metropolitan area (11.6-13.2°E, 41.4-42.4°N)
- Temporal filtering: 2020-2025
- Path/Row filtering: WRS-2 Path 191 (reduces seam lines)

### 2. Cloud Masking
- **Landsat**: QA_PIXEL band (cloud, cloud shadow, snow/ice flags)
- **MODIS**: QC flags (day/night quality assessment)

### 3. Temporal Compositing
- 15-day median composites
- Reduces cloud contamination
- Provides consistent temporal sampling
- Minimizes atmospheric effects

### 4. LST/SST Calculation
- **Algorithm**: Split-window with emissivity correction
- **Emissivity**: 
  - Landsat: Derived from NDVI and fractional vegetation cover
  - MODIS: From land cover classification lookup table
- **LST**: Applied to thermal bands (Landsat ST_B10, MODIS LST_Day/Night)
- **SST**: Same thermal algorithm applied to water bodies

### 5. Spatial Processing
- MODIS resampling: 1km → 30m (bilinear interpolation)
- Co-registration: All products aligned to common grid
- Tiled download: Large areas split into manageable tiles
- Mosaicking: Tiles merged into seamless products

### 6. Integration
- Water mask classification: JRC + ESA sources combined
- Hybrid LSSTSST: `SST_LANDSAT` where water, `LST_LANDSAT` where land
- Dynamic layers: Land cover and emissivity per timestep
- Static layers: DEM and water mask broadcast to all times

---

## Installation

### Requirements
- Python 3.8+
- Google Earth Engine account with service account credentials
- ~50GB disk space
- 8GB RAM minimum (16GB recommended)

### Setup

```bash
# 1. Clone or download pipeline files
cd /path/to/pipeline/

# 2. Install Python dependencies
pip install -r requirements.txt

# 3. Configure GEE authentication
# Edit config_refactored.yaml with your service account details:
nano config_refactored.yaml

# Update these fields:
ee_project: 'your-gee-project'
service_account:
  key_file: '/path/to/your/service-account-key.json'
  email: 'your-service-account@your-project.iam.gserviceaccount.com'
```

---

## Quick Start

### Automated Execution (Recommended)

```bash
# Run complete pipeline
python run_pipeline.py
```

This executes:
1. Download static layers (DEM, water mask)
2. Download Landsat 15-day composites
3. Download MODIS 15-day composites
4. Validate all files (2 cycles with re-download)
5. Create final NetCDF with all products

**Expected runtime**: 12-24 hours for 2020-2025 data

### Manual Execution

```bash
# Step 1: Download static layers (run once)
python download_static_layers.py
# Output: LST_Rome_30m_15day/static/DEM.tif
#         LST_Rome_30m_15day/static/WaterMask.tif

# Step 2: Download Landsat composites
python download_landsat.py
# Output: ~150 files in LST_Rome_30m_15day/landsat/
# Each file: 20 bands (9 original + 6 derived + 2 dynamic + 3 coords)

# Step 3: Download MODIS composites
python download_modis.py
# Output: ~150 files in LST_Rome_30m_15day/modis/
# Each file: 17 bands (12 original + 2 dynamic + 3 coords)

# Step 4: Validate files (optional but recommended)
python validate_tiffs.py --delete
# Checks integrity, removes corrupted files

# Step 5: Create NetCDF
python make_netcdf_final.py
# Output: Rome_MultiSensor_LST_SST_2020_2025.nc
```

---

## Configuration

Edit `config_refactored.yaml` to customize:

### Area of Interest
```yaml
geometry:
  min_lon: 11.6  # Western boundary (includes Tyrrhenian Sea)
  min_lat: 41.4  # Southern boundary
  max_lon: 13.2  # Eastern boundary
  max_lat: 42.4  # Northern boundary
```

### Temporal Coverage
```yaml
temporal:
  start_date: '2020-01-01'
  end_date: '2025-12-31'
  interval_days: 15  # Composite frequency
```

### Spatial Parameters
```yaml
spatial:
  target_resolution: 30  # meters
  tile_size: 0.4         # degrees (for tiled downloads)
```

### Cloud Filtering
```yaml
landsat:
  max_cloud_cover: 80  # percentage
  path_row_filter:
    enabled: true
    paths: [191]  # Rome WRS-2 path (reduces seam lines)
```

### Compositing Method
```yaml
compositing:
  method: 'median'  # Options: 'median', 'mean', 'mosaic'
```

---

## Usage Examples

### Loading Data

```python
import xarray as xr
import matplotlib.pyplot as plt

# Open NetCDF
ds = xr.open_dataset('Rome_MultiSensor_LST_SST_2020_2025.nc')

# Inspect structure
print(ds)
print(ds.dims)  # {'time': 150, 'y': 1000, 'x': 1500}

# View metadata
print(ds.attrs['title'])
print(ds.attrs['creator_name'])
print(ds.attrs['creator_url'])
print(ds.attrs['processing_description'])
```

### Accessing Variables

```python
# Temperature products
lst_landsat = ds['LST_LANDSAT']      # Shape: (time, y, x)
lst_modis = ds['LST_MODIS']          # Shape: (time, y, x)
sst_landsat = ds['SST_LANDSAT']      # Shape: (time, y, x)
hybrid_temp = ds['LSSTSST']          # Shape: (time, y, x) - KEY PRODUCT!

# Landsat spectral bands
blue = ds['LANDSAT_SR_B2']           # Blue band
nir = ds['LANDSAT_SR_B5']            # NIR band
thermal = ds['LANDSAT_ST_B10']       # Thermal band

# Ancillary data
water_mask = ds['WATERMASK']         # Land/water classification
dem = ds['DEM']                      # Elevation
land_cover = ds['LANDCOVER']         # Dynamic land cover
emissivity = ds['EMISSIVITY']        # Dynamic emissivity
```

### Time Series Analysis

```python
# Extract specific date
jan_2020 = ds['LSSTSST'].sel(time='2020-01-15', method='nearest')

# Seasonal average
winter_2020 = ds['LSSTSST'].sel(time=slice('2020-12', '2021-03')).mean(dim='time')
summer_2020 = ds['LSSTSST'].sel(time=slice('2020-06', '2020-09')).mean(dim='time')

# Time series at point
lat, lon = 41.9, 12.5  # Rome city center
point_temp = ds['LSSTSST'].sel(y=lat, x=lon, method='nearest')
point_temp.plot()
plt.title('Temperature Time Series - Rome Center')
plt.ylabel('Temperature (K)')
plt.show()
```

### Spatial Analysis

```python
# Temperature difference land vs sea
jan_lst = ds['LST_LANDSAT'].sel(time='2020-01-15', method='nearest')
jan_sst = ds['SST_LANDSAT'].sel(time='2020-01-15', method='nearest')
land_sea_diff = jan_lst - jan_sst

# Plot hybrid product
import cartopy.crs as ccrs

fig, ax = plt.subplots(figsize=(12, 8), 
                       subplot_kw={'projection': ccrs.PlateCarree()})
jan_2020.plot(ax=ax, cmap='RdYlBu_r', 
              vmin=270, vmax=310,
              cbar_kwargs={'label': 'Temperature (K)'})
ax.coastlines()
ax.gridlines(draw_labels=True)
plt.title('Hybrid LST/SST - January 2020')
plt.show()
```

### Multi-sensor Comparison

```python
# Compare Landsat vs MODIS LST
date = '2020-07-15'
lst_l = ds['LST_LANDSAT'].sel(time=date, method='nearest')
lst_m = ds['LST_MODIS'].sel(time=date, method='nearest')

diff = lst_l - lst_m

fig, axes = plt.subplots(1, 3, figsize=(18, 5))
lst_l.plot(ax=axes[0], cmap='hot', vmin=280, vmax=320)
axes[0].set_title('Landsat LST')
lst_m.plot(ax=axes[1], cmap='hot', vmin=280, vmax=320)
axes[1].set_title('MODIS LST')
diff.plot(ax=axes[2], cmap='RdBu_r', center=0)
axes[2].set_title('Difference (Landsat - MODIS)')
plt.tight_layout()
plt.show()
```

### Land Cover Analysis

```python
# Track land cover changes
lc_2020 = ds['LANDCOVER'].sel(time='2020-06-15', method='nearest')
lc_2024 = ds['LANDCOVER'].sel(time='2024-06-15', method='nearest')

# Classes: 12=Cropland, 13=Urban, 10=Grassland
changes = lc_2024 - lc_2020
print(f"Pixels changed: {(changes != 0).sum().values}")
```

### Water Mask Application

```python
# Separate land and water temperatures
water_mask = ds['WATERMASK'].isel(time=0)  # Same for all times
hybrid = ds['LSSTSST'].sel(time='2020-07-15', method='nearest')

land_temp = hybrid.where(water_mask == 0)  # Land only
water_temp = hybrid.where(water_mask == 1)  # Water only

print(f"Mean land temperature: {land_temp.mean().values:.2f} K")
print(f"Mean sea temperature: {water_temp.mean().values:.2f} K")
```

---

## Output Structure

```
LST_Rome_30m_15day/
├── static/
│   ├── DEM.tif                              # SRTM 30m elevation
│   └── WaterMask.tif                        # Binary land/water mask
│
├── landsat/
│   ├── Landsat_2020-01-01_0000.tif         # 20 bands per file
│   ├── Landsat_2020-01-16_0001.tif
│   └── ... (~150 files)
│
├── modis/
│   ├── MODIS_2020-01-01_0000.tif           # 17 bands per file
│   ├── MODIS_2020-01-16_0001.tif
│   └── ... (~150 files)
│
├── Rome_MultiSensor_LST_SST_2020_2025.nc   # Final NetCDF (~5-10 GB)
└── multisensor_pipeline.log                # Processing log
```

---

## NetCDF Metadata

The final NetCDF file includes comprehensive metadata compliant with CF-1.8 and ACDD-1.3 standards:

### Global Attributes (Selection)
- **Title**: "Multi-Sensor Land and Sea Surface Temperature Dataset - Rome Metropolitan Area"
- **Creator**: Alessandro Sebastianelli
- **Institution**: CMCC Foundation
- **Website**: https://alessandrosebastianelli.github.io/
- **Contact**: alessandro.sebastianelli@cmcc.it
- **Data Sources**: Landsat 8/9, MODIS Terra/Aqua, SRTM, JRC, ESA
- **Processing Description**: Complete 8-step workflow documentation
- **Spatial Coverage**: Rome area (11.6-13.2°E, 41.4-42.4°N)
- **Temporal Coverage**: 2020-2025, 15-day intervals
- **Resolution**: 30 meters
- **Conventions**: CF-1.8, ACDD-1.3
- **References**: Links to all data source documentation

### Variable Attributes (Examples)
```python
LSSTSST:
  long_name: "Hybrid Land and Sea Surface Temperature"
  standard_name: "surface_temperature"
  units: "K"
  description: "Seamless temperature: LST over land, SST over water"
  methodology: "Conditional selection using WATERMASK"
  valid_range: [250.0, 350.0]

LANDSAT_SR_B2:
  long_name: "Landsat Surface Reflectance Band 2 (Blue)"
  wavelength: "0.45-0.51 µm"
  units: "1"
  valid_range: [0.0, 1.0]

WATERMASK:
  long_name: "Water Mask"
  flag_values: [0, 1]
  flag_meanings: "land water"
  source: "JRC Global Surface Water + ESA WorldCover"
```

---

## Validation and Quality Control

### Automated Validation
```bash
# Check file integrity
python validate_tiffs.py

# Remove corrupted files and report
python validate_tiffs.py --delete

# Validate specific directory
python validate_tiffs.py --dir ./LST_Rome_30m_15day/landsat/
```

### Quality Metrics
- Cloud masking reduces atmospheric contamination
- 15-day median compositing minimizes outliers
- Multiple validation cycles ensure data integrity
- Corrupted files automatically removed and re-downloaded

### Known Limitations
- LST/SST accuracy affected by atmospheric conditions
- Viewing geometry impacts (especially MODIS at high angles)
- Emissivity uncertainties (~1-2% typical)
- Cloud masking not perfect in all conditions
- Spatial resolution differences (Landsat 30m native, MODIS resampled)

---

## Performance

### Processing Time (Rome area, 2020-2025)
- Static layers: ~5-10 minutes
- Landsat (~150 composites): ~6-12 hours
- MODIS (~150 composites): ~6-12 hours  
- NetCDF creation: ~10-30 minutes
- **Total**: ~12-24 hours for complete dataset

### Storage Requirements
- Landsat TIFFs: ~10-15 GB
- MODIS TIFFs: ~5-10 GB
- Final NetCDF: ~5-10 GB (compressed)
- **Total**: ~30-40 GB

### Memory Requirements
- Minimum: 8GB RAM
- Recommended: 16GB RAM
- Peak usage: NetCDF creation step

---

## Troubleshooting

### Authentication Issues
```python
# Check GEE credentials
ee.Initialize(project='gee-juno')
# If fails, verify:
# 1. Service account key file exists
# 2. Email address is correct
# 3. Account has Earth Engine permissions
```

### Download Failures
- **Network issues**: Increase delays in config
- **Rate limiting**: Reduce concurrent requests
- **Tile failures**: Lower tile_size in config
- **Timeout errors**: Increase timeout_seconds

### Memory Issues
- Reduce tile_size: `0.4 → 0.2` degrees
- Limit concurrent processing
- Close other applications
- Use smaller date ranges

### Windows Path Issues
✅ **Fixed!** All scripts use short tile names (`img00000_t00.tif`)

---

## Citation

If you use this dataset, please cite:

```
Sebastianelli, A. (2025). Multi-Sensor Land and Sea Surface Temperature Dataset 
for Rome Metropolitan Area. CMCC Foundation. 
Data sources: Landsat 8/9 (USGS), MODIS Terra/Aqua (NASA), SRTM (NASA/JPL).
Processed using Google Earth Engine.
```

### Acknowledgments
- Landsat data courtesy of the U.S. Geological Survey
- MODIS data courtesy of NASA
- SRTM data courtesy of NASA/JPL
- Processed using Google Earth Engine platform
- ESA WorldCover and JRC Global Surface Water datasets

---

## License

Data are provided for scientific research and educational purposes. Users should acknowledge the original data sources (Landsat, MODIS, SRTM) and cite this dataset appropriately.

---

## Contact and Support

**Author**: Alessandro Sebastianelli  
**Email**: alessandro.sebastianelli@cmcc.it  
**Website**: https://alessandrosebastianelli.github.io/  
**Institution**: CMCC Foundation  
**Institution Website**: https://www.cmcc.it/

For questions, issues, or collaborations, please contact the author.

---

## Version History

### v2.0 (December 2025)
- Complete pipeline refactoring
- Multi-sensor integration (Landsat + MODIS)
- Hybrid LST/SST product
- Dynamic land cover and emissivity
- Comprehensive metadata (CF-1.8, ACDD-1.3)
- Service account authentication
- Automated validation cycles
- Single NetCDF output

### v1.0 (Initial)
- Basic Landsat LST download
- Simple processing workflow

---

## Files in This Package

1. **config_refactored.yaml** - Main configuration file
2. **download_static_layers.py** - Downloads DEM and water mask
3. **download_landsat.py** - Downloads Landsat composites
4. **download_modis.py** - Downloads MODIS composites
5. **make_netcdf_final.py** - Creates final NetCDF product
6. **validate_tiffs.py** - Validates and cleans files
7. **run_pipeline.py** - Master orchestrator script
8. **requirements.txt** - Python dependencies
9. **README.md** - This file
10. **FINAL_STRUCTURE.md** - Detailed NetCDF structure documentation

---

**Last Updated**: December 2025  
**Pipeline Version**: 2.0  
**Author**: Alessandro Sebastianelli, CMCC Foundation