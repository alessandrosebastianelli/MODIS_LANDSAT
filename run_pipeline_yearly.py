#!/usr/bin/env python3
# ============================================
# Yearly Pipeline - Process One Year at a Time
# ============================================

import subprocess
import sys
from pathlib import Path
import time
from datetime import datetime
import yaml
import logging
import argparse

# ============================================
# LOGGING SETUP
# ============================================

def setup_logging():
    """Setup logging for yearly pipeline"""
    log_file = 'logs/pipeline_yearly.log'
    
    Path('logs').mkdir(parents=True, exist_ok=True)
    
    logger = logging.getLogger('Pipeline_Yearly')
    logger.setLevel(logging.INFO)
    logger.handlers = []
    logger.propagate = False
    
    # File handler
    fh = logging.FileHandler(log_file, mode='a')
    fh.setLevel(logging.INFO)
    fh.setFormatter(logging.Formatter(
        '%(asctime)s | %(name)-15s | %(levelname)-8s | %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    ))
    logger.addHandler(fh)
    
    # Console handler with colors
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    
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
    logger.addHandler(ch)
    
    return logger

logger = setup_logging()

# ============================================
# CONFIGURATION
# ============================================

with open('config_refactored.yaml', 'r') as f:
    config = yaml.safe_load(f)

# Years to process
START_YEAR = 2007
END_YEAR = 2021
YEARS = list(range(START_YEAR, END_YEAR + 1))

# Output directories
BASE_DIR = Path(config['output']['base_dir'])
LANDSAT_DIR = BASE_DIR / config['output']['subdirs']['landsat']
MODIS_DIR = BASE_DIR / config['output']['subdirs']['modis']
STATIC_DIR = BASE_DIR / config['output']['subdirs']['static']
NETCDF_DIR = BASE_DIR / 'netcdf_yearly'

# Scripts
DOWNLOAD_STATIC = 'download_static_layers.py'
DOWNLOAD_LANDSAT = 'download_landsat.py'
DOWNLOAD_MODIS = 'download_modis.py'
VALIDATE_TIFFS = 'validate_tiffs.py'
MAKE_NETCDF = 'make_netcdf_yearly_v2.py'

# ============================================
# HELPER FUNCTIONS
# ============================================

def run_command(cmd, description):
    """Run a command and return success status"""
    logger.info(f"Running: {description}")
    logger.info(f"Command: {' '.join(cmd)}")
    
    start_time = time.time()
    
    logger.info(f"Executing: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(
            cmd,
            check=True,
            capture_output=False,
            text=True
        )
        
        elapsed = time.time() - start_time
        logger.info(f"✓ Success: {description} ({elapsed:.1f}s)")
        return True
        
    except subprocess.CalledProcessError as e:
        elapsed = time.time() - start_time
        logger.error(f"✗ Failed: {description} ({elapsed:.1f}s)")
        logger.error(f"Error code: {e.returncode}")
        return False
    except FileNotFoundError:
        logger.error(f"✗ Script not found: {cmd[1]}")
        return False

def count_files_for_year(directory, year):
    """Count TIFF files for a specific year"""
    if not directory.exists():
        return 0
    pattern = f"*{year}-*"
    return len(list(directory.glob(pattern)))

def check_netcdf_exists(year):
    """Check if NetCDF for year already exists"""
    nc_file = NETCDF_DIR / f"Rome_LST_SST_{year}_30m_15day.nc"
    return nc_file.exists()

# ============================================
# MAIN PIPELINE
# ============================================

def main():
    parser = argparse.ArgumentParser(description='Process satellite data year by year')
    parser.add_argument('--year', type=int, help='Process specific year only')
    parser.add_argument('--start-year', type=int, default=START_YEAR, help='Start year')
    parser.add_argument('--end-year', type=int, default=END_YEAR, help='End year')
    parser.add_argument('--skip-static', action='store_true', help='Skip static layers download')
    parser.add_argument('--skip-validation', action='store_true', help='Skip validation step')
    parser.add_argument('--skip-netcdf', action='store_true', help='Skip NetCDF creation')
    
    args = parser.parse_args()
    
    # Determine years to process
    if args.year:
        years_to_process = [args.year]
    else:
        years_to_process = list(range(args.start_year, args.end_year + 1))
    
    logger.info("="*70)
    logger.info("YEARLY SATELLITE DATA PIPELINE")
    logger.info("="*70)
    logger.info(f"Years to process: {years_to_process}")
    logger.info(f"Output directory: {BASE_DIR}")
    logger.info("="*70)
    
    start_time_total = datetime.now()
    
    # Statistics
    stats = {
        'years_processed': 0,
        'years_failed': 0,
        'years_skipped': 0,
        'total_landsat': 0,
        'total_modis': 0
    }
    
    # ============================================
    # STEP 0: STATIC LAYERS (once)
    # ============================================
    
    if not args.skip_static:
        logger.info("="*70)
        logger.info("STEP 0: STATIC LAYERS (run once)")
        logger.info("="*70)
        
        dem_file = STATIC_DIR / 'DEM_30m.tif'
        water_file = STATIC_DIR / 'WaterMask_30m.tif'
        
        if dem_file.exists() and water_file.exists():
            logger.info("Static layers already exist, skipping")
        else:
            success = run_command(
                [sys.executable, DOWNLOAD_STATIC],
                "Download static layers"
            )
            
            if not success:
                logger.warning("Static layers download failed, continuing anyway")
    
    # ============================================
    # PROCESS EACH YEAR
    # ============================================
    
    for year in years_to_process:
        logger.info("="*70)
        logger.info(f"PROCESSING YEAR: {year}")
        logger.info("="*70)
        
        year_start = time.time()
        
        # Check if NetCDF already exists
        if check_netcdf_exists(year) and not args.skip_netcdf:
            logger.info(f"NetCDF for {year} already exists, skipping year")
            stats['years_skipped'] += 1
            continue
        
        # ============================================
        # STEP 1: LANDSAT
        # ============================================
        
        logger.info(f"[{year}] Step 1: Download Landsat data")
        
        landsat_before = count_files_for_year(LANDSAT_DIR, year)
        
        success = run_command(
            [sys.executable, DOWNLOAD_LANDSAT, '--year', str(year)],
            f"Landsat {year}"
        )
        
        if not success:
            logger.error(f"Landsat download failed for {year}")
            stats['years_failed'] += 1
            continue
        
        landsat_after = count_files_for_year(LANDSAT_DIR, year)
        landsat_new = landsat_after - landsat_before
        stats['total_landsat'] += landsat_new
        logger.info(f"Landsat files for {year}: {landsat_after} ({landsat_new} new)")
        
        # ============================================
        # STEP 2: MODIS
        # ============================================
        
        logger.info(f"[{year}] Step 2: Download MODIS data")
        
        modis_before = count_files_for_year(MODIS_DIR, year)
        
        success = run_command(
            [sys.executable, DOWNLOAD_MODIS, '--year', str(year)],
            f"MODIS {year}"
        )
        
        if not success:
            logger.error(f"MODIS download failed for {year}")
            stats['years_failed'] += 1
            continue
        
        modis_after = count_files_for_year(MODIS_DIR, year)
        modis_new = modis_after - modis_before
        stats['total_modis'] += modis_new
        logger.info(f"MODIS files for {year}: {modis_after} ({modis_new} new)")
        
        # ============================================
        # STEP 3: VALIDATE (optional)
        # ============================================
        
        if not args.skip_validation:
            logger.info(f"[{year}] Step 3: Validate files")
            
            success = run_command(
                [sys.executable, VALIDATE_TIFFS, '--year', str(year), '--delete'],
                f"Validation {year}"
            )
            
            if not success:
                logger.warning(f"Validation failed for {year}, continuing")
        
        # ============================================
        # STEP 4: CREATE NETCDF
        # ============================================
        
        if not args.skip_netcdf:
            logger.info(f"[{year}] Step 4: Create NetCDF")
            
            success = run_command(
                [sys.executable, MAKE_NETCDF, '--year', str(year)],
                f"NetCDF {year}"
            )
            
            if not success:
                logger.error(f"NetCDF creation failed for {year}")
                stats['years_failed'] += 1
                continue
        
        # ============================================
        # YEAR COMPLETE
        # ============================================
        
        year_elapsed = time.time() - year_start
        stats['years_processed'] += 1
        
        logger.info("="*70)
        logger.info(f"YEAR {year} COMPLETE ({year_elapsed/3600:.1f} hours)")
        logger.info("="*70)
    
    # ============================================
    # FINAL SUMMARY
    # ============================================
    
    end_time_total = datetime.now()
    elapsed_total = end_time_total - start_time_total
    
    logger.info("="*70)
    logger.info("PIPELINE COMPLETE")
    logger.info("="*70)
    logger.info(f"Started: {start_time_total.strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"Finished: {end_time_total.strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"Total time: {elapsed_total}")
    logger.info("")
    logger.info(f"Years processed: {stats['years_processed']}")
    logger.info(f"Years failed: {stats['years_failed']}")
    logger.info(f"Years skipped: {stats['years_skipped']}")
    logger.info(f"Total Landsat files: {stats['total_landsat']}")
    logger.info(f"Total MODIS files: {stats['total_modis']}")
    logger.info("="*70)
    
    # List created NetCDF files
    nc_files = sorted(NETCDF_DIR.glob('*.nc'))
    if nc_files:
        logger.info(f"\nCreated {len(nc_files)} NetCDF files:")
        for nc in nc_files:
            size_mb = nc.stat().st_size / (1024**2)
            logger.info(f"  {nc.name}: {size_mb:.2f} MB")
    
    return 0 if stats['years_failed'] == 0 else 1

# ============================================
# ENTRY POINT
# ============================================

if __name__ == "__main__":
    try:
        exit_code = main()
        sys.exit(exit_code)
    except KeyboardInterrupt:
        logger.error("\nPipeline interrupted by user")
        sys.exit(1)
    except Exception as e:
        logger.error(f"\nUnexpected error: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)