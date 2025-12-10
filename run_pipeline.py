# ============================================
# Master Pipeline Runner
# Orchestrates: Static -> Landsat -> MODIS -> Validate -> NetCDF
# ============================================

import subprocess
import sys
from pathlib import Path
import time
from datetime import datetime
import yaml
import logging

# ============================================
# LOGGING SETUP
# ============================================

def setup_logging(config):
    log_level = getattr(logging, config['logging']['level'])
    log_file = config['logging']['log_file']
    
    logger = logging.getLogger('Pipeline')
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

# ============================================
# PIPELINE STEPS
# ============================================

STEPS = [
    {
        'name': 'Static Layers',
        'script': 'download_static_layers.py',
        'description': 'Download DEM, Land Cover, Water Mask, Emissivity',
        'required': True
    },
    {
        'name': 'Landsat Composites',
        'script': 'download_landsat.py',
        'description': 'Download Landsat 15-day composites with all bands',
        'required': True,
        'validate': True
    },
    {
        'name': 'MODIS Composites',
        'script': 'download_modis.py',
        'description': 'Download MODIS 15-day composites with all bands',
        'required': True,
        'validate': True
    },
    {
        'name': 'Create NetCDF',
        'script': 'make_netcdf.py',
        'description': 'Combine all data into NetCDF with LST/SST hybrid',
        'required': True
    }
]

# ============================================
# HELPER FUNCTIONS
# ============================================

def run_command(cmd, description):
    """Run a command and handle errors"""
    logger.info(f"\n{'='*70}")
    logger.info(f"Running: {description}")
    logger.info(f"{'='*70}")
    logger.debug(f"Command: {' '.join(cmd)}")
    
    start_time = time.time()
    
    try:
        result = subprocess.run(
            cmd,
            check=True,
            capture_output=False,
            text=True
        )
        
        elapsed = time.time() - start_time
        logger.info(f"✓ {description} completed ({elapsed:.1f}s)")
        return True
        
    except subprocess.CalledProcessError as e:
        elapsed = time.time() - start_time
        logger.error(f"✗ {description} failed ({elapsed:.1f}s)")
        logger.error(f"Error code: {e.returncode}")
        return False
    except FileNotFoundError:
        logger.error(f"✗ Script not found: {cmd[1]}")
        return False

def count_files(directory):
    """Count TIFF files in directory"""
    if not directory.exists():
        return 0
    return len(list(directory.glob("*.tif")))

def get_dir_size(directory):
    """Get total size of TIFF files in MB"""
    if not directory.exists():
        return 0.0
    total = sum(f.stat().st_size for f in directory.glob("*.tif"))
    return total / (1024**2)

# ============================================
# MAIN PIPELINE
# ============================================

def main():
    logger.info("\n" + "="*70)
    logger.info("MULTI-SENSOR LST/SST PIPELINE")
    logger.info("="*70)
    logger.info(f"Configuration: config_refactored.yaml")
    logger.info(f"Steps: {len(STEPS)}")
    logger.info(f"Validation cycles: {config['pipeline']['num_validation_cycles']}")
    
    start_time_total = datetime.now()
    logger.info(f"Started: {start_time_total.strftime('%Y-%m-%d %H:%M:%S')}")
    
    base_dir = Path(config['output']['base_dir'])
    landsat_dir = base_dir / config['output']['subdirs']['landsat']
    modis_dir = base_dir / config['output']['subdirs']['modis']
    static_dir = base_dir / config['output']['subdirs']['static']
    
    # Track results
    results = {
        'completed': [],
        'failed': [],
        'skipped': []
    }
    
    # ============================================
    # EXECUTE PIPELINE STEPS
    # ============================================
    
    for i, step in enumerate(STEPS, 1):
        logger.info(f"\n{'='*70}")
        logger.info(f"STEP {i}/{len(STEPS)}: {step['name']}")
        logger.info(f"{'='*70}")
        logger.info(f"Description: {step['description']}")
        
        # Run main script
        success = run_command(
            [sys.executable, step['script']],
            step['name']
        )
        
        if success:
            results['completed'].append(step['name'])
            
            # If step requires validation, run validation cycles
            if step.get('validate', False):
                logger.info(f"\nRunning validation cycles for {step['name']}...")
                
                for cycle in range(config['pipeline']['num_validation_cycles']):
                    logger.info(f"\nValidation cycle {cycle+1}/{config['pipeline']['num_validation_cycles']}")
                    
                    # Validate
                    validate_success = run_command(
                        [sys.executable, 'validate_tiffs.py', '--delete'],
                        f"Validation ({step['name']})"
                    )
                    
                    if not validate_success:
                        logger.warning("Validation failed, but continuing...")
                    
                    # Re-download if needed (only if not last cycle)
                    if cycle < config['pipeline']['num_validation_cycles'] - 1:
                        logger.info(f"Re-running {step['name']} to recover corrupted files...")
                        run_command(
                            [sys.executable, step['script']],
                            f"{step['name']} (re-run {cycle+1})"
                        )
                        
                        time.sleep(config['pipeline']['pause_between_cycles'])
        else:
            results['failed'].append(step['name'])
            if step.get('required', False):
                logger.error(f"Required step failed: {step['name']}")
                logger.error("Stopping pipeline")
                return 1
            else:
                logger.warning(f"Optional step failed, continuing...")
        
        # Pause between steps
        if i < len(STEPS):
            time.sleep(3)
    
    # ============================================
    # FINAL SUMMARY
    # ============================================
    
    logger.info("\n" + "="*70)
    logger.info("PIPELINE COMPLETE")
    logger.info("="*70)
    
    end_time_total = datetime.now()
    elapsed_total = end_time_total - start_time_total
    
    logger.info(f"Started: {start_time_total.strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"Finished: {end_time_total.strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"Total time: {elapsed_total}")
    logger.info("")
    logger.info(f"Completed steps: {len(results['completed'])}/{len(STEPS)}")
    logger.info(f"Failed steps: {len(results['failed'])}/{len(STEPS)}")
    
    if results['completed']:
        logger.info("\nCompleted:")
        for step in results['completed']:
            logger.info(f"  ✓ {step}")
    
    if results['failed']:
        logger.info("\nFailed:")
        for step in results['failed']:
            logger.info(f"  ✗ {step}")
    
    # Check outputs
    logger.info("\n" + "="*70)
    logger.info("OUTPUT SUMMARY")
    logger.info("="*70)
    
    if static_dir.exists():
        count = count_files(static_dir)
        size = get_dir_size(static_dir)
        logger.info(f"Static layers: {count} files, {size:.2f} MB")
    
    if landsat_dir.exists():
        count = count_files(landsat_dir)
        size = get_dir_size(landsat_dir)
        logger.info(f"Landsat: {count} files, {size:.2f} MB")
    
    if modis_dir.exists():
        count = count_files(modis_dir)
        size = get_dir_size(modis_dir)
        logger.info(f"MODIS: {count} files, {size:.2f} MB")
    
    # Check NetCDF
    nc_file = Path(config['output']['netcdf_filename'])
    if nc_file.exists():
        nc_size = nc_file.stat().st_size / (1024**2)
        logger.info(f"\nNetCDF: {nc_file.name}, {nc_size:.2f} MB")
        logger.info("  ✓ Combined dataset ready!")
    else:
        logger.warning("\nNetCDF not created")
    
    logger.info("="*70)
    logger.info("✓ Pipeline execution complete!")
    logger.info("="*70)
    
    return 0 if len(results['failed']) == 0 else 1

# ============================================
# ENTRY POINT
# ============================================

if __name__ == "__main__":
    try:
        exit_code = main()
        sys.exit(exit_code)
    except KeyboardInterrupt:
        logger.warning("\nPipeline interrupted by user")
        sys.exit(1)
    except Exception as e:
        logger.critical(f"\nUnexpected error: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
