# ============================================
# Validate GeoTIFF Files
# Works on landsat/, modis/, static/ subdirectories
# ============================================

from pathlib import Path
import rasterio
import argparse
import yaml
import logging

# ============================================
# LOGGING SETUP
# ============================================

def setup_logging(config):
    log_level = getattr(logging, config['logging']['level'])
    log_file = config['logging']['log_file']
    
    logger = logging.getLogger('Validate')
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
# ARGPARSE
# ============================================

parser = argparse.ArgumentParser(description="Validate GeoTIFF files")
parser.add_argument("--delete", action="store_true",
                    help="Delete corrupted files")
parser.add_argument("--dir", type=str, default=None,
                    help="Specific directory to validate (overrides config)")
args = parser.parse_args()

# ============================================
# LOAD CONFIG
# ============================================

with open('config_refactored.yaml', 'r') as f:
    config = yaml.safe_load(f)

logger = setup_logging(config)

# Determine directories to validate
base_dir = Path(config['output']['base_dir'])

if args.dir:
    dirs_to_validate = [Path(args.dir)]
else:
    # Validate all subdirectories
    dirs_to_validate = [
        base_dir / config['output']['subdirs']['landsat'],
        base_dir / config['output']['subdirs']['modis'],
        base_dir / config['output']['subdirs']['static']
    ]

delete_corrupted = args.delete

logger.info("="*70)
logger.info("GeoTIFF Validation Tool")
logger.info("="*70)
logger.info(f"Mode: {'DELETE corrupted' if delete_corrupted else 'LIST only'}")
logger.info(f"Directories: {len(dirs_to_validate)}")

# ============================================
# VALIDATE FUNCTION
# ============================================

def validate_directory(directory):
    """Validate all TIFFs in a directory"""
    
    if not directory.exists():
        logger.warning(f"Directory does not exist: {directory}")
        return {'valid': 0, 'corrupted': 0, 'size_valid': 0, 'size_corrupted': 0}
    
    tif_files = sorted(directory.glob("*.tif"))
    
    if len(tif_files) == 0:
        logger.info(f"No TIFF files in {directory}")
        return {'valid': 0, 'corrupted': 0, 'size_valid': 0, 'size_corrupted': 0}
    
    logger.info(f"\nValidating {directory.name}/ ({len(tif_files)} files)...")
    
    valid_files = []
    corrupted_files = []
    total_size_valid = 0
    total_size_corrupted = 0
    
    for i, tif in enumerate(tif_files, 1):
        is_valid = True
        error_msg = None
        file_size = tif.stat().st_size
        
        # Check file size
        if file_size < 1000:
            is_valid = False
            error_msg = f"File too small ({file_size} bytes)"
        else:
            # Try to open and read
            try:
                with rasterio.open(tif) as src:
                    if src.count < 1:
                        is_valid = False
                        error_msg = "No bands found"
                    else:
                        # Try reading first band
                        try:
                            _ = src.read(1)
                        except Exception as read_error:
                            is_valid = False
                            error_msg = f"Read error: {str(read_error)[:60]}"
                        
                        if is_valid:
                            valid_files.append(tif)
                            total_size_valid += file_size
            except Exception as e:
                is_valid = False
                error_msg = str(e)[:80]
        
        # Report corrupted files
        if not is_valid:
            corrupted_files.append(tif)
            total_size_corrupted += file_size
            logger.error(f"[{i}/{len(tif_files)}] CORRUPTED: {tif.name}")
            logger.error(f"           {error_msg}")
        else:
            if i % 10 == 0:
                logger.debug(f"[{i}/{len(tif_files)}] Validated {len(valid_files)} files...")
    
    # Summary for this directory
    logger.info(f"\n{directory.name}/ Summary:")
    logger.info(f"  Valid: {len(valid_files)}")
    logger.info(f"  Corrupted: {len(corrupted_files)}")
    
    # Delete if requested
    if corrupted_files and delete_corrupted:
        logger.warning(f"  Deleting {len(corrupted_files)} corrupted files...")
        deleted = 0
        for cf in corrupted_files:
            try:
                cf.unlink()
                deleted += 1
            except Exception as e:
                logger.error(f"  Could not delete {cf.name}: {e}")
        logger.info(f"  ✓ Deleted {deleted}/{len(corrupted_files)} files")
    
    return {
        'valid': len(valid_files),
        'corrupted': len(corrupted_files),
        'size_valid': total_size_valid,
        'size_corrupted': total_size_corrupted,
        'corrupted_list': corrupted_files if not delete_corrupted else []
    }

# ============================================
# VALIDATE ALL DIRECTORIES
# ============================================

total_stats = {
    'valid': 0,
    'corrupted': 0,
    'size_valid': 0,
    'size_corrupted': 0,
    'corrupted_list': []
}

for directory in dirs_to_validate:
    stats = validate_directory(directory)
    total_stats['valid'] += stats['valid']
    total_stats['corrupted'] += stats['corrupted']
    total_stats['size_valid'] += stats['size_valid']
    total_stats['size_corrupted'] += stats['size_corrupted']
    total_stats['corrupted_list'].extend(stats.get('corrupted_list', []))

# ============================================
# FINAL SUMMARY
# ============================================

logger.info("\n" + "="*70)
logger.info("OVERALL VALIDATION SUMMARY")
logger.info("="*70)
logger.info(f"Total files validated: {total_stats['valid'] + total_stats['corrupted']}")
logger.info(f"Valid files: {total_stats['valid']}")
logger.info(f"Corrupted files: {total_stats['corrupted']}")
logger.info(f"Valid size: {total_stats['size_valid']/(1024**2):.2f} MB")
logger.info(f"Corrupted size: {total_stats['size_corrupted']/(1024**2):.2f} MB")

if total_stats['corrupted'] > 0:
    if delete_corrupted:
        logger.info("\n✓ Corrupted files have been deleted")
        logger.info("  Re-run download scripts to fetch them again")
    else:
        logger.info(f"\nTo delete corrupted files, run:")
        logger.info(f"  python validate_tiffs.py --delete")

logger.info("="*70)
