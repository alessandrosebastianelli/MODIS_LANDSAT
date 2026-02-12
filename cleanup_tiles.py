#!/usr/bin/env python3
"""
Cleanup temporary tiles directories
"""

from pathlib import Path
import shutil
import yaml

# Load config
with open('config_refactored.yaml', 'r') as f:
    config = yaml.safe_load(f)

base_dir = Path(config['output']['base_dir'])
landsat_dir = base_dir / config['output']['subdirs']['landsat']
modis_dir = base_dir / config['output']['subdirs']['modis']
static_dir = base_dir / config['output']['subdirs']['static']

# Find and remove all tiles_temp directories
for directory in [landsat_dir, modis_dir, static_dir]:
    if not directory.exists():
        continue
    
    tiles_dirs = list(directory.glob('tiles_temp*'))
    
    for tiles_dir in tiles_dirs:
        if tiles_dir.is_dir():
            try:
                shutil.rmtree(tiles_dir)
                print(f"Removed: {tiles_dir}")
            except Exception as e:
                print(f"Could not remove {tiles_dir}: {e}")

print("Cleanup complete")