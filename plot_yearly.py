#!/usr/bin/env python3
"""
Plot yearly NetCDF composites - unified script
"""

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import argparse

parser = argparse.ArgumentParser(description='Plot yearly composite data')
parser.add_argument('--year', type=int, required=True, help='Year to plot')
parser.add_argument('--mode', type=str, default='map', 
                    choices=['map', 'overview', 'timeseries'],
                    help='Plot mode: map (single variable), overview (4-panel), timeseries')
parser.add_argument('--variable', type=str, default='LANDSAT_LST_K', 
                    help='Variable to plot (e.g., LANDSAT_LST_K, LANDSAT_NDVI)')
parser.add_argument('--timestep', type=int, default=0, 
                    help='Time index to plot (for map/overview modes)')
args = parser.parse_args()

# Load data
nc_file = Path(f'LST_Rome_30m_15day/netcdf_yearly/Rome_LST_SST_{args.year}_30m_15day.nc')

if not nc_file.exists():
    print(f"Error: {nc_file} not found")
    exit(1)

print(f"Loading {nc_file}...")
ds = xr.open_dataset(nc_file)

# ==============================================
# MODE 1: Single variable map
# ==============================================
if args.mode == 'map':
    var_name = args.variable
    
    if var_name not in ds.data_vars:
        print(f"Variable '{var_name}' not found. Available:")
        for v in sorted(ds.data_vars)[:30]:
            print(f"  - {v}")
        exit(1)
    
    data = ds[var_name].isel(time=args.timestep)
    time_val = ds.time.isel(time=args.timestep).values
    
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Handle temperature conversion
    if 'LST' in var_name or 'SST' in var_name:
        data_plot = data - 273.15
        cmap = 'RdYlBu_r'
        label = 'Temperature (°C)'
        vmin, vmax = data_plot.quantile([0.02, 0.98]).values
    elif 'NDVI' in var_name or 'EVI' in var_name or 'SAVI' in var_name:
        data_plot = data
        cmap = 'RdYlGn'
        label = 'Vegetation Index'
        vmin, vmax = -0.2, 0.8
    elif 'NDBI' in var_name or 'UI' in var_name:
        data_plot = data
        cmap = 'YlOrRd'
        label = 'Urban Index'
        vmin, vmax = data_plot.quantile([0.02, 0.98]).values
    elif 'DEM' in var_name:
        data_plot = data
        cmap = 'terrain'
        label = 'Elevation (m)'
        vmin, vmax = None, None
    else:
        data_plot = data
        cmap = 'viridis'
        label = var_name
        vmin, vmax = data_plot.quantile([0.02, 0.98]).values
    
    im = data_plot.plot(ax=ax, cmap=cmap, vmin=vmin, vmax=vmax, 
                         cbar_kwargs={'label': label})
    
    ax.set_title(f'{var_name} - {args.year} - {str(time_val)[:10]}')
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    
    plt.tight_layout()
    out_file = f'plot_{args.year}_{var_name.replace("/", "_")}_t{args.timestep:03d}.pdf'
    plt.savefig(out_file, dpi=150, bbox_inches='tight', format='pdf')
    print(f"Saved: {out_file}")
    
    print(f"\nStatistics:")
    print(f"  Mean: {float(data_plot.mean()):.2f}")
    print(f"  Std:  {float(data_plot.std()):.2f}")
    print(f"  Min:  {float(data_plot.min()):.2f}")
    print(f"  Max:  {float(data_plot.max()):.2f}")

# ==============================================
# MODE 2: Overview (4-panel)
# ==============================================
elif args.mode == 'overview':
    time_val = ds.time.isel(time=args.timestep).values
    
    # Standard variables to plot
    var_config = [
        ('LANDSAT_LST_K', 'LST (°C)', 'RdYlBu_r', True),
        ('LANDSAT_NDVI', 'NDVI', 'RdYlGn', False),
        ('LANDSAT_NDBI', 'NDBI', 'YlOrRd', False),
        ('DEM', 'DEM (m)', 'terrain', False)
    ]
    
    vars_to_plot = [(v, t, c, temp) for v, t, c, temp in var_config if v in ds.data_vars]
    
    if len(vars_to_plot) == 0:
        print("No standard variables found. Available:")
        for v in sorted(ds.data_vars)[:20]:
            print(f"  - {v}")
        exit(1)
    
    n_plots = len(vars_to_plot)
    ncols = 2
    nrows = (n_plots + 1) // 2
    
    fig, axes = plt.subplots(nrows, ncols, figsize=(14, 6*nrows))
    if nrows == 1:
        axes = axes.reshape(1, -1)
    axes = axes.flatten()
    
    for idx, (var, title, cmap, is_temp) in enumerate(vars_to_plot):
        data = ds[var].isel(time=args.timestep)
        
        if is_temp:
            data_plot = data - 273.15
            vmin, vmax = data_plot.quantile([0.02, 0.98]).values
        else:
            data_plot = data
            if 'DEM' in var:
                vmin, vmax = None, None
            elif 'NDVI' in var:
                vmin, vmax = -0.2, 0.8
            else:
                vmin, vmax = data_plot.quantile([0.02, 0.98]).values
        
        data_plot.plot(ax=axes[idx], cmap=cmap, vmin=vmin, vmax=vmax,
                       add_colorbar=True, cbar_kwargs={'label': title})
        
        axes[idx].set_title(title)
        axes[idx].set_xlabel('Longitude')
        axes[idx].set_ylabel('Latitude')
    
    for idx in range(n_plots, len(axes)):
        axes[idx].axis('off')
    
    fig.suptitle(f'Rome Area - {args.year} - {str(time_val)[:10]}', 
                 fontsize=14, y=0.995)
    plt.tight_layout()
    
    out_file = f'overview_{args.year}_t{args.timestep:03d}.pdf'
    plt.savefig(out_file, dpi=150, bbox_inches='tight', format='pdf')
    print(f"Saved: {out_file}")

# ==============================================
# MODE 3: Time series
# ==============================================
elif args.mode == 'timeseries':
    var_name = args.variable
    
    if var_name not in ds.data_vars:
        print(f"Variable not found: {var_name}")
        exit(1)
    
    data = ds[var_name]
    mean_ts = data.mean(dim=['x', 'y'])
    
    if 'LST' in var_name or 'SST' in var_name:
        mean_ts = mean_ts - 273.15
        ylabel = 'Temperature (°C)'
    else:
        ylabel = var_name
    
    fig, ax = plt.subplots(figsize=(12, 5))
    mean_ts.plot(ax=ax, marker='o', linewidth=2)
    ax.set_title(f'{var_name} - {args.year} - Spatial Mean')
    ax.set_xlabel('Date')
    ax.set_ylabel(ylabel)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    out_file = f'timeseries_{args.year}_{var_name.replace("/", "_")}.pdf'
    plt.savefig(out_file, dpi=150, bbox_inches='tight', format='pdf')
    print(f"Saved: {out_file}")

plt.close('all')