import numpy as np 
import xarray as xr
import glob 
import os
import sys 
import h5py
import matplotlib.pyplot as plt
import rioxarray
from scipy import interpolate
import alphashape
import geopandas as gpd
from shapely.geometry import Point
import rasterio
from rasterio.features import geometry_mask
import datetime 
import pandas as pd 
import math 
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pdb 
import argparse
import pyproj 

##########################
def load_crs_from_prj(prj_file_path):
    # Check if the .prj file exists
    if not os.path.exists(prj_file_path):
        raise FileNotFoundError(f"The file {prj_file_path} does not exist.")
    
    # Read the WKT from the .prj file
    with open(prj_file_path, 'r') as prj_file:
        wkt = prj_file.read().strip()

    # Create a CRS object from the WKT
    crs = pyproj.CRS(wkt)

    return crs

##########################
if __name__ == '__main__':
##########################

    parser = argparse.ArgumentParser(description="to plot fci data at fire location")
    parser.add_argument('firename', type=str, help="fire name defined in the code")
    args = parser.parse_args()
    
    if args.firename == 'aveiro':
        ### les Landes
        latf, lonf = 40.689, -8.629
        width = 200.e3 #m
        firename = 'aveiro'
        firename_short = 'averio'
    
    dirin  = '/data/paugam/Data/2024_aveiro/FCI-{:s}-nc/'.format(firename_short)

    ds = xr.open_dataset(dirin+'MTG-'+firename_short+'.nc')
    crs = load_crs_from_prj(dirin+'MTG-'+firename_short+'.prj')
    threshold = 320
  
    
    ir_38_max = ds.ir_38.max(dim='time')
    bmap = np.zeros_like(ir_38_max)
    reference_time = pd.Timestamp("2024-09-15 00:00:00")  

    for j,i in zip(*np.where(ir_38_max > threshold)):

        df = ds.ir_38[:,j,i].to_dataframe()

        # Create a condition for BT > threshold
        df['above_threshold'] = df['ir_38'] > threshold

        # Check for consecutive rows where the condition is True
        df['consecutive'] = df['above_threshold'] & df['above_threshold'].shift(1, fill_value=False)

        # Find the first time this happens
        result = df[df['consecutive']].iloc[0] if not df[df['consecutive']].empty else None

        if result is not None:
            bmap[j,i] = (result.name - reference_time).total_seconds()
        #else: 
        #    bmap[j,i] = np.nan 

    ds_bmap = ir_38_max.copy().rename('bmap') 
    ds_bmap.data = bmap  
    ds_bmap.attrs['units']='seconds since {:s}'.format(str(reference_time))
    ds_bmap.rio.write_crs(crs)
    
    ds_bmap.to_netcdf(dirin+'bmap1.nc')

    #ds_bmap.plot()
    #plt.show()
