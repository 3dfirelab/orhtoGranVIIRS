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

#homebrewed
from fci_plot import *

##########################
if __name__ == '__main__':
##########################
    
    parser = argparse.ArgumentParser(description="to plot fci data at fire location")
    parser.add_argument('firename', type=str, help="fire name defined in the code")
    args = parser.parse_args()

    if args.firename == 'aveiro':
        ### les Landes
        #latf, lonf = 40.689, -8.629
        latf, lonf = 40.689, -8.56
        width = 50.e3 #m
        firename = 'aveiro'
        firename_short = 'averio'

    
    dirin  = '/data/paugam/Data/2024_aveiro/FCI-nc-ortho/'
    
    #define bbox
    bbox = define_bbox_from_center(latf, lonf, width)
    epsgCode = get_utm_epsg(latf, lonf)
    
    dirout  = '/data/paugam/Data/2024_aveiro/FCI-{:s}-nc/'.format(firename_short)
    os.makedirs(dirout, exist_ok=True)
    os.makedirs(dirout+'png/', exist_ok=True)
    
    #select times
    times = []
    bands = []
    sats = []
    fciNCs = glob.glob(dirin+'*.nc')
    for fcinc in fciNCs:
        times.append( fcitime2datetime(*(os.path.basename(fcinc).split('-')[-1].split('.')[:2])) )
        sats.append('MTG')
        bands.append('RGBIR38105')
   
    fcidf = pd.DataFrame({
        'time': times, 
        'band': bands,
        'sat': sats,
        'filename': fciNCs
        })
    fcidf_s = fcidf.sort_values(by=['time', 'band'], ascending=[True, True])

    #loop through times and plot fire img
    ds_list = []
    for idx, fciInfo  in fcidf_s.iterrows(): 
       
        fig, axs = plt.subplots(1, 3, figsize=(16,6), subplot_kw={'projection':  ccrs.epsg(epsgCode)})
        fig.suptitle('{:s}: {:s} | {:s}'.format(firename, fciInfo['sat'], fciInfo['time'].strftime("%Y-%m-%d %H:%M") ) )
        flag_plot = False
        
        #extract subdomain
        #########
       
        ds_subset_utm, boundslatlon = get_UTMsubset_fire(fciInfo,bbox,epsgCode) 
        in_bounds = (boundslatlon['lat_min'] <= latf <= boundslatlon['lat_max']) and \
                    (boundslatlon['lon_min'] <= lonf <= boundslatlon['lon_max'])
        if not(in_bounds): continue
      
        print(firename_short,' extract ' , fciInfo['sat'], fciInfo['time'])
        ds_list.append(ds_subset_utm)
       
        #plot png
        #########

        ax = axs[0]
        pco=ax.pcolormesh(ds_subset_utm.x, ds_subset_utm.y, ds_subset_utm.ir_38.isel(time=0).data, cmap='inferno')
        gridlines = ax.gridlines(draw_labels=True)
        gridlines.top_labels = True
        gridlines.right_labels = False
        gridlines.bottom_labels = False
        gridlines.left_labels = True
        ax.set_aspect(1)
        cbaxes = fig.add_axes([ax.get_position().x0, ax.get_position().y0-0.03, ax.get_position().x1-ax.get_position().x0, 0.05])
        cbar = fig.colorbar(pco ,cax = cbaxes,orientation='horizontal')
        cbar.set_label('MWIR BT (K)')
            
        ax = axs[1]
        pco=ax.pcolormesh(ds_subset_utm.x, ds_subset_utm.y, ds_subset_utm.ir_105.isel(time=0).data, cmap='jet')
        gridlines = ax.gridlines(draw_labels=True)
        gridlines.top_labels = True
        gridlines.right_labels = False
        gridlines.bottom_labels = False
        gridlines.left_labels = True
        ax.set_aspect(1)
        cbaxes = fig.add_axes([ax.get_position().x0, ax.get_position().y0-0.03, ax.get_position().x1-ax.get_position().x0, 0.05])
        cbar = fig.colorbar(pco ,cax = cbaxes,orientation='horizontal')
        cbar.set_label('LWIR BT (K)')
        flag_plot = True
            
        Rmin = ds_subset_utm['R'].min()
        Rmax = ds_subset_utm['R'].max()
        Gmin = ds_subset_utm['G'].min()
        Gmax = ds_subset_utm['G'].max()
        Bmin = ds_subset_utm['B'].min()
        Bmax = ds_subset_utm['B'].max()


        # Stack the bands together
        rgb = xr.concat([(ds_subset_utm['R'].isel(time=0)-Rmin)/Rmax,
                         (ds_subset_utm['G'].isel(time=0)-Gmin)/Gmax,
                         (ds_subset_utm['B'].isel(time=0)-Bmin)/Bmax],
                         dim='bands')

        # Ensure bands are in correct order (Red, Green, Blue)
        # rgb should now have shape (n, m, 3) where 3 represents R, G, B
        rgb = rgb.transpose('y', 'x', 'bands')

        ax = axs[2] 
        pco = ax.pcolormesh(ds_subset_utm.x, ds_subset_utm.y, rgb)
        gridlines = ax.gridlines(draw_labels=True)
        gridlines.top_labels = True
        gridlines.right_labels = True
        gridlines.bottom_labels = False
        gridlines.left_labels = False
        ax.set_aspect(1)
        flag_plot = True
    
        plt.subplots_adjust(bottom=0.2)

        fig.savefig('{:s}/{:s}-{:s}-{:s}.png'.format(dirout+'png/',firename_short,'MTG-fci',fciInfo['time'].strftime("%Y-%m-%d-%H%M")))
        
        plt.close(fig)
    
    ds = xr.concat(ds_list, dim="time")
    ds.rio.write_crs(epsgCode)
    ds.to_netcdf(dirout+'MTG-'+firename_short+'.nc')

    # Extract the CRS information
    crs_wkt = ds.rio.crs.to_wkt()
    # Save the CRS to a .prj file
    prj_file_path = dirout+'MTG-'+firename_short+'.prj'
    with open(prj_file_path, 'w') as prj_file:
        prj_file.write(crs_wkt)
