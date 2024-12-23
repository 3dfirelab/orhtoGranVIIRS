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

##########################
def fcitime2datetime(yearjd,time):
    tmp_ = yearjd[:4]+'-'+yearjd[4:]+'_'+time[:2]+':'+time[2:]
    return datetime.datetime.strptime(tmp_, '%Y-%j_%H:%M')


##########################
def define_bbox_from_center(lat, lon, width_m):
    # Conversion constants
    meters_per_deg_lat = 111320  # Approximate for latitude globally
    meters_per_deg_lon = lambda lat: 111320 * math.cos(math.radians(lat))  # Depends on latitude
    
    # Calculate half-width in degrees
    half_width_deg_lat = width_m / (2 * meters_per_deg_lat)
    half_width_deg_lon = width_m / (2 * meters_per_deg_lon(lat))
    
    # Define bounding box
    bbox = {
        "minx": lon - half_width_deg_lon,
        "miny": lat - half_width_deg_lat,
        "maxx": lon + half_width_deg_lon,
        "maxy": lat + half_width_deg_lat
    }
    return bbox

##########################
def get_utm_epsg(lat, lon):
    # Calculate UTM zone
    zone_number = int((lon + 180) / 6) + 1
    
    # Determine if it's northern or southern hemisphere
    if lat >= 0:
        epsg_code = 32600 + zone_number  # Northern hemisphere
    else:
        epsg_code = 32700 + zone_number  # Southern hemisphere
        
    return epsg_code

##########################
def get_UTMsubset_fire(seviriInfo,bbox,epsgCode):

    ds = xr.open_dataset(seviriInfo.filename)
    ds = ds.rio.write_crs(ds.rio.crs, inplace=True)
    #ds = ds.drop_vars('acq_time')
    ds = ds.rio.reproject('epsg:4326')

    ds_subset = ds.rio.clip_box(
        minx=bbox["minx"],
        miny=bbox["miny"],
        maxx=bbox["maxx"],
        maxy=bbox["maxy"]
    )
    
    boundslatlon = {
            'lat_min' : ds_subset.y.min().values,
            'lat_max' : ds_subset.y.max().values,
            'lon_min' : ds_subset.x.min().values,
            'lon_max' : ds_subset.x.max().values
            }

    return ds_subset.rio.reproject('epsg:{:d}'.format(epsgCode)), boundslatlon


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

    
    dirin  = '/data/paugam/Data/2024_aveiro/FCI-nc-ortho/'
    
    #define bbox
    bbox = define_bbox_from_center(latf, lonf, width)
    epsgCode = get_utm_epsg(latf, lonf)
    
    dirout  = '/data/paugam/Data/2024_aveiro/FCI-{:s}-png/'.format(firename_short)
    os.makedirs(dirout, exist_ok=True)
    
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
    for idx, fciInfo  in fcidf_s.iterrows(): 
       
        fig, axs = plt.subplots(1, 3, figsize=(16,6), subplot_kw={'projection':  ccrs.epsg(epsgCode)})
        fig.suptitle('{:s}: {:s} | {:s}'.format(firename, fciInfo['sat'], fciInfo['time'].strftime("%Y-%m-%d %H:%M") ) )
        flag_plot = False
       
        ds_subset_utm, boundslatlon = get_UTMsubset_fire(fciInfo,bbox,epsgCode) 
        in_bounds = (boundslatlon['lat_min'] <= latf <= boundslatlon['lat_max']) and \
                    (boundslatlon['lon_min'] <= lonf <= boundslatlon['lon_max'])
        if not(in_bounds): continue
      
        print(firename_short, fciInfo['sat'], fciInfo['time'])
        
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

        fig.savefig('{:s}/{:s}-{:s}-{:s}.png'.format(dirout,firename_short,'MTG-fci',fciInfo['time'].strftime("%Y-%m-%d-%H%M")))
        
        plt.close(fig)
    





