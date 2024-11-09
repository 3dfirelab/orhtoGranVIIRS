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

##########################
def viirstime2datetime(yearjd,time):
    tmp_ = yearjd[1:5]+'-'+yearjd[5:]+'_'+time[:2]+':'+time[2:]
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
def get_UTMsubset_fire(viirsInfo):

    ds = xr.open_dataset(viirsInfo.filename)
    ds.rio.write_crs('EPSG:4326', inplace=True)
    ds = ds.rename({'lat':'y'})
    ds = ds.rename({'lon':'x'})

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

    latf, lonf = 41.694, 1.894
    width = 75.e3 #m
    dirin  = '/data/paugam/Data/ElPontDeVilamora/VIIRS/nc-ortho/'
    dirout  = '/data/paugam/Data/ElPontDeVilamora/VIIRS/png/'
    firename = 'El Pont de Vilomara'
    #define bbox
    bbox = define_bbox_from_center(latf, lonf, width)
    epsgCode = get_utm_epsg(latf, lonf)

    #select times
    times = []
    bands = []
    sats = []
    viirsNCs = glob.glob(dirin+'*.nc')
    for viirsnc in viirsNCs:
        times.append( viirstime2datetime(*(os.path.basename(viirsnc).split('.')[1:3])) )
        if 'IMG' in viirsnc:
            bands.append( 'ir')
        elif 'MOD' in viirsnc:
            bands.append( 'rgb')
        else: 
            sys.exit()
        if 'VNP' in viirsnc:
            sats.append( 'Suomi-NPP')
        elif 'VJ1' in viirsnc:
            sats.append( 'NOAA-20' )

    viirsdf = pd.DataFrame({
        'time': times, 
        'band': bands,
        'sat': sats,
        'filename': viirsNCs
        })
    viirsdf_s = viirsdf.sort_values(by=['time', 'band'], ascending=[True, True])

    #loop through times and plot fire img
    for time, group in viirsdf_s.groupby('time'): 
       
        fig, axs = plt.subplots(1, 3, figsize=(16,6), subplot_kw={'projection':  ccrs.epsg(epsgCode)})
        fig.suptitle('{:s}: {:s} | {:s}'.format(firename, group['sat'].iloc[0], time.strftime("%Y-%m-%d %H:%M") ) )
        flag_plot = False
       
        if not(group.band.str.contains('rgb').any()):
            new_row = {'band': 'black', 'time': group.time.iloc[0], 'sat':'na', 'filename':group.filename.iloc[0]}
            # Add the new row
            group.loc[len(group)] = new_row

        for ii, viirsInfo in group.iterrows():
            ds_subset_utm, boundslatlon = get_UTMsubset_fire(viirsInfo) 
            in_bounds = (boundslatlon['lat_min'] <= latf <= boundslatlon['lat_max']) and \
                        (boundslatlon['lon_min'] <= lonf <= boundslatlon['lon_max'])
            if not(in_bounds): continue
          
            print(viirsInfo.filename)
            if viirsInfo.band == 'ir':
                ax = axs[0]
                pco=ax.pcolormesh(ds_subset_utm.x, ds_subset_utm.y, ds_subset_utm.I04.data, cmap='inferno')
                gridlines = ax.gridlines(draw_labels=True)
                gridlines.top_labels = True
                gridlines.right_labels = False
                gridlines.bottom_labels = False
                gridlines.left_labels = True
                ax.set_aspect(1)
                cbaxes = fig.add_axes([ax.get_position().x0, ax.get_position().y0-0.03, ax.get_position().x1-ax.get_position().x0, 0.05])
                cbar = fig.colorbar(pco ,cax = cbaxes,orientation='horizontal')
                cbar.set_label('MWIR Radiance')
                
                ax = axs[1]
                pco=ax.pcolormesh(ds_subset_utm.x, ds_subset_utm.y, ds_subset_utm.I05.data, cmap='jet')
                gridlines = ax.gridlines(draw_labels=True)
                gridlines.top_labels = True
                gridlines.right_labels = False
                gridlines.bottom_labels = False
                gridlines.left_labels = True
                ax.set_aspect(1)
                cbaxes = fig.add_axes([ax.get_position().x0, ax.get_position().y0-0.03, ax.get_position().x1-ax.get_position().x0, 0.05])
                cbar = fig.colorbar(pco ,cax = cbaxes,orientation='horizontal')
                cbar.set_label('LWIR Radiance')
                flag_plot = True
                
            if viirsInfo.band == 'rgb':
                
                # Stack the bands together
                rgb = xr.concat([(ds_subset_utm['M05']-ds_subset_utm['M05'].min())/ds_subset_utm['M05'].max(),
                                 (ds_subset_utm['M04']-ds_subset_utm['M04'].min())/ds_subset_utm['M04'].max(),
                                 (ds_subset_utm['M03']-ds_subset_utm['M03'].min())/ds_subset_utm['M03'].max()],
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
            
            if viirsInfo.band == 'black':
                
                black = np.zeros(ds_subset_utm.I04.shape)

                ax = axs[2] 
                pco = ax.pcolormesh(ds_subset_utm.x, ds_subset_utm.y, black,cmap='Greys_r',vmin=0,vmax=1)
                gridlines = ax.gridlines(draw_labels=True)
                gridlines.top_labels = True
                gridlines.right_labels = True
                gridlines.bottom_labels = False
                gridlines.left_labels = False
                ax.set_aspect(1)
                flag_plot = True
               
        if flag_plot:
            plt.subplots_adjust(bottom=0.2)
            fig.savefig('{:s}/{:s}-{:s}.png'.format(dirout,firename,time.strftime("%Y-%m-%d-%H%M")))
        
        plt.close(fig)
    





