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
import pdb 

# Function to recursively convert HDF5 groups and datasets to a dictionary
def h5_to_dict(h5group):
    result = {}
    variables = []
    for key, item in h5group.items():
        if isinstance(item, h5py.Group):
            result[key], variable = h5_to_dict(item)
        elif isinstance(item, h5py.Dataset):
            if 'scale_factor' in item.attrs.keys():
                scale_factor = item.attrs['scale_factor']
            else: 
                scale_factor = 1
            if 'add_offset' in item.attrs.keys():
                add_offset = item.attrs['add_offset']
            else: 
                add_offset = 0
            if 'radiance_scale_factor' in item.attrs.keys():
                radiance_scale_factor = item.attrs['radiance_scale_factor']
            else: 
                radiance_scale_factor = 1
            if 'radiance_add_offset' in item.attrs.keys():
                radiance_add_offset = item.attrs['radiance_add_offset']
            else:
                radiance_add_offset = 0

            result[key] = ( ( (item[()] * scale_factor) + add_offset ) * radiance_scale_factor + radiance_add_offset)
            #result[key] = item[()] 
            
            variable = key
        variables.append(variable)
    return result, variables

if __name__ == '__main__':

    dirin02 = '/data/paugam/Data/ElPontDeVilamora/VIIRS/L1B-02/'
    dirout  = '/data/paugam/Data/ElPontDeVilamora/VIIRS/nc-ortho/'
    
    viirs02s = glob.glob(dirin02+'*.nc') 

    for viirs02 in viirs02s:
        foutname = dirout + '.'.join(os.path.basename(viirs02).split('.')[:3])+'.nc'

        if not(os.path.isfile(foutname)):
            print('---------------------')
            print(os.path.basename(viirs02))
                
            #asociated L1B03
            try:
                if 'VNP' in viirs02:
                    tmp_ = 'VNP'
                elif 'VJ1' in viirs02:
                    tmp_ = 'VJ1'
                else:
                    print('###########################')
                    print('sat not defined here # stop')
                    sys.exit()
                viirs03 = glob.glob('.'.join(viirs02.replace('L1B-02','L1B-03').replace('{:s}02'.format(tmp_),'{:s}03'.format(tmp_)).split('.')[:3])+'*.nc')[0]
            except:
                print('L1B-03 file not found for: ',  viirs02)
                sys.exit() 
            #load geoloc
            ############
            hdf5_file = viirs03
            with h5py.File(hdf5_file, 'r') as h5file:
                # Convert HDF5 file content to a dictionary
                data_dict, variables = h5_to_dict(h5file)

            nc_dict = {}
            for key in data_dict.keys():
                if key == 'geolocation_data':
                    for key2 in data_dict[key]:
                        if key2 in ['latitude', 'longitude', 'height']: 
                            print(key2)
                            nc_dict[key2] = (['y', 'x'], data_dict[key][key2])
            
            nj,ni = nc_dict['latitude'][1].shape

            dsgeo = xr.Dataset( nc_dict,
                coords={
                    'x': np.arange(ni),
                    'y': np.arange(nj),
                    })
            
            #load radiance/reflectance
            ############
            hdf5_file = viirs02
            
            if 'MOD' in os.path.basename(viirs02): 
                ibands = [5,4,3]
                bandletter = 'M'
                method_interpolation = 'linear'
                resolution = 0.01
            elif 'IMG' in os.path.basename(viirs02):
                ibands = [1,2,3,4,5]
                bandletter = 'I'
                method_interpolation = 'nearest'
                resolution = 0.005
           
            bands = []; bands_flag = []
            for  iband in ibands: 
                bands.append('{:s}{:02d}'.format(bandletter,iband))
                bands_flag.append('{:s}{:02d}_quality_flags'.format(bandletter,iband))

            with h5py.File(hdf5_file, 'r') as h5file:
                # Convert HDF5 file content to a dictionary
                data_dict, variables = h5_to_dict(h5file)

            nc_dict = {}
            for key in data_dict.keys():
                if key == 'observation_data':
                    for key2 in data_dict[key]:
                        if key2 in bands+bands_flag:
                            print(key2)
                            nc_dict[key2] = (['y', 'x'], data_dict[key][key2])
            
            ds = xr.Dataset( nc_dict,
                coords={
                    'latitude' : (('y','x'),dsgeo['latitude'].data),
                    'longitude': (('y','x'),dsgeo['longitude'].data),
                    })

            ds.rio.write_crs('EPSG:4326', inplace=True)
         
            idx = np.where(ds[bands_flag[3]].data != 256)
            # Sample points (you can replace this with your actual points)
            points = [Point(x, y) for x, y in zip(dsgeo['longitude'].data[idx][::100], dsgeo['latitude'].data[idx][::100]) ]
            gdf = gpd.GeoDataFrame(geometry=points)

            # Define alpha parameter based on your dataset
            alpha = 0.3  # Adjust as needed for a tighter fit

            # Generate the concave hull
            concave_hull = alphashape.alphashape(gdf.geometry, alpha)

            # Convert to GeoDataFrame and plot
            concave_hull_gdf = gpd.GeoDataFrame(geometry=[concave_hull])

            # Plot the result
            '''
            fig, ax = plt.subplots()
            gdf.plot(ax=ax, color='blue', marker='o', markersize=5, label='Points')
            concave_hull_gdf.plot(ax=ax, color='lightgreen', alpha=0.5, edgecolor='green', label='Concave Hull')
            plt.legend()
            plt.show()
            '''

            latmin,latmax = dsgeo['latitude'].min(), dsgeo['latitude'].max()+resolution
            lonmin,lonmax = dsgeo['longitude'].min(), dsgeo['longitude'].max()+resolution
            new_lat = np.arange(latmin,latmax,resolution)
            new_lon = np.arange(lonmin,lonmax,resolution)


            # Calculate the number of rows and columns based on the bounding box and resolution
            height =  new_lat.shape[0]
            width  =  new_lon.shape[0]

            # Define transform for the output mask
            transform = rasterio.transform.from_origin(lonmin, new_lat.max(), resolution, resolution)

            # Create the mask array (initialize as all False)
            mask = geometry_mask(
                geometries=concave_hull_gdf.geometry,
                out_shape=(height, width),
                transform=transform,
                invert=True  # Set to True to make features in the gdf True in the mask
            )
            mask = mask[::-1,:].T
           
            ds_list= []
            for band,band_flag in zip(bands,bands_flag):
                if (band not in ds.variables.keys()) & (band not in ['I04','I05', 'I04_quality_flags', 'I05_quality_flags']): 
                    continue
                print('interp ', band)
                idx = np.where(ds[band_flag].data != 256)
                grid_y, grid_x = np.meshgrid(new_lat,new_lon)
                projected_data = interpolate.griddata((dsgeo['longitude'].data[idx], dsgeo['latitude'].data[idx]), 
                                                       ds[band].data[idx], 
                                                      (grid_x, grid_y), method=method_interpolation)

                da = xr.DataArray( projected_data.T,
                    coords={
                        'lat' : (('lat'),new_lat),
                        'lon':  (('lon'),new_lon),
                        })
                da.name = band
                ds_list.append(da.where(mask.T))

            ds = xr.merge(ds_list)

            ds.to_netcdf(foutname)


    

