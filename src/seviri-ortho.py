import numpy as np 
import xarray as xr 
import glob 
from satpy import Scene
import sys
import os 
import zipfile
import math 
import datetime 
import pdb 
import pickle
import shutil

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

def adjust_da_attr(da):
    #expand dict attr 
    # Expand the 'metadata' dictionary into multiple attributes
    for attrname in ['orbital_parameters','time_parameters']:
        metadata = da.attrs.get(attrname, {})  # Get the metadata dictionary

        # Add each key-value pair as a separate attribute
        for key, value in metadata.items():
            da.attrs['{:s}_{:s}'.format(attrname,key)] = value

        # Remove the original 'metadata' attribute
        del da.attrs[attrname]

    #replace datetime attr by string
    todel=[]
    for attr, value in da.attrs.items():
        if isinstance(value, datetime.datetime):
            da.attrs[attr] = value.strftime("%Y-%m-%d %H:%M:%S.%f") 
        if isinstance(value, np.bool):
            da.attrs[attr] = value.astype(np.uint8)
        if isinstance(value, bool):
            da.attrs[attr] = int(value)
        if value is None:
            todel.append(attr)
    
    for attr in todel:
        del da.attrs[attr]
    return da

#############################
if __name__ == '__main__':
#############################
    
    firename = 'pdV'
    latf, lonf = 41.694, 1.894
    width = 75.e3 #m

    dirin = '/data/paugam/Data/ElPontDeVilamora/Severi/NAT/'
    dirout = '/data/paugam/Data/ElPontDeVilamora/Severi/nc-ortho/'
    extraction_path = "/tmp/paugam/SEVIRI/extracted_data/"
    if os.path.isdir: shutil.rmtree(extraction_path)
    zipfiles =sorted( glob.glob(dirin+'*.zip') )
    bbox = define_bbox_from_center(latf, lonf, width) 

    for zipfile_ in zipfiles:
        os.makedirs(extraction_path, exist_ok=True)
        with zipfile.ZipFile(zipfile_, 'r') as zip_ref:
            zip_ref.extractall(extraction_path)

        nat_files = glob.glob(extraction_path+'*.nat')

        for file in nat_files:
        
            # define reader
            reader = "seviri_l1b_native"
            # read the file
            scn = Scene(filenames = {reader:[file]})
           
            #load RGB
            composite_id = ['natural_color']
            scn.load(composite_id, upper_right_corner="NE")

            #load IR
            scn.load(['IR_039','IR_108'], upper_right_corner="NE")

            #crop to zone
            scn_cropped = scn.crop(xy_bbox=(-15E5, 34E5, 15E5, 48E5))

            daR = adjust_da_attr(scn_cropped['natural_color'].sel(bands='R').drop_vars('bands').rename('R'))
            #daR.attrs.clear()
            daG = adjust_da_attr(scn_cropped['natural_color'].sel(bands='G').drop_vars('bands').rename('G'))
            #daG.attrs.clear()
            daB = adjust_da_attr(scn_cropped['natural_color'].sel(bands='B').drop_vars('bands').rename('B'))
            #daB.attrs.clear()
            daIR039 = adjust_da_attr(scn_cropped['IR_039'].rename('IR_039'))
            #daIR039.attrs.clear()
            daIR108 = adjust_da_attr(scn_cropped['IR_108'].rename('IR_108'))
            #daIR108.attrs.clear()

            ds = xr.merge([
                          daR,daG,daB,daIR039,daIR108
                          ], 
                          compat='override'
                          )
            time_img = datetime.datetime.strptime(ds.attrs['time_parameters_nominal_start_time'], "%Y-%m-%d %H:%M:%S.%f").strftime("%Y%j.%H%M") 
            print(time_img)
            crs = ds.attrs["area"].to_cartopy_crs()
            ds = ds.rio.write_crs(crs,inplace=True)
            ds=ds.drop_vars('crs')
           
            listofattrtosaveExt = ['area','_satpy_id','prerequisites']
            with open(dirout+'/seviri-extAttr-S-EU-{:s}.pkl'.format(time_img), 'wb') as f:
                pickle.dump([ds.attrs[xx] for xx in listofattrtosaveExt], f)
            for xx in listofattrtosaveExt:
                del ds.attrs[xx]
            for var in list(ds.data_vars): 
                for xx in listofattrtosaveExt:
                    if xx in ds[var].attrs.keys():
                        del ds[var].attrs[xx]
            
            ds.attrs['crs']=ds.rio.crs.to_string()

            ds.to_netcdf(dirout+'/seviri-S-EU-{:s}.nc'.format(time_img))
        shutil.rmtree(extraction_path)
            
