import numpy as np 
import xarray as xr 
import glob 
import satpy 
from satpy import Scene
from satpy import find_files_and_readers
import sys
import os 
import zipfile
import math 
import datetime 
import pdb 
import pickle
import shutil
import rioxarray
import matplotlib.pyplot as plt 
import contextlib
import warnings 

##########################
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
        if isinstance(value, np.ndarray):
            da.attrs[attr] = ' '.join(['{:d}'.format(xx) for xx in value])
        if isinstance(value, np.uint16):
            da.attrs[attr] = int(value) 
        if isinstance(value, np.float64):
            da.attrs[attr] = float(value) 
        if isinstance(value, satpy.dataset.dataid.WavelengthRange):
            da.attrs[attr] = '{:.2f} {:.2f} {:.2f}'.format(value.min, value.central,value.max)
            
        if value is None:
            todel.append(attr)
    
    for attr in todel:
        del da.attrs[attr]
    return da
#############################
if __name__ == '__main__':
#############################

    dirin = '/data/paugam/Data/2024_aveiro/{:s}'
    dirout = '/data/paugam/Data/2024_aveiro/FCI-nc-ortho/'
    os.makedirs(dirout, exist_ok=True)
    start_time = datetime.datetime(2024,9,15,0,0,0)
    end_time   = datetime.datetime(2024,9,16,23,50,0)
    resolution = datetime.timedelta(minutes=10)

    ds_arr = []
    current_time = start_time 
    while current_time <= end_time:
        print(current_time, end=' ')
        next_time = current_time + resolution
        
        if os.path.isfile(dirout+'/fci-SILEXdomain-{:s}.nc'.format(current_time.strftime("%Y%j.%H%M"))): 
            print('already done')
            current_time = next_time
            continue
        print('')
        files = find_files_and_readers(base_dir=dirin.format(current_time.strftime("%Y%m%d")), reader='fci_l1c_nc', 
                                       start_time=current_time, end_time=next_time)
        
        # read the file
        scn = Scene(filenames=files)
        #load RGB and IR ands resample
        with open(os.devnull, 'w') as fnull:
            with contextlib.redirect_stdout(fnull), contextlib.redirect_stderr(fnull):
                scn.load(['natural_color','ir_38','ir_105'], upper_right_corner="NE")
        scn_resampled = scn.resample("eurol1", resampler='nearest', radius_of_influence=5000)
            
        #composite_id = ['natural_color']
        #scn_resampled.load(['natural_color'], upper_right_corner="NE")

        #get time
        time = scn_resampled['natural_color'].attrs['start_time']
        time_img_ = scn_resampled['natural_color'].attrs['start_time'].strftime("%Y%j.%H%M") 
        if (time-current_time).total_seconds() != 0 : 
            print('pb in time')
            sys.exit()

        #crop to zone
        scn_cropped = scn_resampled.crop(xy_bbox=(-1.2E6, -6.34E6, 1.1E6, -4.78E6))
        
        
        daR = adjust_da_attr(scn_cropped['natural_color'].sel(bands='R').drop_vars('bands').rename('R'))
        #daR.attrs.clear()
        daG = adjust_da_attr(scn_cropped['natural_color'].sel(bands='G').drop_vars('bands').rename('G'))
        #daG.attrs.clear()
        daB = adjust_da_attr(scn_cropped['natural_color'].sel(bands='B').drop_vars('bands').rename('B'))
        #daB.attrs.clear()
        daIR038 = adjust_da_attr(scn_cropped['ir_38'].rename('ir_38'))
        #daIR039.attrs.clear()
        daIR105 = adjust_da_attr(scn_cropped['ir_105'].rename('ir_105'))
        #daIR108.attrs.clear()

        ds = xr.merge([
                      daR,daG,daB,daIR038,daIR105
                      ], 
                      compat='override'
                      )

        crs = ds.attrs["area"].to_cartopy_crs()
        ds = ds.rio.write_crs(crs,inplace=True)
        ds=ds.drop_vars('crs')

        ds = ds.expand_dims({"time": [time]})  # Add time dimension
       
        
        '''
        listofattrtosaveExt = ['area','_satpy_id','prerequisites']
        with open(dirout+'/seviri-extAttr-S-EU-{:s}.pkl'.format(time_img), 'wb') as f:
            pickle.dump([ds.attrs[xx] for xx in listofattrtosaveExt], f)
        for xx in listofattrtosaveExt:
            del ds.attrs[xx]
        for var in list(ds.data_vars): 
            for xx in listofattrtosaveExt:
                if xx in ds[var].attrs.keys():
                    del ds[var].attrs[xx]
        '''
        ds.attrs['crs']=ds.rio.crs.to_string()

        attr2del = ['area','_satpy_id', 'prerequisites', 'optional_prerequisites','ancillary_variables']
        for var in ds.data_vars:
            for attrname in attr2del:
                if attrname in ds[var].attrs:
                    del ds[var].attrs[attrname]
        for attrname in attr2del[:-1]:
            del ds.attrs[attrname]

        #ds_arr.append(ds)
        current_time = next_time
   
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            ds.to_netcdf(dirout+'/fci-SILEXdomain-{:s}.nc'.format(time_img_))
        
        del scn 

    #ds = xr.concat(ds_arr, dim="time")
    #sys.exit()
    
       
