import os
import numpy as np
import xarray as xr

def list_files(dir):
    r = []
    for root, dirs, files in os.walk(dir):
        for name in files:
            r.append(os.path.join(root, name))
    return r

if __name__ == '__main__':
    exp05 = [x for x in list_files('../data/v7_CMIP5_pub/') if x[-8:]=='exp05.nc' and 'mm' in x]
    new_ds = []
    for i, fn in enumerate(exp05):
        ds = xr.open_dataset(fn).sle
        model = fn.split('/')[5]
        if i<2 or i==12:
            model = model+'_'
        print(i, model)
        ds = ds.assign_coords(model=model)
        ds = ds.assign_coords(time=np.arange(2014,2100))
        new_ds.append(ds)
    ismip = xr.concat(new_ds, dim='model')
    ismip.attrs['provenance'] = 'collated ISMIP6 Greenland data in ISMIP.py'
    ismip.to_netcdf('../data/v7_CMIP5_pub/ISMIP6_Greenland.nc')