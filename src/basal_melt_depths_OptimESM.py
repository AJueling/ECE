"""creating netcdf files with minimum and maximum basal melt depths

variables to adjust on your local system:
fn_do       :  file path and name of deptho.nc file
fn_fwf      :  file path and name of areacello.nc file
output_path :  path where files will be saved

OptimESM output files:
`basal_melt_depth{1/2}.nc` with variable `bmdepth`
where {1} is the minimum depth = 220 m (= `zdraft` in sbcfwf.F90)
and   {2} is the maximum depth = min(660 m, bathymetry) (= `zshelf` in sbcfwf.F90)
"""

import os
import xarray as xr
import datetime


fn_do  = '../data/T255-ORCA1/deptho.nc'
fn_fwf = '/ec/res4/hpcperm/nm6/ece3data/nemo/fwf-update/FWF_LRF_y1850_3315Gtpy.nc'
output_path = '.'


def test_input_files(fn):
    """ ensure existence of input files """
    assert os.path.exists(fn), f'file does not exist:  {fn}'


def basal_melt_depths_OptimESM():
    """ where there is forcing, take min(660 m, Depth) """
    do = xr.open_dataset(fn_do).deptho
    fwf = xr.open_dataset(fn_fwf).sorunoff_f.isel(time_counter=0).drop('time_counter')
    for i in [1,2]:
        bm = do.copy()
        bm.name = 'bmdepth'
        if i==1:    # minimum depth = 220 m
            bm.values = xr.where(fwf>0,220,0).values
            sn = 'minimum_depth_of_basal_melt_injection'
            ln = 'minimum depth at which basal melt is injected'
            cm = 'in sbcfwf.F90 subroutine fwf_bm this is the variable `zdraft`'
        elif i==2:   #maximum depth = 600 m
            bm.values = xr.where(fwf>0,660,0).values
            bm.values = xr.where(bm>do, do, bm).values
            sn = 'maximum_depth_of_basal_melt_injection'
            ln = 'maximum depth at which basal melt is injected'
            cm = 'in sbcfwf.F90 subroutine fwf_bm this is the variable `zshelf`'
        bm = bm.astype('float64')
        now_str = datetime.datetime.now().replace(microsecond=0).isoformat()
        bm.attrs = {
            'standard_name': sn,
            'long_name': ln,
            'comment': cm,
            'units': 'm',
            'cell_measures': 'area: areacello',
            'history':f'{now_str} created with basal_melt_depths.py'
        }
        if 'lat' in bm.coords.keys():
            bm = bm.drop(['lat','lon'])
        bm.to_netcdf(f'{output_path}/basal_melt_depth{i}.nc')
    return


def test_output_files(fn):
    """ ensure correct structure of output files
    the files need to have the data in double precision / float64
    the coordinates should be 'i', 'j', 'latitude', and 'longitude'
        (though NEMO reads 'x', 'y', 'nav_lat', and 'nav_lon'
         coordinates without problems as well)
    the 'i'/'j' coordinates should not be of int64 type (not recognized by NEMO)
    """
    ds = xr.open_dataarray(fn)
    assert ds.data.dtype == 'float64'
    assert ds.i.dtype != 'int64'
    assert ds.j.dtype != 'int64'
    coordinates = list(ds.coords.keys())
    coordinates.sort()
    assert coordinates == ['i', 'j', 'latitude', 'longitude']
    return


if __name__=="__main__":
    test_input_files(fn_do)
    test_input_files(fn_fwf)
    basal_melt_depths_OptimESM()
    test_output_files(f'{output_path}/basal_melt_depth1.nc')
    test_output_files(f'{output_path}/basal_melt_depth2.nc')