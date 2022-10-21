""" calculate the JJA 500hPa temperature ('T500') for the 8 regions in Greenland from Lenearts et al (2015) """

import os
import sys
import dask
import intake
import fsspec
import numpy as np
import xarray as xr

from tqdm.autonotebook import tqdm

gm255 = xr.open_dataarray('../data/Lenaerts2015/Lenaerts2015_mask_255.nc')
gm511 = xr.open_dataarray('../data/Lenaerts2015/Lenaerts2015_mask_511.nc')

if __name__ == '__main__':
    store = sys.argv[1]
    print(store)
    
    ext = ''
    if store=='ceda-zarr':
        col_fn = "https://raw.githubusercontent.com/cedadev/cmip6-object-store/master/catalogs/ceda-zarr-cmip6.json"
    elif store=='ceda-nc':
        col_fn = '../../EC-Earth3-data/ceda_nc_highresmip.json'
    elif store=='jasmin-nc':
        col_fn = '../../EC-Earth3-data/jasmin_nc_highresmip.json'
    elif store=='ecmwf-cca-scratch':
        col_fn = '../../EC-Earth3-data/ecmwf_cca_scratch.json'
    ext = '-ext'
    col = intake.open_esm_datastore(col_fn)
    cat = col.search(institution_id="EC-Earth-Consortium", table_id='Amon', variable_id='ta')
    df = cat.df
    
    for i in tqdm(range(len(df))):  # loop over all datasets
#         if i<10:  continue
        print()
        model, member, exp = df.loc[i]['source_id'], df.loc[i]['member_id'], df.loc[i]['experiment_id']
        print(i, model, member, exp+ext)
        fn = f'../results/Lenaerts2015_regions/T500_{model}_{member}_{exp}{ext}.nc'
        print(fn)

        # if os.path.exists(fn):  continue
        if exp not in ['hist-1950','highres-future']:  continue
            
        if model=='EC-Earth3P-HR':
            mask = gm511.astype(int)
            chunks = dict(time=12,plev=19,lat=512,lon=1024)
        else:
            continue
            mask = gm255
            chunks = dict(time=12,plev=19,lat=256,lon=512)
        # print('\nmask\n', mask)
        # print(mask.sum().values)

        if store=='ceda-zarr':
            try:
                fsmap = fsspec.get_mapper(df['zarr_path'][i])    
                ds = xr.open_dataset(fsmap, consolidated=True, engine='zarr', chunks=chunks).ta.sel(plev=50000.)
            except:
                ds = xr.open_mfdataset(df['nc_path'][i], engine='h5netcdf', chunks=chunks).ta.sel(plev=50000.)
        else:
            ds = xr.open_mfdataset(df['nc_path'][i], engine='h5netcdf', chunks=chunks).ta.sel(plev=50000.)
        # print('\nafter loading mfdataset\n', ds)
            
        # create JJA average
        jun = ds.isel(time=slice(5,None,12))
        jul = ds.isel(time=slice(6,None,12))
        aug = ds.isel(time=slice(7,None,12))
        jun = jun.assign_coords(time=jun.time.dt.year)
        jul = jul.assign_coords(time=jul.time.dt.year)
        aug = aug.assign_coords(time=aug.time.dt.year)
        ds = (jun+jul+aug).squeeze()/3
        # print('\nafter time averaging\n', ds)
        # # print(ds.sum(['lat','lon']).values)
        
        if model=='EC-Earth3P-HR':
            ds = ds.assign_coords({'lat':gm511.lat.values, 'lon':gm511.lon.values})
        else:
            ds = ds.assign_coords({'lat':gm255.lat.values, 'lon':gm255.lon.values})
            # print('\n', ds)
        
        # spatial integration
        assert np.all(mask.lon==ds.lon)
        assert np.all(mask.lat==ds.lat)
        out = xr.DataArray(dims=['time','region'], coords=dict(time=ds.time,region=np.arange(1,9)))
        # print('\n', out.shape, '\n', out,'\n')
        for rn in tqdm(np.arange(1,9)):
            # print('\nregions for loop', i, rn, '\n')

            # # print('\nwhere statement\n',mask.where(mask==rn).mean(['lat','lon']).values)
            # # print('\nwhere statement\n',mask.where(mask==rn).sum(['lat','lon']).values)
            # print('\nwhere statement\n',ds.isel(time=0).mean(['lat','lon']).values)
            truth= xr.where(mask==rn, 1, 0)
            # print(truth.lat)
            # print(truth.lon)
            # print(ds.lat)
            # print(ds.lon)
            assert np.all(truth.lat==ds.lat)
            assert np.all(truth.lon==ds.lon)
            # print('\nwhere statement\n',ds.isel(time=0).where(mask==rn).mean(['lat','lon']).values)
            # print('\nwhere statement\n',xr.where(mask==rn, ds.isel(time=0), np.nan).mean(['lat','lon']).values)
            out[:,rn-1] = ds.where(mask==rn).weighted(np.cos(np.deg2rad(ds.lat))).mean(['lat','lon'])
        # print('\noutput after integraton\n', out)
        
        out.to_netcdf(fn)
        print(f'wrote output to {fn}')
