import os
import sys
import glob
# activate_file = os.path.join('/home/users/ajuling/nb-venvs/venv-cmip6-zarr/bin/activate_this.py')
# exec(open(activate_file).read(), dict(__file__=activate_file))

import dask
import numpy as np
import pandas as pd
import xarray as xr

from paths import path_data
from iterate import IterateECE, subselect_sysargv
from tqdm.auto import tqdm


dz = xr.open_dataarray(f'{path_data}/T511-ORCA025/dz.nc')
zub = xr.open_dataarray(f'{path_data}/T511-ORCA025/zub.nc')
zlb = xr.open_dataarray(f'{path_data}/T511-ORCA025/zlb.nc')

# area_255 = xr.open_dataset(f'{path_data}/T255-ORCA1/areas.nc')['O1t0.srf'].rename({'x_3':'i','y_3':'j'})
# area_511 = xr.open_dataset(f'{path_data}/T511-ORCA025/areas.nc')['Ot25.srf'].rename({'x_3':'i','y_3':'j'})
area_255 = xr.open_dataset(f'{path_data}/T255-ORCA1/areacello.nc')['areacello']
area_511 = xr.open_dataset(f'{path_data}/T511-ORCA025/areacello.nc')['areacello']


# Levermann LARMIP regions (their Tbl. 1)
# [lower lat, upper lat, l. lon, u. lon, central depth]
Lregs = dict(rEAIS = [-76,-65,350,175,369],
             rWedd = [-85,-72,295,350,420],
             rAmun = [-76,-70,213,285,305],
             rRoss = [-85,-76,150,215,312],
             rPen1 = [-70,-65,294,300,420], # Erwin has 310 as eastern boundary
             rPen2 = [-75,-70,285,295,420],
)

freg = {'rEAIS': [.50, .50],
        'rRoss': [.68, .32],
        'rAmun': [.29, .71],
        'rWedd': [.62, .38],
        'rPen' : [.27, .73]}


def create_zw(D, dz, verbose=False):
    """ creates (partial) depth array for cells 100 m around depth D, i.e. +/- 50 m """
    zw = xr.zeros_like(dz)                                        # z-weights
    for i in range(75):
        ub, lb = zub.isel(lev=i).values, zlb.isel(lev=i).values # upper and lower bounds
        assert lb>ub
        if (ub<D-50 and lb<D-50) or (ub>D+50 and lb>D+50):    # cell completely above or below
            if verbose:  print(f'{i:2d}, {ub:6.1f}-{lb:6.1f} outside of {D-50,D+50}')
            continue
        elif ub<D-50 and lb>D-50:                               # upper partial cell
            zw[i] = lb-D+50
            if verbose:  print(f'{i:2d}, {ub:6.1f}-{lb:6.1f} partially in {D-50,D+50}')
        elif ub>D-50 and lb>D-50 and ub<D+50 and lb<D+50:     # full cell
            zw[i] = dz[i]
            if verbose:  print(f'{i:2d}, {ub:6.1f}-{lb:6.1f} fully in {D-50,D+50}')
        elif ub<D+50 and lb>D+50:                               # lower partial cell
            zw[i] = D+50-ub
            if verbose:  print(f'{i:2d}, {ub:6.1f}-{lb:6.1f} partially in {D-50,D+50}')
    if verbose:  print(zw.values, zw.sum().values)
    assert(zw.sum().values==100.)
    return zw


def create_Levermann_weights(region, zw, area, thetao):
    """ creates 3D field with volume in Levermann region/depth """
    if area.shape==(292, 362):      # standard resolution
        fn = f'{path_data}/T255-ORCA1/LARMIP_weights_EC-Earth3P_{region}.nc'
    elif area.shape==(1050, 1442):  # high resolution
        fn = f'{path_data}/T511-ORCA025/LARMIP_weights_EC-Earth3P-HR_{region}.nc'
    
    if os.path.exists(fn):
        weights = xr.open_dataarray(fn)
    else:
        region_coords=Lregs[region]
        assert len(region_coords)==5
        llat, ulat = region_coords[0],region_coords[1]
        llon, ulon = region_coords[2],region_coords[3]
        for lat in [llat,ulat]:   assert lat<=90 and lat>=-90
        for lon in [llon,ulon]:   assert lon<=360 and lon>0

        a, b = area, thetao
        if llon>ulon:
            A = a.where(b.latitude>llat).where(b.latitude<ulat)
            A = A.where(b.longitude>llon).fillna(0) + A.where(b.longitude<ulon).fillna(0)
        else:
            A = a.where(b.latitude>llat).where(b.latitude<ulat).where(b.longitude>llon).where(b.longitude<ulon).fillna(0)
        weights = (xr.ones_like(thetao).where(thetao>-5)*A*zw).fillna(0)
        weights.to_netcdf(fn)
    return weights


def calc_Levermann_timeseries(area, thetao):
    """ calculates temperature time series of Levermann regions """
    assert len(area.i)==len(thetao.i) and len(area.j)==len(thetao.j)
    avgs = np.zeros((len(Lregs.keys()),len(thetao.time)))
    for i, reg in enumerate(Lregs.keys()):
        print(reg)
        D = Lregs[reg][-1]
        zw = create_zw(D, dz)
        weights = create_Levermann_weights(region=reg, zw=zw, area=area, thetao=thetao.isel(time=0).squeeze())
        avgs[i,:] = thetao.weighted(weights).mean(['lev','j','i']).values
    ds = xr.DataArray(dims=['region','time'],
                      coords={'time':thetao.time, 'region':list(Lregs.keys())},
                      data=avgs
                     )
    return ds


def read_larmip2_lrf(data_dir, basal_melt):
    """
    copied from https://github.com/KNMI-sealevel/KPZ/blob/main/notebooks/CompLRF.ipynb
    Read LARMIP2 Linear Response Functions downloaded from:
    https://github.com/ALevermann/Larmip2019.
    Basal melt is in m.y-1. BM02, BM04, BM08, BM16 are available. 
    basal_melt = BM08 is used in Levermann et al. 2020.
    """
    
    reg_names = {'R1':'EAIS', 'R2':'Ross', 'R3':'Amundsen', 
                 'R4':'Weddell', 'R5':'Peninsula'}

    for idb, reg in enumerate(reg_names):
        path = f'{data_dir}/RF_*_{basal_melt}_{reg}.dat'
        files = glob.glob(path)

        for idf, f in enumerate(files):
            ds = pd.read_csv(f, names=['RF']).to_xarray()
            f2 = f.split('/')[-1]
            ds = ds.expand_dims({'model': [f2[3:-12]]})
            
            if idf ==0:
                ds2 = ds
            else:
                ds2 = xr.concat([ds2, ds], dim='model')

        ds2 = ds2.expand_dims({'region': [reg_names[reg]]})
        if idb == 0:
            RF = ds2
        else:
            RF = xr.concat([RF, ds2], dim='region')

    RF = RF.rename({'index' : 'time'})
    RF = RF.transpose('region', 'model', 'time')
    
    return RF.RF


def ensemble_fit(da, deg=1):
    """ fits a polynomial to timeseries """
    assert 'exp' in da.coords
    assert 'year' in da.coords
    if len(da.exp.values.shape)==0:
        n = 1
    elif len(da.exp.values.shape)==1:
        n = len(da.exp.values)
    x = np.concatenate(len(da.mem)*n*[da.year.values])
    y = da.values.flatten(order='F')
    idx = np.isfinite(x) & np.isfinite(y)
    x, y = x[idx], y[idx]
    pf = np.polynomial.polynomial.polyfit(x, y, deg)
    return x, y, pf


if __name__ == '__main__':  # calculate the LARMIP temperatures
    source, experiment, member = subselect_sysargv(sys.argv)
    for da, src, exp, mem in IterateECE(var='thetao',
                                        source=source,
                                        experiment=experiment,
                                        member=member,
                                        cat='jasmin-nc',
                                        only_filenames=(source=='EC-Earth3P-HR')
                                       ):
        print(src, exp, mem)
        fn = f'../results/LARMIP/LARMIP_{src}_{mem}_{exp}.nc'
        if os.path.exists(fn):
            continue

        if src=='EC-Earth3P-HR':  # create time series stepwise
            area = area_511
            for i, fn_ in tqdm(enumerate(da)):
                assert type(da)==list and os.path.exists(fn_)
                fn_temp = f'../results/LARMIP/LARMIP_{src}_{mem}_{exp}_{i}.nc'
                if os.path.exists(fn_temp):  continue
                da_ = xr.open_dataset(fn_).thetao
                ds_ = calc_Levermann_timeseries(area=area, thetao=da_)
                ds_.to_netcdf(fn_temp)
            ds = xr.open_mfdataset(f'../results/LARMIP/LARMIP_{src}_{mem}_{exp}_*.nc')
        elif src in ['EC-Earth3P']:  # create time series in one go
            area = area_255
            ds = calc_Levermann_timeseries(area=area, thetao=da)
        else:
            raise ValueError('src not implemented')
        ds.to_netcdf(fn)