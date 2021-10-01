import os
import sys
activate_file = os.path.join('/home/users/ajuling/nb-venvs/venv-cmip6-zarr/bin/activate_this.py')
exec(open(activate_file).read(), dict(__file__=activate_file))

import dask
import numpy as np
import xarray as xr

from iterate import IterateECE
from tqdm.auto import tqdm

path_data = '/home/users/ajuling/ECE/data'

dz = xr.open_dataarray(f'{path_data}/T511-ORCA025/dz.nc')
zub = xr.open_dataarray(f'{path_data}/T511-ORCA025/zub.nc')
zlb = xr.open_dataarray(f'{path_data}/T511-ORCA025/zlb.nc')

area_255 = xr.open_dataset(f'{path_data}/T255-ORCA1/areas.nc')['O1t0.srf'].rename({'x_3':'i','y_3':'j'})
area_511 = xr.open_dataset(f'{path_data}/T511-ORCA025/areas.nc')['Ot25.srf'].rename({'x_3':'i','y_3':'j'})

# Levermann LARMIP regions (their Tbl. 1)
# [lower lat, upper lat, l. lon, u. lon, central depth]
Lregs = dict(rEAIS = [-76,-65,350,173,369],
             rWedd = [-90,-72,295,350,429],
             rAmun = [-90,-70,210,295,305],
             rRoss = [-90,-76,150,210,312],
             rPen1 = [-70,-65,294,310,420],
             rPen2 = [-75,-70,285,295,420],
)


def create_zw(D, dz, verbose=False):
    """ creates (partial) depth array for cells 100 m around depth D """
    zw = xr.zeros_like(dz)                                        # z-weights
    for i in range(75):
        ub, lb = zub.isel(lev=i).values, zlb.isel(lev=i).values # upper and lower bounds
        assert lb>ub
        if (ub<D-100 and lb<D-100) or (ub>D+100 and lb>D+100):    # cell completely above or below
            if verbose:  print(f'{i:2d}, {ub:6.1f}-{lb:6.1f} outside of {D-100,D+100}')
            continue
        elif ub<D-100 and lb>D-100:                               # upper partial cell
            zw[i] = lb-D+100
            if verbose:  print(f'{i:2d}, {ub:6.1f}-{lb:6.1f} partially in {D-100,D+100}')
        elif ub>D-100 and lb>D-100 and ub<D+100 and lb<D+100:     # full cell
            zw[i] = dz[i]
            if verbose:  print(f'{i:2d}, {ub:6.1f}-{lb:6.1f} fully in {D-100,D+100}')
        elif ub<D+100 and lb>D+100:                               # lower partial cell
            zw[i] = D+100-ub
            if verbose:  print(f'{i:2d}, {ub:6.1f}-{lb:6.1f} partially in {D-100,D+100}')
    if verbose:  print(zw.values, zw.sum().values)
    assert(zw.sum().values==200.)
    return zw


def create_Levermann_weights(region, zw, area, thetao):
    """ creates 3D field with volume in Levermann region/depth """
    if area.shape==(292, 362):      # standard resolution
        fn = f'{path_data}/T255-ORCA1/Levermann_weights_{region}-SR.nc'
    elif area.shape==(1050, 1442):  # high resolution
        fn = f'{path_data}/T255-ORCA1/Levermann_weights_{region}-HR.nc'
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

# for da, src, exp, mem in IterateECE(var='thetao', cat='jasmin-nc'):
# for da, src, exp, mem in IterateECE(var='thetao', source='EC-Earth3P', cat='jasmin-nc'):
# for da, src, exp, mem in IterateECE(var='thetao', source='EC-Earth3P-HR', cat='jasmin-nc'):
if __name__ == '__main__':
    if len(list(sys.argv))==4:
        source = sys.argv[1]
        experiment = sys.argv[2]
        member = sys.argv[3]
    else:
        source, experiment, member = None, None, None
        
    for da, src, exp, mem in IterateECE(var='thetao',
                                            source=source,
                                            experiment=experiment,
                                            member=member,
                                            cat='jasmin-nc',
                                            only_filenames=(source=='EC-Earth3P-HR')
                                           ):    
            print(src, exp, mem)
            fn = f'../results/thetao/thetao_Levermann_{src}_{mem}_{exp}.nc'
            if os.path.exists(fn):
                continue
                
            if src=='EC-Earth3P-HR':
                # create time series stepwise
                area = area_511
                for i, fn_ in tqdm(enumerate(da)):
                    assert type(da)==list and os.path.exists(fn_)
                    fn_temp = f'../results/thetao/thetao_Levermann_{src}_{mem}_{exp}_{i}.nc'
                    if os.path.exists(fn_temp):  continue
                    da_ = xr.open_dataset(fn_).thetao
                    ds_ = calc_Levermann_timeseries(area=area, thetao=da_)
                    ds_.to_netcdf(fn_temp)
                ds = xr.open_mfdataset(f'../results/thetao/thetao_Levermann_{src}_{mem}_{exp}_*.nc')
            elif src in ['EC-Earth3P']:
                # create time series in one go
                area = area_255
                ds = calc_Levermann_timeseries(area=area, thetao=da)
            else:
                raise ValueError('src not implemented')
            ds.to_netcdf(fn)