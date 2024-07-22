import os
import xarray as xr
import matplotlib.pyplot as plt

aco = xr.open_dataset('../data/T255-ORCA1/areacello.nc').areacello

pic_dir = '/ec/res4/scratch/nkaj/cmorised-results/piControl'
sdl2_dir = '/ec/res4/scratch/nkaj/cmorised-results/sdl2/v001/CMIP6/SOFIAMIP/KNMI/EC-Earth3/faf-antwater/r1i1p1f1'
sf02_dir = '/ec/res4/scratch/nkaj/cmorized-results/sf02/CMIP6/SOFIAMIP/KNMI/EC-Earth3/faf-antwater/r1i1p1f1'

def Ross_thetao500_timeseries(expn, depth=500):
    """ global mean surface temperature 30-year average """
    
    def preprocess(ds):
        if depth==0:
            lev = 0
        elif depth==500:
            lev = 39
        return ds.isel(lev=lev)
    
    aco_ = aco.where(aco.latitude<-73).where(aco.longitude<180).where(aco.longitude>160).fillna(0)
    
    if expn=='t264':
        kwargs = {'use_cftime':True}
    else:
        kwargs = {}
    
    fn = f'../results/sofiamip-timeseries/Ross_thetao_{depth}_timeseries_{expn}.nc'
    if os.path.exists(fn):
        print(f'loading exp={expn} thetao climatology at {depth} m')
        thetao_clim = xr.open_dataarray(fn)
    else:
        print(f'calculating exp={expn} thetao climatology at {depth} m')
        if expn=='sdl2':
            fns = f'{sdl2_dir}/Omon/thetao/gn/v20240404/*.nc'
            tslice = slice(-12*30,None)
        elif expn=='sf02':
            fns = f'{sf02_dir}/Omon/thetao/gn/v20240404/*.nc'
            tslice = slice(-12*30,None)
        elif expn=='t264':
            fns = f'{pic_dir}/Omon/thetao/*.nc'
            tslice = slice(12*(2300+70-2259),12*(2300+100-2259))
        thetao_clim = xr.open_mfdataset(fns, preprocess=preprocess, **kwargs).thetao.weighted(aco_).mean(['i','j']).compute().chunk(1)
        thetao_clim.to_netcdf(fn)
    return thetao_clim


def tas_timeseries(expn):
    """ global mean surface temperature monthly average time series """
    fn = f'../../results/sofiamip-timeseries/tasga_{expn}.nc'
    fn_60S = f'../../results/sofiamip-timeseries/tas_60S_{expn}.nc'
    
    if expn=='t264':
        kwargs = {'use_cftime':True}
    else:
        kwargs = {}

    if expn=='sdl2':
        fns = f'{sdl2_dir}/Amon/tas/gr/v20240404/*.nc'
    elif expn=='sf02':
        fns = f'{sf02_dir}/Amon/tas/gr/v20240404/*.nc'
    elif expn=='t264':
        fns = f'{pic_dir}/Amon/tas/*.nc'

    if os.path.exists(fn)==False or os.path.exists(fn_60S)==False:
        tas = xr.open_mfdataset(fns, **kwargs).tas

    if os.path.exists(fn)==False:
        print('calculating timeseries')
        tasga = tas.weighted(aca).mean(['lat','lon'])#.compute().chunk(1)
        tasga.to_netcdf(fn)
    else:
        print('loading timeseries')
        tasga = xr.open_dataarray(fn, **kwargs)
    # tasga = None

    if os.path.exists(fn_60S)==False:
        print('calculating timeseries')
        weights = aca.where(aca.lat<=-60).fillna(0)
        print('calculating timeseries')
        tas60S = tas.weighted(weights).mean(['lat','lon']).compute().chunk(1)
        print('calculating timeseries')
        tas60S.to_netcdf(fn_60S)
    else:
        print('loading timeseries')
        tas60S = xr.open_dataarray(fn_60S, **kwargs)

    return tasga, tas60S


if __name__=='__main__':
    r500_pic  = Ross_thetao500_timeseries(expn='t264')
    r500_sf02 = Ross_thetao500_timeseries(expn='sf02')
    r500_sdl2 = Ross_thetao500_timeseries(expn='sdl2')

plt.plot(range(2259,2760), r500_pic )
plt.plot(range(2301,2398), r500_sf02)
plt.plot(range(2300,2400), r500_sdl2)