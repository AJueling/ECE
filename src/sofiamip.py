import os
import copy
import numpy as np
import cftime
import xarray as xr

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from ece_output import open_ece_cmor_output
from ece_output import sofiamip_exps, first_year, table, t264_dir
from ece_output import aco, aca

src_path = os.getcwd()
if src_path[-5:]=='ipynb':
    src_path += '/..'

quants = {
    'tosga'   :'tos',
    'tos60'   :'tos',
    't100so'  :'thetao',
    't500so'  :'thetao',
    't1000so' :'thetao',
    'thetaozm':'thetao',
    'sosga'   :'sos',
    'sos60'   :'sos',
    'amoc'    :'msftyz',
    'siana'   :'siconc',
    'siasa'   :'siconc',
    'hfdsgi'  :'hfds',
    'hfds60'  :'hfds',
    'fwfgi'   :'friver',
    'fwf60'   :'friver',
}

names = {
    'tosga'   :'ocean surface temperature global average',
    'tos60'   :'ocean surface temperature <60°S',
    't100so'  :'ocean potential temperature at 100 m depth <60°S',
    't500so'  :'ocean potential temperature at 500 m depth <60°S',
    't1000so' :'ocean potential temperature at 1000 m depth <60°S',
    'thetaozm':'ocean potential temperature zonal mean',
    'sosga'   :'ocean surface salinity global average',
    'sos60'   :'ocean surface salinity <60°S',
    'amoc'    :'AMOC strength at 30°N, 1000 m',
    'siana'   :'Arctic sea ice area',
    'siasa'   :'Antarctic sea ice area',
    'hfdsgi'  :'ocean surface heat flux global average',
    'hfds60'  :'ocean surface heat flux <60°S',
    'fwfgi'   :'ocean surface freshwater flux global average',
    'fwf60'   :'ocean surface freshwater flux <60°S',
}

t264_quants = {
    'tosga':'sosstsst',
    'sosga':'sosaline',
    'amoc' :'max_amoc_30N',
    'siana':'tot_area_ice_north',
    'siasa':'tot_area_ice_south',
}

grids = {'Ofx':'gn',
         'SImon':'gn',
         'Emon':'gr',
         'Omon':'gn',
         'LImon':'gr',
         'fx':'gr',
         'Lmon':'gr',
         'Amon':'gr'
        }



class SofiaMip():
    """
    input:
    quant
    """
    def __init__(self, quant, plot=False, recalculate=False) -> None:
        print(src_path)
        assert quant in quants.keys()
        self.quant = quant
        self.var = quants[quant]
        self.table = table[self.var]
        self.prep = self.preprocess()

        # piControl data
        self.pifn = f'{src_path}/../results/sofiamip-timeseries/{quant}_pi.nc'
        if os.path.exists(self.pifn) and recalculate==False:
            self.pi = xr.open_dataarray(self.pifn)
        # else:
        #     print('calculating piControl')
        #     pi_fn = 'EC-Earth3_piControl_r1i1p1f1_gn'
        #     t = self.table
        #     v = self.var
        #     fn = f'{t264_dir}/{t}/{v}/{v}_{t}_{pi_fn}_*.nc'
        #     # if os.path.exists(fn):
        #     print(fn)
        #     self.pi = xr.open_mfdataset(fn, decode_times=False, preprocess=self.prep)[v]
        #     self.pi = self.calculations(ts=self.pi)
        #     self.pi.to_netcdf(self.pifn)
        #     # else:
        #     #     self.pi = None

        # t264 data
        if quant in t264_quants.keys():
            t264 = xr.open_dataset(f'{src_path}/../data/t264/ocean/t264_2160_2459_time-series_ocean.nc')
            self.t264 = t264[t264_quants[quant]]
            if quant in ['siana','siasa']:
                self.t264 /= 100

        # experiment data
        self.ts = {}
        for expn in sofiamip_exps:
            print(expn)
            fn = f'{src_path}/../results/sofiamip-timeseries/{quant}_{expn}.nc'
            if os.path.exists(fn) and recalculate==False:
                self.load_timeseries(fn, expn)
            else:
                self.calc_timeseries(fn, expn)

        self.shift_times()

        if plot==True:
            self.plot_timeseries()
        
        return

    def load_timeseries(self, fn, expn):
        """ """
        print(f'loading {self.quant} ...')
        self.ts[expn] = xr.open_dataarray(fn)
        return

    def calc_timeseries(self, fn, expn):
        """ """
        self.ts[expn] = open_ece_cmor_output(expn=expn, var=self.var, preprocess=self.prep)
        self.ts[expn] = self.calculations(self.ts[expn])
        self.ts[expn].to_netcdf(fn)
        return
    
    def preprocess(self):
        """ preprocessing functions for xr.load_mfdataset() """
        def select_500m(ds):
            print('selecting 500m level')
            return ds.sel({'lev':500}, method='nearest')
        def select_100m(ds):
            return ds.sel({'lev':100}, method='nearest')
        def select_1000m(ds):
            return ds.sel({'lev':1000}, method='nearest')
        def zonal_mean(ds):
            return ds.mean('i')

        if self.quant=='t100so':      prep = select_100m
        elif self.quant=='t500so':    prep = select_500m
        elif self.quant=='t1000so':   prep = select_1000m
        elif self.quant=='thetaozm':  prep = zonal_mean
        else:                         prep = None
        return prep

    def calculations(self, ts):
        """ actual xr calculations based in quantity requested """
        print(f'calculating {self.quant} ...')
        if self.quant in ['tosga','sosga']:
            # ocean global average
            ts = ts.weighted(aco).mean(['i','j'])
        elif self.quant in ['tasga','prga']:
            # atmospheric global average
            ts = ts.weighted(aca).mean(['i','j'])
        elif self.quant in ['hfdsgi','fwfgi']:
            # global area integral
            ts = (ts*aco).sum(['i','j'])
        elif self.quant in ['tos60','sos60','t500so']:
            # Southern Ocean average
            ts = ts.where(aco.latitude<-60).weighted(aco).mean(['i','j'])
        elif self.quant in ['hfds60','fwf60']:
            # Southern Ocean area integral
            ts = (ts*aco).where(aco.latitude<-60).sum(['i','j'])
        elif self.quant=='amoc':
            ts = ts.sel({'rlat':30, 'lev':1000}, method='nearest').isel(basin=1)
        elif self.quant=='siana':
            ts = (ts*aco/100).where(aco.latitude>0).sum(['i','j'])
        elif self.quant=='siasa':
            ts = (ts*aco/100).where(aco.latitude<0).sum(['i','j'])
        return ts.load()
    
    def shift_times(self):
        """ shifts t264, pi, and experiments
        
        """
        # if self.table=='Omon':
        t = self.pi.time.values[0]
        # print(f'\nfirst element of self.pi.time:\n{t}\n')
        if type(t)==int:
            # print('convert from integer time in days since 1850 to float year C.E.')
            self.pi = self.pi.assign_coords(time=self.pi.time/365+1850)
        elif type(t)==cftime._cftime.DatetimeProlepticGregorian:
            # print('convert from Gregorian calendar to year')
            times = [t_.year + t_.month/12 - 1/24 for t_ in self.pi.time.values]
            self.pi = self.pi.assign_coords(time=times)
        else:
            print('did not convert')

        self.ts_shifted = copy.deepcopy(self.ts)
        self.pi_shifted = copy.deepcopy(self.ts)

        for expn in sofiamip_exps:
            # shifting experiments onto piControl time  axis
            pi_shift = (first_year[expn]-2259)*12
            ts_len = len(self.ts[expn].values)
            islice = slice(pi_shift, pi_shift+ts_len)
            timeslice = self.pi.time.isel(time=islice)
            self.ts_shifted[expn] = self.ts_shifted[expn].assign_coords(time=timeslice)
            self.pi_shifted[expn].values = self.pi.isel(time=islice)

            # t264 ...
        return


    def plot_timeseries(self, ylabel=None):
        """ """
        fig = plt.figure(figsize=(12,6), tight_layout=True)
        gs = gridspec.GridSpec(2, 2)
        ax0 = fig.add_subplot(gs[0, :])  # full piControl time series
        ax1 = fig.add_subplot(gs[1, 0])  # time series with shifted pi and fits
        ax2 = fig.add_subplot(gs[1, 1])  # anomalies

        ax0.plot(self.pi.time, self.pi.rolling(time=12).mean(), c='grey')
        ax2.axhline(0, c='k', lw=.4)
        for i, expn in enumerate(sofiamip_exps):
            c = f'C{i}'
            # piControl overview
            ts = self.ts_shifted[expn]
            # time = ga.time.dt.year + first_year[expn]-1850 + ga.time.dt.month/12
            ax0.plot(ts.time, ts.rolling(time=12).mean(), c=c)
            # if self.quant in t264_quants:
                # plot additional time series

            # shifted to same starting point
            pi = self.pi_shifted[expn]
            fit = pi.polyfit(dim='time', deg=2).polyfit_coefficients.values
            time = pi.time.astype('float')
            pi_quadfit = fit[0]*time**2 + fit[1]*time + fit[2]*np.ones(len(time))
            ax1.plot(pi.time, pi.rolling(time=12).mean(), c=c, lw=.3)
            ax1.plot(pi.time, pi_quadfit, c=c)
            ax1.plot(self.ts[expn].time, self.ts[expn].rolling(time=12).mean(), c=c, lw=2.)

            # anomalies
            ax2.plot(self.ts[expn].time, (self.ts[expn]-pi_quadfit).rolling(time=12).mean())

        # making plot nicer
        for ax in [ax0, ax1, ax2]:
            ax.set_xlabel('time [years]')
            if ylabel is not None:
                ax.set_ylabel(ylabel)
        for ax in [ax1, ax2]:
            ax.set_xlim((np.datetime64('1845-01-01','ns'),np.datetime64('1955-01-01','ns')))

        fig.suptitle(names[self.quant])
        fig.align_labels()  # same as fig.align_xlabels(); fig.align_ylabels()
        fig.savefig(f'../../results/sofiamip-plots/sofiamip-timeseries-{self.quant}.pdf')
        return

    def plot_maps(self):
        """ """
        return
    
    def plot_zonal_means(self):
        """ """
        return