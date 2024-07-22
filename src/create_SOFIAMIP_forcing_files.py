"""
This class creates the freshwater forcing files that work with the fwf branch of EC-Earth3

Paths might need to be updated

created by André Jüling
last updated 2024-07-18
"""
import os
import sys
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pathlib import Path


spy = 3600*24*365   # [s yr^-1]

path_grid = '/home/nkaj/ECE/data'
path_fwf = '../../ECE/results/FWF'
path_fwfmap_files = '../'

sofiamip_exps = [
    'antwater',
    'hist-antwater-70-01',
    'hist-antwater-70-03',
    'hist-antwater-70-05',
    'hist-antwater-92-11',
    'ssp126-ismip6-water',
    'ssp585-ismip6-water',
    'antwater-60S',
    'antwater-lh',
    'antwater-ambe',
    'antwater-depth',
    'antwater-depth-lh',
    'antwater-ambe-depth-lh',
    ]

class FwfMaps:
    """ create freshwater forcing files for EC-Earth experiments 
    output:
    netcdf file with 
    - runoff flux [kg/m^2/s]
    - basal melt flux [kg/m^2/s]
    - time 12 months (important that time axis is unlimited)
    - horizontal coordinates (i,j) as the model resolution requires
    
    """
    def __init__(self, exp, create=True, test=True):
        """
        - reads areacello.nc for respective grid
        - initialize with 
        exp     experiment names
        """
        self.exp = exp
        assert exp in sofiamip_exps

        self.gnem, self.gifs = '1', '255'
        fn_areacello = f'{path_grid}/T{self.gifs}-ORCA{self.gnem}/areacello.nc'
        if os.path.exists(fn_areacello):    
            self.areacello = xr.open_dataset(fn_areacello).areacello
        else:
            raise ValueError(f'areacello filename: {fn_areacello} does not exist')
        
        if self.exp.startswith('ssp') or self.exp.startswith('hist'):
            self.fn_out = f'/ec/res4/hpcperm/nkaj/forcing/SOFIAMIP/{exp}/FWF_{exp}_ORCA{self.gnem}'
        else:
            self.fn_out = f'/ec/res4/hpcperm/nkaj/forcing/SOFIAMIP/{exp}/FWF_{exp}_ORCA{self.gnem}.nc'

        if create==True:
            self.create_forcing_files()
        if test==True:
            self.test_output()
        return

    def create_forcing_files(self, year=None):
        """ """

        Path(f'/ec/res4/hpcperm/nkaj/forcing/SOFIAMIP/{self.exp}').mkdir(parents=True, exist_ok=True)
        if self.exp=='antwater':
            self.create_antwater_file()
        elif self.exp=='antwater-lh':
            self.create_antwater_lh_file()
        elif self.exp=='antwater-depth-lh':
            self.create_antwater_depth_lh_files()
        elif self.exp.startswith('hist') or self.exp.startswith('ssp'):
            self.create_hist_ssp_files()
        return

    def create_Antarctic_coast(self):
        ac = self.areacello
        ac = ac.fillna(0)
        rmf = xr.open_dataarray(f'{path_grid}/T{self.gifs}-ORCA{self.gnem}/runoff_maps_AIS_GrIS_ORCA{self.gnem}.nc')

        # select coastal cells around Antarctica
        coast = xr.where(ac-ac.shift(i= 1)==ac, ac, 0) + \
                xr.where(ac-ac.shift(i=-1)==ac, ac, 0) + \
                xr.where(ac-ac.shift(j= 1)==ac, ac, 0) + \
                xr.where(ac-ac.shift(j=-1)==ac, ac, 0)
        coast = coast/coast*ac

        ac_antwater = coast.where(rmf==66).where(coast>0)
        ac_antwater = xr.where(ac_antwater,ac_antwater,0)
        print(f'integrated area of circum-Antarctic cells: {ac_antwater.sum().values/1e12} million km^2')
        print(f'length circum-Antarctic cells at 100 km width: {ac_antwater.sum().values/1e5/1e3} km')
        return ac_antwater

    def create_dataset_structure(self, year=1850):
        # create "scaffold" dataset
        aco = self.areacello
        time = np.array([np.datetime64(f'{year}-{m:02d}-01 12:00:00.000000000') for m in range(1,13)])
        ds = xr.Dataset(coords={'time_counter':time,
                                'i':aco.i.values,
                                'j':aco.j.values,
                                'longitude': (('j','i'), aco.longitude.values),
                                'latitude' : (('j','i'), aco.latitude.values),
                                },
                        data_vars={'sofwfrnf':(('time_counter','j','i'),np.zeros((12,len(aco.j),len(aco.i)))),
                                   'sofwfcal':(('time_counter','j','i'),np.zeros((12,len(aco.j),len(aco.i)))),
                                   'sofwfisf':(('time_counter','j','i'),np.zeros((12,len(aco.j),len(aco.i)))),
                                   }
                        )
        for i, quant in enumerate(['sofwfrnf','sofwfcal','sofwfisf']):
            if quant=='sofwfrnf': 
                long_name = 'forced runoff flux'
            elif quant=='sofwfcal':
                long_name = 'forced calving flux'
            elif quant=='sofwfisf':
                long_name = 'forced basal melt flux'
            ds[quant] = ds[quant].fillna(0)
            ds[quant].attrs = {'long_name':long_name, 'units':'kg/m^2/s'}
        return ds

    def create_antwater_forcing(self):
        ac_antwater = self.create_Antarctic_coast()
        # scale by zonal extent only, i.e. divide by meridional extent of cells
        # distribute 0.1 Sv = 1e5 m^3/s over the area to create a flux density map [m^3/s/m^2] = [m/s]
        latdiff = np.cos(np.deg2rad(ac_antwater.latitude))
        # latdiff = (ac_antwater.latitude.shift(j=-1)-ac_antwater.latitude.shift(j=1))/2
        fwf = 1e5/(ac_antwater/latdiff).sum(['i','j'])*(ac_antwater/latdiff)  # rate in  [m^3/s]
        fwf = fwf/ac_antwater*1e3  # [m^3/s] -> [m/s] -> [kg/m^2/s] @ 1000 kg/m^3
        return fwf

    def create_antwater_file(self):
        """
        - climatological file without yearly variations
        - uniformly distributes 0.1 Sv in Antarctica-adjacent NEMO cells
        - assign 0.1 Sv to runoff forcing, calving & basal melt remain 0
        """
        ds = self.create_dataset_structure()
        fwf = self.create_antwater_forcing()
        for t in range(12):
            ds['sofwfrnf'][t,:,:] = fwf
            ds['sofwfrnf'] = ds['sofwfrnf'].fillna(0)
        ds.to_netcdf(self.fn_out, unlimited_dims=['time_counter'])
        return

    def create_antwater_lh_file(self):
        """ assign 0.1 Sv to calving forcing, runoff & basal melt remain 0 """
        ds = self.create_dataset_structure()
        fwf = self.create_antwater_forcing()
        for t in range(12):
            ds['sofwfcal'][t,:,:] = fwf
            ds['sofwfcal'] = ds['sofwfcal'].fillna(0)
        ds.to_netcdf(self.fn_out, unlimited_dims=['time_counter'])
        return

    def create_antwater_depth_lh_files(self):
        """ assign 0.1 Sv to calving forcing, runoff & basal melt remain 0 """
        ds = self.create_dataset_structure()
        fwf = self.create_antwater_forcing()
        for t in range(12):
            ds['sofwfcal'][t,:,:] = fwf
            ds['sofwfcal'] = ds['sofwfcal'].fillna(0)
        ds.to_netcdf(self.fn_out, unlimited_dims=['time_counter'])
        return

    def time_dependent_factor(self, year, month):
        if self.exp.startswith('hist'):
            # the 10 in there because the antwater forcing it is multiplied with is 0.1 Sv
            # so now the factor*antwater adds up to the required amount
            if self.exp.endswith('70-01'):
                t = (year-1970) + month/12
                rate = 0.1e-3  # Sv/yr
            elif self.exp.endswith('70-03'):
                t = (year-1970) + month/12
                rate = 0.3e-3
            elif self.exp.endswith('70-05'):
                t = (year-1970) + month/12
                rate = 0.5e-3
            elif self.exp.endswith('92-11'):
                t = (year-1992) + month/12
                rate = 1.1e-3
            factor = 10*rate*t
        elif self.exp=='ssp126-ismip6-water':
            factor = 0.15
            # this multiplied with 0.1 Sv results in the 0.015 Sv
        elif self.exp=='ssp585-ismip6-water':
            # see eqn A1 and Table A2 in Swart et al.
            # fit in Gt/yr
            t = (year-2015) + month/12
            A = 5.18e2
            K = 3.14e3
            B = 0.21
            v = 3.85e-5
            Q = 1.48e1
            factor = (A+(K-A)/(1+Q*np.exp(-B*t)))/.55  # in Gt/yr
            Sv_per_Gtpy = 3.17e-5
            factor *= Sv_per_Gtpy  # [Gt/yr] -> [Sv]
        else:
            raise ValueError(f'exp should be hist or ssp, but exp={self.exp}')
        return factor
    
    def create_hist_ssp_files(self):
        fwf = self.create_antwater_forcing()
        if self.exp.startswith('hist-antwater-70'):
            start_year = 1970
            end_year = 2020
        if self.exp.startswith('hist-antwater-92'):
            start_year = 1992
            end_year = 2020
        if self.exp.startswith('ssp'):
            start_year = 2015
            end_year = 2100
        for y in range(start_year,end_year+1):
            ds = self.create_dataset_structure(year=y)
            for t in range(12):
                factor = self.time_dependent_factor(year=y,month=t+1)
                ds['sofwfrnf'][t,:,:] = factor*fwf
                ds['sofwfrnf'] = ds['sofwfrnf'].fillna(0)
            ds.to_netcdf(f'{self.fn_out}_y{y}.nc', unlimited_dims=['time_counter'])
        return

    def test_output(self):
        """ calculates the annual average """
        ac = self.areacello
        if self.exp.startswith('ssp') or self.exp.startswith('hist'):
            ds = xr.open_mfdataset(self.fn_out+'*')
            rnf = (ac*ds.sofwfrnf).sum(['i','j']).values
            cal = (ac*ds.sofwfcal).sum(['i','j']).values
            isf = (ac*ds.sofwfisf).sum(['i','j']).values
            tot = rnf + cal + isf
            plt.figure()
            time = ds['time_counter']
            plt.plot(time, rnf/1e9, label='rnf')
            plt.plot(time, cal/1e9, label='cal')
            plt.plot(time, isf/1e9, label='isf')
            plt.plot(time, tot/1e9, label='tot', ls='--')
            plt.legend()
            plt.title(self.exp)
            plt.ylabel('time  [year]')
            plt.ylabel('freshwater flux  [Sv]')
            print(f'The average runoff     freshwater flux is: {rnf.mean()/1e9:.4f} Sv')  # [kg/m^2/s]*[m^2] -> [Sv]
            print(f'The average calving    freshwater flux is: {cal.mean()/1e9:.4f} Sv')  # [kg/m^2/s]*[m^2] -> [Sv]
            print(f'The average basal melt freshwater flux is: {isf.mean()/1e9:.4f} Sv')  # [kg/m^2/s]*[m^2] -> [Sv]
            print(f'The average total      freshwater flux is: {tot.mean()/1e9:.4f} Sv')  # [kg/m^2/s]*[m^2] -> [Sv]

        else:  # climatological files
            assert os.path.exists(self.fn_out)
            ds = xr.open_dataset(self.fn_out)
            mean_rate_rnf = (ac*ds.sofwfrnf).mean('time_counter').sum(['i','j']).values
            mean_rate_cal = (ac*ds.sofwfcal).mean('time_counter').sum(['i','j']).values
            mean_rate_isf = (ac*ds.sofwfisf).mean('time_counter').sum(['i','j']).values
            mean_rate_tot = mean_rate_rnf + mean_rate_cal + mean_rate_isf
            print(f'The average runoff     freshwater flux is: {mean_rate_rnf /1e9} Sv')  # [kg/m^2/s]*[m^2] -> [Sv]
            print(f'The average calving    freshwater flux is: {mean_rate_cal/1e9} Sv')  # [kg/m^2/s]*[m^2] -> [Sv]
            print(f'The average basal melt freshwater flux is: {mean_rate_isf/1e9} Sv')  # [kg/m^2/s]*[m^2] -> [Sv]
            print(f'The average total      freshwater flux is: {mean_rate_tot/1e9} Sv')  # [kg/m^2/s]*[m^2] -> [Sv]
        return