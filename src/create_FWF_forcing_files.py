""" create forcing files

PRIMAVERA
- read time series, created in `create_FWF_timeseries.ipynb`
- ORCA1 and ORCA025
- Greenland only complex
- Greenland and Antarctica simple
- Greenland and Antarctica complex
- complex basal melt forcing saved as ...

simple: only via sbcfwf module; variable names sofwfrnf, socalving_f
complex: also via isf module, variable name (also)


SOFIAMIP
- antwater: coast surface, fixed 0.1 Sv
- antwater-hist: coastal surface, increasing
- antwater-depth: coast surface, fixed 0.1 Sv
- antwater-60S
- antwater-ambe
- antwater-ambe-depth

"""

import datetime

from testing import test_file_existence

output_path = '/ec/res4/hpcperm/nkaj/forcing'

"""


"""
import os
import sys
import numpy as np
import xarray as xr

spy = 3600*24*365   # [s yr^-1]

path_fwf = '../results/FWF'
path_data = '../data'
path_output = '/ec/res4/hpcperm/nkaj/forcing'


class FwfMaps:
    """ create freshwater forcing files for EC-Earth experiments 
    input:
    time series of forcing

    output:
    netcdf file with 
    - runoff flux [kg/m^2/s]
    - calving flux [kg/m^2/s]
    - basal melt flux [kg/m^2/s]
    - time 12 months (important that time axis is unlimited)
    - horizontal coordinates (i,j) as the model resolution requires
    """

    def __init__(self, res, exp):
        """
        - initialize with resolution ('SR' or 'HR')and experiment names
        - reads areacello
        """
        self.res = res
        self.exp = exp
        assert exp in ['simple_GrIS_AIS', 'complex_GrIS_AIS', 'complex_GrIS_only', 'SOFIAMIP']

        if res=='SR':    self.gnem, self.gifs = '1', '255'
        elif res=='HR':  self.gnem, self.gifs = '025', '511'

        fn_areacello = f'{path_data}/T{self.gifs}-ORCA{self.gnem}/areacello.nc'
        if os.path.exists(fn_areacello):    
            self.areacello = xr.open_dataset(fn_areacello).areacello
        else:
            raise ValueError(f'areacello filename: {fn_areacello} does not exist')

        return

    def create_map(self, year=None):
        """ """
        if self.exp in ['simple_GrIS_AIS', 'complex_GrIS_AIS', 'complex_GrIS_only']:
            assert year is not None
            dG = xr.open_dataset(f'{path_fwf}/GrIS_FWF.nc')
            if self.exp in ['simple_GrIS_AIS', 'complex_GrIS_AIS']:
                dA = xr.open_dataset(f'{path_fwf}/AIS_FWF.nc')

            if year=='all':
                years = np.arange(1985,2101)
            elif type(year)==int:
                assert year in dG.year
                years = np.arange(year,year+1)
            else:
                raise ValueError('year must be "all" or integer in range(1985,2101)')

            for yr in years:
                # Greenland discharge and runoff
                dG_ = dG.sel(year=yr)
                self.GD = dG_.D*dG.D_seasonality
                self.GR = dG_.R*dG.R_seasonality

                # Antarctica discharge and basal melt
                if self.exp in ['simple_GrIS_AIS', 'complex_GrIS_AIS']:
                    dA_ = dA.sel(year=yr)
                    self.AD = dA_.D*dA.D_seasonality
                    self.AB = dA_.B*dA.B_seasonality

                self.fn_out = f'{path_output}/{self.exp}_forcing/FWF_{self.exp}_ORCA{self.gnem}_y{yr}.nc'

                if self.exp=='simple_GrIS_AIS':
                    self.create_simple_map(year=yr)
                elif self.exp in ['complex_GrIS_only','complex_GrIS_AIS']:
                    self.create_complex_maps(year=yr)
                
                self.test_output()

        elif self.exp=='SOFIAMIP':
            self.fn_out = f'../fafmip_forcing/FWF_fafmip_antwater_ORCA{self.gnem}.nc'
            self.create_fafmip_map()
            self.test_output()


    def create_simple_map(self, year):
        """
        - calculates total Greenland and Antarctic surface areas
        - uniformly distributes meltwater in
        # output: creates 19 GB of files for ORCA025 (for 2 variables only)
        # zipped them with `tar -cvzf simple_forcing.tar.gz simple_forcing`
        """
        rmf = xr.open_dataarray(f'{path_data}/T{self.gifs}-ORCA{self.gnem}/runoff_maps_AIS_GrIS_ORCA{self.gnem}.nc')
        Aa = self.areacello.where(rmf==66).sum()
        Ga = self.areacello.where(rmf==1).sum()

        for i, quant in enumerate(['sofwfrnf','sofwfcal']):
            if quant=='sofwfrnf': 
                # GrIS: Runoff
                Q = xr.where(rmf==1,self.GR/Ga,0)
                long_name = 'forced runoff flux'
            elif quant=='sofwfcal':
                # GrIS & AIS: Discharge
                Q = (xr.where(rmf==1,self.GD/Ga,0) + xr.where(rmf==66,(self.AD+self.AB)/Aa,0))
                long_name = 'forced calving flux'
            Q = Q*1e12/spy  # [Gt/yr] -> [kg/s]
            Q.attrs = {'long_name':long_name, 'units':'kg/m^2/s'}
            Q = Q.rename({'lat':'latitude','lon':'longitude'})
            t = np.array([np.datetime64(f'{year}-{m:02d}-01 12:00:00.000000000') for m in Q.month.values])
            Q = Q.assign_coords({'time_counter':('month', t)})
            Q = Q.drop(['month','year']).rename({'month':'time_counter'})
            if i==0:
                Q.name = quant
                FWF = Q.to_dataset().transpose('time_counter',...)
            else:
                FWF[quant] = Q.transpose('time_counter',...)
        FWF.to_netcdf(self.fn_out, unlimited_dims=['time_counter'])
        return

    def create_complex_maps(self, year):
        """
        runoff: Antarctic = negligible 
        calving: iceberg melt = Merino and MArsh melt pattern,
        basal melt: using isf module to vertically distribute along ice shelf fronts
        """
        # spatial distribution
        # Greenland melt from Bamber2018
        # Greenlans calving iceberg melt from Marsh2015
        # Antarctic calving from Merino...
        # Antarctic basal melt from Mathiot2017

        rmf = xr.open_dataarray(f'{path_data}/T{self.gifs}-ORCA{self.gnem}/runoff_maps_AIS_GrIS_ORCA{self.gnem}.nc')
        Aa = self.areacello.where(rmf==66).sum()
        Ga = self.areacello.where(rmf==1).sum()

        for i, quant in enumerate(['sofwfrnf','sofwfcal','sofwfisf']):
            if quant=='sofwfrnf': 
                # GrIS: Runoff
                Q = xr.where(rmf==1,self.GR/Ga,0)
                long_name = 'forced runoff flux'
            elif quant=='sofwfcal':
                # GrIS & AIS: Discharge
                Q = (xr.where(rmf==1,self.GD/Ga,0)
                if self.exp=='complex_GrIS_AIS':
                    Q += xr.where(rmf==66,self.AD/Aa,0))
                long_name = 'forced calving flux'
            elif quant=='sofwfisf' and self.exp=='complex_GrIS_AIS':
                # AIS:  Basal Melt
                Q = (xr.where(rmf==66,self.AB/Aa,0))
                long_name = 'forced basal melt flux'

            Q = Q*1e12/spy  # [Gt/yr] -> [kg/s]
            Q.attrs = {'long_name':long_name, 'units':'kg/m^2/s'}
            Q = Q.rename({'lat':'latitude','lon':'longitude'})
            t = np.array([np.datetime64(f'{year}-{m:02d}-01 12:00:00.000000000') for m in Q.month.values])
            Q = Q.assign_coords({'time_counter':('month', t)})
            Q = Q.drop(['month','year']).rename({'month':'time_counter'})
            if i==0:
                Q.name = quant
                FWF = Q.to_dataset().transpose('time_counter',...)
            else:
                FWF[quant] = Q.transpose('time_counter',...)
        FWF.to_netcdf(self.fn_out, unlimited_dims=['time_counter'])
        pass

    # def create_fafmip_map(self):
    #     """
    #     - climatological file without yearly variations
    #     - uniformly distributes 0.1 Sv in Antarctica-adjacent NEMO cells
    #     """
    #     ac = self.areacello
    #     ac = ac.fillna(0)
    #     rmf = xr.open_dataarray(f'{path_data}/T{self.gifs}-ORCA{self.gnem}/runoff_maps_AIS_GrIS_ORCA{self.gnem}.nc')

    #     # select coastal cells around Antarctica
    #     coast = xr.where(ac-ac.shift(i= 1)==ac, ac, 0) + \
    #             xr.where(ac-ac.shift(i=-1)==ac, ac, 0) + \
    #             xr.where(ac-ac.shift(j= 1)==ac, ac, 0) + \
    #             xr.where(ac-ac.shift(j=-1)==ac, ac, 0)
    #     coast = coast/coast*ac

    #     ac_fafmip = coast.where(rmf==66).where(coast>0)
    #     ac_fafmip = xr.where(ac_fafmip,ac_fafmip,0)
    #     print(f'integrated area of circum-Antarctic cells: {ac_fafmip.sum().values/1e12} million km^2')
    #     print(f'length circum-Antarctic cells at 100 km width: {ac_fafmip.sum().values/1e5/1e3} km')

    #     # scale by zonal extent only, i.e. divide by meridional extent of cells
    #     # distribute 0.1 Sv = 1e5 m^3/s over the area to create a flux density map [m^3/s/m^2] = [m/s]
    #     latdiff = np.cos(np.deg2rad(ac_fafmip.latitude))
    #     # latdiff = (ac_fafmip.latitude.shift(j=-1)-ac_fafmip.latitude.shift(j=1))/2
    #     fwf = 1e5/(ac_fafmip/latdiff).sum(['i','j'])*(ac_fafmip/latdiff)  # rate in  [m^3/s]
    #     fwf = fwf/ac_fafmip*1e3  # [m^3/s] -> [m/s] -> [kg/m^2/s] @ 1000 kg/m^3

    #     # create "scaffold" dataset
    #     time = np.array([np.datetime64(f'1850-{m:02d}-01 12:00:00.000000000') for m in range(1,13)])
    #     ds = xr.Dataset(coords={'time_counter':time,
    #                             'i':ac.i.values,
    #                             'j':ac.j.values,
    #                             'longitude': (('j','i'), ac.longitude.values),
    #                             'latitude' : (('j','i'),  ac.latitude.values),
    #                             },
    #                     data_vars={'sofwfrnf' :(('time_counter','j','i'),np.zeros((12,len(ac.j),len(ac.i)))),
    #                                'socalving_f':(('time_counter','j','i'),np.zeros((12,len(ac.j),len(ac.i))))}
    #                     )
    #     # assign 0.1 Sv to runoff forcing, calving remains 0
    #     for t in range(12):
    #         ds['sofwfrnf'][t,:,:] = fwf
        
    #     for i, quant in enumerate(['sofwfrnf','socalving_f']):
    #         if quant=='sofwfrnf': 
    #             long_name = 'forced runoff flux'
    #         elif quant=='socalving_f':
    #             long_name = 'forced calving flux'
    #         ds[quant] = ds[quant].fillna(0)
    #         ds[quant].attrs = {'long_name':long_name, 'units':'kg/m^2/s'}

    #     ds.to_netcdf(self.fn_out, unlimited_dims=['time_counter'])
    #     return

    # def create_eveline_map(self):
    #     pass
    
    # def test_output(self):
        """ """
        ds = xr.open_dataset(self.fn_out)
        ac = self.areacello
        mean_rate_runoff  = (ac*ds.sofwfrnf ).mean('time_counter').sum(['i','j']).values
        mean_rate_calving = (ac*ds.socalving_f).mean('time_counter').sum(['i','j']).values
        mean_rate_total   = mean_rate_runoff + mean_rate_calving
        print(f'The average runoff  freshwater flux is: {mean_rate_runoff /1e9} Sv')  # [kg/m^2/s]*[m^2] -> [Sv]
        print(f'The average calving freshwater flux is: {mean_rate_calving/1e9} Sv')  # [kg/m^2/s]*[m^2] -> [Sv]
        print(f'The average total   freshwater flux is: {mean_rate_total  /1e9} Sv')  # [kg/m^2/s]*[m^2] -> [Sv]
        
        return


if __name__=='__main__':
    """ just for eveline's on-the-fly creation of the maps """
    print(sys.argv)
    resolution = sys.argv[1]
    expname = sys.argv[2]
    variable = sys.argv[3]
    print(f'hello world: {variable}')
    for res in ['SR','HR']:
        for scenario in ['simple_GrIS_AIS','complex_GrIS_only','complex_GrIS_AIS']
    FwfMaps()



#     FwfMaps(res=resolution, exp=expname).create_simple_map(variable)


# def create_PRIMAVERA_forcing():
#     '''
#     freshwater_forcing_yxxxx.nc: sofwfisf, sofwfrnf, sofwfcal

    
#     '''
#     # load time series -> how much when
#     '../results/FWF/GrIS_FWF.nc'
#     '../results/FWF/GrIS_AIS.nc'
#     # load weighted distribution maps -> where

#     for scn in ['GrIS_complex','GrIS_AIS_complex','GrIS_AIS_simple']:
#         for res in ['1','025']:
#             for yr in range(1985,2101):
#                 fn_out = f'{output_path}/PRIMAVERA/FWF_{scn}_ORCA{res}_y{yr:04d}.nc'
#                 if os.path.exists(fn_out):
#                     continue
#     return

# def create_SOFIAMIP_forcing():
#     '''
#     freshwater_forcing_yxxxx.nc: sofwfisf, sofwfrnf, sofwfcal
#     '''
#     for exp in ['antwater']:
#         fn_out = f'{output_path}/SOFIAMIP/FWF_{exp}_ORCA1.nc'
#         if os.path.exists(fn_out):
#             continue
#     return


# if __name__=="__main__":
#     # PRIMAVERA
#     test_file_existence()
#     create_PRIMAVERA_forcing()
#     # create_SOFIAMIP_forcing()