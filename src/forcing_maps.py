"""


"""
import os
import sys
import numpy as np
import xarray as xr

spy = 3600*24*365   # [s yr^-1]

path_grid = '../../ECE/data'
path_fwf = '../../ECE/results/FWF'
path_fwfmap_files = '../'

class FwfMaps:
    """ create freshwater forcing files for EC-Earth experiments 
    output:
    netcdf file with 
    - runoff flux [kg/m^2/s]
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
        assert exp in ['simple', 'complex', 'fafmip', 'eveline']

        if res=='SR':    self.gnem, self.gifs = '1', '255'
        elif res=='HR':  self.gnem, self.gifs = '025', '511'

        fn_areacello = f'{path_grid}/T{self.gifs}-ORCA{self.gnem}/areacello.nc'
        if os.path.exists(fn_areacello):    
            self.areacello = xr.open_dataset(fn_areacello).areacello
        else:
            raise ValueError(f'areacello filename: {fn_areacello} does not exist')

        return

    def create_map(self, year=None):
        """ """
        if self.exp in ['simple', 'complex']:
            assert year is not None
            dA = xr.open_dataset(f'{path_fwf}/AIS_FWF.nc')
            dG = xr.open_dataset(f'{path_fwf}/GrIS_FWF.nc')

            if year=='all':
                years = np.arange(1985,2101)
            elif type(year)==int:
                assert year in dA.year
                years = np.arange(year,year+1)
            else:
                raise ValueError('year must be "all" or integer in range(1985,2101)')

            for yr in years:
                dA_ = dA.sel(year=yr)
                dG_ = dG.sel(year=yr)

                # Antarctica discharge and basal melt
                self.AD = dA_.D*dA.D_seasonality
                self.AB = dA_.B*dA.B_seasonality

                # Greenland discharge and runoff
                self.GD = dG_.D*dG.D_seasonality
                self.GR = dG_.R*dG.R_seasonality

                self.fn_out = f'{path_fwfmap_files}/{self.exp}_forcing/FWF_{self.exp}_ORCA{self.gnem}_y{yr}.nc'

                if self.exp=='simple':
                    self.create_simple_map(year=yr)
                elif self.exp=='complex':
                    self.create_complex_maps(year=yr)
                
                self.test_output()

        elif self.exp=='fafmip':
            self.fn_out = f'../fafmip_forcing/FWF_fafmip_antwater_ORCA{self.gnem}.nc'
            self.create_fafmip_map()
            self.test_output()

        elif self.exp=='eveline':
            pass

    def create_simple_map(self, year):
        """
        - calculates total Greenland and Antarctic surface areas
        - uniformly distributes meltwater in
        # output: creates 19 GB of files for ORCA025
        # zipped them with `tar -cvzf simple_forcing.tar.gz simple_forcing`
        """
        rmf = xr.open_dataarray(f'{path_grid}/T{self.gifs}-ORCA{self.gnem}/runoff_maps_AIS_GrIS_ORCA{self.gnem}.nc')
        Aa = self.areacello.where(rmf==66).sum()
        Ga = self.areacello.where(rmf==1).sum()

        for i, quant in enumerate(['sorunoff_f','socalving_f']):
            if quant=='sorunoff_f': 
                # GrIS: Runoff
                Q = xr.where(rmf==1,self.GR/Ga,0)
                long_name = 'forced runoff flux'
            elif quant=='socalving_f':
                # GrIS: Discharge;  AIS:  Discharge + Basal Melt
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
        """ iceberg melt = Merino pattern, basal melt = vertically distributed along ice shelf fronts, runoff = negligible """
        # iceberg melt

        # basal melt

        FWF.to_netcdf(self.fn_out, unlimited_dims=['time_counter'])
        pass

    def create_fafmip_map(self):
        """
        - climatological file without yearly variations
        - uniformly distributes 0.1 Sv in Antarctica-adjacent NEMO cells
        """
        ac = self.areacello
        ac = ac.fillna(0)
        rmf = xr.open_dataarray(f'{path_grid}/T{self.gifs}-ORCA{self.gnem}/runoff_maps_AIS_GrIS_ORCA{self.gnem}.nc')

        # select coastal cells around Antarctica
        coast = xr.where(ac-ac.shift(i= 1)==ac, ac, 0) + \
                xr.where(ac-ac.shift(i=-1)==ac, ac, 0) + \
                xr.where(ac-ac.shift(j= 1)==ac, ac, 0) + \
                xr.where(ac-ac.shift(j=-1)==ac, ac, 0)
        coast = coast/coast*ac

        ac_fafmip = coast.where(rmf==66).where(coast>0)
        ac_fafmip = xr.where(ac_fafmip,ac_fafmip,0)
        print(f'integrated area of circum-Antarctic cells: {ac_fafmip.sum().values/1e12} million km^2')
        print(f'length circum-Antarctic cells at 100 km width: {ac_fafmip.sum().values/1e5/1e3} km')

        # scale by zonal extent only, i.e. divide by meridional extent of cells
        # distribute 0.1 Sv = 1e5 m^3/s over the area to create a flux density map [m^3/s/m^2] = [m/s]
        latdiff = np.cos(np.deg2rad(ac_fafmip.latitude))
        # latdiff = (ac_fafmip.latitude.shift(j=-1)-ac_fafmip.latitude.shift(j=1))/2
        fwf = 1e5/(ac_fafmip/latdiff).sum(['i','j'])*(ac_fafmip/latdiff)  # rate in  [m^3/s]
        fwf = fwf/ac_fafmip*1e3  # [m^3/s] -> [m/s] -> [kg/m^2/s] @ 1000 kg/m^3

        # create "scaffold" dataset
        time = np.array([np.datetime64(f'1850-{m:02d}-01 12:00:00.000000000') for m in range(1,13)])
        ds = xr.Dataset(coords={'time_counter':time,
                                'i':ac.i.values,
                                'j':ac.j.values,
                                'longitude': (('j','i'), ac.longitude.values),
                                'latitude' : (('j','i'),  ac.latitude.values),
                                },
                        data_vars={'sorunoff_f' :(('time_counter','j','i'),np.zeros((12,len(ac.j),len(ac.i)))),
                                   'socalving_f':(('time_counter','j','i'),np.zeros((12,len(ac.j),len(ac.i))))}
                        )
        # assign 0.1 Sv to runoff forcing, calving remains 0
        for t in range(12):
            ds['sorunoff_f'][t,:,:] = fwf
        
        for i, quant in enumerate(['sorunoff_f','socalving_f']):
            if quant=='sorunoff_f': 
                long_name = 'forced runoff flux'
            elif quant=='socalving_f':
                long_name = 'forced calving flux'
            ds[quant] = ds[quant].fillna(0)
            ds[quant].attrs = {'long_name':long_name, 'units':'kg/m^2/s'}

        ds.to_netcdf(self.fn_out, unlimited_dims=['time_counter'])
        return

    def create_eveline_map(self):
        pass
    
    def test_output(self):
        """ """
        ds = xr.open_dataset(self.fn_out)
        ac = self.areacello
        mean_rate_runoff  = (ac*ds.sorunoff_f ).mean('time_counter').sum(['i','j']).values
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
    exopname = sys.argv[2]
    variable = sys.argv[3]
    print(f'hello world: {variable}')

    FwfMaps(res=resolution, exp=expname).create_simple_map(variable)