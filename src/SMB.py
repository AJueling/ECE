# calculating monthly precip, sublimation, and runoff for GrIS and AIS
# need to activate ECE conda environment: `source activate ECE`

import os
import sys
import time
activate_file = os.path.join('/home/users/ajuling/nb-venvs/venv-cmip6-zarr/bin/activate_this.py')
exec(open(activate_file).read(), dict(__file__=activate_file))

# import dask
import numpy as np
import xarray as xr

from iterate import IterateECE, subselect_sysargv
from tqdm.auto import tqdm

if __name__ == '__main__':
    source, experiment, member = subselect_sysargv(sys.argv)

    for var in ['pr','mrro','evspsbl']:
        start = time.time()
        print(var, start)
        try:
            for da, src, exp, mem in IterateECE(var=var,
                                                source=source,
                                                experiment=experiment,
                                                cat='jasmin-nc',
                                                member=member,
                                                verbose=True,
                                               ):  # will cycle through experiments and ensemble members
                if src=='EC-Earth3P-HR':
                    res = 'T511-ORCA025'
                elif src=='EC-Earth3P':
                    res = 'T255-ORCA1'
                area = xr.open_dataset(f'../data/{res}/areacella.nc').areacella
                Amask = xr.open_dataarray(f'../data/{res}/Antarctica_sftlf_{res}.nc')
                Gmask = xr.open_dataarray(f'../data/{res}/Greenland_sftlf_{res}.nc')
                
                start_ = time.time()
                
                # GrIS
                fn = f'../results/{var}/{var}_GrIS_{src}_{mem}_{exp}.nc'
                if os.path.exists(fn)==False:
                    da_ = (da.assign_coords(dict(lat=area.lat,lon=area.lon))*area*Gmask/100).sum(['lon','lat']).compute()                    
                    da_.to_netcdf(fn)
                
                # AIS
                fn = f'../results/{var}/{var}_AIS_{src}_{mem}_{exp}.nc'
                if os.path.exists(fn)==False:
                    da_ = (da.assign_coords(dict(lat=area.lat,lon=area.lon))*area*Amask/100).sum(['lon','lat']).compute()                    
                    da_.to_netcdf(fn)
                    
                print(f'{src:10} {exp:12} {mem:12} {var:7} {time.time()-start:.1f} {time.time()-start_:.1f}')
        except:
            print(f'does not exist:   {src} {var}')