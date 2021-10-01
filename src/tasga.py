import os
import sys
activate_file = os.path.join('/home/users/ajuling/nb-venvs/venv-cmip6-zarr/bin/activate_this.py')
exec(open(activate_file).read(), dict(__file__=activate_file))

import numpy as np

from iterate import IterateECE

for da, src, exp, mem in IterateECE(var='tas', source='EC-Earth3P-HR', experiment='hist-1950'):
    print(src, exp, mem)
    fn = f'../results/tas/tasga_{src}_{mem}_{exp}.nc'
    gmst = da.weighted(np.cos(np.deg2rad(da.lat))).mean(['lat','lon'])
    gmst.to_netcdf(fn)