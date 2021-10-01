import os
import sys
activate_file = os.path.join('/home/users/ajuling/nb-venvs/venv-cmip6-zarr/bin/activate_this.py')
exec(open(activate_file).read(), dict(__file__=activate_file))

import xarray as xr

from iterate import IterateECE

area_255 = xr.open_dataset('../data/T255/areas.nc')['O1t0.srf'].rename({'x_3':'i','y_3':'j'})
area_511 = xr.open_dataset('../data/T511/areas.nc')['Ot25.srf'].rename({'x_3':'i','y_3':'j'})

for da, src, exp, mem in IterateECE(var='thetaoga', cat='jasmin-nc'):
    print(src, exp, mem)
    print(da)
    da.to_netcdf(f'../results/thetao/thetaoga_{src}_{mem}_{exp}.nc')