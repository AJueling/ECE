import os
import sys
activate_file = os.path.join('/home/users/ajuling/nb-venvs/venv-cmip6-zarr/bin/activate_this.py')
exec(open(activate_file).read(), dict(__file__=activate_file))

import xarray as xr

from iterate import IterateECE

area_255 = xr.open_dataset('../data/T255/areas.nc')['O1t0.srf'].rename({'x_3':'i','y_3':'j'})
area_511 = xr.open_dataset('../data/T511/areas.nc')['Ot25.srf'].rename({'x_3':'i','y_3':'j'})

for da, src, exp, mem in IterateECE(var='siconc', cat='jasmin-nc'):  # , source='EC-Earth3P-HR'
    if src in ['EC-Earth3P','EC-Earth3','EC-Earth3-Veg']:
        area = area_255
    elif src=='EC-Earth3P-HR':
        area = area_511
    else:
        raise ValueError(f'src={src} not known')
    print(src, exp, mem)
    fn_NH = f'../results/siconc/siarea_NH_{src}_{mem}_{exp}.nc'
    fn_SH = f'../results/siconc/siarea_SH_{src}_{mem}_{exp}.nc'
    (da*area).where(da.j>len(da.j)/2).sum(['i','j']).to_netcdf(fn_NH)
    (da*area).where(da.j<len(da.j)/2).sum(['i','j']).to_netcdf(fn_SH)