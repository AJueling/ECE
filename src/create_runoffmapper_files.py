""" adapts standard runoffmapper files
- exclude Greenland and/or Antarctica
- for PRIMAVERA-Greenland, Canadian Arctic (Ellesmere-AxelHeiberg-Devon) needs to be separated
- PRIMAVERA simpler runoffmapper file or CMIP6 more complex runoffmapper file
- runoffmapper files valid for both ORCA1 and ORCA025?
"""

import os
import xarray as xr
import datetime

from testing import test_file_existence

output_path = '/ec/res4/hpcperm/nkaj/forcing'


def create_PRIMAVERA_runoffmaps():
    """ removes either only GrIS or GrIS & AIS
    from the PRIMAVERA runoffmaps.nc
    see create_new_runoff_maps.ipynb
    """
    for reg in ['GrIS_only','GrIS_AIS']:
        for res in ['1','025']:
            fn_out = f'{output_path}/runoffmaps_{reg}_ORCA{res}.nc'
            if os.path.exists(fn_out):
                continue
    return


def create_OptimESM_runoffmaps():
    """ removes the Antarctic regions
    from the CMIP6 runoffmaps.nc
    """
    fn_out = f'{output_path}/'
    if os.path.exists(fn_out):
        pass
    return


if __name__=="__main__":
    create_PRIMAVERA_runoffmaps()
    # create_OptimESM_runoffmaps()