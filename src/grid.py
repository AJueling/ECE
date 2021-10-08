import numpy as np
import xarray as xr

from constants import R_earth

def spher_surf_element(r, dtheta, lat_N, lat_S):
    """ surface area of element of sphere """
    return r**2 * dtheta * (np.sin(lat_N)-np.sin(lat_S))

def lat_lon_area(da):
    """ generates surface area [m^2] of a rectangular grid with 
    where the input `da` is a dataraay with `lat` and `lon` coordinates
    useful for IFS grids
    """
    assert 'lat' in da.coords and 'lon' in da.coords
    assert len(da.lat.shape)==1 and len(da.lon.shape)==1
    lats, lons = da.lat, da.lon
    nx, ny = len(da.lon), len(da.lat)
    area = xr.DataArray(dims=['lat','lon'],
                        coords=dict(lat=lats, lon=lons))
#     area = xr.zeros_like(da)
    area = area.squeeze()
    for j, l in enumerate(da.lat[:-1]):
        area[j,:] = spher_surf_element(R_earth, 2*np.pi/nx,
                                       da.lat[j+1]*np.pi/180,
                                       da.lat[j]  *np.pi/180)
    area.name = 'grid cell surface area'
    area.attrs['long_name'] = 'grid cell surface area'
    area.attrs['standard_name'] = 'areacella'
    area.attrs['units'] = 'm^2'
    return area


