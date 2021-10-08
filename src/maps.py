import numpy as np


def rect_polygon(extent):
    assert type(extent)==tuple
    (lonmin,lonmax,latmin,latmax) = extent
    n=50
    xs = [np.linspace(lonmin,lonmax,n), np.linspace(lonmax,lonmax,n),
          np.linspace(lonmax,lonmin,n), np.linspace(lonmin,lonmin,n)]
    ys = [np.linspace(latmin,latmin,n), np.linspace(latmin,latmax,n),
          np.linspace(latmax,latmax,n), np.linspace(latmax,latmin,n)]
    xs = [item for sublist in xs for item in sublist]
    ys = [item for sublist in ys for item in sublist]
    poly_coords = np.swapaxes(np.array([xs, ys]),0,1)
    return poly_coords