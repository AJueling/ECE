import numpy as np
import matplotlib.path as mpath

proj_GrIS = {'central_latitude':72,
             'central_longitude':-42,
             'standard_parallels':[66,78],
             'cutoff':55}

# for round circular plot
# ax.set_boundary(circle, transform=ax.transAxes)
theta = np.linspace(0, 2*np.pi, 100)
center, radius = [0.5, 0.5], 0.5
verts = np.vstack([np.sin(theta), np.cos(theta)]).T
circle = mpath.Path(verts * radius + center)


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