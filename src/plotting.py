import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from constants import Gt_per_mm, Sv_per_Gtpy

index = np.arange(1970,2101,1)

def create_secondary_axes(ax):
    """ creates a set of secondary axes for an (m, (rate [Gt/yr], cumulative [Gt])) plot"""
    assert type(ax)==np.ndarray
    ax2 = ax.copy()

    if len(ax.shape)==1: # (1,2)
        (n,) = ax.shape
        assert n==2
        ax2[0] = ax[0].secondary_yaxis('right', functions=(lambda x : x*Sv_per_Gtpy, lambda x : x/Sv_per_Gtpy))
        ax2[1] = ax[1].secondary_yaxis('right', functions=(lambda x : x/Gt_per_mm  , lambda x : x*Gt_per_mm  ))
    elif len(ax.shape)==2: # (m,2)
        (m,n) = ax.shape
        assert n==2
        for i in range(m):
            ax2[i,0] = ax[i,0].secondary_yaxis('right', functions=(lambda x : x*Sv_per_Gtpy, lambda x : x/Sv_per_Gtpy))
            ax2[i,1] = ax[i,1].secondary_yaxis('right', functions=(lambda x : x/Gt_per_mm  , lambda x : x*Gt_per_mm  ))
    else:
        raise ValueError('not the right shape')
    
    return ax2


def MB_units(axs):
    """ adds rate and cumulative units to a given pair of axis objects """
    assert len(axs)==2
    axs[0].text(0, 1.05, '[Gt/yr]'    , ha='right', transform=axs[0].transAxes)
    axs[0].text(1, 1.05, '[Sv]'    , ha='left' , transform=axs[0].transAxes)
    axs[1].text(0, 1.05,r'[$10^3$ Gt]', ha='right', transform=axs[1].transAxes)
    axs[1].text(1, 1.05, '[m]'        , ha='left' , transform=axs[1].transAxes)
    return


def MBint(ts):
    """ converts irregular time series to [1970,2100,1] time series """
    assert type(ts)==pd.core.series.Series
    
    return pd.Series(index=index, data=np.interp(x=index, xp=ts.index, fp=ts.values))


def MB2rate(ts):
    """ takes a (non-regular) timeseries of cumulatively summed and converts it into rate via gradient """
    assert type(ts)==pd.core.series.Series
    return pd.Series(data=np.gradient(MBint(ts)), index=index)


def MB2cumsum(ts):
    """ """
    assert type(ts)==pd.core.series.Series
    return pd.Series(data=np.cumsum(MBint(ts)), index=index)