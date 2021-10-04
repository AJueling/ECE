import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

from paths import path_data
from constants import Gt_per_mm, Sv_per_Gtpy

index = np.arange(1970,2101,1)

""" figure elements """

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


""" Antarctic components """

def plot_IMBIE_AIS(ax_MB, ax_SMB, c='k'):
    # Shepherd2018: IMBIE AIS, MB only, c='k'
    A = pd.read_csv(f'{path_data}/IMBIE2018_Antarctica/imbie_dataset-2018_07_23_AIS.csv').set_index('Year')
    I_MB = A['Cumulative ice mass change (Gt)']
    line, = ax_MB[1].plot(I_MB.index, I_MB/1e3, c='k',label='IMBIE Antarctica')
    ax_MB[0].plot(I_MB.index, np.gradient(I_MB, 1/12), c='k')
    
    E = pd.read_csv(f'{path_data}/IMBIE2018_Antarctica/ExtData1_AIS_SMB_1992-2017.csv')
    for i in range(5):
        ax_SMB[0].plot([1992,2017],2*[E.loc[i]['SMB (Gt/yr)']], c='grey', label=E.loc[i]['Model'])
    return

def plot_SROCC_AIS(ax_MB):
    """ SROCC """
    S = pd.read_csv(f'{path_data}/SROCC/SROCC_T4.4.csv')
    fac = -1000*Gt_per_mm  # given in meters SLR
    for i, rcp in enumerate(['RCP2.6','RCP4.5','RCP8.5']):
        c = ['green', 'orange', 'red'][i]
        for [s,e] in [[2031,2050], [2046,2065], [2081,2100]]:
            value = fac*float(S[S['term']==f'Antarctica {s}-{e}'][rcp].values[0].split()[0])/1e3
            ax_MB[1].plot([s,e], 2*[value], c=c, label=f'SROCC {rcp}')
    return

def plot_Golledge2019_AIS(ax_MB, ax_SMB, ax_D, c='C8'):
    # Golledge2019: MB, SMB, D, c='C8'
    for i, q in enumerate(['MB','SMB','D']):
        ax = [ax_MB, ax_SMB, ax_D][i]
        G = pd.read_csv(f'{path_data}/Golledge2019/A{q}.csv', sep=';', decimal=',', names=[q], header=None, index_col=0)
        ax[0].plot(G.index, G[q].values, c='C8', label='Golledge2019')
        ts = MB2cumsum(G[q])/1000
        Gl, = ax[1].plot(ts.index, ts-ts.iloc[45], c='C8', label='Golledge2019')
    return

def plot_vdBerk2014_AIS(ax_D, ax_FW, c='C0'):
    t = np.arange(1970,2100,1)-2000
    D1 = 237+np.array([237*(1.013**t_-1) if t_<=30 else 177*7 for t_ in t])
    D2 = 388+388*np.array([1.013**t_-1 if t_<=30 else (1.013**30-1)*np.exp(0.0385*(t_-30)) for t_ in t])
    D3 = 107+107*np.array([1.013**t_-1 if t_<=30 else (1.013**30-1)*np.exp(0.0375*(t_-30)) for t_ in t])
    SD = D1+D2+D3
    ax_D[0] .plot(t+2000, SD                                   , c=c, label=r'vdBerk14')
    ax_D[1] .plot(t+2000, (np.cumsum(SD)-np.cumsum(SD)[45])/1e3, c=c, label=r'vdB14 $\Sigma D_{j}$')
    ax_FW[0].plot(t+2000, SD                                   , c=c, label=r'vdB14 $\Sigma D_{j}$')
    ax_FW[1].plot(t+2000, (np.cumsum(SD)-np.cumsum(SD)[45])/1e3, c=c, label=r'vdB14 $\Sigma D_{j}$')
    return

def plot_Bintanja2017(ax_FW, c='C9'):
    """ Bintanja 2017"""
    for i, fwf in enumerate([10,20,60,120]):
        l = ['Bintanja17',None,None,None][i]
        t = np.arange(2006,2046)
        f = fwf*np.arange(1,41)
        ax_FW[0].plot(t, f, c=c, label=l)
        ax_FW[1].plot(t, (np.cumsum(f)-np.cumsum(f)[10])/1e3, c=c)
    return

def plot_AIS_scenario(ax_FW, c='C1'):
    """ 
    (1) scenario = LARMIP + CORDEX
    (2) model-internal FWF = rerouted FW in excess of control
    """
    da = xr.open_dataarray(f'../../results/SMB/SMB_components_AIS.nc')
    da = da.sel(src='EC-Earth3P-HR').mean('mem', skipna=True)
    da = da.sel(exp='highres-future') \
       - da.sel(exp='control-1950').mean('year')
    fwf = da.sel(var='pr')+da.sel(var='evspsbl')
    ax_FW[0].plot(da.year, 10*fwf, c='r', ls='--', label='model-internal')
    ax_FW[1].plot(da.year, 10*np.cumsum(fwf)/1e3, c='r', ls='--')
    return


""" Greenland components """

def plot_IMBIE_GrIS(ax_MB, ax_SMB, ax_D, c='k'):
    IG = pd.read_csv('../../data/IMBIE2019_Greenland/imbie_dataset_greenland_dynamics-2020_02_28_mass.csv').set_index('Year')
    for k, form in enumerate(['Rate of', 'Cumulative']):
        fac = [1,1/1000][k]
        ax_MB[k].set_title(['Rate', 'Cumulative'][k])
        units = ['Gt/yr','Gt'][k]
        I20, = ax_MB[k].plot(IG.index, IG[f'{form} ice sheet mass change ({units})']*fac, c='k', label='IMBIE2020')
        ax_SMB[k].plot(IG.index, IG[f'{form} surface mass balance anomaly ({units})']*fac, c='k', label='IMBIE2020 anomaly')
        ax_D[k].plot(IG.index, IG[f'{form} ice dynamics anomaly ({units})']        *fac, c='k', label='IMBIE2020 anomaly')
    return

def plot_ISMIP6_GrIS(ax_MB, c='c'):
    ismip = xr.open_dataarray('../../data/v7_CMIP5_pub/ISMIP6_Greenland.nc')
    fac = Gt_per_mm
    ISMIP, = ax_MB[1].plot(ismip.time, ismip.mean(dim='model')*fac, c=c, label='ISMIP6')
    ax_MB[1].fill_between(ismip.time, ismip.min(dim='model')*fac, ismip.max(dim='model')*fac, color=c, alpha=.3)
    return

def plot_GRACE(ax_MB):
    GR = pd.read_fwf('../../data/GRACE/greenland_mass_200204_202102.txt', names=['time','mass anomaly [Gt]','uncertainty'] ,skiprows=31)
    GRACE, = ax_MB[1].plot(GR['time'], GR['mass anomaly [Gt]']/1000, c='grey', label='GRACE')
    return
    
def plot_SROCC_GrIS(ax_MB, ax_SMB, ax_D):
    S = pd.read_csv('../../data/SROCC/SROCC_T4.4.csv')
    fac = 1000*Gt_per_mm
    for i, rcp in enumerate(['RCP2.6','RCP4.5','RCP8.5']):
        c = ['green', 'orange', 'red'][i]
        smb = -fac*float(S[S['term']=='Greenland SMB'][rcp].values[0].split()[0])
        d = fac*float(S[S['term']=='Greenland DYN'][rcp].values[0].split()[0])
        ax_MB[1].plot([2081,2100], 2*[smb/1000]     , c=c, label=f'SROCC {rcp}')
        ax_SMB[1].plot([2081,2100], 2*[smb/1000]     , c=c, label=f'SROCC {rcp}')
        ax_D[1].plot([2081,2100], 2*[(-smb-d)/1000], c=c, label=f'SROCC {rcp}')
    return

def plot_Golledge2019_GrIS(ax_MB, ax_SMB, ax_D, c='C8'):
    for i, q in enumerate(['MB','SMB','D']):
        ax = [ax_MB, ax_SMB, ax_D][i]
        G = pd.read_csv(f'../../data/Golledge2019/G{q}.csv', sep=';', decimal=',', names=[q], header=None, index_col=0)
        G19, = ax[0].plot(G.index, G[q].values, c=c, label='Golledge2019')
        G_ = MB2cumsum(G[q])/1e3
        ax[1].plot(G_.index, G_-G_.iloc[45], c=c)
    return

def plot_Choi2021(ax_MB, ax_SMB, ax_D, c='C4'):
    for i, cmip in enumerate([5,6]):
        for j, comp in enumerate(['MB','SMB','D']):
            
            ax = [ax_MB, ax_SMB, ax_D][j]
            df = pd.read_csv(f'../../data/Choi2021/CMIP{cmip}_{comp}.csv', names=['time','values']).set_index('time')
            if comp=='D':  # we define discharge positive to ocean
                df *= -1
            ts = MBint(df.squeeze())
            C21, = ax[1].plot(ts.index, ts-ts.iloc[45], ls=['-','--'][i], color=c, label=f'Choi21 CMIP{cmip} {comp}')
            ts = MB2rate(df.squeeze()*1000)
            ax[0].plot(ts.index, ts, ls=['-','--'][i], color=c)
    return

def plot_vdBerk2014_GrIS(ax_D, ax_R, ax_FW, c='C0'):
    t = np.arange(1970,2100,1)-2000
    R = (0.013+(2.96e-4*t))/Sv_per_Gtpy
    D1 = 69.5*(3/104*(t+4)+1)
    D2 = 81.7*np.array([((t_+4)/54+1) if t_<=50 else 1 for t_ in t])
    D3 = 36 + (4/100*115*t)
    SD = D1+D2+D3
    ax_D[0].plot(t+2000, SD  , c='C0', label=r'vdB14 $\Sigma D_{j}$')
    ax_R[0].plot(t+2000, R                                        , c=c, label=r'vdB14 $R$')
    ax_FW[0].plot(t+2000, R+SD                                     , c=c, label=r'vdB14 $R + \Sigma D_{j}$')
    ax_D[1].plot(t+2000, (np.cumsum(SD  )-np.cumsum(SD  )[45])/1e3, c=c, label=r'vdB14 $\Sigma D_{j}$')
    ax_R[1].plot(t+2000, (np.cumsum(R   )-np.cumsum(R   )[45])/1e3, c=c, label=r'vdB14 $R$')
    ax_FW[1].plot(t+2000, (np.cumsum(R+SD)-np.cumsum(R+SD)[45])/1e3, c=c, label=r'vdB14 $R + \Sigma D_{j}$')
    ax_D[0].annotate('retreat to land', xy=(2050,550), xytext=(2060,500), arrowprops={'arrowstyle':'->'}, color=c)
    return

def plot_GrIS_scenario(ax_R, ax_FW, c='C1'):
    """ 
    (1) scenario = Lennaerts15 runoff paramaterization + Choi21 discharge
    (2) model-internal FWF = rerouted FW in excess of control
    """
    RECE = xr.open_dataarray('../../results/Lenaerts2015_regions/param_runoff_EC-Earth3P-HR.nc').isel(time=slice(-131,None))
    L15, = ax_R[0].plot(RECE.time, RECE             , c=c, label='ECE3-HR L15 param.')
    ax_R[1].plot(RECE.time, RECE.cumsum('time')/1000, c=c, label='ECE3-HR L15 param.')
    
    df = -pd.read_csv(f'../../data/Choi2021/CMIP5_D.csv', names=['time','values']).set_index('time')
    ts = MBint(df.squeeze())
    ns, = ax_FW[1].plot(RECE.time, RECE.cumsum('time')/1000+(ts-ts.iloc[45]), c='r', label='our scenario')
    da = xr.open_dataarray(f'../../results/SMB/SMB_components_GrIS.nc')
    da = da.sel(src='EC-Earth3P-HR').mean('mem', skipna=True)
    da = da.sel(exp='highres-future') \
       - da.sel(exp='control-1950').mean('year')
    fwf = da.sel(var='pr')+da.sel(var='evspsbl')
    ax_FW[0].plot(da.year, 10*fwf, c='r', ls='--', label='model-internal')
    ax_FW[1].plot(da.year, 10*np.cumsum(fwf)/1e3, c='r', ls='--')
    
    return