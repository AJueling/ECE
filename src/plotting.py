import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

from paths import path_data
from constants import Gt_per_mm, Sv_per_Gtpy
from LARMIP import freg

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
        for j, [s,e] in enumerate([[2031,2050], [2046,2065], [2081,2100]]):
            label = [f'SROCC {rcp}',None,None][j]
            value = fac*float(S[S['term']==f'Antarctica {s}-{e}'][rcp].values[0].split()[0])/1e3
            ax_MB[1].plot([s,e], 2*[value], c=c, label=label)
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
    ax_D[0] .plot(t+2000, SD                                   , c=c)
    ax_D[1] .plot(t+2000, (np.cumsum(SD)-np.cumsum(SD)[45])/1e3, c=c, label=r'vdBerk2014 $\Sigma D_{j}$')
    ax_FW[0].plot(t+2000, SD                                   , c=c)
    ax_FW[1].plot(t+2000, (np.cumsum(SD)-np.cumsum(SD)[45])/1e3, c=c, label=r'vdBerk2014 ($=\Sigma D_{j}$)')
    return

def plot_Bintanja2015(ax_FW, c='C9'):
    """ Bintanja 2015"""
    for i, fwf in enumerate([10,20,60,120]):
        if i==0:  label = 'Bintanja2015'
        else :    label = None
        l = ['Bintanja17',None,None,None][i]
        t = np.arange(2006,2046)
        f = fwf*np.arange(1,41)
        ax_FW[0].plot(t, f, c=c, label=l)
        ax_FW[1].plot(t, (np.cumsum(f)-np.cumsum(f)[10])/1e3, c=c, label=label)
    return

def plot_CORDEX(ax_R, c='C7'):
    areacella = xr.open_dataset('../../data/CORDEX-RACMO/areacella_ANT-44_ECMWF-ERAINT_evaluation_r1i1p1_KNMI-RACMO21P_v1_fx.nc').areacella  # [m^2]
    ds = xr.open_mfdataset('../../data/CORDEX-RACMO/mrro_ANT-44_MOHC-HadGEM2-ES_rcp85*.nc')
    
    da = (ds['mrro']*areacella*3600*24*30/1e12).sum(['rlat','rlon'])
    ax_R[0].plot(da.time.dt.year[::12], 
                 da.groupby('time.year').sum(), c=c)
    ax_R[1].plot(da.time.dt.year[::12],
                 da.groupby('time.year').sum().cumsum('year')/1e3, c=c, label=r'CORDEX $R$')
    return

def plot_AIS_scenario(ax_D, ax_B, ax_FW, c='C1'):
    """ 
    (1) scenario = LARMIP + CORDEX
    (2) model-internal FWF = rerouted FW in excess of control
    """
    # total mass LRFs
    da = xr.open_dataarray('../../results/LARMIP/LARMIP_response.nc').sel(model='FETISH_ULB')  # in 10^3 Gt cumulative
    regions = ['EAIS', 'Ross', 'Amundsen', 'Weddell', 'Peninsula']
    cumsum  = np.zeros(len(da.year))
    cumsum0 = np.zeros(len(da.year))
    cumsum1 = np.zeros(len(da.year))
    for i, reg in enumerate(regions):
        if reg=='Weddell':      reg_ = 'rWedd'
        elif reg=='Amundsen':   reg_ = 'rAmun'
        elif reg=='Peninsula':  reg_ = 'rPen'
        else:                   reg_ = 'r'+reg
        cumsum  += da.sel(reg=reg_).values                # D+B
        cumsum0 += da.sel(reg=reg_).values*freg[reg_][0]  # D
        cumsum1 += da.sel(reg=reg_).values*freg[reg_][1]  # B
    ax_D [1].plot(da.year, cumsum0, color='r', label=r'LRF ${D}$'  )  # D
    ax_B [1].plot(da.year, cumsum1, color='r', label=r'LRF ${B}$'  )  # B
    ax_FW[1].plot(da.year, cumsum , color='r', label=r'LRF ${D+B}$')  # D+B
    
    ax_D [0].plot(da.year, np.gradient(cumsum0)*1e3, color='r')  # D
    ax_B [0].plot(da.year, np.gradient(cumsum1)*1e3, color='r')  # B
    ax_FW[0].plot(da.year, np.gradient(cumsum )*1e3, color='r')  # D+B
    
    # model-internal
    da = xr.open_dataarray(f'../../results/SMB/SMB_components_AIS.nc')
    da = da.sel(src='EC-Earth3P-HR').mean('mem', skipna=True)
    da = da.sel(exp='highres-future') \
       - da.sel(exp='control-1950').mean('year')
    fwf = da.sel(var='pr')+da.sel(var='evspsbl')
    ax_FW[0].plot(da.year, fwf, c='r', ls='--')
    ax_FW[1].plot(da.year, np.cumsum(fwf)/1e3, c='r', ls='--', label='model-internal')
    return


""" Greenland components """

def plot_IMBIE_GrIS(ax_MB, ax_SMB, ax_D, c='k'):
    IG = pd.read_csv('../../data/IMBIE2019_Greenland/imbie_dataset_greenland_dynamics-2020_02_28_mass.csv').set_index('Year')
    for k, form in enumerate(['Rate of', 'Cumulative']):
        l1 = [None,'IMBIE2020'][k]
        l2 = [None,'IMBIE2020 anomaly'][k]
        fac = [1,1/1000][k]
        ax_MB[k].set_title(['Rate', 'Cumulative'][k])
        units = ['Gt/yr','Gt'][k]
        ax_MB[k].plot(IG.index, IG[f'{form} ice sheet mass change ({units})']*fac, c='k', label=l1)
        ax_SMB[k].plot(IG.index, IG[f'{form} surface mass balance anomaly ({units})']*fac, c='k', label=l2)
        ax_D[k].plot(IG.index, IG[f'{form} ice dynamics anomaly ({units})']        *fac, c='k', label=l2)
    return

def plot_ISMIP6_GrIS(ax_MB, c='c'):
    ismip = xr.open_dataarray('../../data/v7_CMIP5_pub/ISMIP6_Greenland.nc')
    fac = Gt_per_mm
    ax_MB[1].plot(ismip.time, ismip.mean(dim='model')*fac, c=c, label='ISMIP6')
    ax_MB[1].fill_between(ismip.time, ismip.min(dim='model')*fac, ismip.max(dim='model')*fac, color=c, alpha=.3)
    return

def plot_GRACE(ax_MB):
    GR = pd.read_fwf('../../data/GRACE/greenland_mass_200204_202102.txt', names=['time','mass anomaly [Gt]','uncertainty'] ,skiprows=31)
    ax_MB[0].plot(GR['time'], np.gradient(GR['mass anomaly [Gt]'].rolling(12).sum()), c='grey')
    ax_MB[1].plot(GR['time'], GR['mass anomaly [Gt]']/1000, c='grey', label='GRACE')
    return
    
def plot_SROCC_GrIS(ax_MB, ax_SMB, ax_D):
    S = pd.read_csv('../../data/SROCC/SROCC_T4.4.csv')
    fac = 1000*Gt_per_mm
    for i, rcp in enumerate(['RCP2.6','RCP4.5','RCP8.5']):
        c = ['green', 'orange', 'red'][i]
        smb = -fac*float(S[S['term']=='Greenland SMB'][rcp].values[0].split()[0])
        d = fac*float(S[S['term']=='Greenland DYN'][rcp].values[0].split()[0])
        ax_MB[1] .plot([2081,2100], 2*[smb/1000]     , c=c, label=f'SROCC {rcp}')
        ax_SMB[1].plot([2081,2100], 2*[smb/1000]     , c=c, label=f'SROCC {rcp}')
        ax_D[1]  .plot([2081,2100], 2*[(-smb-d)/1000], c=c, label=f'SROCC {rcp}')
    return

def plot_Bamber2018(ax_D, ax_R, ax_FW, c='C6'):
    """ """
    da = xr.open_dataarray('../../data/Bamber2018/Bamber2018_L15_regions.nc')
    da = da.sel(region='GrIS').groupby(da.time.dt.year).sum()
    ax_D[0] .plot(da.year, da.sel(component='solid_ice')                                   , c=c, label='Bamaber2018')
    ax_R[0] .plot(da.year, da.sel(component='runoff_ice')+da.sel(component='runoff_tundra'), c=c, label='Bamaber2018')
    ax_FW[0].plot(da.year, da.sum('component')                                             , c=c, label='Bamaber2018')
    return

def plot_Mankoff2021(ax_MB, ax_SMB, ax_D, ax_R, ax_FW, c='C7'):
    """ """
    df = pd.read_csv('../../data/Mankoff2021/MB_SMB_D_BMB_ann.csv')
    ax_MB [0].plot(df.time, df['MB'] , c=c, label='Mankoff2021')
    ax_SMB[0].plot(df.time, df['SMB'], c=c, label='Mankoff2021')
    ax_D  [0].plot(df.time, df['D']  , c=c, label='Mankoff2021')
#     ax_R  .plot(df.time, df['R'])
#     ax_FW .plot(df.time, df['D']+df['R']+df['BMB'])
    return

def plot_Golledge2019_GrIS(ax_MB, ax_SMB, ax_D, c='C8'):
    for i, q in enumerate(['MB','SMB','D']):
        ax = [ax_MB, ax_SMB, ax_D][i]
        G = pd.read_csv(f'../../data/Golledge2019/G{q}.csv', sep=';', decimal=',', names=[q], header=None, index_col=0)
        ax[0].plot(G.index, G[q].values, c=c)
        G_ = MB2cumsum(G[q])/1e3
        ax[1].plot(G_.index, G_-G_.iloc[45], c=c, label='Golledge2019')
    return

def plot_Choi2021(ax_MB, ax_SMB, ax_D, c='C4'):
    for i, cmip in enumerate([5,6]):
        for j, comp in enumerate(['MB','SMB','D']):
            
            ax = [ax_MB, ax_SMB, ax_D][j]
            df = pd.read_csv(f'../../data/Choi2021/CMIP{cmip}_{comp}.csv', names=['time','values']).set_index('time')
            if comp=='D':  # we define discharge positive to ocean
                df *= -1
            ts = MBint(df.squeeze())
            ax[1].plot(ts.index, ts-ts.iloc[45], ls=['-','--'][i], color=c, label=f'Choi21 CMIP{cmip} {comp}')
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
    ax_D[0] .plot(t+2000, SD  , c=c)
    ax_R[0] .plot(t+2000, R   , c=c)
    ax_FW[0].plot(t+2000, R+SD, c=c)
    ax_D[1] .plot(t+2000, (np.cumsum(SD  )-np.cumsum(SD  )[45])/1e3, c=c, label=r'vdBerk2014 $\Sigma D_{j}$')
    ax_R[1] .plot(t+2000, (np.cumsum(R   )-np.cumsum(R   )[45])/1e3, c=c, label=r'vdBerk2014 $R$')
    ax_FW[1].plot(t+2000, (np.cumsum(R+SD)-np.cumsum(R+SD)[45])/1e3, c=c, label=r'vdBerk2014 $R + \Sigma D_{j}$')
    ax_D[0].annotate('retreat to land', xy=(2050,550), xytext=(2060,500), arrowprops={'arrowstyle':'->'}, color=c)
    return

def plot_GrIS_scenario(ax_R, ax_FW, c='C1'):
    """ 
    (1) scenario = Lennaerts15 runoff paramaterization + Choi21 discharge
    (2) model-internal FWF = rerouted FW in excess of control
    """
    # runoff parameterization & Choi21 discharge
    RECE = xr.open_dataarray('../../results/Lenaerts2015_regions/param_runoff_EC-Earth3P-HR.nc').isel(time=slice(-131,None))
    ax_R[0].plot(RECE.time, RECE             , c=c)
    ax_R[1].plot(RECE.time, RECE.cumsum('time')/1000, c=c, label='ECE3-HR L15 param.')
    df = -pd.read_csv(f'../../data/Choi2021/CMIP5_D.csv', names=['time','values']).set_index('time')
    ts = MBint(df.squeeze())
    ax_FW[1].plot(RECE.time, RECE.cumsum('time')/1000+(ts-ts.iloc[45]), c='r', label='our scenario')
    
    # model internal
    da = xr.open_dataarray(f'../../results/SMB/SMB_components_GrIS.nc')
    da = da.sel(src='EC-Earth3P-HR').mean('mem', skipna=True)
    da = da.sel(exp='highres-future') \
       - da.sel(exp='control-1950').mean('year')
    fwf = da.sel(var='pr')+da.sel(var='evspsbl')
    ax_FW[0].plot(da.year, fwf, c='r', ls='--')
    ax_FW[1].plot(da.year, np.cumsum(fwf)/1e3, c='r', ls='--', label='model-internal')
    
    return