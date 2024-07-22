import xarray as xr

t264_dir = '/ec/res4/hpcperm/nkaj/cmip6-output/piControl/r1i1p1f1'
data_dir = '/ec/res4/scratch/nkaj/ecearth3/'  # on hpc2020
# data_dir = "../../../Downloads/"  # on local computer

aca = xr.open_dataset('../../data/T255-ORCA1/areacella.nc').areacella
aco = xr.open_dataset('../../data/T255-ORCA1/areacello.nc').areacello.fillna(0)

sofiamip_exps = ['sf02','sf04','sf06','sf08'] # 'sf01'
first_year = {
    'sf01':2260,
    'sf02':2300,
    'sf04':2380,
    'sf06':2460,
    'sf08':2540,
}

file_dict = {
    'tos'   :'grid_T',
    'sos'   :'grid_T',
    'siconc':'icemod',
}

table = {
    'tos'   :'Omon',
    'sos'   :'Omon',
    'so'    :'Omon',
    'thetao':'Omon',
    'hfds'  :'Omon',
    'friver':'Omon',
    'msftyz':'Omon',
    'tas'   :'Amon',
    'siconc':'SImon',
}

# pi = xr.open_mfdataset(f'{t264_dir}/Omon/tos/tos_Omon_EC-Earth3_piControl_r1i1p1f1_gn_*.nc', decode_times=False)


def open_ece_nemo_output(expn, var=None, reassign_time=False):
    """ return xarray dataset of raw EC-Earth NEMO/LIM output
    """
    if var is None:
        grid = 'grid_T'
    else:
        try:
            grid = file_dict[var]
        except:
            raise ValueError('var must be in file_dict')
    
    ds = xr.open_mfdataset(f'{data_dir}/{expn}/output/nemo/0*/{expn}_1m_*0101_*1231_{grid}.nc', decode_times=False)
    if grid=='grid_T':
        rename_dict = {'x':'i',
                       'y':'j',
                       'time_counter':'time',
                       'nav_lat':'latitude',
                       'nav_lon':'longitude'}
    elif grid=='icemod':
        rename_dict = {'x_grid_T':'i',
                       'y_grid_T':'j',
                       'time_counter':'time',
                       'nav_lat_grid_T':'latitude',
                       'nav_lon_grid_T':'longitude'}
    ds = ds.rename(rename_dict)
    ds = ds.drop('time_centered')
    if var is not None:
        ds = ds[var]

    if reassign_time==True:
        assert expn in list(first_year.keys())
        pi_time = (first_year[expn]-2259)*12
        pi_slice = slice(pi_time, pi_time+len(ds.time))
        ds['time'] = pi.time.isel(time=pi_slice).values
    return ds


def open_ece_cmor_output(expn, var, control=False, preprocess=None):
    """ open cmorized model output 
    """
    if expn in ['sf02','sf04','sf06','sf08']:
        cmor_dir = f'/ec/res4/scratch/nkaj/cmorized-results/{expn}/CMIP6/SOFIAMIP/KNMI/EC-Earth3/faf-antwater/r1i1p1f1'
    elif expn==['sdl2']:
        cmor_dir = f'/ec/res4/scratch/nkaj/cmorised-results/{expn}/v001/CMIP6/SOFIAMIP/KNMI/EC-Earth3/faf-antwater/r1i1p1f1'
    else:
        raise ValueError(f'the experiment name expn={expn} is unknown')
    var_dir = f'{table[var]}/{var}/gn/*/'
    if preprocess is None:
        print(f'no preprocessing for {var}')
        ds = xr.open_mfdataset(f'{cmor_dir}/{var_dir}/*.nc', concat_dim='time')
    else:
        print(f'preprocessing {var}')
        ds = xr.open_mfdataset(f'{cmor_dir}/{var_dir}/*.nc', concat_dim='time', preprocess=preprocess)
    return ds[var]


def open_ece_cmip6_output(var, control=False):
    """ open cmip6 output 
    """
    assert var in ['friver']
    cmip6_dir = f'/ec/res4/hpcperm/nkaj/cmip6-output/piControl/r1i1p1f1/Omon/{var}/'
    da = xr.open_mfdataset(f'{cmip6_dir}/*.nc', concat_dim='time')[var]
    return da