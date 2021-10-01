
# example files

# gws: group workspace
# badc: holds way more data; name of old archive (British Atmopsheric Data Centre)

# ceda_ECE3P = '/badc/cmip6/data/CMIP6/HighResMIP/EC-Earth-Consortium/EC-Earth3P'
ceda_ECE3P = '/gws/nopw/j04/primavera2/stream1/CMIP6/HighResMIP/EC-Earth-Consortium/EC-Earth3P'
# ceda_ECE3P_HR = '/badc/cmip6/data/CMIP6/HighResMIP/EC-Earth-Consortium/EC-Earth3P-HR'
ceda_ECE3P_HR = '/gws/nopw/j04/primavera2/stream1/CMIP6/HighResMIP/EC-Earth-Consortium/EC-Earth3P-HR'

# JASMIN
exf_ORCA1_tos = f'{ceda_ECE3P}/hist-1950/r3i1p2f1/Omon/tos/gn/v20190215/*.nc'
exf_ORCA025_tos = f'{ceda_ECE3P_HR}/hist-1950/r1i1p2f1/Omon/tos/gn/v20181212/*.nc'

exf_T255_mrro = f'{ceda_ECE3P}/hist-1950/r2i1p2f1/Lmon/mrro/gr/v20190812/*.nc'
exf_T511_mrro = f'{ceda_ECE3P_HR}/hist-1950/r1i1p2f1/Lmon/mrro/gr/v20181212/*.nc'
exf_T255_pr = f'{ceda_ECE3P}/hist-1950/r2i1p2f1/Amon/pr/gr/v20190215/*.nc'
exf_T511_pr = f'{ceda_ECE3P_HR}/hist-1950/r2i1p2f1/Amon/pr/gr/v20190625/*.nc'
