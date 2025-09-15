from importlib import reload
import sys
sys.path.insert(0, '/home/dpe000/Class4_Currents/python')

import xarray as xr
import numpy as np
import time

import find_value_at_point
import read_grid

file='/fs/site5/eccc/mrd/rpnenv/dpe000/SynObs/cutoff_no_us_PR_PF/is_tmp/nrt_global_j3_phy_l3_20211022_20211102.nc'

nav_lon, nav_lat, grid_area = read_grid.read_coord_mesh(hall='hall5')
nx, ny = nav_lon.shape

SLA_DS = xr.open_dataset(file)
ISVAL = np.where(~np.isnan(SLA_DS['sla_filtered']))
SLA_QC = SLA_DS.isel(time=ISVAL)

lat = SLA_QC['latitude']
lon = SLA_QC['longitude']
nobs = len(lat)

dimn='time'

time0 = time.time()
func = lambda lon, lat : find_value_at_point.find_nearest_point_list(lon, lat, nav_lon, nav_lat)
IJPXS = xr.apply_ufunc(func, lon, lat, output_core_dims=[['NPTS']], dask="parallalized")
print('TIME', time.time() - time0)

time0 = time.time()
fund = lambda lon, lat : find_value_at_point.find_nearest_gcircle_list(lon[:100], lat[:100], nav_lon, nav_lat)
IJPYS = xr.apply_ufunc(fund, lon, lat, output_core_dims=[['NPTS']], dask="parallalized")
print('TIME', time.time() - time0)

IJPTS = np.empty((nx,ny), dtype=object)
for ii in range(nx):
   for jj in range(ny):
       IJPTS[ii, jj] = []
       
NNPTS = np.zeros((nx, ny))
       
time0 = time.time()
for iobs in range(nobs): 
    (ii, jj) = IJPXS[iobs,:].data
    NNPTS[ii, jj] = NNPTS[ii, jj] + 1
    IJPTS[ii, jj].append(iobs)

print('TIME', time.time() - time0)

IVAL = np.where(NNPTS > 0 )

NNPTS_OBS = NNPTS[IVAL]
IJPTS_OBS = IJPTS[IVAL]
NOobs = len(NNPTS_OBS)    

time0 = time.time()
SLA_LS = []   
for iobs in range(NOobs):
    ## I ASSUME WE USE THE FILTERED OBS?
    sla_obs = SLA_QC.isel(time=IJPTS_OBS[iobs])['sla_filtered'].values
    if ( len(sla_obs)%2 == 0 ):
        med_obs = np.median(np.sort(sla_obs[1:]))
    else:
        med_obs = np.median(sla_obs)
    imobs = np.where(sla_obs == med_obs)[0][0]
    #print(iobs, imobs)
    SLA_EL = SLA_QC.isel(time=IJPTS_OBS[iobs][imobs])
    SLA_EL['sla_std'] = np.std(sla_obs)
    SLA_LS.append(SLA_QC.isel(time=IJPTS_OBS[iobs][imobs]))

SLA_NU = xr.concat(SLA_LS, dim='time')
print('TIME', time.time() - time0)

SLA_NU.to_netcdf(file.replace('.nc', '_sobb.nc')

