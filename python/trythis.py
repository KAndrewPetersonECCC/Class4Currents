from importlib import reload
import sys
sys.path.insert(0, '/home/dpe000/Class4_Currents/python')

import xarray as xr
import numpy as np
import time

import find_value_at_point
import read_grid


SAL_DS = xr.open_dataset('/fs/site5/eccc/mrd/rpnenv/dpe000/maestro_archives/dev-4.0.0.ensembles_v3.0/SAM2/20211027/DIAG/trial/20211027_OLA_IS_SLA.nc')

nav_lon, nav_lat, grid_area = read_grid.read_coord_mesh()

lat = SAL_DS['LATITUDE']
lon = SAL_DS['LATITUDE']
lat_pt = lat[0].values
lon_pt = lon[0].values

ipt, jpt = find_value_at_point.find_nearest_point(lon_pt, lat_pt, nav_lon, nav_lat)
ipc, jpc = find_value_at_point.find_nearest_gcircle(lon_pt, lat_pt, nav_lon, nav_lat)

lon_list = list(lon[0:100].values)
lat_list = list(lat[0:100].values)

lon_llist = lon.values.astype(list)
lat_llist = lat.values.astype(list)

time0 = time.time()
IJPTS = find_value_at_point.find_nearest_point_list(lon_list, lat_list, nav_lon, nav_lat)
print(time.time()-time0)

time0 = time.time()
IJPCS = find_value_at_point.find_nearest_gcircle_list(lon_list, lat_list, nav_lon, nav_lat)
print(time.time()-time0)

time0 = time.time()
IJPTS = find_value_at_point.find_nearest_point_list(lon_list, lat_list, nav_lon, nav_lat, mp=True)
print(time.time()-time0)

time0 = time.time()
IJPCS = find_value_at_point.find_nearest_gcircle_list(lon_list, lat_list, nav_lon, nav_lat, mp=True)
print(time.time()-time0)

func = lambda lon, lat : find_value_at_point.find_nearest_point_list(lon, lat, nav_lon, nav_lat)
for PTS in [1000, 2000, 10000, 100000, 1000000, -1]:
    time0 = time.time()
    IJPXS = xr.apply_ufunc(func, lon[0:PTS], lat[0:PTS], output_core_dims=[['NPTS']], dask="parallalized")
    print(len(lon[0:PTS]), time.time()-time0, IJPXS.shape)
