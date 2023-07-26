import geopy.distance
import numpy as np
import multiprocessing
import itertools
import os

lon_list = [156, 165, -170, -155, -140, -140, -140, -140, -71.76, 137]
lat_list = [  0,   0,    0,    0,   -1,    0,    1,    2, -30.32,  13]
NCPUS = len(os.sched_getaffinity(0))

def find_nearest_geopt(lon_pt, lat_pt, lon_grid, lat_grid):
    ## lon: 0 -> 360
    ## lat: -90 -> 90

    min_distance=[]
    icc = []
    jcc = []
    distance_grid, distance_min, imin, jmin = grid_geopy_distance(lon_pt, lat_pt, lon_grid, lat_grid)
    return imin, jmin

def grid_geopy_distance(lon_pt, lat_pt, lon_grid, lat_grid):
    nx, ny = lon_grid.shape
    distance_grid = np.zeros((nx,ny))
    distance_min = 40e6 # circumference of earth in metres ( minimum distance SHOULD be smaller than this!) 
    for ix in range(nx):
        for iy in range(ny):
            distance_pt = geopy.distance.distance((lat_grid[ix,iy],lon_grid[ix,iy]),(lat_pt, lon_pt)).m
            #distance_pt = 0  #THIS OBVIOUSLY DOES NOT WORK
            if ( distance_pt < distance_min ):
                distance_min = distance_pt
                imin = ix
                jmin = iy
            distance_grid[ix,iy] = distance_pt
    return distance_grid, distance_min, imin, jmin

def find_nearest_point(lon_pt, lat_pt, lon_grid, lat_grid):
    ## lon: 0 -> 360
    ## lat: -90 -> 90
    
    min_distance=[]
    icc = []
    jcc = []
    for cc in [0, -1, 1]:
        lon_cc = lon_pt + cc*360.0 
        distance = np.square(lat_grid-lat_pt)+np.square(lon_grid-lon_cc)
        if ( distance.ndim > 1 ):
            ipt, jpt = np.unravel_index(np.argmin(distance), distance.shape)
            min_distance.append(distance[ipt, jpt])
        else:
            ipt = np.argmin(distance)
            jpt = 1
            min_distance.append(distance[ipt])
        icc.append(ipt)
        jcc.append(jpt)
    ##print('min', min_distance)
    val, idx = min((val, idx) for (idx, val) in enumerate(min_distance))
    ##print(idx, val)
    ipt=icc[idx]
    jpt=jcc[idx]
    return ipt, jpt

def find_nearest_glbpt(lon_pt, lat_pt, lon_grid, lat_grid):
    ## lon: 0 -> 360
    ## lat: -90 -> 90
    
    min_distance=[]
    icc = []
    jcc = []
    for cc in [0, -1, 1]:
        lon_cc = lon_pt + cc*360.0 
        distance = np.square(lat_grid-lat_pt)+np.square(np.cos((np.pi/180.0)*(0.5*lat_grid+0.5*lat_pt)))*np.square(lon_grid-lon_cc)
        ipt, jpt = np.unravel_index(np.argmin(distance), distance.shape)
        min_distance.append(distance[ipt, jpt])
        icc.append(ipt)
        jcc.append(jpt)
    #print(min_distance)
    val, idx = min((val, idx) for (idx, val) in enumerate(min_distance))
    #print(idx, val)
    ipt=icc[idx]
    jpt=jcc[idx]
    return ipt, jpt

## THERE ARE ERRORS IN THIS CALCULATION AT EXTREME SOUTHERN LATITUDES
def find_nearest_glcpt(lon_pt, lat_pt, lon_grid, lat_grid):
    ## lon: 0 -> 360
    ## lat: -90 -> 90
    
    min_distance=[]
    icc = []
    jcc = []
    for cc in [0, -1, 1]:
        lon_cc = lon_pt + cc*360.0 
        converged_lon_cc = lon_cc * np.cos( ( np.pi / 180.0 ) * lat_pt )
        converged_lon_gd = lon_grid * np.cos( ( np.pi / 180.0 ) *  lat_grid )
        distance = np.square(lat_grid-lat_pt)+np.square(converged_lon_gd - converged_lon_cc )
        ipt, jpt = np.unravel_index(np.argmin(distance), distance.shape)
        min_distance.append(distance[ipt, jpt])
        icc.append(ipt)
        jcc.append(jpt)
    ##print('min_distance', min_distance)
    val, idx = min((val, idx) for (idx, val) in enumerate(min_distance))
    ##print(idx, val)
    ipt=icc[idx]
    jpt=jcc[idx]
    return ipt, jpt

def find_nearest_twopt(lon, lat, lon_grid, lat_grid):
    ipt, jpt = find_nearest_point(lon, lat, lon_grid, lat_grid)
    ipg, jpg = find_nearest_geopt(lon, lat, lon_grid, lat_grid)
    return ipt, jpt, ipg, jpg
        
def find_nearest_point_list(lon_list, lat_list, lon_grid, lat_grid, mp=False):
    if ( len(lon_list) != len(lat_list) ):
        print('NEED equal length lists')
        return None
    npts = len(lon_list)
    IJPTS=[]
    if ( mp ):
      nproc = np.max([NCPUS, npts])
      PPOOL = multiprocessing.Pool(nproc)
      IZIP = zip(lon_list, lat_list, itertools.repeat(lon_grid), itertools.repeat(lat_grid))
      IJPTS = PPOOL.starmap(find_nearest_point, IZIP)
      PPOOL.close()
      PPOOL.join()
    else:
      for ipt in range(npts):
        lon_pt=lon_list[ipt]
        lat_pt=lat_list[ipt]
        IJPTS.append( find_nearest_point(lon_pt, lat_pt, lon_grid, lat_grid) )
    return IJPTS

def find_nearest_geopt_list(lon_list, lat_list, lon_grid, lat_grid, mp=True):
    if ( len(lon_list) != len(lat_list) ):
        print('NEED equal length lists')
        return None
    npts = len(lon_list)
    IJPTS=[]
    if ( mp ):
      nproc = np.max([NCPUS, npts])
      PPOOL = multiprocessing.Pool(nproc)
      IZIP = zip(lon_list, lat_list, itertools.repeat(lon_grid), itertools.repeat(lat_grid))
      IJPTS = PPOOL.starmap(find_nearest_geopt, IZIP)
      PPOOL.close()
      PPOOL.join()
    else:
      for ipt in range(npts):
        lon_pt=lon_list[ipt]
        lat_pt=lat_list[ipt]
        IJPTS.append( find_nearest_geopt(lon_pt, lat_pt, lon_grid, lat_grid) )
    return IJPTS

def find_nearest_twopt_list(lon_list, lat_list, lon_grid, lat_grid, mp=True):
    if ( len(lon_list) != len(lat_list) ):
        print('NEED equal length lists')
        return None
    npts = len(lon_list)
    IJPTS=[]
    if ( mp ):
      nproc = np.max([NCPUS, npts])
      PPOOL = multiprocessing.Pool(nproc)
      IZIP = zip(lon_list, lat_list, itertools.repeat(lon_grid), itertools.repeat(lat_grid))
      IJPTS = PPOOL.starmap(find_nearest_twopt, IZIP)
      PPOOL.close()
      PPOOL.join()
    else:
      for ipt in range(npts):
        lon_pt=lon_list[ipt]
        lat_pt=lat_list[ipt]
        IJPTS.append( find_nearest_twopt(lon_pt, lat_pt, lon_grid, lat_grid) )
    return IJPTS
