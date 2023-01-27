import numpy as np
        
import scipy.interpolate
#import geopy.distance
import rpnpy.librmn.all as rmn

missing=-999

def interpolate_to_point(TM_grid, lon_grid, lat_grid, lon_pt, lat_pt, method='linear', convlon=False):
    if ( isinstance(TM_grid, np.ma.core.MaskedArray) ):
        ivalid = np.where(TM_grid.mask == False)
        TM_flat = TM_grid[ivalid]
        lon_flat = lon_grid[ivalid]%360
        lat_flat = lat_grid[ivalid]
    else:
        TM_flat = TM_grid.flatten()
        lon_flat = lon_grid.flatten()%360
        lat_flat = lat_grid.flatten()

    lon_in = lon_pt%360
    if ( convlon ):
        lon_flat = lon_flat * np.cos(lat_flat * np.pi / 180.0 )
        lon_in = lon_in * np.cos(lat_pt * np.pi / 180.0 )
    TM_pt = scipy.interpolate.griddata( (lon_flat,lat_flat), TM_flat, (lon_in, lat_pt), method=method, fill_value=missing)     
    if ( isinstance(lat_pt, float) ):
        TM_pt = float(TM_pt)
    return TM_pt
    
def ezinterpolate_to_point(lon_pt, lat_pt, TM_grid, file, field, src='A', grid_lat=None, grid_lon=None):
    if ( src == '0' ):
        TM_pt = np.NaN
    if ( ( src == 'A' ) or ( src == 'S' ) or ( src == 'W' ) ):  # Yin Yan grids
        gid_tuple, flds = stfd.read_fstd_gid(file, field)
        (gid, gid0, gid1) = gid_tuple


    if ( ( src == 'A' ) or ( src == 'S' ) ): # Scalar fields.  Direct call to gdllsval should work
        TM_pt  = float( rmn.gdllsval(gid,  lat_pt,  lon_pt,  TM_grid)  )  

    if ( src == 'W' ):  # Vector points.  Need gdllvval -- but also need to work on subgrids.
        if ( isinstance(grid_lon, type(None) ) ):  # Can either pass grid_lon or read from file.
            (__, __, grid_lon, grid_lat) , __, __ =  stfd.read_yy_fstd_multi(file, field)
        UU_grid = TM_grid[0]
        VV_grid = TM_grid[1]

        # Need to break these up into the 2 sub-grids 
        
        ( UU_grid0, UU_grid1 ) = stfd.make_subgrid_fields(UU_grid, gid_tuple)
        ( VV_grid0, VV_grid1 ) = stfd.make_subgrid_fields(VV_grid, gid_tuple)
        ( lon0, lon1 ) = stfd.make_subgrid_fields(grid_lon, gid_tuple)
        ( lat0, lat1 ) = stfd.make_subgrid_fields(grid_lat, gid_tuple)

        (uupt0,vvpt0)  = rmn.gdllvval(gid0,  lat_pt,  lon_pt,  UU_grid0, VV_grid0)
        (uupt1,vvpt1)  = rmn.gdllvval(gid1,  lat_pt,  lon_pt,  UU_grid1, VV_grid1)

        spdpt0 = np.zeros(1).astype(np.float32)
        dirpt0 = np.zeros(1).astype(np.float32)
        spdpt1 = np.zeros(1).astype(np.float32)
        dirpt1 = np.zeros(1).astype(np.float32)
        # This should now give the meteorological direction relative to true north (plus speed).
        ier = rmn.c_gdwdfuv(gid0, spdpt0, dirpt0, uupt0, vvpt0, np.array(lat_pt).astype(np.float32), np.array(lon_pt).astype(np.float32), 1)
        ier = rmn.c_gdwdfuv(gid1, spdpt1, dirpt1, uupt1, vvpt1, np.array(lat_pt).astype(np.float32), np.array(lon_pt).astype(np.float32), 1)
        
        uept0 = -1.0 * spdpt0 * np.sin(dirpt0 * np.pi / 180.0) 
        vnpt0 = -1.0 * spdpt0 * np.cos(dirpt0 * np.pi / 180.0) 
        uept1 = -1.0 * spdpt1 * np.sin(dirpt1 * np.pi / 180.0) 
        vnpt1 = -1.0 * spdpt1 * np.cos(dirpt1 * np.pi / 180.0) 
        print('COMPONENTS', uept0, vnpt0, uept1, vnpt1)
        uept0, vnpt0 = wind_components_from_direction( spdpt0, dirpt0)
        uept1, vnpt1 = wind_components_from_direction( spdpt1, dirpt1)
        print('COMPONENTS', uept0, vnpt0, uept1, vnpt1)

        # Only one of these is correct
        ipt, jpt = find_nearest_point(lon_pt, lat_pt, grid_lon, grid_lat)
        
        in0 = len( np.where( ( lat0 == grid_lat[ipt, jpt] ) & ( lon0 == grid_lon[ipt, jpt] ) )[0]) > 0
        in1 = len( np.where( ( lat1 == grid_lat[ipt, jpt] ) & ( lon1 == grid_lon[ipt, jpt] ) )[0]) > 0
        
        if ( in0 ): 
            uept = uept0
            vnpt = vnpt0
            spdpt = spdpt0
            dirpt = dirpt0
        if ( in1): 
            uept = uept1
            vnpt = vnpt1
            spdpt = spdpt1
            dirpt = dirpt1
        TM_pt = ( float(uept), float(vnpt), float(spdpt), float(dirpt))
    return TM_pt



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
    #print(min_distance)
    val, idx = min((val, idx) for (idx, val) in enumerate(min_distance))
    #print(idx, val)
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
    #print(min_distance)
    val, idx = min((val, idx) for (idx, val) in enumerate(min_distance))
    #print(idx, val)
    ipt=icc[idx]
    jpt=jcc[idx]
    return ipt, jpt

def grid_geopy_distance(lon_pt, lat_pt, lon_grid, lat_grid):
    nx, ny = lon_grid.shape
    distance_grid = np.zeros((nx,ny))
    distance_min = 40e6 # circumference of earth in metres ( minimum distance SHOULD be smaller than this!) 
    for ix in range(nx):
        for iy in range(ny):
            #distance_pt = geopy.distance.distance((lat_grid[ix,iy],lon_grid[ix,iy]),(lat_pt, lon_pt)).m
            distance_pt = 0  #THIS OBVIOUSLY DOES NOT WORK
            if ( distance_pt < distance_min ):
                distance_min = distance_pt
                imin = ix
                jmin = iy
            distance_grid[ix,iy] = distance_pt
    return distance_grid, distance_min, imin, jmin

def find_nearest_geopt(lon_pt, lat_pt, lon_grid, lat_grid):
    ## lon: 0 -> 360
    ## lat: -90 -> 90

    min_distance=[]
    icc = []
    jcc = []
    distance_grid, distance_min, imin, jmin = grid_geopy_distance(lon_pt, lat_pt, lon_grid, lat_grid)
    return imin, jmin

