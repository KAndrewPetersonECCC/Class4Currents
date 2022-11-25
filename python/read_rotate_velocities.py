import sys
sys.path.insert(0, '/home/dpe000/Class4_Currents/python')
import datetime
import numpy as np

import stfd
import grdfld
import cst_interp

test_file = '/fs/site5/eccc/prod/ops/suites/geps_20220621/e1/gridpt/prog/ens.glboce/2022110700_384_000'
oper_dir = '/fs/site5/eccc/prod/ops/suites/geps_20220621/e1/gridpt/prog/ens.glboce/'
geps_tmp='/fs/site5/eccc/mrd/rpnenv/dpe000/tmpdir/GEPS_TMP'

def UU_ORCA_to_NE(UV):
    hall='hall5'
    fila='/space/'+hall+'/sitestore/eccc/mrd/rpnenv/socn000/env_rhel-8-icelake-64/datafiles/constants/oce/repository/master/CONCEPTS/orca025/grids/orca025grid_new.std'
    LONA, LATA, TH = stfd.read_fstd_var(fila, 'LAAN')
    UO, VO = UV
    UE = UO * np.cos(TH) - VO * np.sin(TH)
    VN = UO * np.sin(TH) + VO * np.cos(TH)
    return [UE, VN]

def UU_list_ORCA_to_NE(UW_list, VW_list):
    neu = len(UW_list)
    nev = len(VW_list)
    if ( neu != nev ):
        print('U and V ensemble not same length', len(UW_list), len(VW_list))
        return None
    else:
        ne = neu
    UE_list = []
    VN_list =[]
    for ie in range(ne):
        UW = UW_list[ie]
        VW = VW_list[ie]
        UV_ROT = UU_ORCA_to_NE((UW, VW))
        UE, VN = UV_ROT
        UE_list.append(UE)
        VN_list.append(VN)
    return UE_list, VN_list

def speed_and_angle(u, v):
    speed =  np.sqrt(u**2+v**2)
    angle = 0*u + 0*v
    nzero = np.where( (u!=0) | (v!=0) ) 
    angle[nzero] = np.arctan2(v[nzero], u[nzero])
    return speed, angle
    
def speed_and_angle_easy(U):
    S = 0.0*U.copy()
    if ( U.ndim == 3 ):  S[:,0,:], S[:,1,:] = speed_and_angle(U[:,0,:], U[:,1,:])
    if ( U.ndim == 4 ):  S[:,0,:,:], S[:,1,:,:] = speed_and_angle(U[:,0,:,:], U[:,1,:,:])
    if ( U.ndim == 5 ):  S[:,0,:,:,:], S[:,1,:,:,:] = speed_and_angle(U[:,0,:,:,:], U[:,1,:,:,:])
    return S

def read_velocities_file(file=test_file, ip1_str='sfc', uvar='UUW', vvar='VVW' ):
    if ( isinstance(ip1_str, int) ): 
        ip1 = ip1_str
    else:
        try:
            ip1 = int(ip1_str)
        except:
            pass
    
    if ( ip1_str == 'sfc' ):
        ip1=15728640  
        ip1a=10979785 

    if ( ip1_str == '15m' ):
        ip1=26314400  
        ip1a=26314400 

    if ( ip1_str == '10ma' ):
        if ( uvar == 'UUW' ):  uvar='UU2W'
        if ( vvar == 'VVW' ):  vvar='VV2W'
        if ( uvar == 'TM' ):   uvar='HCW'
        if ( vvar == 'SALW' ): vvar='SCW'
        ip1=75597472  
        ip1a=ip1

    if ( ip1_str == '20ma' ):
        if ( uvar == 'UUW' ):  uvar='UU2W'
        if ( vvar == 'VVW' ):  vvar='VV2W'
        if ( uvar == 'TM' ):   uvar='HCW'
        if ( vvar == 'SALW' ): vvar='SCW'
        ip1=75697472
        ip1a=ip1

    if ( ip1_str == '30ma' ):
        if ( uvar == 'UUW' ):  uvar='UU2W'
        if ( vvar == 'VVW' ):  vvar='VV2W'
        if ( uvar == 'TM' ):   uvar='HCW'
        if ( vvar == 'SALW' ): vvar='SCW'
        ip1=75797472
        ip1a=ip1

    if ( ip1_str == '10-20ma'  or ip1_str=='10-20mi' ):
        if ( uvar == 'UUW' ):  uvar='UU2W'
        if ( vvar == 'VVW' ):  vvar='VV2W'
        if ( uvar == 'TM' ):   uvar='HCW'
        if ( vvar == 'SALW' ): vvar='SCW'
        ip1=75697472
        ip1a=ip1
        ip1b=75597472


    try:
        print(ip1)
        lon, lat, UW = stfd.read_fstd_var(file, uvar, typvar='P@', ip1=ip1)
        lon, lat, VW = stfd.read_fstd_var(file, vvar, typvar='P@', ip1=ip1)
    except:
        print(ip1a)
        lon, lat, UW = stfd.read_fstd_var(file, uvar, typvar='P@', ip1=ip1a)
        lon, lat, VW = stfd.read_fstd_var(file, vvar, typvar='P@', ip1=ip1a)


    if ( ip1_str=='10-20ma'):  ## VALID FOR DATES AFTER 2021120106
        lon, lat, U1W = stfd.read_fstd_var(file, uvar, typvar='P@', ip1=ip1b)
        lon, lat, V1W = stfd.read_fstd_var(file, vvar, typvar='P@', ip1=ip1b)
        UW = two_depth_average( UW, U1W, D2=20, D1=10, newnorm=True)
        VW = two_depth_average( VW, V1W, D2=20, D1=10, newnorm=True)
        if ( np.max(np.absolute(UW)) > 20.0 or np.max(np.absolute(VW)) > 20.0):  
            print("LOOKS LIKE NORMALIZATION WRONG")
            print("TRY USING ip1_str=10-20mi")
            print("file", file)
    if ( ip1_str=='10-20mi'):  ## VALID FOR DATES BEFORE 2021120106
        lon, lat, U1W = stfd.read_fstd_var(file, uvar, typvar='P@', ip1=ip1b)
        lon, lat, V1W = stfd.read_fstd_var(file, vvar, typvar='P@', ip1=ip1b)
        UW = two_depth_average( UW, U1W, D2=20, D1=10, newnorm=False)
        VW = two_depth_average( VW, V1W, D2=20, D1=10, newnorm=False)
        if ( np.max(np.absolute(UW)) < 0.05 or np.max(np.absolute(VW)) ):  
            print("LOOKS LIKE NORMALIZATION WRONG")
            print("TRY USING ip1_str=10-20ma")
            print("file", file)
        
    return lon, lat, UW, VW
    
def two_depth_average(U2, U1, D2=20, D1=10, newnorm=True):
    if ( newnorm ):
        UA = (D2*U2 - D1*U1) / (D2-D1)
    else:
        UA = ( U2 - U1 ) / (D2-D1)
    return UA 

ref_grid_local='/fs/site6/eccc/mrd/rpnenv/dpe000/GDPS/gdps_Class4/dest_grid.std.25'
def create_interp_files(date, lhr, enslist=list(range(21)), dir=oper_dir, tmpdir=geps_tmp):
    if ( isinstance(date, datetime.datetime) or isinstance(date, datetime.date) ):  datestr=date.strftime('%Y%m%d%H')
    if ( isinstance(date, int) ): datestr=str(date)
    if ( isinstance(date, str) ): datestr=date
    if ( len(datestr) == 8 ): datestr=datestr+'00'  # default to 00Z in necessary

    lhrstr = str(lhr).zfill(3)
    for ens in enslist:
        ensstr = str(ens).zfill(3)  
        base_filename=datestr+'_'+lhrstr+'_'+ensstr
        file=dir+'/'+base_filename
        tmpfile=tmpdir+'/'+base_filename+'.z25'
        cst_interp.cst_interpolation(file, tmpfile, ref_grid=ref_grid_local)
    return

def read_velocities_timestamp(date, lhr, ens, dir=oper_dir, ip1_str='sfc', uvar='UUW', vvar='VVW', suffix=''):
    if ( isinstance(date, datetime.datetime) or isinstance(date, datetime.date) ):  datestr=date.strftime('%Y%m%d%H')
    if ( isinstance(date, int) ): datestr=str(date)
    if ( isinstance(date, str) ): datestr=date
    if ( len(datestr) == 8 ): datestr=datestr+'00'  # default to 00Z in necessary

    lhrstr = str(lhr).zfill(3)
    ensstr = str(ens).zfill(3)    
    file=dir+'/'+datestr+'_'+lhrstr+'_'+ensstr+suffix

    print(file)
    lon, lat, UW, VW = read_velocities_file(file=file, ip1_str=ip1_str, uvar=uvar, vvar=vvar)
    
    return lon, lat, UW, VW

def read_ensemble_velocities(date, lhr, enslist=list(range(21)), dir=oper_dir, ip1_str='sfc', uvar='UUW', vvar='VVW', suffix='' ):
    UW_list = []
    VW_list = []
    for ens in enslist:
        lon, lat, UW, VW = read_velocities_timestamp(date, lhr, ens, dir=dir, ip1_str=ip1_str, uvar=uvar, vvar=vvar, suffix=suffix)
        if ( ens == enslist[0] ):
            lon0=lon
            lat0=lat
        UW_list.append(UW)
        VW_list.append(VW)
    return lon0, lat0, UW_list, VW_list
    
def ens_mean(FLD_list):
    ne = len(FLD_list)
    FLD_ENM = sum(FLD_list) / ne
    return FLD_ENM
     
def ens_mean_variance(FLD_list):
    ne = len(FLD_list)
    FLD_ENM = ens_mean(FLD_list)
    FLD_sqlist = []
    for ie in range(ne):
        FLD_ANOM = FLD_list[ie] - FLD_ENM
        FLD_sqlist.append(np.square(FLD_ANOM))
    FLD_VAR = ens_mean(FLD_sqlist)
    return FLD_ENM, FLD_VAR

def grd_fld(lon, lat, field, ddeg=0.25, method='nearest', central_longitude=0):
    grid_lon, grid_lat, grid_fld=grdfld.grdfld(lon, lat, field, ddeg=ddeg, method=method, central_longitude=central_longitude)
    return grid_lon, grid_lat, grid_fld
    
def grd_list(lon, lat, UV_list):
    ne=len(UV_list)
    grUV_list = []
    for ie in range(ne):
        lon_grid, lat_grid, UV_grid = grdfld.grdfld(lon, lat, UV_list[ie], ddeg=0.25, method='2sweep')        
        grUV_list.append(UV_grid)
        if ( ie == 0 ):
            lon_grid0 = lon_grid
            lat_grid0 = lat_grid
    return lon_grid0, lat_grid0, grUV_list
