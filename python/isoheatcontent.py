import numpy as np
import socket
import time
import sys
        
import stfd
import read_grid
import datetime

import tmpcp
import find_hall

geops_dir='/fs/homeu2/eccc/mrd/ords/rpnenv/dpe000/GEOPS/GEOPS/python'
sys.path.insert(0, geops_dir)
import read_class4
del sys.path[0]

hall = find_hall.find_hall()

rcp = 3991.86795711963
rau0 = 1026.0
rau0_rcp =    rcp * rau0
cmap='RdYlBu_r'
KCONV=273.16

def f_rcp():  return rcp
def f_rau0(): return rau0
def f_rau0_rcp(): return rau0_rcp

def Temperature_in_C(T):
    T_C = T - KCONV
    return T_C
    
def Temperature_in_K(T):
    T_K = T + KCONV
    return T_K
        
def isotherm_from_below(T_3D, mask, lev_1d, bottom, Tlevel=20):
    # T_3D/Tlevel MUST BOTH be in degrees Celsius OR Kelvin
    # should depth be depth below surface?
    #  add ssh to e3t[0,:,:] for non-VVL case
    
    nz, nx, ny = T_3D.shape  # note the odd shape of standard file 3D files.
    D_iso = bottom
    # search upwards from bottom to 2nd most top layer
    
    is_3dlev = False
    if ( isinstance(lev_1d, np.ndarray) ):
        if ( lev_1d.ndim == 3 ):
            is_3dlev = True

    for iz in sorted(range(1,nz), reverse=True):
       Tk = T_3D[iz,:,:]
       Mk = mask[iz,:,:]
       dT = T_3D[iz-1,:,:] - T_3D[iz,:,:]   ## assumed positive
       if ( not is_3dlev ):
           dz = lev_1d[iz] - lev_1d[iz-1]
       else:
           dz = lev_1d[iz+1,:,:] - lev_1d[iz,:,:]
       dTdz = dT / dz
       iocean = np.where( Mk > 0 )
       D_ocean = D_iso[iocean]
       ibelow = np.where( Tk[iocean] <= Tlevel ) 
       
       if ( not is_3dlev ):
           D_ocean[ibelow] = lev_1d[iz] - (Tlevel - Tk[iocean][ibelow])/ dTdz[iocean][ibelow]
       else:
           D_ocean[ibelow] = lev_1d[iz,:,:] - (Tlevel - Tk[iocean][ibelow])/ dTdz[iocean][ibelow]
       D_iso[iocean] = D_ocean
       
    # Finally, if top layer is less than isotherm, or is land, isotherm depth is zero    
    ibelow = np.where(T_3D[0,:,:] <= Tlevel )
    D_iso[ibelow] = 0.0
    # Finally, if top layer is land, isotherm depth is zero    
    iland = np.where(mask[0,:,:] == 0 )
    D_iso[iland] = 0.0
    #print('Final below', np.min(D_iso), np.max(D_iso))
    return D_iso
    
def isotherm_from_above(T_3D, mask, lev_1d, bottom, Tlevel=20):
    # T_3D MUST be in degrees Celsius (or Tlevel must be changed)
    # should depth be depth below surface?
    #  add ssh to e3t[0,:,:] for non-VVL case
    
    nz, nx, ny = T_3D.shape  # note the odd shape of standard file 3D files.
    #print('nz, nx, ny ', nz, nx, ny)
    D_iso = np.zeros((nx,ny))
    D_thr = np.zeros((nx,ny))  ## in order not to get lower layer in inversion, need to know when threshold first meet.
    
    is_3dlev = False
    if ( isinstance(lev_1d, np.ndarray) ):
        if ( lev_1d.ndim == 3 ):
            is_3dlev = True

    for iz in range(nz-1):
       Tk = T_3D[iz,:,:]
       Mk = mask[iz+1,:,:]
       dT = T_3D[iz,:,:] - T_3D[iz+1,:,:]   ## assumed positive
       if ( not is_3dlev ):
           dz = lev_1d[iz+1] - lev_1d[iz]
       else:
           dz = lev_1d[iz+1,:,:] - lev_1d[iz,:,:]
       dTdz = dT / dz
       iabove = np.where( ( Tk > Tlevel ) & (D_thr == 0) )
       ibelow = np.where( Tk < Tlevel )
       D_above = D_iso[iabove]
       D_thr[ibelow] = 1

       iocean = np.where( Mk[iabove] > 0 )
       ibottom = np.where( Mk[iabove] == 0 )

       if ( not is_3dlev ):
           D_above[iocean] = lev_1d[iz] + (Tk[iabove][iocean] - Tlevel ) / dTdz[iabove][iocean] 
       else:
           D_above[iocean] = lev_1d[iz,:,:] + (Tk[iabove][iocean] - Tlevel ) / dTdz[iabove][iocean] 
       D_above[ibottom] = bottom[iabove][ibottom]
       
       D_iso[iabove] = D_above
       #if (len(D_iso[iabove][iocean])>0): print('ocean level =', iz, np.max(D_iso[iabove][iocean]))
       #if (len(D_iso[iabove][ibottom])>0): print('bottom level =', iz, np.max(D_iso[iabove][ibottom]))

    # Finally, if bottom layer is greater than isotherm, isotherm depth is bottom    
    Tk = T_3D[nz-1,:,:]
    Mk = mask[iz,:,:]
    iabove = np.where( Tk > Tlevel )
    D_above = D_iso[iabove]
    iocean = np.where( Mk[iabove] > 0 )
    D_above[iocean] = bottom[iabove][iocean]
    D_iso[iabove] = D_above
    #print('Final above', np.min(D_iso), np.max(D_iso))
    
    return D_iso

def isotherm(T_3D, mask, lev_1d, bottom, Tlevel=20, from_above=True):
    # from below is faster
    if ( from_above ): 
        D_iso = isotherm_from_above(T_3D, mask, lev_1d, bottom, Tlevel=Tlevel)
    else:
        D_iso = isotherm_from_below(T_3D, mask, lev_1d, bottom, Tlevel=Tlevel)
    return D_iso

def isotherm_class4_slow(OBS, DEPTH,Tlevel=20, from_above=True):
    # CAN USE ISOTHERM ROUTINE.  Not Necessarily efficient: Reverse array order
    if ( DEPTH.ndim > 1 ):
        nobs, ndepths = DEPTH.shape
    else: 
        nobs=1
        ndepths = len(DEPTH)
        print(nobs, ndepths)
        DEPTH = np.reshape(DEPTH, (nobs, ndepths) )

    if ( OBS.ndim == 3 ): 
        T_OBS = OBS[:,0,:]
    elif (OBS.ndim == 1 ):
        T_OBS = np.reshape(OBS, (nobs, ndepths))
    else:
        T_OBS = OBS[:,:]
        
    MASK = 1 - DEPTH.mask     
    BOTM = np.max(DEPTH,1)
    D_iso = []
    for iobs in range(nobs):
        #print('iobs = ', iobs, nobs)
        T_3D = np.reshape(T_OBS[iobs,:], (ndepths, 1, 1))
        
        M_3D = np.reshape(MASK[iobs,:], (ndepths, 1, 1))
        B_3D = np.reshape(BOTM[iobs], (1,1))
        L_1D = DEPTH[iobs,:]
        D_tmp = isotherm(T_3D, M_3D, L_1D, B_3D, Tlevel=Tlevel, from_above=from_above)
        #print('iobs = ', iobs, nobs, D_tmp.shape)
        D_iso.append(D_tmp[0,0])
    D_iso = np.array(D_iso)
    return D_iso

def rearrange_class4_fld( FLD ):
    nobs, ndepths = FLD.shape  
    F_3D = np.reshape( np.transpose(FLD), (ndepths, nobs, 1) )
    return F_3D
        
def mk1D_class4_fld( FLD ):
    nobs, none = FLD.shape  
    F_1D = np.reshape(FLD, nobs)
    return F_1D
        
def make_mask(FLD):
    MASK = (1 - FLD.mask.astype(int))
    return MASK

def isotherm_class4(OBS, DEPTH, Tlevel=20, from_above=True):
    # CAN USE ISOTHERM ROUTINE.  Not Necessarily efficient: Reverse array order
    nobs, ndepths = DEPTH.shape  
    if ( OBS.ndim == 3 ): 
        T_OBS = OBS[:,0,:]
    else:
        T_OBS = OBS[:,:]
        
    MASK = read_class4.make_mask(T_OBS) 
    BOTM = np.max(DEPTH,1)
    
    T_3D = rearrange_class4_fld(T_OBS)
    M_3D = rearrange_class4_fld(MASK)
    L_3D = rearrange_class4_fld(DEPTH)
    B_3D = np.reshape( BOTM, (nobs, 1) )
    D_iso = isotherm(T_3D, M_3D, L_3D, B_3D, Tlevel=Tlevel, from_above=from_above)
    D_iso = mk1D_class4_fld(D_iso)
    return D_iso
    
def depth_integral(F_3D, e3t, mask, depth=10):
    # NOTE: E3t[0] = E3T[0]+SSH for non-VVL cases.

    nz, nx, ny = F_3D.shape
    zthick = np.zeros((nx, ny))
    depth_fld = np.zeros((nx, ny))
    if ( ( type(depth) is int ) or ( type(depth) is float ) ): 
        depth_2d = depth*np.ones((nx,ny))
    else:
        depth_2d = depth
    
    for iz in range(nz):
        dzthick = np.minimum( np.maximum(depth_2d-zthick, 0) , e3t[iz, :, :]*mask[iz,:,:] )
        zthick = zthick + dzthick*mask[iz,:,:]
        depth_fld = depth_fld + F_3D[iz,:,:] * dzthick* mask[iz,:,:]

    return depth_fld, zthick

def heat_content(T_3D, e3t, mask, depth=10):
    # NOTE T_3d needs to be in K for proper heat content
    depth_fld, zthick = depth_integral(T_3D, e3t, mask, depth=depth)
    htc = rau0_rcp * depth_fld
    return htc

def heat_content_class4(OBS, DEPTH, depth=10):
    if ( OBS.ndim == 3 ): 
        T_OBS = OBS[:,0,:]
    else:
        T_OBS = OBS[:,:]
    nobs, ndepths = T_OBS.shape
    MASK = read_class4.make_mask(T_OBS)
    E3T, zedge = read_class4.calc_e3t_inc(DEPTH)
    #print(nobs, ndepths)    
    T_3D = rearrange_class4_fld(T_OBS.data)
    E_3D = rearrange_class4_fld(E3T)
    M_3D = rearrange_class4_fld(MASK)
    #print(T_3D.shape, E_3D.shape, M_3D.shape)
    htc = heat_content(T_3D, E_3D, M_3D, depth=depth)
    htc = mk1D_class4_fld(htc)
    #print(htc.shape)
    return htc   
    
def heat_content_diff(T_3D, e3t, mask, depth=[10, 100]):
    # NOTE T_3d needs to be in K for proper heat content
    depth_fld_1, zthick = depth_integral(T_3D, e3t, mask, depth=depth[0])
    depth_fld_2, zthick = depth_integral(T_3D, e3t, mask, depth=depth[1])
    htc = rau0_rcp * (depth_fld_2 - depth_fld_1)
    return htc
    
def heat_to_degC(H, depth=10, anomaly=False):
    TH = H / rau0_rcp / depth
    if (not anomaly):
        TH = Temperature_in_C(TH)
    return TH
    
def salt_content(S_3D, e3t, mask, depth=10):
    depth_fld, zthick = depth_integral(S_3D, e3t, mask, depth=depth)
    stc = rau0 * depth_fld
    return stc
    
def salt_content_class4(OBS, DEPTH, depth=10):
    if ( OBS.ndim == 3 ): 
        S_OBS = OBS[:,1,:]
    else:
        S_OBS = OBS[:,:]
    nobs, ndepths = DEPTH.shape
    MASK = make_mask(S_OBS)
    E3T, zedge = read_class4.calc_e3t_inc(DEPTH)
    #print(nobs, ndepths)    
    S_3D = rearrange_class4_fld(S_OBS.data)
    E_3D = rearrange_class4_fld(E3T)
    M_3D = rearrange_class4_fld(MASK)
    stc = salt_content(S_3D, E_3D, M_3D, depth=depth)
    stc = mk1D_class4_fld(stc)
    return stc   

def salt_content_diff(S_3D, e3t, mask, depth=[10, 100]):
    depth_fld_1, zthick = depth_integral(S_3D, e3t, mask, depth=depth[0])
    depth_fld_2, zthick = depth_integral(S_3D, e3t, mask, depth=depth[1])
    stc = rau0 * (depth_fld_2 - depth_fld_1)
    return stc

def salt_to_PSU(S, depth=10):
    SH = S / rau0 / depth
    return SH
    
def depth_mean_velocity(U_3D, e3u, mask, depth=30):
    nz, nx, ny = U_3D.shape
    depth_fld, zthick = depth_integral(U_3D, e3u, mask, depth=depth)
    inonzero = np.where(zthick > 0.0 )
    utc = np.zeros((nx, ny))
    utc[inonzero] = depth_fld[inonzero] / zthick[inonzero]
    return utc
    
def hurricane_heat_potential(T_3D, lev, e3t, mask, T_threshold=26):
    nz, nx, ny = T_3D.shape
    T_C = Temperature_in_C(T_3D)
    bottom = read_grid.bottom_depth_from_e3t(e3t, mask)
    depth_iso = isotherm(T_C, mask, lev, bottom, Tlevel=T_threshold, from_above=True)
    mask_iso = np.zeros((nx, ny))
    ispositive = np.where( depth_iso > 0 )
    mask_iso[ispositive] = 1
    mask_tot = mask * mask_iso
    T_wrt_threshold = T_C - T_threshold
    depth_fld, zthick = depth_integral(T_wrt_threshold, e3t, mask_tot, depth=depth_iso)
    hhp = depth_fld * rau0_rcp
    return hhp

def hurricane_heat_potential_class4(OBS, DEPTH, T_threshold=26):
    if ( OBS.ndim == 3 ): 
        T_OBS = OBS[:,0,:]
    else:
        T_OBS = OBS[:,:]
    nobs, ndepths = T_OBS.shape
    depth_iso = isotherm_class4(T_OBS, DEPTH, Tlevel=T_threshold)
    mask = read_class4.make_mask(T_OBS)
    E3T, zedge = read_class4.calc_e3t_inc(DEPTH)
    mask_iso = np.zeros(nobs)
    ispositive = np.where( depth_iso > 0 )
    mask_iso[ispositive] = 1
    mask_tot = np.transpose(np.transpose(mask) * mask_iso)
    T_wrt_threshold = T_OBS.data - T_threshold
    T_3D = rearrange_class4_fld(T_wrt_threshold)
    E3T.shape
    E_3D = rearrange_class4_fld(E3T)
    M_3D = rearrange_class4_fld(mask_tot)
    D_2D = np.reshape(depth_iso, (nobs, 1))
    depth_fld, zthick = depth_integral(T_3D, E_3D, M_3D, depth=D_2D)
    hhp = depth_fld * rau0_rcp
    hhp = mk1D_class4_fld(hhp)
    return hhp



def calculate_from_file(file, var='H10', hall=hall):
    if ( len(file.split(":")) == 2 ): file=tmpcp.tmpcp(file)
    # Need U/V grid nav_lon/lat too, or instead.
    # Ditto umask and vmask
    for ii in range(len(var)):
        print(var[ii])
    bad=True
    if ( var == 'TM' ): 
        fvar = 'TM' 
        evar = 'e3t_0'
        mvar = 'tmask'
        gvar = 'T'
        bad = False
    if ( var[0] == 'H' ): 
        print(var[0])
        fvar = 'TM' 
        evar = 'e3t_0'
        mvar = 'tmask'
        gvar = 'T'
        bad = False
    try: 
        var2 = var[2]
    except:
        var2 = 'X'
    if ( var2 == 'D' ): 
        #print(var2)
        fvar = 'TM' 
        evar = 'e3t_0'
        mvar = 'tmask'
        gvar = 'T'
        bad = False
    if ( var[0] == 'S' ): 
        print(var[0])
        fvar = 'SALW'
        evar = 'e3t_0'
        mvar = 'tmask'
        gvar = 'T'
        bad = False
    if ( var[0] == 'U' ): 
        print(var[0])
        fvar = 'UUW'
        #evar = 'e3u_0'
        #mvar = 'umask'
        #gvar = 'U'
        # In standard file, everything is on T-grid
        evar = 'e3t_0'
        mvar = 'tmask'
        gvar = 'T'
        bad = False
    if ( var[0] == 'V' ): 
        print(var[0])
        fvar = 'VVW'
        #evar = 'e3v_0'
        #mvar = 'vmask'
        #gvar = 'V'
        # In standard file, everything is on T-grid
        evar = 'e3t_0'
        mvar = 'tmask'
        gvar = 'T'
        bad = False
    if ( var == 'MLW' ): 
        print(var[0])
        fvar = 'MLW'
        evar = 'e3t_0'
        mvar = 'tmask'
        gvar = 'T'
        bad = False
        
    if ( bad ):
        raise Exception(var+' variable did not register.  Check variable name.')

    nav_lon, nav_lat, SSH = stfd.read_fstd_var(file, 'SSH', typvar='P@')
    nav_lon, nav_lat, area = read_grid.read_coord(hall=hall,grid=gvar)
    e3t = read_grid.read_e3t_mesh(hall=hall,var=evar)
    mask = read_grid.read_tmask_mesh(hall=hall,var=mvar)
    lev , T_3D = stfd.read_fstd_multi_lev(file, fvar, typvar='P@')
    if ( ( T_3D.shape[0] != mask.shape[0] ) and ( T_3D.shape[0] != 1 ) ): #  Mis-match between mask and T_3D.
        # May have read in higher frequency variables
        lev, T_3D = stfd.read_fstd_multi_lev(file, fvar, vfreq=24, typvar='P@')
    if ( ( T_3D.shape[0] != mask.shape[0] ) and ( T_3D.shape[0] != 1 ) ): #  Mis-match between mask and T_3D STILL A PROBLEM
        # May need to remove bottom temperature from list of levels
        if ( 1.0 in lev ):
          ibottom = lev.index(1.0)
          lev.pop(ibottom)
          T_3D = np.delete(T_3D, ibottom, 0)
    if ( ( T_3D.shape[0] != mask.shape[0] ) and ( T_3D.shape[0] != 1 ) ): #  Mis-match between mask and T_3D STILL A PROBLEM
        raise Exception(fvar+' has wrong shape '+str(T_3D.shape))        
    bottom = read_grid.bottom_depth_from_e3t(e3t, mask)
    e3t[0,:,:] = e3t[0,:,:] + SSH
    botton = read_grid.bottom_depth_from_e3t(e3t, mask)

    try: 
        var2 = var[2]
    except:
        var2 = 'X'
    if ( var2 == 'D' ):  # D20, D26, D28
        Tlevel = float(var[0:2])
        FLD = isotherm(T_3D, mask, lev, bottom, Tlevel=Temperature_in_K(Tlevel), from_above=True)
        #FLD = isotherm(Temperature_in_C(T_3D), mask, lev, bottom, Tlevel=Tlevel, from_above=True)
        
    if ( ( var[0] == 'H' ) or ( var[0] == 'S' ) or ( var[0] == 'U' ) or ( var[0] == 'V' ) ):
        try:
            level=float(var[1:3])
        except:
            try: 
                var2 = var[2]
            except:
                var2 = 'X'
            try:
                if ( var2 == 'C' ):
                    level=float(var[1])*100.0
                if ( var2 == 'K' ):
                    level=float(var[1])*1000.0 
            except:
                level=0.0

    if ( ( var[0] == 'H' ) and ( var != 'HHP' ) ):                    
        FLD = heat_content(T_3D, e3t, mask, depth=level)    
    if ( ( var[0] == 'S' ) ):                    
        FLD = salt_content(T_3D, e3t, mask, depth=level)    
    if ( ( var[0] == 'U' ) or ( var[0] == 'V' ) ):                    
        FLD = depth_mean_velocity(T_3D, e3t, mask, depth=level)    
    if ( ( var == 'HHP' ) ):                    
        FLD = hurricane_heat_potential(T_3D, lev, e3t, mask, T_threshold=26)                                 
    if ( ( var == 'MLW' ) ):
        FLD = np.squeeze(T_3D)
    if ( ( var == 'TM' ) ):
        FLD = np.squeeze(T_3D[0])
    return nav_lon, nav_lat, FLD

def calculate_from_ncfile(file, fssh, var='H10', hall=hall):
    if ( len(file.split(":")) == 2 ): file=tmpcp.tmpcp(file)
    # Need U/V grid nav_lon/lat too, or instead.
    # Ditto umask and vmask
    for ii in range(len(var)):
        print(var[ii])
    bad=True
    if ( var == 'TM' ): 
        fvar = 'thetao' 
        evar = 'e3t_0'
        mvar = 'tmask'
        gvar = 'T'
        bad = False
    if ( var[0] == 'H' ): 
        print(var[0])
        fvar = 'thetao' 
        evar = 'e3t_0'
        mvar = 'tmask'
        gvar = 'T'
        bad = False
    try: 
        var2 = var[2]
    except:
        var2 = 'X'
    if ( var2 == 'D' ): 
        #print(var2)
        fvar = 'thetao' 
        evar = 'e3t_0'
        mvar = 'tmask'
        gvar = 'T'
        bad = False
    if ( var[0] == 'S' ): 
        print(var[0])
        fvar = 'so'
        evar = 'e3t_0'
        mvar = 'tmask'
        gvar = 'T'
        bad = False
    if ( var[0] == 'U' ): 
        print(var[0])
        fvar = 'uo'
        evar = 'e3u_0'
        mvar = 'umask'
        gvar = 'U'
        # In standard file, everything is on T-grid
        #evar = 'e3t_0'
        #mvar = 'tmask'
        #gvar = 'T'
        bad = True
    if ( var[0] == 'V' ): 
        print(var[0])
        fvar = 'vo'
        evar = 'e3v_0'
        mvar = 'vmask'
        gvar = 'V'
        # In standard file, everything is on T-grid
        #evar = 'e3t_0'
        #mvar = 'tmask'
        #gvar = 'T'
        bad = True
    if ( var == 'MLW' ): 
        print(var[0])
        fvar = 'MLW'
        evar = 'e3t_0'
        mvar = 'tmask'
        gvar = 'T'
        bad = True
        
    if ( bad ):
        raise Exception(var+' variable did not register.  Check variable name.')

    SSH_NC, LON_NC, LAT_NC, DEP, DAT = read_grid.read_netcdf(fssh, 'zos',depname=None, tunits='seconds', DAY0=datetime.datetime(1950, 1, 1, 0) )
    # NOW NEED TO AVERAGE IN TIME AND TRANSPOSE
    SSH = np.transpose(np.mean(SSH_NC, axis=0) )
    #nav_lon, nav_lat, SSH = stfd.read_fstd_var(file, 'SSH', typvar='P@')
    nav_lon, nav_lat, area = read_grid.read_coord(hall=hall,grid=gvar)
    e3t = read_grid.read_e3t_mesh(hall=hall,var=evar)
    mask = read_grid.read_tmask_mesh(hall=hall,var=mvar)
    #lev , T_3D = stfd.read_fstd_multi_lev(file, fvar, typvar='P@')
    T_NC, LON_NC, LAT_NC, lev, DAT = read_grid.read_netcdf(file, fvar,tunits='seconds', DAY0=datetime.datetime(1950, 1, 1, 0) )
    T_3D = np.transpose(np.squeeze(T_NC), [0, 2, 1])   # PROBLEM FOR 2D FIELD # NEED TO CONVERT TO KELVIN
    if ( fvar == 'thetao' ): # Convert from Celsius to Kelvin
        T_3D = T_3D + KCONV 

    if ( ( T_3D.shape[0] != mask.shape[0] ) and ( T_3D.shape[0] != 1 ) ): #  Mis-match between mask and T_3D STILL A PROBLEM
        raise Exception(fvar+' has wrong shape '+str(T_3D.shape))        
    bottom = read_grid.bottom_depth_from_e3t(e3t, mask)
    e3t[0,:,:] = e3t[0,:,:] + SSH
    botton = read_grid.bottom_depth_from_e3t(e3t, mask)

    try: 
        var2 = var[2]
    except:
        var2 = 'X'
    if ( var2 == 'D' ):  # D20, D26, D28
        Tlevel = float(var[0:2])
        FLD = isotherm(T_3D, mask, lev, bottom, Tlevel=Temperature_in_K(Tlevel), from_above=True)
        #FLD = isotherm(Temperature_in_C(T_3D), mask, lev, bottom, Tlevel=Tlevel, from_above=True)
        
    if ( ( var[0] == 'H' ) or ( var[0] == 'S' ) or ( var[0] == 'U' ) or ( var[0] == 'V' ) ):
        try:
            level=float(var[1:3])
        except:
            try: 
                var2 = var[2]
            except:
                var2 = 'X'
            try:
                if ( var2 == 'C' ):
                    level=float(var[1])*100.0
                if ( var2 == 'K' ):
                    level=float(var[1])*1000.0 
            except:
                level=0.0

    if ( ( var[0] == 'H' ) and ( var != 'HHP' ) ):                    
        FLD = heat_content(T_3D, e3t, mask, depth=level)    
    if ( ( var[0] == 'S' ) ):                    
        FLD = salt_content(T_3D, e3t, mask, depth=level)    
    if ( ( var[0] == 'U' ) or ( var[0] == 'V' ) ):                    
        FLD = depth_mean_velocity(T_3D, e3t, mask, depth=level)    
    if ( ( var == 'HHP' ) ):                    
        FLD = hurricane_heat_potential(T_3D, lev, e3t, mask, T_threshold=26)                                 
    if ( ( var == 'MLW' ) ):
        FLD = np.squeeze(T_3D)
    if ( ( var == 'TM' ) ):
        FLD = np.squeeze(T_3D[0])
    return nav_lon, nav_lat, FLD
