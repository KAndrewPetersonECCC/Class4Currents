import os
import sys
import traceback
#sys.path.insert(0, '/home/dpe000/Class4_Currents/python')
#from importlib import reload
       
import netCDF4
import shutil 
import numpy as np
import math
import datetime
import time
import subprocess
import glob
import copy
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import scipy.interpolate as si

import multiprocessing
import itertools
from functools import partial
num_cpus = len(os.sched_getaffinity(0))

import find_fcst_file
import get_archive
import stfd
import read_grid
import isoheatcontent
import find_value_at_point
import cplot
import find_hall
import datadatefile
import taylor_diagram
sys.path.insert(0, '/home/dpe000/python/properscoring-0.1')
import properscoring as ps
del sys.path[0]

KCONV=273.16

cmap_full_field='gist_stern_r'
cmap_anom_field='seismic'
cmap_anom_fmask='RdYlBu_r'
cmap_anom_field = copy.copy(cm.seismic)
cmap_anom_field.set_bad('g', 1.0)
cmap_anom_rield = copy.copy(cm.seismic_r)
cmap_anom_rield.set_bad('g', 1.0)
cmap_full_field = copy.copy(cm.gist_stern_r)
cmap_full_field.set_bad('r', 1.0)
cmap_anom_fmask = copy.copy(cm.RdYlBu_r)
cmap_anom_fmask.set_bad('g', 1.0)


mask = read_grid.read_mask(var='tmask')
maskt = read_grid.read_mask(var='tmask')
masku = read_grid.read_mask(var='umask')
maskv = read_grid.read_mask(var='vmask')
mask0 = np.squeeze(mask[0,:,:])
e3t = read_grid.read_e3t_mesh(var='e3t_0')
e1t = read_grid.read_mesh_var('e1t')
e2t = read_grid.read_mesh_var('e2t')
nav_lon, nav_lat, grid_area = read_grid.read_coord(grid='T')

psyfile='CLASS4_currents/class4_20211031_PSY4V3R1_orca12_currents.nc'
obsfile='CLASS4_currents_CCMEP/class4_20211031_GIOPS_orca025_currents.nc'
#shutil.copy(psyfile, obsfile)
psyfile='CLASS4_currents_CHARLY_MARCH/class4_20190101_PSY4V3R1_orca12_currents.nc'
obsfile='CLASS4_currents_CCMEP_FILT/class4_20190101_GIOPS_orca025_currents.nc'
psyffile='CLASS4_currents_CHARLY_MARCH/class4_20190101_PSY4V3R1_orca12_currents-filtr.nc'
obsffile='CLASS4_currents_CCMEP_FILT/class4_20190101_GIOPS_orca025_currents-filtr.nc'

file_best = '/fs/site6/eccc/mrd/rpnenv/dpe000//tmpdir/gu/ORCA025-CMC-ANAL_1d_grid_T_2021110100.nc'
file_besu = '/fs/site6/eccc/mrd/rpnenv/dpe000//tmpdir/gu/ORCA025-CMC-ANAL_1d_grid_U_2021110100.nc'
file_besv = '/fs/site6/eccc/mrd/rpnenv/dpe000//tmpdir/gu/ORCA025-CMC-ANAL_1d_grid_V_2021110100.nc'
file_fcts = ' tmpdir value'

site=find_hall.get_site()
datadir='/fs/'+site+'/eccc/mrd/rpnenv/dpe000/'
tempdir=datadir+'/'+'tmpdir'

missing=-999

date_ic2 = datetime.datetime(2019,7,3)
date_ic3 = datetime.datetime(2021,12,1,12)
date_mg1 = datetime.datetime(2020,1,21,12)

def UU_ORCA_to_NE(UV):
    hall=find_hall.find_hall()
    fila='/space/'+hall+'/sitestore/eccc/mrd/rpnenv/socn000/env_rhel-8-icelake-64/datafiles/constants/oce/repository/master/CONCEPTS/orca025/grids/orca025grid_new.std'
    LONA, LATA, TH = stfd.read_fstd_var(fila, 'LAAN')
    UO, VO = UV
    UE = UO * np.cos(TH) - VO * np.sin(TH)
    VN = UO * np.sin(TH) + VO * np.cos(TH)
    return [UE, VN]

def U_Local_to_NE(UVT):
    U0, V0, TH = UVT
    UE = UO * np.cos(TH) - VO * np.sin(TH)
    VN = UO * np.sin(TH) + VO * np.cos(TH)
    return [UE, VN]
    
def calc_m15(FDin, e3tin, maskin):
    if ( isinstance(FDin, list) ):
        F15 = []
        for FDii in FDin:
            F15.append(calc_m15(FDii, e3tin, maskin))
    else:
        F_20, __ = isoheatcontent.depth_integral(FDin, e3tin, maskin, depth=20)
        F_10, __ = isoheatcontent.depth_integral(FDin, e3tin, maskin, depth=10)
        F15 = ( F_20 - F_10 ) / 10
    return F15
    
def calc_u15(UVin, e3tin, maskin):
    UUin, VVin = UVin
    [u15, v15] = calc_m15(UVin, e3tin, maskin)
    return [u15, v15]

def calc_u15_from_3hr(file_3hrs, vfreq=3, newnorm=True):   ## ALSO SHOULD WORK FOR daily files if vfreq=24
    U10_list = []
    U20_list = []
    V10_list = []
    V20_list = []
    for file_3hr in file_3hrs:
        ldu, UD = stfd.read_fstd_multi_lev(file_3hr, 'UU2W', vfreq=vfreq, typvar='P@')
        #print(ldu)
        idu=[ldu.index(10.0), ldu.index(20.0)]
        U10_list.append(UD[idu[0]])
        U20_list.append(UD[idu[1]])
        ldv, VD = stfd.read_fstd_multi_lev(file_3hr, 'VV2W', vfreq=vfreq, typvar='P@')
        #print(ldv)
        idv=[ldv.index(10), ldv.index(20)]
        V10_list.append(VD[idv[0]])
        V20_list.append(VD[idv[1]])
    U10 = sum(U10_list)/len(U10_list)
    U20 = sum(U20_list)/len(U20_list) 
    V10 = sum(V10_list)/len(V10_list)
    V20 = sum(V20_list)/len(V20_list)
    #NOTE: SIMPLICATION OF U15D = (20.0*U20 - 10.0*U10) / (20.0-10.0)
    if ( newnorm ):
        U15D = 2.0*U20 - U10
        V15D = 2.0*V20 - V10
    else:
        U15D = ( U20 - U10 ) / 10.0
        V15D = ( V20 - V10 ) / 10.0
    return U15D, V15D
    
def calc_u00_from_3hr(file_3hrs, depth):
    U00_list = []
    V00_list = []
    T00_list = []
    for file_3hr in file_3hrs:
        ldt, TD = stfd.read_fstd_multi_lev(file_3hr, 'TM', vfreq=3, typvar='P@')
        try:
            idt=ldt.index(depth)
        except:
            idt=0
        T00_list.append(TD[idt])
        ldu, UD = stfd.read_fstd_multi_lev(file_3hr, 'UUW', vfreq=3, typvar='P@')
        idu=ldu.index(depth)
        U00_list.append(UD[idu])
        ldv, VD = stfd.read_fstd_multi_lev(file_3hr, 'VVW', vfreq=3, typvar='P@')
        idv=ldv.index(depth)
        V00_list.append(VD[idv])
    U00 = sum(U00_list)/len(U00_list)
    V00 = sum(V00_list)/len(V00_list)
    T00 = sum(T00_list)/len(T00_list)
    return U00, V00, T00
    
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
     
def speed_and_angle_list(Ulist):
    Slist=[]
    for U in Ulist:
        Slist.append(speed_and_angle_easy(U))
    return Slist
    
     
       
def fine_grid_mask(FD, e3tin):
    maskf=np.ones(FD.shape)

    e3t0=np.max(e3tin, axis=(1,2))
    e3tf=np.zeros(FD.shape)

    nzf, nxf, nyf = e3tf.shape
    print('DIM',nzf, nxf, nyf,e3t0.shape)

    for iz in range(nzf):
        e3tf[iz,:,:] = e3t0[iz] * np.ones((nxf,nyf))
            
    return maskf, e3tf

def find_obs_value(FD_LIST, NEMO_LONLAT, OBS_LONLAT):

    LONPT, LATPT = OBS_LONLAT
    LONN,  LATN  = NEMO_LONLAT

    ipt, jpt = find_value_at_point.find_nearest_point(LONPT, LATPT, LONN, LATN)
    
    FDpt = []
    FDpi = []
    for ifld, FDN in enumerate(FD_LIST):
        #ivalid = np.where(FDN.mask == False)
        #print('PRINT', ifld, FDN.shape, LONN.shape, LATN.shape, LONPT, LATPT)
        if ( LONN.ndim == 1 ):
          FDpt.append(FDN[ipt])
        else:
          FDpt.append(FDN[ipt, jpt])
        FDpi.append(find_value_at_point.interpolate_to_point(FDN, LONN, LATN, LONPT, LATPT))
    print('POINT VALUES', FDpt, FDpi)
    return FDpt, FDpi

def find_nearest_points(LONLATPT, LONLAT_LIST):
    LONPT, LATPT = LONLATPT
    IJPTS=[]
    for LONLATin in LONLAT_LIST:
        LONin, LATin = LONLATin
        ipto, jpto = find_value_at_point.find_nearest_point(LONPT, LATPT, LONin, LATin)
        IJPTS.append((ipto,jpto))
    return IJPTS

def msk_flds(FLDSin, LONLATin, maskin=mask0):
    LONin, LATin = LONLATin
    FLDSout=[]
    for FLD in FLDSin:
        LONM, LATM, FLDm = cplot.mask_field(LONin, LATin, FLD, maskin)
        FLDSout.append(FLDm)
    return FLDSout, (LONM, LATM)
        
def put_flds_latlon(FLDSin, LONLATin, ddeg=0.1, method='2sweeplinear'):
    LONin, LATin = LONLATin
    FLDSout=[]
    for FLD in FLDSin:
        LONf, LATf, FLDf = cplot.grdfld(LONin, LATin, FLD, ddeg=ddeg, method=method)
        FLDSout.append(FLDf)
    return FLDSout, (LONf, LATf)
        
def read_anal_file(file, grid='T', dates=[]):
    if ( grid=='T'): 
      fvar='thetao'
      depname='deptht'
      maskf = maskt
    if ( grid=='U'): 
      fvar='uo'
      depname='depthu'
      maskf = masku
    if ( grid=='V'):
      fvar='vo'
      depname='depthv'
      maskf = maskv
    FLD, LON, LAT, lev, DAT = read_grid.read_netcdf(file, fvar, depname=depname, tunits='seconds', DAY0=datetime.datetime(1950, 1, 1, 0) )
    FLD = np.squeeze(FLD)
    if ( FLD.ndim == 3 ):
        FLD = np.transpose(FLD, [0, 2, 1])   # PROBLEM FOR 2D FIELD # NEED TO CONVERT TO KELVIN
    elif ( FLD.ndim == 4 ):
        FLD = np.transpose(FLD, [0, 1, 3, 2])   # PROBLEM FOR 2D FIELD # NEED TO CONVERT TO KELVIN
    if ( fvar == 'thetao' ): # Convert from Celsius to Kelvin
        FLD = FLD + KCONV 
    if ( FLD.ndim == 3 ): 
        FLD = np.ma.masked_array(FLD, mask=1-maskf)
    elif ( FLD.ndim == 4 ):
        for iitime in range(len(FLD)):
            FLD[iitime, :, :, :] = np.ma.masked_array(FLD[iitime, :, :, :], mask=1-maskf)
    if ( len(dates) > 0 ):
        new_DAT = []
        indices = []
        for date in dates:
            if ( date in DAT ):
                index = DAT.index(date)
                new_DAT.append(DAT[index])
                indices.append(index)
        DAT=new_DAT
        if ( FLD.ndim == 4 ):
            FLD=FLD[indices,:,:,:]
    LON=np.transpose(LON)
    LAT=np.transpose(LAT)
    return FLD, LON,  LAT, lev, DAT
    
def map_to_A_grid(FLD, LON, LAT, grid='T'):
    if ( FLD.ndim==2 ):
        nx, ny = FLD.shape
        FLD=reshape(FLD, (1, nx, ny))
    nz, nx, ny = FLD.shape
    FLDN=np.zeros(FLD.shape)
    LONN=np.zeros(LON.shape)
    LATN=np.zeros(LAT.shape)
    if ( grid == 'T' ): 
        FLDN=FLD
        LONN=LON
        LATN=LAT
    elif ( grid == 'U' ):
        for ii in range(1,nx-1):
            FLDN[:, ii, :] = 0.5*FLD[:,ii-1,:]*masku[:,ii-1,:] + 0.5*FLD[:, ii,:]*masku[:, ii,:]
            LONN[   ii, :] = 0.5*LON[  ii-1,:]                 + 0.5*LON[   ii,:]
            LATN[   ii, :] = 0.5*LAT[  ii-1,:]                 + 0.5*LAT[   ii,:]
        FLDN[:, 0, :] = FLDN[:, nx-2, :]
        LONN[   0, :] = LONN[   nx-2, :]
        LATN[   0, :] = LATN[   nx-2, :]
        FLDN[:, nx-1, :] = FLDN[:, 1, :]
        LONN[   nx-1, :] = LONN[   1, :]
        LATN[   nx-1, :] = LATN[   1, :]
    elif ( grid=='V' ):
       for jj in range(1, ny-1):
           FLDN[:,:,jj] = 0.5*FLD[:,:,jj-1]*maskv[:,:,jj-1] + 0.5*FLD[:,:,jj]*maskv[:,:,jj]
           LONN[  :,jj] = 0.5*LON[  :,jj-1]                 + 0.5*LON[  :,jj]
           LATN[  :,jj] = 0.5*LAT[  :,jj-1] + 0.5*LAT[  :,jj]
       FLDN[:,:,0] = FLD[:,:,1]*maskv[:,:,1]  # Free Slip condition
       LONN[  :,0] = LON[  :,1]  # ALMOST?
       LATN[:,0] = 1.5*LAT[:,1] - 0.5*LAT[:,2]  # ?? 
       for ii in range(nx):
           FLDN[:,ii,ny-1] = FLDN[:, nx-2-ii, ny-2]
           LONN[  ii,ny-1] = LONN[   nx-2-ii, ny-2]
           LATN[  ii,ny-1] = LATN[   nx-2-ii, ny-2]
    FLDN = np.ma.array(FLDN * maskt, mask=1-maskt)
    
    return np.squeeze(FLDN), LONN, LATN

grid_dir='/fs/site6/eccc/mrd/rpnenv/dpe000/GDPS/gdps_Class4/'
script_dir='/fs/homeu2/eccc/mrd/ords/rpnenv/dpe000/Class4_Currents/jobscripts/'
template_dir='/fs/homeu2/eccc/mrd/ords/rpnenv/dpe000/Class4_Currents/templates/'
grid_ref='dest_grid.std.2'

def cst_interpolation(source_file, destination_file, ref_grid=grid_dir+grid_ref):
    script=script_dir+'perform_interpolation.sh'
    rc = subprocess.call(['bash', script, '-fs='+source_file, '-fd='+destination_file])
    return rc

def filter_standard_file(filter, source_file, destination_file):
    script=script_dir+'perform_filter.sh'
    clobber=''
    if ( destination_file == source_file ): 
        clobber="--clobber"
        destination_file = source_file+'.tmp'
    rc=subprocess.call(['bash' , script, '-f='+filter, '-i='+source_file, '-o='+destination_file, clobber])
    return rc
    
def test_interpolation(date, lead, ref_grid=grid_dir+grid_ref):
    subprocess.call(['echo', 'Good Morning Drew'])
    fcst_file=find_fcst_file.find_fcst_file('OPD', date, lead, 0, src='ocn', execute=True)
    # Filter OUT TM, UUW and VVW
    filter_standard_file(template_dir+'/flt.ksh', fcst_file, fcst_file)
    intr_file=fcst_file+'.Z01'
    print( fcst_file, intr_file )
    rc = cst_interpolation(fcst_file, intr_file, ref_grid=ref_grid)
    if ( rc != 0 ):
       print("I'm afraid I can't do that Drew", rc)
    LONI, LATI, TI = stfd.read_fstd_var(intr_file, 'TM', typvar='P@')
    leu, UI = stfd.read_fstd_multi_lev(intr_file, 'UUW',vfreq=24, typvar='P@')
    lev, VI = stfd.read_fstd_multi_lev(intr_file, 'VVW',vfreq=24, typvar='P@')
    let, TI = stfd.read_fstd_multi_lev(intr_file, 'TM',vfreq=24, typvar='P@')
    return [UI, VI, TI], [LONI, LATI], [leu, lev, let]

def plot_fields(FLDS, LONLATFLD, suptitle=None, grid=False, outfile_prefix='PLOT/'):
    (LONFLD, LATFLD) = LONLATFLD
    # [TFLD, UFLD, VFLD, SPEED, ANGLE] = FLDS  ## DON'T UNCOMMENT -- FAILS if change number of fields and still want to work.
    for ifld, FLD in enumerate(FLDS):
        CLEV = np.arange(-5, 5.5, 0.5)
        if ( ifld == 0 ): CLEV = np.arange(-2, 34, 2)
        if ( ifld == 3 ): CLEV = np.arange(0, 5.2, 0.2)
        if ( ifld == 4 ): CLEV = np.pi * np.arange(-1.0,1.1,0.1)
        if ( ifld == 0 ): FLD = FLD - KCONV
        if ( ifld == 0 ): title='T'
        if ( ifld == 1 ): title='U'
        if ( ifld == 2 ): title='V'
        if ( ifld == 3 ): title='S'
        if ( ifld == 4 ): title='A'
        cmap_use=cmap_anom_field
        if ( ifld == 0 ): cmap_use=cmap_full_field
        if ( ifld == 3 ): cmap_use=cmap_full_field

        outfile=outfile_prefix+title+'.png'
        if ( not grid ):
            cplot.pcolormesh(LONFLD, LATFLD, FLD, levels=CLEV, cmap=cmap_use, project='PlateCarree', 
               outfile=outfile, make_global=True, title=title, suptitle=suptitle, 
               cbar=True, obar='vertical', fontsizes=None, add_gridlines=True)
        else:
            cplot.grd_pcolormesh(LONFLD, LATFLD, FLD, levels=CLEV, cmap=cmap_use, project='PlateCarree', 
               outfile=outfile, ddeg=0.2, make_global=True, title=title, suptitle=suptitle, 
               cbar=True, obar='vertical', fontsizes=None, add_gridlines=True)
        outfile=outfile_prefix+title+'.NP.png'
        cplot.pcolormesh(LONFLD, LATFLD, FLD, levels=CLEV, cmap=cmap_use, project='NorthPolarStereo', 
               outfile=outfile, title=title, suptitle=suptitle, 
               cbar=True, obar='vertical', fontsizes=None, box=[-180, 180, 65, 90])
              
    return

def calc_error(obser, model, isangle=-1, etype='mean'):

    zero_obs = 1
    zero_fld = 1
    if ( etype == 'obs' ): 
        zero_fld=0
        zero_obs=1
    if ( etype == 'fld' ): 
        zero_obs=0
        zero_fld=-1    
    if ( ( isinstance(model, list) ) or ( isinstance(model, tuple) ) ):
        ERROR = []
        for imodel in model:
            iERROR = calc_error(obser, imodel, isangle=isangle, etype=etype)
            ERROR.append(iERROR)
        return ERROR
        
    nobss, nvars, ndeps = obser.shape
    mobss, mvars, mfcst, mensm, mdeps = (0, 0, 0, 0, 0)
    if ( model.ndim == 3 ):
        mobss, mvars, mdeps = model.shape
    elif ( model.ndim == 4 ): 
        mobss, mvars, mfcst, mdeps = model.shape
    elif ( model.ndim == 5 ): 
        mobss, mvars, mfcst, mensm, mdeps = model.shape

    if ( (nobss != mobss) or (nvars != mvars) or (ndeps != mdeps) ):
        print('SHAPE MISMATCH')
        print((nobss,mobss),(nvars, mvars), (ndeps, mdeps))
        return None
        
    if ( mfcst == 0 ):
        ERROR = zero_obs*obser - zero_fld*model
        if ( isangle >= 0 ):
            for iobss in range(mobss):
              for ideps in range(mdeps):
                if ( ERROR[iobss,isangle,ideps] >=      np.pi ): ERROR[iobss,isangle,ideps] = ERROR[iobss,isangle,ideps] - 2*np.pi
                if ( ERROR[iobss,isangle,ideps] <= -1.0*np.pi ): ERROR[iobss,isangle,ideps] = 2*np.pi + ERROR[iobss,isangle,ideps]
    if ( mfcst > 0 and mensm == 0 ):
        ERROR = 0.0*model.copy()
        for ifcst in range(mfcst):
            ERROR[:,:,ifcst,:] = zero_obs*obser[:,:,:] - zero_fld*model[:,:,ifcst,:]
            if ( isangle >= 0 ):
                for iobss in range(mobss):
                  for ideps in range(mdeps):
                    if ( ERROR[iobss,isangle,ifcst,ideps] >=      np.pi ): ERROR[iobss,isangle,ifcst,ideps] = ERROR[iobss,isangle,ifcst,ideps] - 2*np.pi
                    if ( ERROR[iobss,isangle,ifcst,ideps] <= -1.0*np.pi ): ERROR[iobss,isangle,ifcst,ideps] = 2*np.pi + ERROR[iobss,isangle,ifcst,ideps]
    if ( mensm > 0 ):
        ERROR = 0.0*model.copy()
        for iensm in range(mensm):
          for ifcst in range(mfcst):
            ERROR[:,:,ifcst,iensm,:] = zero_obs*obser[:,:,:] - zero_fld*model[:,:,ifcst, iensm, :]
            if ( isangle >= 0 ):
                for iobss in range(mobss):
                  for ideps in range(mdeps):
                    if ( ERROR[iobss,isangle,ifcst,iensm,ideps] >=      np.pi ): ERROR[iobss,isangle,ifcst,iensm,ideps] = ERROR[iobss,isangle,ifcst,iensm,ideps] - 2*np.pi
                    if ( ERROR[iobss,isangle,ifcst,iensm,ideps] <= -1.0*np.pi ): ERROR[iobss,isangle,ifcst,iensm,ideps] = 2*np.pi + ERROR[iobss,isangle,ifcst,iensm,ideps]

    if ( ( etype == 'mean' ) or ( etype == 'obs' ) or ( etype == 'fld' ) ):
        pass
    elif ( etype == 'square' ):
        ERROR = np.square(ERROR)
    elif ( etype == 'absolute' ):
        ERROR = np.absolute(ERROR)
        

    if ( etype == 'cross' ):
        MOD = calc_error(obser, model, etype='obs', isangle=isangle)
        OBS = calc_error(obser, model, etype='fld', isangle=isangle)
        ERROR = np.multiply(MOD, OBS)
            
    return ERROR

def calc_mean_error(obser, model, isangle=-1):
    if ( ( isinstance(model, list) ) or ( isinstance(model, tuple) ) ):
        mean_error = []
        for imodel in model:
            iERROR = calc_mean_error(obser, imodel, isangle=isangle)
            mean_error.append(iERROR)
        return mean_error

    nobss, nvars, ndeps = obser.shape
    mobss, mvars, mfcst, mdeps = (0, 0, 0, 0)
    if ( model.ndim == 3 ):
        mobss, mvars, mdeps = model.shape
    elif ( model.ndim == 4 ): 
        mobss, mvars, mfcst, mdeps = model.shape

    if ( (nobss != mobss) or (nvars != mvars) or (ndeps != mdeps) ):
        print('SHAPE MISMATCH')
        return None

    error = calc_error(obser, model, isangle=isangle)
    mean_error = np.mean(error, axis=0)
    rmse_error = np.mean(np.square(error), axis=0)
    abse_error = np.mean(np.absolute(error), axis=0)
    
    return mean_error, rmse_error, abse_error
        
def calc_crps(obser, model, ens_axis=2):  # WHAT COULD BE DONE FOR AN ANGLE?
    if ( ( isinstance(model, list) ) or ( isinstance(model, tuple) ) ):
        crps_error = []
        for imodel in model:
            iERROR = calc_crps(obser, imodel)
            mean_error.append(iERROR)
        return crps_error
    nobss, nvars, ndeps = obser.shape
    nensm, mobss, mvars, nfcst, mdeps = (1, 0, 0, 1, 0)
    if ( model.ndim == 3 ):
        nensm = 1
        mobss, mvars, mdeps = model.shape
    if ( model.ndim == 4 and ens_axis == 2 ):
        mobss, mvars, nensm, mdeps = model.shape
    if ( model.ndim == 4 and ens_axis == -1 ):
        nensm = 1
        mobss, mvars, nfcst, mdeps = model.shape
    if ( model.ndim == 5 and ens_axis == 3 ):
        mobss, mvars, nfcst, nensm, mdeps = model.shape
    #print( 'nensm, nobss, nvars, ndeps, nfcst', nensm, nobss, nvars, ndeps, nfcst)   
    #print( 'nensm, mobss, mvars, mdeps, nfcst', nensm, mobss, mvars, mdeps, nfcst)   
    if ( (nobss != mobss) or (nvars != mvars) or (ndeps != mdeps) ):
        print('SHAPE MISMATCH')
        return None
     
    if ( ens_axis == -1 ):
        crps_error = np.ma.array(np.zeros(model.shape), mask=False)
    elif ( ens_axis == 2 ):
        crps_error = np.ma.array(np.zeros(model[:,:,0,:].shape), mask=False)
    elif ( ens_axis == 3 ):
       crps_error = np.ma.array(np.zeros(model[:,:,:,0,:].shape), mask=False) 
    #print('SHAPE', crps_error.shape)

    for iobss in range(nobss):
       for ivars in range(nvars):
           for ideps in range(ndeps):
               obs=obser[iobss, ivars, ideps]
               if ( nfcst == 1 ):
                   if ( nensm == 1 ):
                       flaa = np.reshape(model, [nobss, nvars, 1, ndeps])[iobss, ivars, :, ideps]
                   else:
                       flaa = model[iobss, ivars, :, ideps]
                   flda = flaa.tolist()
                   if ( np.all(flaa.mask) ):
                       crps_error.mask[iobss, ivars, ideps] = True
                   else:
                       crps_error[iobss, ivars, ideps] = ps.crps_ensemble(obs, flda)
               elif ( nfcst > 1 ):
                   for ifcst in range(nfcst):
                       if ( nensm == 1 ):
                           flaa=np.reshape(model, [nobss, nvars, nfcst, 1, ndeps])[iobss, ivars, ifcst, :, ideps]
                       else:
                           flaa = model[iobss, ivars, ifcst, :, ideps]
                       flda = flaa.tolist()
                       if ( np.all(flaa.mask) ):
                           crps_error.mask[iobss, ivars, ifcst, ideps] = True
                       else:
                           crps_error[iobss, ivars, ifcst, ideps] = ps.crps_ensemble(obs, flda)
    return crps_error


tate=datetime.datetime(2021, 10, 31)
def process_obsfile(date=tate, TEST_SINGLE=False, Plot=False, CHARLY=True, filter=True, nPersist=True):
    datestr=date.strftime('%Y%m%d')   
    time0=time.time()
    if ( not CHARLY ):
        psyfile='CLASS4_currents/class4_'+datestr+'_PSY4V3R1_orca12_currents.nc'
        obsfile1='CLASS4_currents_CCMEP/class4_'+datestr+'_GIOPS_orca025_currents.1.nc'   # Nearest neighbour
        obsfile2='CLASS4_currents_CCMEP/class4_'+datestr+'_GIOPS_orca025_currents.2.nc'   # Nearest neighbour on 0.2 lat/lon grid
        obsfile3='CLASS4_currents_CCMEP/class4_'+datestr+'_GIOPS_orca025_currents.3.nc'   # Nearest neighbour on CST interp 0.2 lat/lon grid
        obsfile4='CLASS4_currents_CCMEP/class4_'+datestr+'_GIOPS_orca025_currents.4.nc'   # Nearest neighbour on 0.2 lat/lon grid from U10/20 outputs.  
        obsfile5='CLASS4_currents_CCMEP/class4_'+datestr+'_GIOPS_orca025_currents.5.nc'   # Interp to Obs point 
        obsfile6='CLASS4_currents_CCMEP/class4_'+datestr+'_GIOPS_orca025_currents.6.nc'   # Interp to Obs pont from U10/20 outputs
    elif ( CHARLY ):
      if ( not filter ): 
        psyfile='CLASS4_currents_CHARLY/class4_'+datestr+'_PSY4V3R1_orca12_currents.nc'
        obsfile1='CLASS4_currents_CCMEP_UFIL/class4_'+datestr+'_GIOPS_orca025_currents.1.nc'
        obsfile2='CLASS4_currents_CCMEP_UFIL/class4_'+datestr+'_GIOPS_orca025_currents.2.nc'
        obsfile3='CLASS4_currents_CCMEP_UFIL/class4_'+datestr+'_GIOPS_orca025_currents.3.nc'
        obsfile4='CLASS4_currents_CCMEP_UFIL/class4_'+datestr+'_GIOPS_orca025_currents.4.nc'
        obsfile5='CLASS4_currents_CCMEP_UFIL/class4_'+datestr+'_GIOPS_orca025_currents.5.nc'  
        obsfile6='CLASS4_currents_CCMEP_UFIL/class4_'+datestr+'_GIOPS_orca025_currents.6.nc'   
      elif ( filter ): 
        psyfile='CLASS4_currents_CHARLY/class4_'+datestr+'_PSY4V3R1_orca12_currents-filtr.nc'
        obsfile1='CLASS4_currents_CCMEP_FILT/class4_'+datestr+'_GIOPS_orca025_currents.f1.nc'
        obsfile2='CLASS4_currents_CCMEP_FILT/class4_'+datestr+'_GIOPS_orca025_currents.f2.nc'
        obsfile3='CLASS4_currents_CCMEP_FILT/class4_'+datestr+'_GIOPS_orca025_currents.f3.nc'
        obsfile4='CLASS4_currents_CCMEP_FILT/class4_'+datestr+'_GIOPS_orca025_currents.f4.nc'
        obsfile5='CLASS4_currents_CCMEP_FILT/class4_'+datestr+'_GIOPS_orca025_currents.f5.nc'   
        obsfile6='CLASS4_currents_CCMEP_FILT/class4_'+datestr+'_GIOPS_orca025_currents.f6.nc' 
        
    print('obsfile', obsfile1, obsfile2, obsfile3)
    # READ IN OBS FILE
    (LONO, LATO, depth), (obser, beste, fcstv, persi) = read_obsfile(psyfile)

    inite = beste.copy()
    nersi = persi.copy()

    Mobser = speed_and_angle_easy(obser)
    Mbeste = speed_and_angle_easy(beste)
    Minite = speed_and_angle_easy(inite)
    Mfcstv = speed_and_angle_easy(fcstv)
    Mpersi = speed_and_angle_easy(persi)
    Mnersi = speed_and_angle_easy(nersi)
    

    nobss, nvars, nfcst, ndeps = fcstv.shape
    if ( TEST_SINGLE ):
        nobss_loop = 1
    else:
        nobss_loop = nobss
        
    ## COPY OF OBS FILE VARIABLES FOR EACH OUTPUT
    tmp_beste = 0.0*beste.copy()
    tmp_inite = 0.0*beste.copy()
    tmp_fcstv = 0.0*fcstv.copy()
    tmp_persi = 0.0*persi.copy()
    tmp_nersi = 0.0*persi.copy()
    tmp_beste.mask = False
    tmp_fcstv.mask = False
    tmp_persi.mask = False
    BESTE_LIST = [tmp_beste.copy(), tmp_beste.copy(), tmp_beste.copy(), tmp_beste.copy(), tmp_beste.copy(), tmp_beste.copy()]
    INITE_LIST = [tmp_inite.copy(), tmp_inite.copy(), tmp_inite.copy(), tmp_inite.copy(), tmp_inite.copy(), tmp_inite.copy()]
    FCSTV_LIST = [tmp_fcstv.copy(), tmp_fcstv.copy(), tmp_fcstv.copy(), tmp_fcstv.copy(), tmp_fcstv.copy(), tmp_fcstv.copy()]
    PERSI_LIST = [tmp_persi.copy(), tmp_persi.copy(), tmp_persi.copy(), tmp_persi.copy(), tmp_persi.copy(), tmp_persi.copy()]
    NERSI_LIST = [tmp_nersi.copy(), tmp_nersi.copy(), tmp_nersi.copy(), tmp_nersi.copy(), tmp_nersi.copy(), tmp_nersi.copy()]
    
    tempb=np.zeros((nobss, ndeps))
    tempi=np.zeros((nobss, ndeps))
    tempf=np.zeros((nobss, nfcst, ndeps))
    tempp=np.zeros((nobss, nfcst, ndeps))
    tempn=np.zeros((nobss, nfcst, ndeps))

    TEMPB_LIST = [tempb.copy(), tempb.copy(), tempb.copy(), tempb.copy(), tempb.copy(), tempb.copy()]
    TEMPI_LIST = [tempi.copy(), tempi.copy(), tempi.copy(), tempi.copy(), tempi.copy(), tempi.copy()]
    TEMPF_LIST = [tempf.copy(), tempf.copy(), tempf.copy(), tempf.copy(), tempf.copy(), tempf.copy()]
    TEMPP_LIST = [tempp.copy() ,tempp.copy(), tempp.copy(), tempp.copy(), tempp.copy(), tempp.copy()]
    TEMPN_LIST = [tempn.copy(), tempn.copy(), tempn.copy(), tempn.copy(), tempn.copy(), tempn.copy()]

    MERCATOR_MEAN_ERRORS =  calc_mean_error(obser, (beste, inite, fcstv, persi, nersi))   
    mr_ERRMb, mr_ERRMi, mr_ERRMf, mr_ERRMp, mr_ERRMn = MERCATOR_MEAN_ERRORS
    mean_ERRMb, rmse_ERRMb, abse_ERRMb = mr_ERRMb
    mean_ERRMi, rmse_ERRMi, abse_ERRMi = mr_ERRMi
    mean_ERRMf, rmse_ERRMf, abse_ERRMf = mr_ERRMf 
    mean_ERRMp, rmse_ERRMp, abse_ERRMp = mr_ERRMp
    mean_ERRMn, rmse_ERRMn, abse_ERRMn = mr_ERRMn
    MERCATOR_MMEAN_ERRORS = calc_mean_error(Mobser, (Mbeste, Minite, Mfcstv, Mpersi, Mnersi),isangle=1)
    mr_ERMMb, mr_ERMMi, mr_ERMMf, mr_ERMMp, mr_ERMMn = MERCATOR_MMEAN_ERRORS
    mean_ERMMb, rmse_ERMMb, abse_ERMMb = mr_ERMMb
    mean_ERMMi, rmse_ERMMi, abse_ERMMi = mr_ERMMi
    mean_ERMMf, rmse_ERMMf, abse_ERMMf = mr_ERMMf 
    mean_ERMMp, rmse_ERMMp, abse_ERMMp = mr_ERMMp
    mean_ERMMn, rmse_ERMMn, abse_ERMMn = mr_ERMMn

    if ( not TEST_SINGLE ):  ## DON"T DO FOR TEST  
        
        #nobss, nvars, nfcst, ndeps = fcstv.shape
        ## Keep 15m 
        outpre='ERRORS/PSY4'
        tmppre='PSY4'
        outdep=1
        if ( CHARLY ):
          outpre='ERRORS_UFIL/PSY4'
          tmppre='uPSY4'
          outdep=0
          if ( filter ): 
            outpre='ERRORS_FILT/PSY4'
            tmppre='fPSY4'

        write_mean_errors(outpre, int(datestr),  MERCATOR_MEAN_ERRORS, vec=['u','v'], tmp_prefix='site6/Class4_Currents/TMP/'+tmppre, depth_index=outdep)
        write_mean_errors(outpre, int(datestr), MERCATOR_MMEAN_ERRORS, vec=['s','a'], tmp_prefix='site6/Class4_Currents/TMP/'+tmppre, depth_index=outdep)
        #nobss, nvars, nfcst, ndeps = fcstv.shape
    timen=time.time()
    timep=timen-time0
    print("TIMING :: Total Observation Preparation Time: "+str(timep))

    time0=timen
    # READ IN MODEL standard FILE
    (LONN, LATN) = (nav_lon, nav_lat)
    bedate=date + datetime.timedelta(days=1)
    #find days to next gd analysis
    andiff = ( 2-bedate.weekday() ) % 7
    andate = bedate + datetime.timedelta(days=andiff)
    anal=False   # retrieve TRIAL
    if ( andiff == 0 ): anal=True   # retrieve anal

    GD='gd'
    GU='gu' 
    if ( andate <= date_ic2-datetime.timedelta(days=7) ): GD='pd'
    if ( bedate <= date_ic2 ): GU='pu'
    file_best, __ =find_fcst_file.find_anal_file(andate, system=GD, anal=anal, var='T',execute=True)
    file_besu, __ =find_fcst_file.find_anal_file(andate, system=GD, anal=anal, var='U',execute=False)
    file_besv, __ =find_fcst_file.find_anal_file(andate, system=GD, anal=anal, var='V',execute=False)
    file_init, __ =find_fcst_file.find_anal_file(bedate, system=GU, anal=True, var='T',execute=True)
    file_iniu, __ =find_fcst_file.find_anal_file(bedate, system=GU, anal=True, var='U',execute=False)
    file_iniv, __ =find_fcst_file.find_anal_file(bedate, system=GU, anal=True, var='V',execute=False)
    print('file_best', file_best)

    timen=time.time()
    timep=timen-time0
    print("TIMING :: Best Estimate Retieval Time: "+str(timep))
    time0=timen

    ## BEST ESTIMATE    
    TM, LONT, LATT, LEVT, DAT = read_anal_file(file_best, grid='T', dates=[bedate]); TM=np.squeeze(TM)
    UW, LONU, LATU, LEVU, DAT = read_anal_file(file_besu, grid='U', dates=[bedate]); UW=np.squeeze(UW)
    VW, LONV, LATV, LEVV, DAT = read_anal_file(file_besv, grid='V', dates=[bedate]); VW=np.squeeze(VW)
    UU, LONUU, LATUU = map_to_A_grid(np.squeeze(UW), LONU, LATU, grid='U')
    VV, LONVV, LATVV = map_to_A_grid(np.squeeze(VW), LONV, LATV, grid='V')
    # Important for U15 calculations that arrays not be masked.
    print(type(TM), type(UU), type(VV))
    if ( isinstance(TM, np.ma.core.MaskedArray) ): TM = TM.data
    if ( isinstance(UU, np.ma.core.MaskedArray) ): UU = UU.data
    if ( isinstance(VV, np.ma.core.MaskedArray) ): VV = VV.data
    print(type(TM), type(UU), type(VV))
    [T00, U00, V00] = [TM[0,:,:], UU[0,:,:], VV[0,:,:]]
    S00, A00 = speed_and_angle(U00, V00)
    [T15, U15, V15] = calc_m15([TM, UU, VV], e3t, mask)
    S15, A15 = speed_and_angle(U15, V15)
    if ( Plot ):
        suptitle='0m'+datestr+'Native Grid Best Estimate'
        prefix='PLOTS/00_BENG_'+datestr+'_'
        plot_fields((T00, U00, V00, S00, A00), (LONN, LATN), suptitle=suptitle, grid=True, outfile_prefix=prefix)
        suptitle='15m'+datestr+'Native Grid Best Estimate'
        prefix='PLOTS/15_BENG_'+datestr+'_'
        plot_fields((T15, U15, V15, S15, A15), (LONN, LATN), suptitle=suptitle, grid=True, outfile_prefix=prefix)
    [U00, V00] = UU_ORCA_to_NE([U00, V00])
    [U15, V15] = UU_ORCA_to_NE([U15, V15])

    # this actually removes non-ocean points
    (T00m, U00m, V00m, T15m, U15m, V15m), (LONM, LATM) = msk_flds([T00, U00, V00, T15, U15, V15], [LONN, LATN], mask0)
    # now grid to lat long
    (T00f, U00f, V00f, T15f, U15f, V15f), (LONF, LATF) = put_flds_latlon([T00m, U00m, V00m, T15m, U15m, V15m], [LONM, LATM], ddeg=0.2, method='2sweeplinear')
    # now grid to observatoin locations
    (TOOo, UOOo, VOOo, T15o, U15o, V15o) = interp_direct_obs( (LONM, LATM), [T00m, U00m, V00m, T15m, U15m, V15m], [LONO, LATO], method='2sweeplinear')
    S00f, A00f = speed_and_angle(U00f, V00f)
    S15f, A15f = speed_and_angle(U15f, V15f)
    if ( Plot == True ):
        suptitle='0m'+datestr+'Best Estimate'
        prefix='PLOTS/00_BEST_'+datestr+'_'
        plot_fields((T00f, U00f, V00f, S00f, A00f), (LONF, LATF), suptitle=suptitle, grid=False, outfile_prefix=prefix)
        suptitle='15m'+datestr+'Best Estimate'
        prefix='PLOTS/15_BEST_'+datestr+'_'
        plot_fields((T15f, U15f, V15f, S15f, A15f), (LONF, LATF), suptitle=suptitle, grid=False, outfile_prefix=prefix)

    if ( ndeps == 2 ):
        BEST=[
          [[T00m, U00m, V00m], [T15m, U15m, V15m]],
          [[T00f, U00f, V00f], [T15f, U15f, V15f]],
          [[T00o, U00o, V00o], [T15o, U15o, V15o]]
             ]
    elif ( ndeps == 1 ):
        BEST=[
          [[T15m, U15m, V15m]],
          [[T15f, U15f, V15f]],
          [[T15o, U15o, V15o]]
             ]

    ## INIT ESTIMATE    
    TM, LONT, LATT, LEVT, DAT = read_anal_file(file_init, grid='T', dates=[bedate]); TM=np.squeeze(TM)
    UW, LONU, LATU, LEVU, DAT = read_anal_file(file_iniu, grid='U', dates=[bedate]); UW=np.squeeze(UW)
    VW, LONV, LATV, LEVV, DAT = read_anal_file(file_iniv, grid='V', dates=[bedate]); VW=np.squeeze(VW)
    UU, LONUU, LATUU = map_to_A_grid(np.squeeze(UW), LONU, LATU, grid='U')
    VV, LONVV, LATVV = map_to_A_grid(np.squeeze(VW), LONV, LATV, grid='V')
    # Important for U15 calculations that arrays not be masked.
    print(type(TM), type(UU), type(VV))
    if ( isinstance(TM, np.ma.core.MaskedArray) ): TM = TM.data
    if ( isinstance(UU, np.ma.core.MaskedArray) ): UU = UU.data
    if ( isinstance(VV, np.ma.core.MaskedArray) ): VV = VV.data
    print(type(TM), type(UU), type(VV))
    [T00, U00, V00] = [TM[0,:,:], UU[0,:,:], VV[0,:,:]]
    S00, A00 = speed_and_angle(U00, V00)
    [T15, U15, V15] = calc_m15([TM, UU, VV], e3t, mask)
    S15, A15 = speed_and_angle(U15, V15)
    if ( Plot ):
        suptitle='0m'+datestr+'Native Grid Init Estimate'
        prefix='PLOTS/00_INNG_'+datestr+'_'
        plot_fields((T00, U00, V00, S00, A00), (LONN, LATN), suptitle=suptitle, grid=True, outfile_prefix=prefix)
        suptitle='15m'+datestr+'Native Grid Init Estimate'
        prefix='PLOTS/15_INNG_'+datestr+'_'
        plot_fields((T15, U15, V15, S15, A15), (LONN, LATN), suptitle=suptitle, grid=True, outfile_prefix=prefix)
    [U00, V00] = UU_ORCA_to_NE([U00, V00])
    [U15, V15] = UU_ORCA_to_NE([U15, V15])

    # this actually removes non-ocean points
    (T00m, U00m, V00m, T15m, U15m, V15m), (LONM, LATM) = msk_flds([T00, U00, V00, T15, U15, V15], [LONN, LATN], mask0)
    # now grid to lat long
    (T00f, U00f, V00f, T15f, U15f, V15f), (LONF, LATF) = put_flds_latlon([T00m, U00m, V00m, T15m, U15m, V15m], [LONM, LATM], ddeg=0.2, method='2sweeplinear')
    # now grid to observation locations
    (TOOo, UOOo, VOOo, T15o, U15o, V15o) = interp_direct_obs( (LONM, LATM), [T00m, U00m, V00m, T15m, U15m, V15m], [LONO, LATO], method='2sweeplinear')
    S00f, A00f = speed_and_angle(U00f, V00f)
    S15f, A15f = speed_and_angle(U15f, V15f)
    if ( Plot == True ):
        suptitle='0m'+datestr+'Best Estimate'
        prefix='PLOTS/00_BEST_'+datestr+'_'
        plot_fields((T00f, U00f, V00f, S00f, A00f), (LONF, LATF), suptitle=suptitle, grid=False, outfile_prefix=prefix)
        suptitle='15m'+datestr+'Best Estimate'
        prefix='PLOTS/15_BEST_'+datestr+'_'
        plot_fields((T15f, U15f, V15f, S15f, A15f), (LONF, LATF), suptitle=suptitle, grid=False, outfile_prefix=prefix)

    if ( ndeps == 2 ):
        INIT=[
          [[T00m, U00m, V00m], [T15m, U15m, V15m]],
          [[T00f, U00f, V00f], [T15f, U15f, V15f]],
          [[T00o, U00o, V00o], [T15o, U15o, V15o]]

             ]
    elif ( ndeps == 1 ):
        INIT=[
          [[T15m, U15m, V15m]],
          [[T15f, U15f, V15f]],
          [[T15o, U15o, V15o]]
             ]

    timen=time.time()
    timep=timen-time0
    print("TIMING :: Best Estimate Processing Time: "+str(timep))
    time0=timen
    timel=time0
    
     
    FCST=[]
    PERS=[]
    NERS=[]

    # GET A WHOLE BUNCH OF FILES AHEAD OF REQUIRED TIME:  LESS CALLS TO RARC?  
    SYS='NUL'    
    if ( bedate-datetime.timedelta(days=10) > date_ic2 ): 
        SYS='OPD'
        tmpdir=tempdir+'/'+SYS+'/'
        branch='operation.forecasts.giops.prog.glboce'
    if ( bedate < date_ic2 ): 
        SYS='PSD'
        tmpdir=tempdir+'/'+SYS+'/'
        branch='parallel.forecasts.giops.prog.glboce'
    if ( SYS != 'NUL'):
        rc = get_archive.get_archive_leads(tmpdir, branch, [bedate-datetime.timedelta(days=10), bedate], (np.arange(0,240,3)+3).tolist(), ensnum=None, execute=True)

    for ifcst in range(nfcst):
        jfcst=ifcst+1
        fchour=jfcst*24
        fchours=(np.arange(fchour-24, fchour, 3)+3).tolist()
        fcdate = bedate - datetime.timedelta(days=jfcst)
        andiff = (2 - fcdate.weekday() ) % 7
        andate = fcdate + datetime.timedelta(days=andiff) 
        anal=False
        if ( andiff == 0 ): anal=True
        npdate = bedate + datetime.timedelta(days=jfcst)
        nndiff = (2 - npdate.weekday() ) % 7
        nndate = npdate + datetime.timedelta(days=nndiff)
        nnal=False
        if ( nndiff == 0 ): nnal=True
        
        print( 'Persistence files', fcdate, andate)
        print( 'Nersistence files', npdate, nndate)
        if ( not nPersist ): 
            npdate=fcdate
            nndate=andate
            print( 'Nersistence files', npdate, nndate)

        CHOOSE='best'
        if ( CHOOSE == 'init' ):
            GU='gu'
            if ( fcdate <= date_ic2 ):  GU='pu'
            file_best, __ = find_fcst_file.find_anal_file(fcdate, system=GU, var='T',execute=True)
            file_besu, __ = find_fcst_file.find_anal_file(fcdate, system=GU, var='U',execute=True)
            file_besv, __ = find_fcst_file.find_anal_file(fcdate, system=GU, var='V',execute=True)
            GU='gu'
            if ( npdate <= date_ic2 ):  GU='pu'
            file_nest, __ = find_fcst_file.find_anal_file(npdate, system=GU, var='T',execute=True)
            file_nesu, __ = find_fcst_file.find_anal_file(npdate, system=GU, var='U',execute=True)
            file_nesv, __ = find_fcst_file.find_anal_file(npdate, system=GU, var='V',execute=True)
        elif ( CHOOSE == 'best' ):
            GD='gd'
            if ( andate <= date_ic2 - datetime.timedelta(days=7) ):  GD='pd'
            file_best, __ = find_fcst_file.find_anal_file(andate, system=GD, anal=anal, var='T',execute=True)
            file_besu, __ = find_fcst_file.find_anal_file(andate, system=GD, anal=anal, var='U',execute=True)
            file_besv, __ = find_fcst_file.find_anal_file(andate, system=GD, anal=anal, var='V',execute=True)
            GD='gd'
            if ( nndate <= date_ic2 - datetime.timedelta(days=7) ):  GD='pd'
            if ( nndate == andate ):
                file_nest = file_best
                file_nesu = file_besu
                file_nesv = file_besv
            else:
                file_nest, __ = find_fcst_file.find_anal_file(nndate, system=GD, anal=nnal, var='T',execute=True)
                file_nesu, __ = find_fcst_file.find_anal_file(nndate, system=GD, anal=nnal, var='U',execute=True)
                file_nesv, __ = find_fcst_file.find_anal_file(nndate, system=GD, anal=nnal, var='V',execute=True)
        # NEED ALL LEADS FOR DAY -- MIGHT AS WELL RETRIEVE AHEAD OF TIME.
        SYS='OPD'
        tmpdir=tempdir+'/'+SYS+'/'
        branch='operation.forecasts.giops.prog.glboce'
        if ( fcdate <= date_ic2 ): 
            SYS='PSD'
            tmpdir=tempdir+'/'+SYS+'/'
            branch='parallel.forecasts.giops.prog.glboce'
        rc = get_archive.get_archive_leads(tmpdir, branch, fcdate, fchours, ensnum=None, execute=True)
        # THIS SHOULD NOW BE FINDING A FILE ALREADY ON THE TMPDIR.
        file_fcst     = find_fcst_file.find_fcst_file(SYS, fcdate, fchour, 0, src='ocn', execute=True)
        filter_standard_file(template_dir+'/flt.ksh', file_fcst, file_fcst)
        # And SO SHOULD THIS NOW BE FINDING A FILE ALREADY ON THE TMPDIR.
        file_3hrs=[]
        for fc3hour in fchours:
            file_3hr=find_fcst_file.find_fcst_file(SYS, fcdate, fc3hour, 0, src='ocn', execute=True)
            file_3hrs.append(file_3hr)
        # Produce a cst interpolated grid
        ## Can't afford memory for 0.1deg grids -- try 0.2deg grids (should be 4 times less).
        file_intr = file_fcst+'.Z02'
        if ( os.path.isfile(file_intr) ):
            try:
                var, varf = stfd.file_query(file_intr, nomvar='UUW')
                if ( len(var) == 0 ): 
                    print('rm '+file_intr+' NO UUW')
                    subprocess.call('RM', file_intr)
            except:
                print('RM '+file_intr+' NOT STFD')
                subprocess.call(['rm', file_intr])
        if ( not os.path.isfile(file_intr) ):
            rc = cst_interpolation(file_fcst, file_intr, ref_grid=grid_dir+'dest_grid.std.2')
        print('file_best', file_best, os.path.isfile(file_best))
        print('file_nest', file_nest, os.path.isfile(file_nest))
        print('file_fcst', file_fcst, os.path.isfile(file_fcst))
        print('file_intr', file_intr, os.path.isfile(file_intr))
        for file_3hr in file_3hrs:
           print('file_3hr', file_3hr, os.path.isfile(file_3hr))  
    
        TM, LONT, LATT, LEVT, DAT = read_anal_file(file_best, grid='T',dates=[fcdate]); TM=np.squeeze(TM)
        UW, LONU, LATU, LEVU, DAT = read_anal_file(file_besu, grid='U',dates=[fcdate]); UW=np.squeeze(UW)
        VW, LONV, LATV, LEVV, DAT = read_anal_file(file_besv, grid='V',dates=[fcdate]); VW=np.squeeze(VW)
        UU, LONUU, LATUU = map_to_A_grid(np.squeeze(UW), LONU, LATU, grid='U')
        VV, LONVV, LATVV = map_to_A_grid(np.squeeze(VW), LONV, LATV, grid='V')
        # Important for U15 calculations that arrays not be masked.
        if ( isinstance(TM, np.ma.core.MaskedArray) ): TM = TM.data
        if ( isinstance(UU, np.ma.core.MaskedArray) ): UU = UU.data
        if ( isinstance(VV, np.ma.core.MaskedArray) ): VV = VV.data
    
        [T00, U00, V00] = [TM[0,:,:], UU[0,:,:], VV[0,:,:]]
        [T15, U15, V15] = calc_m15([TM, UU, VV], e3t, mask)
        [U00, V00] = UU_ORCA_to_NE([U00, V00])
        [U15, V15] = UU_ORCA_to_NE([U15, V15])

        # this actually removes non-ocean points
        (T00m, U00m, V00m, T15m, U15m, V15m), (LONM, LATM) = msk_flds([T00, U00, V00, T15, U15, V15], [LONN, LATN], mask0)
        # now grid to lat long
        (T00f, U00f, V00f, T15f, U15f, V15f), (LONF, LATF) = put_flds_latlon([T00m, U00m, V00m, T15m, U15m, V15m], [LONM, LATM], ddeg=0.2, method='2sweeplinear')
        # now grid to observation locations
        (TOOo, UOOo, VOOo, T15o, U15o, V15o) = interp_direct_obs( (LONM, LATM), [T00m, U00m, V00m, T15m, U15m, V15m], [LONO, LATO], method='2sweeplinear')

        if ( ndeps == 2 ):
            PERS.append([
                     [[T00m, U00m, V00m],[T15m, U15m, V15m]],
                     [[T00f, U00f, V00f],[T15f, U15f, V15f]],
                     [[T00o, U00o, V00o],[T15o, U15o, V15o]]
                        ])
        elif ( ndeps == 1 ):
            PERS.append([
                     [[T15m, U15m, V15m]],
                     [[T15f, U15f, V15f]],
                     [[T15o, U15o, V15o]]
                        ])
 
        TM, LONT, LATT, LEVT, DAT = read_anal_file(file_nest, grid='T',dates=[npdate]); TM=np.squeeze(TM)
        UW, LONU, LATU, LEVU, DAT = read_anal_file(file_nesu, grid='U',dates=[npdate]); UW=np.squeeze(UW)
        VW, LONV, LATV, LEVV, DAT = read_anal_file(file_nesv, grid='V',dates=[npdate]); VW=np.squeeze(VW)
        UU, LONUU, LATUU = map_to_A_grid(np.squeeze(UW), LONU, LATU, grid='U')
        VV, LONVV, LATVV = map_to_A_grid(np.squeeze(VW), LONV, LATV, grid='V')
        # Important for U15 calculations that arrays not be masked.
        if ( isinstance(TM, np.ma.core.MaskedArray) ): TM = TM.data
        if ( isinstance(UU, np.ma.core.MaskedArray) ): UU = UU.data
        if ( isinstance(VV, np.ma.core.MaskedArray) ): VV = VV.data
    
        [T00, U00, V00] = [TM[0,:,:], UU[0,:,:], VV[0,:,:]]
        [T15, U15, V15] = calc_m15([TM, UU, VV], e3t, mask)
        [U00, V00] = UU_ORCA_to_NE([U00, V00])
        [U15, V15] = UU_ORCA_to_NE([U15, V15])

        # this actually removes non-ocean points
        (T00m, U00m, V00m, T15m, U15m, V15m), (LONM, LATM) = msk_flds([T00, U00, V00, T15, U15, V15], [LONN, LATN], mask0)
        # now grid to lat long
        (T00f, U00f, V00f, T15f, U15f, V15f), (LONF, LATF) = put_flds_latlon([T00m, U00m, V00m, T15m, U15m, V15m], [LONM, LATM], ddeg=0.2, method='2sweeplinear')
        # now grid to observation locations
        (TOOo, UOOo, VOOo, T15o, U15o, V15o) = interp_direct_obs( (LONM, LATM), [T00m, U00m, V00m, T15m, U15m, V15m], [LONO, LATO], method='2sweeplinear')

        if ( ndeps == 2 ):
            NERS.append([
                     [[T00m, U00m, V00m],[T15m, U15m, V15m]],
                     [[T00f, U00f, V00f],[T15f, U15f, V15f]],
                     [[T00o, U00o, V00o],[T15o, U15o, V15o]]
                        ])
        elif ( ndeps == 1 ):
            NERS.append([
                     [[T15m, U15m, V15m]],
                     [[T15f, U15f, V15f]],
                     [[T15o, U15o, V15o]]
                        ])

        LONS, LATS, TF = stfd.read_fstd_var(file_fcst, 'TM', typvar='P@')
        leu, UF = stfd.read_fstd_multi_lev(file_fcst, 'UUW',vfreq=24, typvar='P@')
        lev, VF = stfd.read_fstd_multi_lev(file_fcst, 'VVW',vfreq=24, typvar='P@')
        let, TF = stfd.read_fstd_multi_lev(file_fcst, 'TM',vfreq=24, typvar='P@')
        # READ IN INTERPOLATED FILES TOO
        LONI, LATI, TI = stfd.read_fstd_var(file_intr, 'TM', typvar='P@')
        LATG, LONG = np.meshgrid(LATI, LONI)
        liu, UI = stfd.read_fstd_multi_lev(file_intr, 'UUW',vfreq=24, typvar='P@')
        liv, VI = stfd.read_fstd_multi_lev(file_intr, 'VVW',vfreq=24, typvar='P@')
        lit, TI = stfd.read_fstd_multi_lev(file_intr, 'TM',vfreq=24, typvar='P@')
        # Important for U15 calculations that arrays not be masked.
        ## BUT NEED THE INTERPOLATED MASK -- which isn't actually masked.
        if ( isinstance(TF, np.ma.core.MaskedArray) ): TF = TF.data
        if ( isinstance(UF, np.ma.core.MaskedArray) ): UF = UF.data
        if ( isinstance(VF, np.ma.core.MaskedArray) ): VF = VF.data
        MI, e3ti = fine_grid_mask(TI, e3t)  ## MI is a np.ones field for TF
        if ( isinstance(TI, np.ma.core.MaskedArray) ): TI = TI.data
        if ( isinstance(UI, np.ma.core.MaskedArray) ): UI = UI.data
        if ( isinstance(VI, np.ma.core.MaskedArray) ): VI = VI.data

        [T00, U00, V00] = [TF[0,:,:], UF[0,:,:], VF[0,:,:]]
        [T15, U15, V15] = calc_m15([TF, UF, VF], e3t, mask)
        # READ IN 3hr U15/V15 files too
        try:
            U15D, V15D = calc_u15_from_3hr(file_3hrs)
            print("SUCCESSFUL retrieval U15 from U10/20")
            T15D = T15
        except:
            U15D, V15D, T15D = U15, V15, T15
            print("UNSUCCESSFUL retrieval U15 from U10/20")
        try:
            U00D, V00D, T00D = calc_u00_from_3hr(file_3hrs, 0.0)
        except:
            U00D, V00D, T00D = U00, V00, T00
        try:
            U00d, V00d. T00d = calc_u00_from_3hr(file_3hrs, 15.0)
        except:
            U00d, V00d, T00d = U15, V15, T15
        S00, A00 = speed_and_angle(U00, V00)
        S15, A15 = speed_and_angle(U15, V15)

        [T00i, U00i, V00i] = [TI[0,:,:], UI[0,:,:], VI[0,:,:]]
        [T15i, U15i, V15i] = calc_m15([TI, UI, VI], e3ti, MI)
        S00i, A00i = speed_and_angle(U00i, V00i)
        S15i, A15i = speed_and_angle(U15i, V15i)
        print("Interpolated Shapes", T00i.shape, T15i.shape, LONG.shape, LATG.shape)
        
        if ( Plot == True ):
            suptitle='0m'+datestr+'Native Grid Lead '+str(ifcst+1)
            prefix='PLOTS/00_NG'+str(ifcst+1)+'_'+datestr+'_'
            plot_fields((T00f, U00f, V00f, S00f, A00f), (LONF, LATF), suptitle=suptitle, grid=False, outfile_prefix=prefix)
            suptitle='15m'+datestr+'Native Grid Lead '+str(ifcst+1)
            prefix='PLOTS/15_NG'+str(ifcst+1)+'_'+datestr+'_'
            plot_fields((T15f, U15f, V15f, S15f, A15f), (LONF, LATF), suptitle=suptitle, grid=False, outfile_prefix=prefix)

        ## ROTATE NATIVE GRID VECTORS TO EAST and NORTH
        [U00, V00] = UU_ORCA_to_NE([U00, V00])
        [U15, V15] = UU_ORCA_to_NE([U15, V15])
        [U00D, V00D] = UU_ORCA_to_NE([U00D, V00D])
        [U15D, V15D] = UU_ORCA_to_NE([U15D, V15D])

        # this actually removes non-ocean points
        (T00m, U00m, V00m, T15m, U15m, V15m), (LONM, LATM) = msk_flds([T00, U00, V00, T15, U15, V15], [LONN, LATN], mask0)
        (T00n, U00n, V00n, T15n, U15n, V15n), (LONM, LATM) = msk_flds([T00D, U00D, V00D, T15D, U15D, V15D], [LONN, LATN], mask0)
        # now grid to lat long
        (T00f, U00f, V00f, T15f, U15f, V15f), (LONF, LATF) = put_flds_latlon([T00m, U00m, V00m, T15m, U15m, V15m], [LONM, LATM], ddeg=0.2, method='2sweeplinear')
        (T00g, U00g, V00g, T15g, U15g, V15g), (LONF, LATF) = put_flds_latlon([T00n, U00n, V00n, T15n, U15n, V15n], [LONM, LATM], ddeg=0.2, method='2sweeplinear')
        # now grid to observation locations
        (TOOo, UOOo, VOOo, T15o, U15o, V15o) = interp_direct_obs( (LONM, LATM), [T00m, U00m, V00m, T15m, U15m, V15m], [LONO, LATO], method='2sweeplinear')
        (TOOp, UOOp, VOOp, T15p, U15p, V15p) = interp_direct_obs( (LONM, LATM), [T00n, U00n, V00n, T15n, U15n, V15n], [LONO, LATO], method='2sweeplinear')
        print('OBSERVATION GRID SIZES', T15p.shape, U15p.shape, V15p.shape)
        S00f, A00f = speed_and_angle(U00f, V00f)
        S15f, A15f = speed_and_angle(U15f, V15f)
        S00g, A00g = speed_and_angle(U00g, V00g)
        S15g, A15g = speed_and_angle(U15g, V15g)

        if ( Plot == True ):
            suptitle='0m'+datestr+'Forecast Lead '+str(ifcst+1)
            prefix='PLOTS/00_FC'+str(ifcst+1)+'_'+datestr+'_'
            plot_fields((T00f, U00f, V00f, S00f, A00f), (LONF, LATF), suptitle=suptitle, grid=False, outfile_prefix=prefix)
            suptitle='15m'+datestr+'Forecast Lead '+str(ifcst+1)
            prefix='PLOTS/15_FC'+str(ifcst+1)+'_'+datestr+'_'
            plot_fields((T15f, U15f, V15f, S15f, A15f), (LONF, LATF), suptitle=suptitle, grid=False, outfile_prefix=prefix)
            suptitle='0m'+datestr+'cstintr Forecast Lead '+str(ifcst+1)
            prefix='PLOTS/00_FI'+str(ifcst+1)+'_'+datestr+'_'
            plot_fields((T00i, U00i, V00i, S00i, A00i), (LONG, LATG), suptitle=suptitle, grid=False, outfile_prefix=prefix)
            suptitle='15m'+datestr+'cstintr Forecast Lead '+str(ifcst+1)
            prefix='PLOTS/15_FI'+str(ifcst+1)+'_'+datestr+'_'
            plot_fields((T15i, U15i, V15i, S15i, A15i), (LONG, LATG), suptitle=suptitle, grid=False, outfile_prefix=prefix)

        if ( ndeps == 2 ):
            FCST.append([
                     [[T00m, U00m, V00m],[T15,  U15,  V15]],
                     [[T00f, U00f, V00f],[T15f, U15f, V15f]],
                     [[T00i, U00i, V00i],[T15i, U15i, V15i]],
                     [[T00g, U00g, V00g],[T15g, U15g, V15g]],
                     [[T00o, UOOo, VOOo],[T15o, U15o, V15o]],
                     [[T00p, U00p, V00p],[T15p, U15p, V15p]]
                        ])
        elif ( ndeps == 1 ):
            FCST.append([
                     [[T15m, U15m, V15m]],
                     [[T15f, U15f, V15f]],
                     [[T15i, U15i, V15i]],
                     [[T15g, U15g, V15g]],
                     [[T15o, U15o, V15o]],
                     [[T15p, U15p, V15p]]
                     
                        ])
        timen=time.time()
        timep=timen-timel
        timet=timen-time0
        print("TIMING :: LEAD "+str(ifcst+1)+" Processing Time: "+str(timep))
        print("TIMING :: Cumalative Forecast Processing and Retrieval Time "+str(timet))
        timel=timen

    (LONN, LATN) = (nav_lon, nav_lat)
    time0=time.time()
    timen=time0
    for iobs in range(nobss_loop):

        [LONP, LATP] = LONO[iobs], LATO[iobs]
        IJPTS = find_nearest_points((LONP, LATP), [(LONM, LATM), (LONF, LATF), (LONG, LATG), (LONF, LATF) ])
        IJPT, IJPF, IJPI, IJPG = IJPTS
        print('IJPTS', IJPTS)

        for igrid, GRID in enumerate(BEST[:2]):
            beste = BESTE_LIST[igrid]
            tempb = TEMPB_LIST[igrid]
            IJPT = IJPTS[igrid]
            
            for kk, KLEV in enumerate(GRID):
                bestl = []
                for ifld, FLD in enumerate(KLEV):
                    if ( IJPT[1] != 1 ):
                         FLP=FLD[IJPT]
                    else: 
                        if ( FLD.ndim == 1 ):
                          FLP=FLD[IJPT[0]]
                        else:
                          print("Error:  Mismatch in points and field shapes")
                          print( igrid, kk, ifld, FLD.shape, IJPT ) 
                          FLP=missing
                    bestl.append(FLP)
                  
                Tpi, Upi, Vpi = bestl   ## These are the best values for level kk
                ## Repeatedly fill in the best values.  Initially Orca grid, then orca grid masked, then lat lon grid.
                beste[iobs, 0, kk] = Upi
                beste[iobs, 1, kk] = Vpi
                tempb[iobs, kk] = Tpi
                print (iobs, igrid, kk, 'Best Values', Upi, Vpi, Tpi, obser[iobs, :, kk] )
            BESTE_LIST[igrid] = beste
            TEMPB_LIST[igrid] = tempb 
        # no 3rd IGRID.  COPY LAST
        BESTE_LIST[2] = BESTE_LIST[1]
        TEMPB_LIST[2] = TEMPB_LIST[1] 
        # no 4th IGRID.  COPY LAST
        BESTE_LIST[3] = BESTE_LIST[1]
        TEMPB_LIST[3] = TEMPB_LIST[1] 
        
        for igrid, GRID in enumerate(INIT[:2]):
            inite = INITE_LIST[igrid]
            tempi = TEMPI_LIST[igrid]
            IJPT = IJPTS[igrid]
            
            for kk, KLEV in enumerate(GRID):
                initl = []
                for ifld, FLD in enumerate(KLEV):
                    if ( IJPT[1] != 1 ):
                        FLP=FLD[IJPT]
                    else: 
                        if ( FLD.ndim == 1 ):
                          FLP=FLD[IJPT[0]]
                        else:
                          print ("Error:  Mismatch in points and field shapes")
                          print( igrid, kk, ifld, FLD.shape, IJPT ) 
                          FLP=missing
                    initl.append(FLP)
                  
                Tpi, Upi, Vpi = initl   ## These are the best values for level kk
                ## Repeatedly fill in the best values.  Initially Orca grid, then orca grid masked, then lat lon grid.
                inite[iobs, 0, kk] = Upi
                inite[iobs, 1, kk] = Vpi
                tempi[iobs, kk] = Tpi
                print (iobs, igrid, kk, 'Init Values', Upi, Vpi, Tpi, obser[iobs, :, kk] )
            INITE_LIST[igrid] = inite
            TEMPI_LIST[igrid] = tempi 
        
        for ip, PRED in enumerate([FCST, PERS, NERS]): 
          for ld, LEAD in enumerate(PRED):
            #print('LEN LEAD', len(LEAD), ld, ip)
            if ( ip == 0 ): LEAD=LEAD[:4]
            if ( ip == 1 ): LEAD=LEAD[:2]
            if ( ip == 2 ): LEAD=LEAD[:2]
            for ig, GRID in enumerate(LEAD): 
              IJPT = IJPTS[ig]
              if ( ip == 0 ):
                fcstv = FCSTV_LIST[ig]
                tempf = TEMPF_LIST[ig]
              if ( ip == 1 ):
                persi = PERSI_LIST[ig]
                tempp = TEMPP_LIST[ig]
              if ( ip == 2 ):
                nersi = NERSI_LIST[ig]
                tempn = TEMPN_LIST[ig]
              for kk, KLEV in enumerate(GRID):
                bestl = []
                for ifld, FLD in enumerate(KLEV):
                  if (  IJPT[1] != 1 ):
                    FLP=FLD[IJPT]
                  else: 
                    if ( FLD.ndim == 1 ):
                      FLP=FLD[IJPT[0]]
                    else:
                      print ("Error:  Mismatch in points and field shapes")
                      print( ip, ld, ig, kk, ifld, FLD.shape, IJPT ) 
                      FLP=missing
                  bestl.append(FLP)
                  
                Tpi, Upi, Vpi = bestl   ## These are the best values for level kk
                ## Repeatedly fill in the best values.  Initially Orca grid, then orca grid masked, then lat lon grid.
                #print (ip, ld, ig, kk, 'Best Values', Upi, Vpi, Tpi)
                if ( ip == 0 ):
                    fcstv[iobs, 0, ld, kk] = Upi
                    fcstv[iobs, 1, ld, kk] = Vpi
                    tempf[iobs, ld, kk] = Tpi
                elif ( ip == 1 ):
                    persi[iobs, 0, ld, kk] = Upi
                    persi[iobs, 1, ld, kk] = Vpi
                    tempp[iobs, ld, kk] = Tpi
                elif ( ip == 2 ):
                    nersi[iobs, 0, ld, kk] = Upi
                    nersi[iobs, 1, ld, kk] = Vpi
                    tempn[iobs, ld, kk] = Tpi
              if ( ip == 0 ):
                FCSTV_LIST[ig] = fcstv
                TEMPF_LIST[ig] = tempf
              if ( ip == 1 ):
                PERSI_LIST[ig] = persi
                TEMPP_LIST[ig] = tempp
              if ( ip == 2 ):
                NERSI_LIST[ig] = nersi
                TEMPN_LIST[ig] = tempn
 
    # no 3rd IGRID.  COPY LAST
    INITE_LIST[2] = INITE_LIST[1]
    TEMPI_LIST[2] = TEMPI_LIST[1] 
    # no 4th IGRID.  COPY LAST
    INITE_LIST[3] = INITE_LIST[1]
    TEMPI_LIST[3] = TEMPI_LIST[1] 
    # no 3rd IGRID for persi/nersi.  COPY LAST
    PERSI_LIST[2] = PERSI_LIST[1]
    TEMPP_LIST[2] = TEMPP_LIST[1]
    NERSI_LIST[2] = NERSI_LIST[1]
    TEMPN_LIST[2] = TEMPN_LIST[1]
    # no 4th IGRID for persi/nersi.  COPY LAST
    PERSI_LIST[3] = PERSI_LIST[1]
    TEMPP_LIST[3] = TEMPP_LIST[1]
    NERSI_LIST[3] = NERSI_LIST[1]
    TEMPN_LIST[3] = TEMPN_LIST[1]
    timen = time.time()
    timet = timen-time0
    print('TIMING :: Total Observation Processing Time: '+str(timet)+' seconds')
    for ig in range(6):
        if ( np.all(BESTE_LIST[ig] == BESTE_LIST[0] )):
            print('SAME I', ig, 0)

    time0 = time.time()

    tempb = TEMPB_LIST[4]
    beste = BESTE_LIST[4]
    for kk, KLEV in enumerate(BEST[2]):
        tempb[:, kk] = BEST[2][kk][0]
        beste[:, 0, kk] = BEST[2][kk][1]
        beste[:, 1, kk] = BEST[2][kk][2]
    TEMPB_LIST[4] = tempb
    TEMPB_LIST[5] = tempb
    BESTE_LIST[4] = beste
    BESTE_LIST[5] = beste

    tempi = TEMPI_LIST[4]
    inite = INITE_LIST[4]

    for kk, KLEI in enumerate(INIT[2]):
        tempi[:, kk] = INIT[2][kk][0]
        inite[:, 0, kk] = INIT[2][kk][1]
        inite[:, 1, kk] = INIT[2][kk][2]
    TEMPI_LIST[4] = tempi
    TEMPI_LIST[5] = tempi
    INITE_LIST[4] = beste
    INITE_LIST[5] = beste
            
    for ig in range(4,6):
        print('FAST', ig)
        fcstv = FCSTV_LIST[ig]
        tempf = TEMPF_LIST[ig]
        persv = PERSI_LIST[ig]
        tempp = TEMPP_LIST[ig]
        nersv = NERSI_LIST[ig]
        tempn = TEMPN_LIST[ig]
        for ld in range(nfcst):
            for kk in range(ndeps):
                print('SIZES', ig, ld, kk, tempf[:,ld,kk].shape, FCST[ld][ig][kk][0].shape)
                print('SIZES', ig, ld, kk, fcstv[:,: ld,kk].shape, FCST[ld][ig][kk][1].shape, FCST[ld][ig][kk][2].shape)
                tempf[:, ld, kk] = FCST[ld][ig][kk][0]
                fcstv[:, 0, ld, kk] = FCST[ld][ig][kk][1]
                fcstv[:, 1, ld, kk] = FCST[ld][ig][kk][2]
                print('SIZES', ig, ld, kk, tempp[:,ld,kk].shape, len(PERS[ld][2][kk][0]))
                print('SIZES', ig, ld, kk, persv[:,: ld,kk].shape, len(PERS[ld][2][kk][1]), len(PERS[ld][2][kk][2]))
                tempp[:, ld, kk] = PERS[ld][2][kk][0]
                persv[:, 0, ld, kk] = PERS[ld][2][kk][1]
                persv[:, 1, ld, kk] = PERS[ld][2][kk][2]
                print('SIZES', ig, ld, kk, tempn[:,ld,kk].shape, len(NERS[ld][2][kk][0]))
                print('SIZES', ig, ld, kk, nersv[:,: ld,kk].shape, len(NERS[ld][2][kk][1]), len(NERS[ld][2][kk][2]))
                tempn[:, ld, kk] = NERS[ld][2][kk][0]
                nersv[:, 0, ld, kk] = NERS[ld][2][kk][1]
                nersv[:, 1, ld, kk] = NERS[ld][2][kk][2]
        TEMPF_LIST[ig] = tempf
        FCSTV_LIST[ig] = fcstv
        TEMPP_LIST[ig] = tempn
        PERSI_LIST[ig] = persv
        TEMPN_LIST[ig] = tempn
        NERSI_LIST[ig] = nersv
        
    timen = time.time()
    timet = timen-time0
    print('TIMING  :: Quick Observation Processing time: '+str(timet)+' seconds')
    for ig in range(6):
        if ( np.all(beste == BESTE_LIST[ig]) ):
            print('NO CHANGE -- WHY', ig)
        if ( np.all(BESTE_LIST[ig] == BESTE_LIST[0] )):
            print('SAME E', ig, 0)
        beste = BESTE_LIST[ig]
        inite = INITE_LIST[ig]
        fcstv = FCSTV_LIST[ig]
        persi = PERSI_LIST[ig]
        nersi = NERSI_LIST[ig]
        tempb = TEMPB_LIST[ig]
        tempi = TEMPI_LIST[ig]
        tempf = TEMPF_LIST[ig]
        tempp = TEMPP_LIST[ig]
        tempn = TEMPN_LIST[ig]
        
        print(ig, '0 best', beste[0, :, 0])
        print(ig, 'E best', beste[-1, :, 0])
                
        Mobser = speed_and_angle_easy(obser)
        Mbeste = speed_and_angle_easy(beste)
        Minite = speed_and_angle_easy(inite)
        Mfcstv = speed_and_angle_easy(fcstv)
        Mpersi = speed_and_angle_easy(persi)
        Mnersi = speed_and_angle_easy(nersi)
        
        CCMEP_MEAN_ERRORS = calc_mean_error(obser, (beste, inite, fcstv, persi, nersi))
        mr_ERRCb, mr_ERRCi, mr_ERRCf, mr_ERRCp, mr_ERRCn = CCMEP_MEAN_ERRORS
        mean_ERRCb, rmse_ERRCb, abse_ERRCb = mr_ERRCb
        mean_ERRCi, rmse_ERRCi, abse_ERRCi = mr_ERRCi
        mean_ERRCf, rmse_ERRCf, abse_ERRCf = mr_ERRCf
        mean_ERRCp, rmse_ERRCp, abse_ERRCp = mr_ERRCp
        mean_ERRCn, rmse_ERRCn, abse_ERRCn = mr_ERRCn
        CCMEP_MMEAN_ERRORS = calc_mean_error(Mobser, (Mbeste, Minite, Mfcstv, Mpersi, Mnersi), isangle=1)
        mr_ERMCb, mr_ERMCi, mr_ERMCf, mr_ERMCp, mr_ERMCn = CCMEP_MMEAN_ERRORS
        mean_ERMCb, rmse_ERMCb, abse_ERMCb = mr_ERMCb
        mean_ERMCi, rmse_ERMCi, abse_ERMCi = mr_ERMCi
        mean_ERMCf, rmse_ERMCf, abse_ERMCf = mr_ERMCf
        mean_ERMCp, rmse_ERMCp, abse_ERMCp = mr_ERMCp
        mean_ERMCn, rmse_ERMCn, abse_ERMCn = mr_ERMCn

        print(datestr+' RESULTS')
    
        #print('Best Error Mercator', mean_ERRMb, rmse_ERRMb, mean_ERMMb, rmse_ERMMb)
        print(ig, 'Best Error CCMEP', mean_ERRCb, rmse_ERRCb, mean_ERMCb, rmse_ERMCb)
        print(ig, 'Init Error CCMEP', mean_ERRCi, rmse_ERRCi, mean_ERMCi, rmse_ERMCi)
    
        #print('Forecast Error Mercator', mean_ERRMf, rmse_ERRMf, mean_ERMMf, rmse_ERMMf)
        #print(ig, 'Forecast Error CCMEP', mean_ERRCf, rmse_ERRCf, mean_ERMCf, rmse_ERMCf)

        #print('Persistence Error Mercator', mean_ERRMp, rmse_ERRMp, mean_ERMMp, rmse_ERMMp)
        #print(ig, 'Persistence Error CCMEP', mean_ERRCp, rmse_ERRCp, mean_ERMCp, rmse_ERMCp)
    
        #print(ig, 'Negative Persistence Error CCMEP', mean_ERRCn, rmse_ERRCn, mean_ERMCn, rmse_ERMCn)

        if ( not TEST_SINGLE ):  ## DON"T DO FOR TEST   
            if ( ig == 0 ) : obsfile_ig=obsfile1
            if ( ig == 1 ) : obsfile_ig=obsfile2
            if ( ig == 2 ) : obsfile_ig=obsfile3
            if ( ig == 3 ) : obsfile_ig=obsfile4
            if ( ig == 4 ) : obsfile_ig=obsfile5
            if ( ig == 5 ) : obsfile_ig=obsfile6
            write_model_obsfile_plus(obsfile_ig, psyfile, (beste, inite, fcstv, persi, nersi))
    
        if ( not TEST_SINGLE ):  ## DON"T DO FOR TEST  

            if ( ig == 0 ) : add='.1'
            if ( ig == 1 ) : add='.2'
            if ( ig == 2 ) : add='.3'
            if ( ig == 3 ) : add='.4'
            if ( ig == 4 ) : add='.5'
            if ( ig == 5 ) : add='.6'

            outpre='ERRORS/GIOPS'
            tmppre='GIOPS'
            outdep=1
            if ( CHARLY ):
              outpre='ERRORS_UFIL/GIOPS'
              tmppre='uGIOPS'
              outdep=0
              if ( filter ): 
                outpre='ERRORS_FILT/GIOPS'
                tmppre='fGIOPS'

            write_mean_errors(outpre+add, int(datestr),  CCMEP_MEAN_ERRORS, vec=['u','v'], tmp_prefix='site6/Class4_Currents/TMP/'+tmppre+add, depth_index=outdep)
            write_mean_errors(outpre+add, int(datestr), CCMEP_MMEAN_ERRORS, vec=['s','a'], tmp_prefix='site6/Class4_Currents/TMP/'+tmppre+add, depth_index=outdep)
    
    return        

def get_geps_files(date):
    date=check_date(date, outtype=datetime.datetime)
    nfcst=16
    nenss=21
    ens_list=range(nenss)
    bedate=date + datetime.timedelta(days=1)
    SYS='OP'
    branch='operation.ensemble.prog.ens.glboce'
    tmpdir=tempdir+'/'+SYS+'/'+'O'+'/'
    missing=[]
    for ifcst in range(nfcst):
        jfcst=ifcst+1
        fchour=jfcst*24
        fcdate = bedate - datetime.timedelta(days=jfcst)
        target=tmpdir+'/'+branch+'/'+fcdate.strftime('%Y%m%d%H')+'_'+str(fchour).zfill(3)+'_???'
        if ( len(glob.glob(target)) < nenss ): 
            rc = get_archive.get_archive(tmpdir, branch, fcdate, fchour, ensnum=list(range(nenss)), execute=True) 
        for iensm in ens_list:
            # THIS SHOULD NOW BE FINDING A FILE ALREADY ON THE TMPDIR.
            file_fcst     = find_fcst_file.find_fcst_file(SYS, fcdate, fchour, iensm, src='ocn', execute=False)
            if ( not os.path.isfile(file_fcst) ): missing.append(file_fcst)
    print("MISSING", missing)
    return
                
def process_geps_obs(date=tate, ens_list=[0], Plot=False, filter=True, DO_CLEAN=True):

    nfcst=16
    #nenss=21
    nenss=len(ens_list)
    ens_str=str(ens_list[0]).zfill(3)
    
    datestr=date.strftime('%Y%m%d')   

    if ( not filter ): 
        obsfile1='CLASS4_currents_CCMEP_UFIL/class4_'+datestr+'_GIOPS_orca025_currents.1.nc'
        obsfile2='CLASS4_currents_CCMEP_UFIL/class4_'+datestr+'_GIOPS_orca025_currents.2.nc'
        oesfile1='CLASS4_currents_GEPS_UFIL/class4_'+datestr+'_GEPS_orca025_currents.1.'+ens_str+'.nc'
        oesfile2='CLASS4_currents_GEPS_UFIL/class4_'+datestr+'_GEPS_orca025_currents.2.'+ens_str+'.nc'
    elif ( filter ): 
        obsfile1='CLASS4_currents_CCMEP_FILT/class4_'+datestr+'_GIOPS_orca025_currents.f1.nc'
        obsfile2='CLASS4_currents_CCMEP_FILT/class4_'+datestr+'_GIOPS_orca025_currents.f2.nc'
        oesfile1='CLASS4_currents_GEPS_FILT/class4_'+datestr+'_GEPS_orca025_currents.f1.'+ens_str+'.nc'
        oesfile2='CLASS4_currents_GEPS_FILT/class4_'+datestr+'_GEPS_orca025_currents.f2.'+ens_str+'.nc'
    OBSFILE_LIST = [obsfile1, obsfile2]
    OESFILE_LIST = [oesfile1, oesfile2] 

    FCSTV_LIST = []
    print(obsfile1)
    (LONO, LATO, depth), (obser, beste, inite, fcstv, persi, nersi) = read_obsfile_plus(obsfile1)
    nobss, nvars, ofcst, ndeps = fcstv.shape
    fcstv = np.zeros((nobss, nvars, nfcst, nenss, ndeps))
    fcstv = np.ma.array(fcstv)
    
    for obsfile in [obsfile1, obsfile2]:
        FCSTV_LIST.append(fcstv.copy())

    (LONN, LATN) = (nav_lon, nav_lat)
    bedate=date + datetime.timedelta(days=1)
    FCST=[]

    SYS='OP'
    tmpdir=tempdir+'/'+SYS+'/'+'O'+'/'
    branch='operation.ensemble.prog.ens.glboce'
    # GET FILES NEEDED FOR TOMORROWS PROCESSING
    time0a = time.time()
    target=tmpdir+'/'+branch+'/'+bedate.strftime('%Y%m%d%H')+'_???_???'
    if ( ens_list[0] == 0 ):
      if ( len(glob.glob(target)) < (nenss*nfcst) ): 
        rc = get_archive.get_archive_leads(tmpdir, branch, bedate, (np.arange(0,384,24)+24).tolist(), ensnum=list(range(21)), execute=True)  
    timea = time.time() - time0a
    print("DEARCHIVING TIME ", bedate, timea)

    SYS='OP'
    missing=[]
    for ifcst in range(nfcst):
        jfcst=ifcst+1
        fchour=jfcst*24
        fcdate = bedate - datetime.timedelta(days=jfcst)
        for iensm in ens_list:
            file_fcst     = find_fcst_file.find_fcst_file(SYS, fcdate, fchour, iensm, src='ocn', execute=False)
            if ( not os.path.isfile(file_fcst) ): missing.append(file_fcst)
    print("MISSING", missing)
        
    time0f = time.time()
    for ifcst in range(nfcst):
      time0f = time.time()
      jfcst=ifcst+1
      fchour=jfcst*24
      fcdate = bedate - datetime.timedelta(days=jfcst)
      SYS='OP'
      tmpdir=tempdir+'/'+SYS+'/'+'O'+'/'
      branch='operation.ensemble.prog.ens.glboce'
      for iensm in ens_list:
          target=tmpdir+'/'+branch+'/'+fcdate.strftime('%Y%m%d%H')+'_'+str(fchour).zfill(3)+'_'+str(iensm).zfill(3)
          if ( not os.path.isfile(target) ):
              rc = get_archive.get_archive(tmpdir, branch, fcdate, fchour, ensnum=iensm, execute=True) 
      LEAD=[]
      time0e = time.time()
      for iensm in ens_list:
        # THIS SHOULD NOW BE FINDING A FILE ALREADY ON THE TMPDIR.
        file_fcst     = find_fcst_file.find_fcst_file(SYS, fcdate, fchour, iensm, src='ocn', execute=True)
        print('file_fcst', file_fcst, os.path.isfile(file_fcst))
        # Read in Forecast File
        LONS, LATS, TF = stfd.read_fstd_var(file_fcst, 'TM', typvar='P@')
        #leu, UF = stfd.read_fstd_multi_lev(file_fcst, 'UU2W',vfreq=24, typvar='P@')
        #lev, VF = stfd.read_fstd_multi_lev(file_fcst, 'VV2W',vfreq=24, typvar='P@')
        if ( fcdate < date_mg1 ): 
            print('FIELDS NOT AVAILABLE BEFORE', date_mg1)
            return 99
        newnorm=True
        if ( fcdate < date_ic3 ): newnorm=False
        (U15, V15) = calc_u15_from_3hr([file_fcst], vfreq=24, newnorm=newnorm)
        if ( np.max(U15) > 20.0 ):  
            print("LOOKS LIKE NORMALIZATION WRONG, newnorm")
            newnorm=False
            print("USE newnorm = ", newnorm)
            print("DATE", fcdate)
            (U15, V15) = calc_u15_from_3hr([file_fcst], vfreq=24, newnorm=False)
        if ( np.max(U15) < 0.05 ):  
            print("LOOKS LIKE NORMALIZATION WRONG, newnorm")
            newnorm=True
            print("USE newnorm = ", newnorm)
            print("DATE", fcdate)
            (U15, V15) = calc_u15_from_3hr([file_fcst], vfreq=24, newnorm=True)
        [U15F, V15F] = UU_ORCA_to_NE([U15, V15])
        
        # this actually removes non-ocean points
        (U15m, V15m), (LONM, LATM) = msk_flds([U15F, V15F], [LONN, LATN], mask0)
        # now grid to lat long
        (U15f, V15f), (LONF, LATF) = put_flds_latlon([U15m, V15m], [LONM, LATM], ddeg=0.2, method='2sweeplinear')
        ENSM = [ [[U15F, V15F]], [[U15f, V15f]] ]
        LEAD.append(ENSM)   ## EXTRA BRACKETS ARE HISTORICAL BAGGAGE.
        timee = time.time() - time0e
        print("PROCESSING TIME for ENSM ", iensm, timee)
      FCST.append(LEAD)   ## INDICES:  nfcst=16, nens=21, ngrid=2, ideps=1, nflds=2
      timef = time.time() - time0f
      print("PROCESSING TIME for LEAD ", ifcst, timef)
        
    time0o = time.time()
    for iobs in range(nobss):
        time0o = time.time()
        [LONP, LATP] = LONO[iobs], LATO[iobs]
        IJPTS = find_nearest_points((LONP, LATP), [(LONN, LATN), (LONF, LATF)])
        IJPT, IJPF = IJPTS
        print('IJPTS', IJPTS)

        for ld, LEAD in enumerate(FCST): 
          for ie, ENSM in enumerate(LEAD):
            for ig, GRID in enumerate(ENSM):
              IJPT = IJPTS[ig]
              fcstv = FCSTV_LIST[ig].copy()
              for kk, KLEV in enumerate(GRID):
                fcstl = []
                for ifld, FLD in enumerate(KLEV):
                    FLP=FLD[IJPT]
                    fcstl.append(FLP)
                Upi, Vpi = fcstl   ## These are the best values for level kk, grid ig, ensemble ie, lead ld
                fcstv[iobs, 0, ld, ie, kk] = Upi
                fcstv[iobs, 1, ld, ie, kk] = Vpi
              FCSTV_LIST[ig] = fcstv.copy()
    timeo = time.time() - time0o
    print("PROCESSING TIME per observation ",timeo, nobss, timeo/nobss)

    for ig in range(len(FCSTV_LIST)):
        fcstv = FCSTV_LIST[ig].copy()
        obsfile = OBSFILE_LIST[ig]
        oesfile = OESFILE_LIST[ig]
        print(obsfile, oesfile, fcstv.shape)
        rc = write_model_obsfile_ensemble(oesfile, obsfile, fcstv, clobber=True)

    for ifcst in range(nfcst):
      jfcst=ifcst+1
      fchour=jfcst*24
      fcdate = bedate - datetime.timedelta(days=jfcst)
      SYS='OP'
      tmpdir=tempdir+'/'+SYS+'/'+'O'+'/'
      branch='operation.ensemble.prog.ens.glboce'
      LEAD=[]
      for iensm in enumerate(ens_list):
        # THIS FILES SHOULD NO LONGER BE NEEDED
        file_fcst     = find_fcst_file.find_fcst_file(SYS, fcdate, fchour, iensm, src='ocn', execute=False)
        if ( os.path.isfile(file_fcst) and DO_CLEAN ):
            subprocess.call('rm', file_fcst)

    return 0
        
def read_obsfile_variable(obsfile, variable):
    # READ IN OBS FILE
    obsset = netCDF4.Dataset(obsfile,mode='r')
    LONO=obsset.variables['longitude'][:]
    LATO=obsset.variables['latitude'][:]
    depth=obsset['depth'][:]

    field=obsset[variable][:]
    obsset.close()

    return (LONO, LATO, depth), field


def read_obsfile(obsfile):
    # READ IN OBS FILE
    obsset = netCDF4.Dataset(obsfile,mode='r')
    LONO=obsset.variables['longitude'][:]
    LATO=obsset.variables['latitude'][:]

    depth=obsset['depth'][:]
    obser=obsset['observation'][:]
    beste=obsset['best_estimate'][:]
    persi=obsset['persistence'][:]
    fcstv=obsset['forecast'][:]

    obsset.close()

    return (LONO, LATO, depth), (obser, beste, fcstv, persi)    

def read_obsfile_variables(obsfile, variables=['observation', 'best_estimate']):

    # READ IN OBS FILE
    obsset = netCDF4.Dataset(obsfile,mode='r')
    LONO=obsset.variables['longitude'][:]
    LATO=obsset.variables['latitude'][:]
    depth=obsset['depth'][:]
    
    FLDLIST=[]
    for variable in variables:
        field=obsset[variable][:]
        FLDLIST.append(field)
    obsset.close()

    return (LONO, LATO, depth), FLDLIST
        
def read_obsfile_plus(obsfile):
    # READ IN OBS FILE
    obsset = netCDF4.Dataset(obsfile,mode='r')
    LONO=obsset.variables['longitude'][:]
    LATO=obsset.variables['latitude'][:]

    depth=obsset['depth'][:]
    obser=obsset['observation'][:]
    beste=obsset['best_estimate'][:]
    fcstv=obsset['forecast'][:]
    if ( 'persistence'in obsset.variables.keys() ):
        persi=obsset['persistence'][:]
    else:
        persi=fcstv.copy()
    if ( 'negative_persistence' in obsset.variables.keys() ):
       nersi = obsset['negative_persistence'][:]
    else:
       nersi = persi.copy()
    if ( 'init_estimate' in obsset.variables.keys() ):
       inite = obsset['init_estimate'][:]
    else:
       inite = beste.copy()

    obsset.close()

    return (LONO, LATO, depth), (obser, beste, inite, fcstv, persi, nersi)    

stoke_variables=['tides', 'stokes_best']

def read_stokes_variables(obsfile):
    (LONO, LATO, depth), (tidev, stokv) = read_obsfile_variables(obsfile, variables=stoke_variables)
    return tidev, stokv
    
def read_obsfile_allt(obsfile):
    # READ IN OBS FILE
    obsset = netCDF4.Dataset(obsfile,mode='r')
    LONO=obsset.variables['longitude'][:]
    LATO=obsset.variables['latitude'][:]

    depth=obsset['depth'][:]
    obser=obsset['observation'][:]
    beste=obsset['best_estimate'][:]
    if ( 'init_estimate' in obsset.variables.keys() ):
        inite=obsset['init_estimate'][:]
    else:
        inite=beste.copy()
    fcstv=obsset['forecast'][:]
    if ( 'persistence' in obsset.variables.keys() ):
        persi=obsset['persistence'][:]
    else:
        persi=fcstv.copy()
    (nobs, nvars, nfcst, ndeps) = fcstv.shape
    if ( 'smoc_drift_best' in obsset.variables.keys() ):
        smocv=obsset['smoc_drift_best'][:]
    else:
        smocv=obsset['smoc_drift'][:]
    if ( 'tides' in obsset.variables.keys() ):
        tidev=obsset['tides'][:]
    else:
        tidev= obsset['tide_current'][:]
    if ( 'stokes_best' in obsset.variables.keys() ):
        stokv=obsset['stokes_best'][:]
    else:
        stokv=obsset['stokes_drift'][:]
    if ( 'bathymetry' in obsset.variables.keys() ):
        bathy=obsset['bathymetry'][:]
    else:
        bathy = np.ma.array(np.zeros((nobs, ndeps)), np.ones((nobs, ndeps)))
    
 
    depth_qc = obsset['depth_qc'][:]
    otime_qc = obsset['time_qc'][:]
    posit_qc = obsset['position_qc'][:]
 
    if ( 'observation_qc' in obsset.variables.keys() ):
        obser_qc = obsset['observation_qc'][:]
    else:
        obser_qc = obsset['obs_qc']

    if ( 'plateform_code' in obsset.variables.keys() ):        
        plf_code = obsset['plateform_code'][:]
    else:
        plf_code = obsset['dc_reference'][:]

    obs_time = obsset['obs_time'][:]
    
    obsset.close()

    return (LONO, LATO, depth, plf_code), (obser, beste, inite, fcstv, persi, smocv, stokv, tidev, bathy), (depth_qc, otime_qc, posit_qc, obser_qc)

def write_model_obsfile(obsfile, tplfile, fields):
    (beste, fcstv, persi) = fields
    # cp template file to new obsfile.  Then enter model data.
    shutil.copy(tplfile, obsfile)
    obsset = netCDF4.Dataset(obsfile,mode='r+')
    obsset['best_estimate'][:]=beste
    obsset['forecast'][:]=fcstv
    obsset['persistence'][:]=persi
    obsset.close()
    return

def write_model_obsfile_plus(obsfile, tplfile, fields):
    (beste, inite, fcstv, persi, nersi)=fields
    # cp template file to new obsfile.  Then enter model data.
    shutil.copy(tplfile, obsfile)
    obsset = netCDF4.Dataset(obsfile,mode='r+')
    obsset['best_estimate'][:]=beste
    obsset['forecast'][:]=fcstv
    obsset['persistence'][:]=persi
    inite_var = obsset.createVariable('init_estimate', np.float32,
                                      ('numobs', 'numvars', 'numdeps'),
                                      fill_value=obsset['best_estimate']._FillValue)
    nersi_var = obsset.createVariable('negative_persistence', np.float32, 
                                      ('numobs', 'numvars','numfcsts', 'numdeps'),
                                      fill_value=obsset['persistence']._FillValue)
    inite_var.units = "m s-1"
    nersi_var.units = "m s-1"
    inite_var.long_name = "Model Initial Estimate"
    nersi_var.long_name = "Model negative persistence: Best estimate n-days after observation." ;
    inite_var[:] = inite
    nersi_var[:] = nersi
    obsset.close()
    return

#bash jobscripts/GEPSdate_Class4Currents.sh --date=20210105 -e=${ie} -f -s
def fill_missing_GEPS_dates(date_range, filter=True, it=2, fill=False, submit=True, fcst='' , mean=False):
    if ( fcst == 'GEPS_STO2X' ): fcst='ENAN'
    if ( fcst == 'STO2X' ): fcst='ENAN'
    date_start = check_date(date_range[0], outtype=datetime.datetime)
    date_final = check_date(date_range[1], outtype=datetime.datetime)
    print(date_start, date_final)
    date_now = date_start
    MISS=0
    while ( date_now <= date_final ):
        doFILL='LIST'
        if ( fill ): doFILL="FILL"
        print(doFILL, date_now)
        if ( fill ):
            if ( MISS >= 300 ):
                print("ABORT", MISS)
                fill=False
            missing_list = fill_missing_GEPS(date_now, filter=filter, it=it, fcst=fcst, submit=submit, mean=mean)
        else:
            missing_list = find_missing_GEPS(date_now, filter=filter, it=it, fcst=fcst, mean=mean)
            print(missing_list)
        MISS = MISS + len(missing_list)
        date_now = date_now + datetime.timedelta(days=1)
    return MISS
 
def fill_missing_GEPS(date, filter=True, it=2, fcst='', submit=True, mean=False):
    if ( fcst == 'GEPS_STO2X' ): fcst='ENAN'
    if ( fcst == 'STO2X' ): fcst='ENAN'
    if ( isinstance(date, list) ):
        for idate in date:
            fill_missing_GEPS(idate, filter=filter, it=it, fcst=fcst, submit=submit)
        return
    datestr=check_date(date)
    missing_list = find_missing_GEPS(date, filter=filter, it=it, fcst=fcst, mean=mean)
    print(missing_list)
    base_command='/home/dpe000/Class4_Currents/jobscripts/GEPSdate_Class4Currents.sh'
    if ( filter ):
        filstr='-f'
    else:
        filstr='-u'
    EA=''
    pass_fcst=''
    if ( len(fcst) > 0 ):
        pass_fcst='--fcst=GEPS_'+fcst
        EA='--EA'
    if ( fcst == 'ENAN' ):
        pass_fcst='--fcst=GEPS_STO2X'
    SUBMIT = ''
    if ( submit ): SUBMIT='-s'
        
    for missing in missing_list:
        estr=str(missing)
        if ( estr == 'enm' ): 
            estr='A'
            pass_fcst='--fcst='+fcst
        command=['bash', base_command, '-d='+datestr, '-e='+estr, filstr, EA, pass_fcst, SUBMIT]
        while '' in command: command.remove('')
        print(command)
        subprocess.call(command)
    return missing_list
    
def find_missing_GEPS(date, filter=True, it=2, fcst='', nens=21, mean=False):
    if ( fcst == 'GEPS_STO2X' ): fcst='ENAN'
    if ( fcst == 'STO2X' ): fcst='ENAN'
    if ( isinstance(date, list) ): 
        missing_listlist = []
        for idate in date:
            missing_list = find_missing_GEPS(idate, filter=filter, it=it, fcst=fcst)
            missing_listlist.append(missing_list)
        return missing_listlist
        
    insuffix='GEPS_orca025_currents' 
    if ( len(fcst) > 0 ):  insuffix=fcst+'_orca025_currents'
    if ( filter ):
        indir='CLASS4_currents_GEPS_FILT'
        if ( len(fcst) > 0 ):  indir='CLASS4_currents_'+fcst+'_FILT'
        insuffix=insuffix+'.f'
    else:
        indir='CLASS4_currents_GEPS_UFIL'
        if ( len(fcst) > 0 ):  indir='CLASS4_currents_'+fcst+'_UFIL'
        insuffix=insuffix+'.'
    nens_in = nens
    if ( mean ):  nens_in=0
    missing_list = find_missing_members(date, indir=indir, insuffix=insuffix, it=it, nens=nens_in)
    return missing_list   
        
def find_missing_members(date, indir='CLASS4_currents_GEPS_FILT', insuffix='GEPS_orca025_currents.f', it=2, nens=21):
    missing_list = []
    datestr=check_date(date)
    bile=indir+'/class4_'+datestr+'_'+insuffix+str(it)
    if ( nens == 0 ):
      ensstr='enm'
      qile=bile+'.'+ensstr+'.nc'
      exist=os.path.isfile(qile)
      #print(qile, exist)
      if ( not exist ): missing_list.append('enm')
    else:
      for ie in range(nens):
        ensstr=str(ie).zfill(3)
        qile=bile+'.'+ensstr+'.nc'
        exist=os.path.isfile(qile)
        #print(qile, exist)
        if ( not exist ): missing_list.append(ie)
    return missing_list
    
def cat_date_ensemble(date, indir='CLASS4_currents_GEPS_FILT', insuffix='GEPS_orca025_currents.f2', clobber=True, nens=21):
    if ( isinstance(date, list) ):
        for idate in date:
            cat_date_ensemble(idate, indir=indir, insuffix=insuffix, clobber=clobber, nens=nens)
    datestr=check_date(date)
    obsfile=indir+'/class4_'+datestr+'_'+insuffix+'.nc'
    tplfile=indir+'/class4_'+datestr+'_'+insuffix+'.000.nc'
    cat_obsfile_ensemble(obsfile, tplfile, nens=nens, clobber=clobber)
    return

def cat_obsfile_ensemble(obsfile, tplfile, nens=21, clobber=True):
    for iens in range(nens):
        ensstr=str(iens).zfill(3)
        ensfile=obsfile.replace('.nc', '.'+ensstr+'.nc')
        (LONO, LATO, depth), fcste = read_obsfile_variable(ensfile, 'forecast')
        if ( iens == 0 ):
            nobs, nvars, nfcsts, nenss, ndeps = np.shape(fcste)
            fcstv=np.zeros((nobs, nvars, nfcsts, nens, ndeps))
        fcstv[:,:,:,iens,:] = fcste[:,:,:,0,:]
    fcste=np.mean(fcstv, axis=3)
        
    write_model_obsfile_ensemble(obsfile, tplfile, fcstv, clobber=clobber)
    write_model_obsfile_ensemble(obsfile.replace('.nc', '.enm.nc'), tplfile, fcste, clobber=clobber)
    return

def assemble_ensemble_date(date, obspre='CLASS4_currents_GEPS_FILT/class4', obssuf='GEPS_orca025_currents', iters=['f1','f2'], nens=21, clobber=True):
    datestr=check_date(date)
    for itstr in iters:
        obsfile=obspre+'_'+datestr+'_'+obssuf+'.'+itstr+'.nc'
        print(obsfile)
        tplfile=obsfile.replace('nc','000.nc')
        cat_obsfile_ensemble(obsfile, tplfile, nens=nens, clobber=clobber)
    
def write_model_obsfile_ensemble(obsfile, tplfile, fcstv, clobber=True, fullpath='/fs/site6/eccc/mrd/rpnenv/dpe000/Class4_Currents'):
    #print('fcstv SHAPE', np.shape(fcstv))
    if ( fcstv.ndim == 5 ):
      nobs, nvars, nfcsts, nenss, ndeps = np.shape(fcstv)
    elif (fcstv.ndim == 4 ):
      nobs, nvars, nfcsts, ndeps = np.shape(fcstv)
    # ncks template file to new obsfile -- but remove forecast -- which we will need to re-create.
    # removing all three numfcsts length variables SHOULD remove numfcst dimension as well.  [THIS MAY FAIL?]
    # Using same module, we should be able to recreate persistence with a longer time line too.
    # Fullpath added for use with ssh.  Not needed for local host execution.
    if ( clobber ):
      rc=subprocess.call(['ncks','-O','-x','-v','forecast,persistence,negative_persistence',fullpath+'/'+tplfile, fullpath+'/'+obsfile], stdout=subprocess.PIPE, stderr=subprocess.PIPE) 
      print('Job RC = ', rc)
      #rc=subprocess.Popen(['ssh', 'ppp6', '/usr/bin/ncks','-O','-x','-v','forecast,persistence,negative_persistence',fullpath+'/'+tplfile, fullpath+'/'+obsfile], stdout=subprocess.PIPE, stderr=subprocess.PIPE) 
      #print('Job RC = ', rc.returncode)
      #print(rc.stdout.read(), rc.stderr.read())
      #print(''.join(['ssh', 'ppp6', '/usr/bin/ncks','-O','-x','-v','forecast,persistence,negative_persistence',fullpath+'/'+tplfile, fullpath+'/'+obsfile]))

    obsset = netCDF4.Dataset(obsfile,mode='r+')
    nfcsts_d = obsset.createDimension('numfcsts', nfcsts)
    if ( fcstv.ndim == 5 ):
        nens_d = obsset.createDimension('numens', nenss)
        fcstv_var = obsset.createVariable('forecast', np.float32,
                                      ('numobs', 'numvars', 'numfcsts', 'numens', 'numdeps'),
                                      fill_value=obsset['best_estimate']._FillValue)
    if ( fcstv.ndim == 4 ):
        fcstv_var = obsset.createVariable('forecast', np.float32,
                                      ('numobs', 'numvars', 'numfcsts', 'numdeps'),
                                      fill_value=obsset['best_estimate']._FillValue)
   
    fcstv_var.units = "m s-1"
    fcstv_var.long_name = "Model forecast counterpart of obs. value" ;
    fcstv_var[:] = fcstv
    obsset.close()
    return

def test_plot_fields(date=datetime.datetime(2021,10,30),lead=24):
    file_fcst=find_fcst_file.find_fcst_file('OPD', date, lead, 0, src='ocn', execute=True)
    LONS, LATS, TF = stfd.read_fstd_var(file_fcst, 'TM', typvar='P@')
    leu, UF = stfd.read_fstd_multi_lev(file_fcst, 'UUW',vfreq=24, typvar='P@')
    lev, VF = stfd.read_fstd_multi_lev(file_fcst, 'VVW',vfreq=24, typvar='P@')
    let, TF = stfd.read_fstd_multi_lev(file_fcst, 'TM',vfreq=24, typvar='P@')
    
    (T00, U00, V00) = (TF[0,:,:], UF[0,:,:], VF[0,:,:])
    plot_fields((T00, U00, V00), (LONS, LATS), suptitle=None, grid=True, outfile_prefix='PLOT')
    return
    


def write_mean_errors(out_prefix, dateint, MEAN_ERRORS, vec=['u','v'], tmp_prefix='NULL',depth_index=1):

    if ( isinstance(dateint, datetime.datetime) or isinstance(dateint, datetime.date) ): dateint=date.strftime("%Y%m%d")
    if ( isinstance(dateint, str) ): dateint=int(dateint)    
    #(mean_beste, rmse_beste, abse_beste), (mean_inite, rmse_inite, abse_inite), (mean_fcstv, rmse_fcstv, abse_fcstv), (mean_persi, rmse_persi, abse_persi), (mean_nersi, rmse_nersi, abse_nersi) = MEAN_ERRORS
    #KEEP 15m errors
    MEANU_ERRORS, MEANV_ERRORS, RMSEU_ERRORS, RMSEV_ERRORS, ABSEU_ERRORS, ABSEV_ERRORS = data_structure_errors(MEAN_ERRORS, depth_index=depth_index)
    file1=out_prefix+'_mean'+vec[0]+'.dat'
    file2=out_prefix+'_mean'+vec[1]+'.dat'   
    file3=out_prefix+'_rmse'+vec[0]+'.dat'
    file4=out_prefix+'_rmse'+vec[1]+'.dat'   
    file5=out_prefix+'_abse'+vec[0]+'.dat'
    file6=out_prefix+'_abse'+vec[1]+'.dat'   
    
    if ( tmp_prefix == 'NULL' ): tmp_prefix=out_prefix
    tile1=tmp_prefix+'_mean'+vec[0]+'.tmp'
    tile2=tmp_prefix+'_mean'+vec[1]+'.tmp'
    tile3=tmp_prefix+'_rmse'+vec[0]+'.tmp'
    tile4=tmp_prefix+'_rmse'+vec[1]+'.tmp'
    tile5=tmp_prefix+'_abse'+vec[0]+'.tmp'
    tile6=tmp_prefix+'_abse'+vec[1]+'.tmp'
     
    datadatefile.add_to_file(dateint, MEANU_ERRORS, file=file1, tmpfile=tile1)
    datadatefile.add_to_file(dateint, MEANV_ERRORS, file=file2, tmpfile=tile2)
    datadatefile.add_to_file(dateint, RMSEU_ERRORS, file=file3, tmpfile=tile3)
    datadatefile.add_to_file(dateint, RMSEV_ERRORS, file=file4, tmpfile=tile4)
    datadatefile.add_to_file(dateint, ABSEU_ERRORS, file=file5, tmpfile=tile5)
    datadatefile.add_to_file(dateint, ABSEV_ERRORS, file=file6, tmpfile=tile6)
    
    return
   
psyfile='CLASS4_currents/class4_20211031_PSY4V3R1_orca12_currents.nc'
obsfile='CLASS4_currents_CCMEP/class4_20211031_GIOPS_orca025_currents.nc'
def write_mean_errors_from_obsfile(date, indir, insuffix, out_prefix, tmp_prefix='NULL',depth_index=1,verbose=False):
    datestr=check_date(date)
    MEAN_ERRORS, MMEAN_ERRORS = load_mean_errors_from_obsfile(date, indir, insuffix)
    if ( verbose ):
        print( MEAN_ERRORS[0], MMEAN_ERRORS[0] )
        print('MEAN_ERRORS, MMEAN_ERRORS = Class4Current.load_mean_errors_from_obsfile("'+datestr+'",'+'"'+indir+'","'+insuffix+'")')
    write_mean_errors(out_prefix, datestr,  MEAN_ERRORS, vec=['u','v'], tmp_prefix=tmp_prefix,depth_index=depth_index)
    write_mean_errors(out_prefix, datestr, MMEAN_ERRORS, vec=['s','a'], tmp_prefix=tmp_prefix,depth_index=depth_index)
    return

def write_3mean_errors_from_obsfile(date, indir, insuffix, out_prefix, tmp_prefix='NULL',depth_index=1,verbose=False):
    datestr=check_date(date)
    MEAN_ERRORS_LIST, MMEAN_ERRORS_LIST = load_3mean_errors_from_obsfile(date, indir, insuffix)
    if ( verbose ):
        print( MEAN_ERRORS_LIST[0][0], MMEAN_ERRORS_LIST[0][0] )
        print('MEAN_ERRORS_LIST, MMEAN_ERRORS_LIST = Class4Current.load_3mean_errors_from_obsfile("'+datestr+'",'+'"'+indir+'","'+insuffix+'")')
    for ierr, MEAN_ERRORS in enumerate(MEAN_ERRORS_LIST):
        if ( ierr == 0 ): add='uncorr'
        if ( ierr == 1 ): add='stokes'
        if ( ierr == 2 ): add='stides'
        write_mean_errors(out_prefix+'_'+add, datestr,  MEAN_ERRORS_LIST[ierr], vec=['u','v'], tmp_prefix=tmp_prefix+'_'+add,depth_index=depth_index)
        write_mean_errors(out_prefix+'_'+add, datestr, MMEAN_ERRORS_LIST[ierr], vec=['s','a'], tmp_prefix=tmp_prefix+'_'+add,depth_index=depth_index)
    return

def find_date_in_obsfile(obsfile):
    date_index_start = obsfile.find('class4_')+7
    date_index_final = date_index_start+8
    datestr=obsfile[date_index_start:date_index_final]
    return datestr

def check_date(date, outtype=str):
    if ( (outtype==str) or (outtype==int) ):
      if ( isinstance(date, datetime.datetime) or isinstance(date, datetime.date) ):  datestr=date.strftime("%Y%m%d")
      if ( isinstance(date, int) ): datestr=str(date)
      if ( isinstance(date, str) ): datestr=date
    if ( outtype ==int ): datestr=int(datestr)
    if ( outtype== datetime.datetime ):
      if ( isinstance(date, int) ): date=str(date)
      if ( isinstance(date, str) and ( len(date) == 8 ) ):
        datestr=datetime.datetime.strptime(date, '%Y%m%d')  
      elif ( isinstance(date, str) and ( len(date) == 10 ) ):
        datestr=datetime.datetime.strptime(date, '%Y%m%d%H')  
      if ( isinstance(date, datetime.datetime) ): datestr=date
      
    return datestr

def write_mean_errors_model(dates, model, IIT=4, filter=True, verbose=False, Three=False):
    if ( ( model == 'PSY' ) or ( model == 'MERCATOR' ) ):
        indir='CLASS4_currents_CHARLY'
        if ( filter ): 
            insuffix='PSY4V3R1_orca12_currents-filtr'
            out_prefix='ERRORS_FILT/PSY4'
            tmp_prefix='site6/Class4_Currents/TMP/fPSY'
        else:
            insuffix='PSY4V3R1_orca12_currents'
            out_prefix='ERRORS_UFIL/PSY4'
            tmp_prefix='site6/Class4_Currents/TMP/uPSY'
    if ( ( model == 'GIOPS' ) or ( model=='CCMEP' ) ):
        if ( filter ):
            indir='CLASS4_currents_CCMEP_FILT'
            insuffix='GIOPS_orca025_currents.f'+str(IIT)
            out_prefix='ERRORS_FILT/GIOPS.'+str(IIT)
            tmp_prefix='site6/Class4_Currents/TMP/fGIOPS.'+str(IIT)
        else:
            indir='CLASS4_currents_CCMEP_UFIL'
            insuffix='GIOPS_orca025_currents.'+str(IIT)
            out_prefix='ERRORS_UFIL/GIOPS.'+str(IIT)
            tmp_prefix='site6/Class4_Currents/TMP/uGIOPS.'+str(IIT)
    if ( ( model == 'GEPS' ) or ( model=='ENSEMBLE' ) ):
        if ( filter ):
            indir='CLASS4_currents_GEPS_FILT'
            insuffix='GEPS_orca025_currents.f'+str(IIT)+'.enm'
            out_prefix='ERRORS_FILT/GEPS.'+str(IIT)
            tmp_prefix='site6/Class4_Currents/TMP/fGEPS.'+str(IIT)
        else:
            indir='CLASS4_currents_GEPS_UFIL'
            insuffix='GEPS_orca025_currents.'+str(IIT)+'.enm'
            out_prefix='ERRORS_UFIL/GEPS.'+str(IIT)
            tmp_prefix='site6/Class4_Currents/TMP/uGEPS.'+str(IIT)

    write_mean_errors_from_obsfiles(dates, indir, insuffix, out_prefix, tmp_prefix=tmp_prefix, depth_index=0, verbose=verbose, Three=Three)
    return
                
def write_mean_errors_from_obsfiles(dates, indir, insuffix, out_prefix, tmp_prefix='NULL', depth_index=1, verbose=False, Three=False):
    if ( isinstance(dates, str) ):
        if ( ( dates == 'all') or (dates == 'ALL') or ( dates == 'ls' ) ):
            datestr='????????'
        else:
            datestr=dates
            if ( len(datestr) != 8 ):
                print('Warning: Date String not CCYYMMDD', datestr)
            if ( len(datestr) > 8):
                print('CRITICAL WARNING:  Date String too long.  May not glob')
        obsfiles = sorted(glob.glob(indir+'/class4_'+datestr+'_'+insuffix+'.nc'))
        date_list = []
        for obsfile in obsfiles:
            date_list.append(find_date_in_obsfile(obsfile))
    elif ( isinstance(dates, list) ):
        date_list=dates
        obsfiles=[]
        for date in date_list:  
           datestr=check_date(date)
           obsfiles.append(indir+'/class4_'+datestr+'_'+insuffix+'.nc')
    else:
        date_list=[date]
          
    for date in date_list: 
        if ( not Three ):
            write_mean_errors_from_obsfile(date, indir, insuffix, out_prefix, tmp_prefix=tmp_prefix,depth_index=depth_index, verbose=verbose)
        else:
            write_3mean_errors_from_obsfile(date, indir, insuffix, out_prefix, tmp_prefix=tmp_prefix,depth_index=depth_index, verbose=verbose)
    return

def load_mean_errors_from_obsfile(date, indir, insuffix):
    datestr=check_date(date)
    obsfile=indir+'/class4_'+datestr+'_'+insuffix+'.nc'
    (LONO, LATO, depth), (obser, beste, inite, fcstv, persi, nersi) = read_obsfile_plus(obsfile)
    #nersi = persi.copy()  # nersi is actually the best estimate from the day AFTER the observation.  
                          # Not (yet?) in the observation file.

    Mobser = 0.0*obser.copy()
    Mbeste = 0.0*beste.copy()
    Minite = 0.0*inite.copy()
    Mfcstv = 0.0*fcstv.copy()
    Mpersi = 0.0*persi.copy()
    Mnersi = 0.0*nersi.copy()

    Mobser[:,0,:], Mobser[:,1,:] = speed_and_angle(obser[:,0,:], obser[:,1,:])
    Mbeste[:,0,:], Mbeste[:,1,:] = speed_and_angle(beste[:,0,:], beste[:,1,:])
    Minite[:,0,:], Minite[:,1,:] = speed_and_angle(inite[:,0,:], inite[:,1,:])
    Mfcstv[:,0,:,:], Mfcstv[:,1,:,:] = speed_and_angle(fcstv[:,0,:,:], fcstv[:,1,:,:])
    Mpersi[:,0,:,:], Mpersi[:,1,:,:] = speed_and_angle(persi[:,0,:,:], persi[:,1,:,:])
    Mnersi[:,0,:,:], Mnersi[:,1,:,:] = speed_and_angle(nersi[:,0,:,:], nersi[:,1,:,:])

    MEAN_ERRORS = calc_mean_error(obser, (beste, inite, fcstv, persi, nersi))
    MMEAN_ERRORS = calc_mean_error(Mobser, (Mbeste, Minite, Mfcstv, Mpersi, Mnersi), isangle=1)
    
    return MEAN_ERRORS, MMEAN_ERRORS
 
def load_3mean_errors_from_obsfile(date, indir, insuffix):
    datestr=check_date(date)
    obsfile=indir+'/class4_'+datestr+'_'+insuffix+'.nc'
    #print(obsfile)
    (LONO, LATO, depth, plf_code), (obser, beste, inite, fcstv, persi, smocv, stokv, tidev, bathy), (depth_qc, otime_qc, posit_qc, obser_qc) = read_obsfile_allt(obsfile)
    
    numobs, numvars, numfcsts, numdeps = fcstv.shape

    MEAN_ERRORS_LIST = []
    MMEAN_ERRORS_LIST = []    
    for ierror in range(3):
        if ( ierror == 0 ):
           pass # non adjustment of velocity
        if ( ierror == 1 ):  # add 0.55 * stokes drift
            beste = beste + 0.55*stokv
            inite = inite + 0.55*stokv
            for ifcst in range(numfcsts):
                fcstv[:,:,ifcst,:] = fcstv[:,:,ifcst,:] + 0.55*stokv
                persi[:,:,ifcst,:] = persi[:,:,ifcst,:] + 0.55*stokv
        if ( ierror == 2 ):  # add tides
            beste = beste + tidev
            inite = inite + tidev
            for ifcst in range(numfcsts):
                fcstv[:,:,ifcst,:] = fcstv[:,:,ifcst,:] + tidev
                persi[:,:,ifcst,:] = persi[:,:,ifcst,:] + tidev
    
        Mobser = 0.0*obser.copy()
        Mbeste = 0.0*beste.copy()
        Minite = 0.0*inite.copy()
        Mfcstv = 0.0*fcstv.copy()
        Mpersi = 0.0*persi.copy()

        Mobser[:,0,:], Mobser[:,1,:] = speed_and_angle(obser[:,0,:], obser[:,1,:])
        Mbeste[:,0,:], Mbeste[:,1,:] = speed_and_angle(beste[:,0,:], beste[:,1,:])
        Minite[:,0,:], Minite[:,1,:] = speed_and_angle(inite[:,0,:], inite[:,1,:])
        Mfcstv[:,0,:,:], Mfcstv[:,1,:,:] = speed_and_angle(fcstv[:,0,:,:], fcstv[:,1,:,:])
        Mpersi[:,0,:,:], Mpersi[:,1,:,:] = speed_and_angle(persi[:,0,:,:], persi[:,1,:,:])

        MEAN_ERRORS = calc_mean_error(obser, (beste, inite, fcstv, persi))
        MMEAN_ERRORS = calc_mean_error(Mobser, (Mbeste, Minite, Mfcstv, Mpersi), isangle=1)
        
        MEAN_ERRORS_LIST.append(MEAN_ERRORS)
        MMEAN_ERRORS_LIST.append(MMEAN_ERRORS)
    
    return MEAN_ERRORS_LIST, MMEAN_ERRORS_LIST

def load_mean_errors_from_obsfiles(dates, indir, insuffix):
    if ( isinstance(dates, str) ):
        if ( ( dates == 'all') or (dates == 'ALL') or ( dates == 'ls' ) ):
            datestr='????????'
            obsfiles = sorted(glob.glob(indir+'/class4_'+datestr+'_'+insuffix+'.nc'))
            date_list = []
            for obsfile in obsfiles:
                date_list.append(find_date_in_obsfile(obsfile))
    elif ( isinstance(dates, list) ):
        date_list=dates
        obsfiles=[]
        for datestr in dates:  
           obsfiles.append(indir+'/'+indir+'/class4_'+datestr+'_'+insuffix+'.nc')
          
    MEAN_ERRORS_LIST = []
    MMEAN_ERRORS_LIST = []
    dates_good = []
    for date in date_list: 
        try:
            MEAN_ERRORS, MMEAN_ERRORS = load_mean_errors_from_obsfile(date, indir, insuffix)
        except:
            print('load from obsfile failed', date)
        else:
            dates_good.append(date)
            MEAN_ERRORS_LIST.append(MEAN_ERRORS)
            MMEAN_ERRORS_LIST.append(MMEAN_ERRORS)
    return dates_good, MEAN_ERRORS_LIST, MMEAN_ERRORS_LIST
    
def load_3mean_errors_from_obsfiles(dates, indir, insuffix):
    if ( isinstance(dates, str) ):
        if ( ( dates == 'all') or (dates == 'ALL') or ( dates == 'ls' ) ):
            datestr='????????'
            obsfiles = sorted(glob.glob(indir+'/class4_'+datestr+'_'+insuffix+'.nc'))
            date_list = []
            for obsfile in obsfiles:
                date_list.append(find_date_in_obsfile(obsfile))
    elif ( isinstance(dates, list) ):
        date_list=dates
        obsfiles=[]
        for this_date in dates:
           datestr=check_date(this_date)  
           obsfiles.append(indir+'/'+indir+'/class4_'+datestr+'_'+insuffix+'.nc')
          
    SUPER_MEAN_ERRORS_LIST = [ [], [], [] ]
    SUPER_MMEAN_ERRORS_LIST = [ [], [], [] ]
    dates_good = []
    for date in date_list: 
        try:
            MEAN_ERRORS, MMEAN_ERRORS = load_3mean_errors_from_obsfile(date, indir, insuffix)
        except:
            print('load from obsfile failed', date)
        else:
            dates_good.append(date)
            for ierr in range(3):
              SUPER_MEAN_ERRORS_LIST[ierr].append(MEAN_ERRORS[ierr])
              SUPER_MMEAN_ERRORS_LIST[ierr].append(MMEAN_ERRORS[ierr])
    return dates_good, SUPER_MEAN_ERRORS_LIST, SUPER_MMEAN_ERRORS_LIST
    
def data_structure_errors(MEAN_ERRORS, depth_index=1):

    if ( len(MEAN_ERRORS) == 5 ):    
        (mean_beste, rmse_beste, abse_beste), (mean_inite, rmse_inite, abse_inite), (mean_fcstv, rmse_fcstv, abse_fcstv), (mean_persi, rmse_persi, abse_persi), (mean_nersi, rmse_nersi, abse_nersi) = MEAN_ERRORS
        #KEEP 15m errors (depth_index=1)
        MEANU_ERRORS=np.concatenate((np.array([mean_beste[0,depth_index]]), np.array([mean_inite[0,depth_index]]), mean_fcstv[0,:,depth_index].flatten(), mean_persi[0,:,depth_index].flatten(), mean_nersi[0,:,depth_index].flatten()))
        MEANV_ERRORS=np.concatenate((np.array([mean_beste[1,depth_index]]), np.array([mean_inite[1,depth_index]]), mean_fcstv[1,:,depth_index].flatten(), mean_persi[1,:,depth_index].flatten(), mean_nersi[1,:,depth_index].flatten()))
        RMSEU_ERRORS=np.concatenate((np.array([rmse_beste[0,depth_index]]), np.array([rmse_inite[0,depth_index]]), rmse_fcstv[0,:,depth_index].flatten(), rmse_persi[0,:,depth_index].flatten(), rmse_nersi[0,:,depth_index].flatten()))
        RMSEV_ERRORS=np.concatenate((np.array([rmse_beste[1,depth_index]]), np.array([rmse_inite[1,depth_index]]), rmse_fcstv[1,:,depth_index].flatten(), rmse_persi[1,:,depth_index].flatten(), rmse_nersi[1,:,depth_index].flatten()))
        ABSEU_ERRORS=np.concatenate((np.array([abse_beste[0,depth_index]]), np.array([abse_inite[0,depth_index]]), abse_fcstv[0,:,depth_index].flatten(), abse_persi[0,:,depth_index].flatten(), abse_nersi[0,:,depth_index].flatten()))
        ABSEV_ERRORS=np.concatenate((np.array([abse_beste[1,depth_index]]), np.array([abse_inite[1,depth_index]]), abse_fcstv[1,:,depth_index].flatten(), abse_persi[1,:,depth_index].flatten(), abse_nersi[1,:,depth_index].flatten()))
    if ( len(MEAN_ERRORS) == 4 ):    
        (mean_beste, rmse_beste, abse_beste), (mean_inite, rmse_inite, abse_inite), (mean_fcstv, rmse_fcstv, abse_fcstv), (mean_persi, rmse_persi, abse_persi) = MEAN_ERRORS
        #KEEP 15m errors (depth_index=1)
        MEANU_ERRORS=np.concatenate((np.array([mean_beste[0,depth_index]]), np.array([mean_inite[0,depth_index]]), mean_fcstv[0,:,depth_index].flatten(), mean_persi[0,:,depth_index].flatten()))
        MEANV_ERRORS=np.concatenate((np.array([mean_beste[1,depth_index]]), np.array([mean_inite[1,depth_index]]), mean_fcstv[1,:,depth_index].flatten(), mean_persi[1,:,depth_index].flatten()))
        RMSEU_ERRORS=np.concatenate((np.array([rmse_beste[0,depth_index]]), np.array([rmse_inite[0,depth_index]]), rmse_fcstv[0,:,depth_index].flatten(), rmse_persi[0,:,depth_index].flatten()))
        RMSEV_ERRORS=np.concatenate((np.array([rmse_beste[1,depth_index]]), np.array([rmse_inite[1,depth_index]]), rmse_fcstv[1,:,depth_index].flatten(), rmse_persi[1,:,depth_index].flatten()))
        ABSEU_ERRORS=np.concatenate((np.array([abse_beste[0,depth_index]]), np.array([abse_inite[0,depth_index]]), abse_fcstv[0,:,depth_index].flatten(), abse_persi[0,:,depth_index].flatten()))
        ABSEV_ERRORS=np.concatenate((np.array([abse_beste[1,depth_index]]), np.array([abse_inite[1,depth_index]]), abse_fcstv[1,:,depth_index].flatten(), abse_persi[1,:,depth_index].flatten()))

    return MEANU_ERRORS, MEANV_ERRORS, RMSEU_ERRORS, RMSEV_ERRORS, ABSEU_ERRORS, ABSEV_ERRORS
    
def data_structure_list_of_errors(MEAN_ERRORS_LIST, depth_index=1):

    MEANU_ERRORS_LIST = []
    RMSEU_ERRORS_LIST = []
    ABSEU_ERRORS_LIST = []
    MEANV_ERRORS_LIST = []
    RMSEV_ERRORS_LIST = []
    ABSEV_ERRORS_LIST = []
    for MEAN_ERRORS in MEAN_ERRORS_LIST:
        MEANU_ERRORS, MEANV_ERRORS, RMSEU_ERRORS, RMSEV_ERRORS, ABSEU_ERRORS, ABSEV_ERRORS = data_structure_errors(MEAN_ERRORS, depth_index=depth_index)
        MEANU_ERRORS_LIST.append(MEANU_ERRORS)
        RMSEU_ERRORS_LIST.append(RMSEU_ERRORS)
        MEANV_ERRORS_LIST.append(MEANV_ERRORS)
        RMSEV_ERRORS_LIST.append(RMSEV_ERRORS)

    MEANU_ERRORS_DATA = np.transpose(np.array(MEANU_ERRORS_LIST))
    RMSEU_ERRORS_DATA = np.transpose(np.array(RMSEU_ERRORS_LIST))
    ABSEU_ERRORS_DATA = np.transpose(np.array(ABSEU_ERRORS_LIST))
    MEANV_ERRORS_DATA = np.transpose(np.array(MEANV_ERRORS_LIST))
    RMSEV_ERRORS_DATA = np.transpose(np.array(RMSEV_ERRORS_LIST))
    ABSEV_ERRORS_DATA = np.transpose(np.array(ABSEV_ERRORS_LIST))
   
    return MEANU_ERRORS_DATA, RMSEU_ERRORS_DATA, ABSEU_ERRORS_DATA, MEANV_ERRORS_DATA, RMSEV_ERRORS_DATA, ABSEV_ERRORS_DATA

date_start = datetime.datetime(2019, 4, 1)
date_final = datetime.datetime(2022, 3, 1)

ERR_list={'best': [0,1],'init' : [1,2] ,'fcst' : [2,12], 'pers' : [12,22], 'ners' : [22,32] }
#MEANU_ERRORS=np.concatenate((np.array([mean_beste[0,depth_index]]), mean_persi[0,:,depth_index].flatten(), mean_fcstv[0,:,depth_index].flatten(), mean_nersi[0,:,depth_index].flatten()))
def plot_errors(EXPTS, LABELS=None, in_prefix='ERRORS', out_prefix='PLOTS/', date_range=(date_start, date_final) ):
    if ( isinstance(LABELS, type(None)) ): LABELS=EXPTS
    (date_min, date_max) = date_range
    plt.rc('font', family='serif')
    #plt.rc('text', usetex=True)
    clrs = ['r','b','g']
    clls = ['m','c','y']    
    for typerror in ['rmseu', 'rmsev', 'rmses', 'rmsea', 'abseu', 'absev', 'abses', 'absea', 'meanu', 'meanv', 'means', 'meana']:
        tfig, taxe = plt.subplots()
        lfig, laxe = plt.subplots()

        for ie, EXPT in enumerate(EXPTS):
            LABE=LABELS[ie]
            clr = clrs[ie%3]
            cll = clls[ie%3]
            errfile=in_prefix+'/'+EXPT+'_'+typerror+'.dat'
            print(errfile)
            intdate, errors = datadatefile.read_file(errfile)
            print(len(intdate))
            dates = datadatefile.convert_strint_datelist(intdate)
    
            new_dates=[]
            new_errors=[]
    
            for idate,date in enumerate(dates):
                if ( ( date >= date_min ) and (date <= date_max ) ):
                    new_dates.append(date)
                    new_errors.append(errors[:,idate])
          
            dates = new_dates
            errors = np.transpose( np.array(new_errors) )  # np.array will flip the time, variable indices.
            print(len(dates))
    
            errors_list = errors.tolist()
            fcst = []
            pers = []
            ners = []
            for ierr, error in enumerate(errors_list):
                skip = False
                if   ( ierr == ERR_list['best'][0] ):  # BEST ESTIMATE
                    best=np.mean(error)
                    lclr=clr
                elif ( ierr == ERR_list['init'][0] ) :   # INIT ESTIMATE
                   lclr=clr
                   lsty='-'
                   init=np.mean(error)
                   skip = True
                elif ( ierr < ERR_list['fcst'][1] ): # FORECAST
                   lclr=clr
                   lsty='--'
                   fcst.append(np.mean(error))
                   skip = True
                   if ( ierr == ERR_list['fcst'][0] ): skip = False
                elif ( ierr < ERR_list['pers'][1] ): # PERSISTENCE
                    lclr=cll
                    pers.append(np.mean(error))
                    lsty=':'        
                    skip = True
                elif ( ierr < ERR_list['ners'][1] ): # Negative Persist
                    lclr=cll
                    ners.append(np.mean(error))
                    inow=len(ners)
                    if ( ners[inow-1] == pers[inow-1] ): skip = True        
                    lsty='-.'
                    skip = True
                if ( ierr == 0 ):
                     taxe.plot(dates, error, color=lclr, linewidth=0.5, label=LABE+' best estimate')
                elif ( not skip ):
                    taxe.plot(dates, error, color=lclr, linestyle=lsty, linewidth=0.5)

            fcst.insert(0, init)
            pers.insert(0, best)
            ners.insert(0, best)
            laxe.plot(range(len(fcst)), fcst, color=clr, label=LABE+' forecast')
            laxe.plot(range(len(pers)), pers, color=cll, linestyle='--', label=LABE+' persistence')
            if ( ( np.mean(pers) != np.mean(ners) ) and ( not skip) ):
                laxe.plot(range(len(ners)), ners, color=cll, linestyle=':', label=LABE+' negative persistence')

        taxe.legend()
        laxe.legend()

        tfig.savefig(out_prefix+typerror+'_time.png')
        lfig.savefig(out_prefix+typerror+'_lead.png')
    
        plt.close(tfig)
        plt.close(lfig)
    
    
    return

def create_dates(date_start, date_final, date_inc=1):
    dates_list = []
    this_date = check_date(date_start, outtype=datetime.datetime)
    stop_date = check_date(date_final, outtype=datetime.datetime)
    while ( this_date <= stop_date ):
        dates_list.append(this_date)
        this_date = this_date + datetime.timedelta(days=date_inc)
    return dates_list
       
def load_3errors_from_obsfile(date, indir, insuffix, ierror=0, etype='mean',magnitude=False):
    datestr=check_date(date)
    obsfile=indir+'/class4_'+datestr+'_'+insuffix+'.nc'
    #print(obsfile)
    (LONO, LATO, depth, plf_code), (obser, beste, inite, fcstv, persi, smocv, stokv, tidev, bathy), (depth_qc, otime_qc, posit_qc, obser_qc) = read_obsfile_allt(obsfile)
    #print( np.mean(np.square(obser-beste)))         
    
    numobs, numvars, numfcsts, numdeps = fcstv.shape

    if ( ierror == 0 ):
       pass # non adjustment of velocity
    if ( ierror == 1 ):  # add 0.55 * stokes drift
        beste = beste + 0.55*stokv
        inite = inite + 0.55*stokv
        for ifcst in range(numfcsts):
            fcstv[:,:,ifcst,:] = fcstv[:,:,ifcst,:] + 0.55*stokv
            persi[:,:,ifcst,:] = persi[:,:,ifcst,:] + 0.55*stokv
    if ( ierror == 2 ):  # add tides
        beste = beste + tidev
        inite = inite + tidev
        for ifcst in range(numfcsts):
            fcstv[:,:,ifcst,:] = fcstv[:,:,ifcst,:] + tidev
            persi[:,:,ifcst,:] = persi[:,:,ifcst,:] + tidev

    if ( magnitude ):  
        (obser, beste, inite, fcstv, persi) = speed_and_angle_list((obser, beste, inite, fcstv, persi))
        beste, inite, fcstv, persi = calc_error(obser, (beste, inite, fcstv, persi), etype=etype, isangle=1)
    else:
        beste, inite, fcstv, persi = calc_error(obser, (beste, inite, fcstv, persi), etype=etype)
    #print( np.mean(np.square(beste)))         

    return (LONO, LATO, depth, plf_code), (obser, beste, inite, fcstv, persi)    

def load_3errors_from_obsfile_single(date, indir, insuffix, variable, ierror=0, etype='mean',magnitude=False):
    (LONO, LATO, depth), (obser, field) = load_3values_from_obsfile_single(date, indir, insuffix, variable, ierror=ierror, etype=etype,magnitude=magnitude)

    if ( magnitude ):  
        error = calc_error(obser, field, etype=etype, isangle=1)
    else:
        error = calc_error(obser, field, etype=etype)

    return (LONO, LATO, depth), error    

def adjust_velocities(vfield, stokes_vars, ierror=0):
    tidev, stokv = stokes_vars
    if ( vfield.ndim == 3 ):
        numobs, numvars, numdeps = vfield.shape
    if ( vfield.ndim == 4 ):
        numobs, numvars, numfcsts, numdeps = vfield.shape
    if ( vfield.ndim == 5 ):
        numobs, numvars, numfcsts, numensm, numdeps = vfield.shape

    nfield = vfield.copy()    

    if ( ierror == 0 ):
       pass
    if ( ierror == 1 ):  # add 0.55 * stokes drift
        if ( vfield.ndim == 3 ):
            nfield = nfield + 0.55*stokv
        if ( vfield.ndim == 4 ):
            for ifcst in range(numfcsts):
                nfield[:,:,ifcst,:] = nfield[:,:,ifcst,:] + 0.55*stokv
        if ( vfield.ndim == 5 ):
            for ifcst in range(numfcsts):
              for iensm in range(numensm):
                nfield[:,:,ifcst, iensm, :] = nfield[:,:,ifcst,iensm, :] + 0.55*stokv
    if ( ierror == 2 ):  # add tides
        if ( vfield.ndim == 3 ):
            nfield = nfield + tidev
        if ( vfield.ndim == 4 ):
            for ifcst in range(numfcsts):
                nfield[:,:,ifcst,:] = nfield[:,:,ifcst,:] + tidev
        if ( vfield.ndim == 5 ):
            for ifcst in range(numfcsts):
              for iensm in range(numensm):
                nfield[:,:,ifcst, iensm, :] = nfield[:,:,ifcst,iensm, :] + tidev

    return nfield
    
def load_3values_from_obsfile_single(date, indir, insuffix, variable, ierror=0, magnitude=False):
    datestr=check_date(date)
    obsfile=indir+'/class4_'+datestr+'_'+insuffix+'.nc'
    #print(obsfile)
    tidev, stokv = read_stokes_variables(obsfile)
    (LONO, LATO, depth), (obser, field) = read_obsfile_variables(obsfile, ['observation', variable])
    #print( np.mean(np.square(obser-beste)))         

    field = adjust_velocities( field, (tidev, stokv), ierror=ierror)    

    if ( magnitude ):  
        obser, field = speed_and_angle_list((obser, field))

    return (LONO, LATO, depth), (obser, field)    

def calc_field_errors( obser, field, isangle=-1):
    merror = calc_error(obser, field, etype='mean', isangle=isangle)
    obspad = calc_error(obser, field, etype='obs')
    sqrerr = np.square(merror)
    abserr = np.absolute(merror)
    covari = np.multiply(obspad, field)
    varfld = np.square(field)
    varobs = np.square(obspad)
    return merror, sqrerr, abserr, covari, varfld, varobs
    
    
def csum_field( field, MN_field=None, axis=0):
    if ( isinstance(field, list) or isinstance(field, tuple) ):
        NMN_field = []
        for ii, ifield in enumerate(field):
            if ( isinstance(MN_field, type(None)) ):
                iMN = None
            else:
                iMN = MN_field[ii]
            iNN, iGC = csum_field(ifield, iMN)
            NMN_field.append(iNN)
        return NMN_field, iGC   #  Assumes all fields have same length?
    if ( isinstance(MN_field, type(None)) ):
        NMN_field = np.sum(field, axis=axis)
    else:
        NMN_field = MN_field + np.sum(field, axis=axis) 
    GCN = np.ma.size(field, axis=axis)
    return NMN_field, GCN

def calc_errors_date_analysis(date, expt, variable, insuffix='GIOPS_orca025_currents.f2', ierror=1, ddeg=4.0):
        (LONO, LATO, depth), (obser, field) = load_3values_from_obsfile_single(date, expt, insuffix, variable, ierror=1)
        (LONM, LATM, depthM), (Mobser, Mfield) = load_3values_from_obsfile_single(date, expt, insuffix, variable, ierror=1, magnitude=True)
        obser = np.concatenate([obser, Mobser[:,[0]]], axis=1)  ## add speed to u,v variables -- if want angle don't [0] index.
        field = np.concatenate([field, Mfield[:,[0]]], axis=1)
        merror, sqrerr, abserr, covari, varfld, varobs = calc_field_errors(obser, field)
        # THIS IS NOT NECESSARY IF THE FIELD AND OBSERVATION ARE SAME DIMENSION -- BUT IT's EASILY CORRECTED
        obser_padded =  calc_error(obser, field, etype='obs')
        obser_madded =  calc_error(Mobser, Mfield, etype='obs')
        #print(obser_padded.shape, field.shape, varobs.shape, varfld.shape, covari.shape, merror.shape, sqrerr.shape, abserr.shape)
        FIELDS=(obser_padded, field, varobs, varfld, covari, merror, sqrerr, abserr)
        # GLOBAL SUMS
        MN_FIELDS, GCN = csum_field( FIELDS, None)
        GD_FIELDS = []
        for fld in FIELDS:
            (LONG, LATG), grid, CNT = make_error_grid(fld.shape, ddeg=ddeg)
            GD_FIELDS.append(grid)
        (LONG, LATG), GD_FIELDS, CNT = add_error_grid_list( (LONO,LATO), FIELDS, GD_FIELDS, CNT, ddeg=ddeg)
        return (LONG, LATG), MN_FIELDS, GD_FIELDS, GCN, CNT, date       

def calc_errors_date_analysis_ensemble(date, expt, variable, insuffix='GIOPS_orca025_currents.f2', ierror=1, ddeg=4.0, ens_axis=2):
        (LONO, LATO, depth), (obser, field) = load_3values_from_obsfile_single(date, expt, insuffix, variable, ierror=ierror)
        (LONM, LATM, depthM), (Mobser, Mfield) = load_3values_from_obsfile_single(date, expt, insuffix, variable, ierror=ierror, magnitude=True)
        obser = np.concatenate([obser, Mobser[:,[0]]], axis=1)  ## add speed to u,v variables -- if want angle don't [0] index.
        field = np.concatenate([field, Mfield[:,[0]]], axis=1)
        if ( ens_axis > 0 ):
            ensm_field = np.nanmean(field, axis=ens_axis)
        else:
            ensm_field = field
        merror, sqrerr, abserr, covari, varfld, varobs = calc_field_errors(obser, ensm_field)
        crps_error = calc_crps(obser, field, ens_axis=ens_axis)
        # THIS IS NOT NECESSARY IF THE FIELD AND OBSERVATION ARE SAME DIMENSION -- BUT IT's EASILY CORRECTED
        obser_padded =  calc_error(obser, ensm_field, etype='obs')
        #print(obser_padded.shape, field.shape, varobs.shape, varfld.shape, covari.shape, merror.shape, sqrerr.shape, abserr.shape, crps_error.shape)
        FIELDS=(obser_padded, ensm_field, varobs, varfld, covari, merror, sqrerr, abserr, crps_error)
        # GLOBAL SUMS
        MN_FIELDS, GCN = csum_field( FIELDS, None)
        GD_FIELDS = []
        for fld in FIELDS:
            (LONG, LATG), grid, CNT = make_error_grid(fld.shape, ddeg=ddeg)
            GD_FIELDS.append(grid)
        (LONG, LATG), GD_FIELDS, CNT = add_error_grid_list( (LONO,LATO), FIELDS, GD_FIELDS, CNT, ddeg=ddeg)
        return (LONG, LATG), MN_FIELDS, GD_FIELDS, GCN, CNT, date       

def add_numpy_list( LISTOFTIMES ):
    # First list is in time
    # Second list is FLDS
    ntimes = len( LISTOFTIMES)
    nfield = len(LISTOFTIMES[0])
    LISTOFSUMS=[]
    for ifield in range(nfield):
        FLDSINTIME = [ LISTOFTIMES[itime][ifield] for itime in range(ntimes) ]
        SUMSINTIME = sum(FLDSINTIME)
        LISTOFSUMS.append(SUMSINTIME)
    return LISTOFSUMS
    
def loop_analysis_dates(dates, expt, variable, insuffix='GIOPS_orca025_currents.f2', ierror=1, ddeg=4.0, mp=True, ens_axis=0):

    ndates=len(dates)
    MN_LIST = []
    GD_LIST = []
    GC_LIST = []
    CN_LIST = []
    DT_LIST = []
    time0 = time.time()
    if  ( not mp ):
        for date in dates:
            print(date)
            # MN_FIELDS=Global / GD_FIELDS=Binned / GCN=Global Count / CNT=Binned Count
            if ( ens_axis == 0 ): 
                (LONG, LATG), MN_FIELDS, GD_FIELDS, GCN, CNT, DT = calc_errors_date_analysis(date, expt, variable, insuffix=insuffix, ierror=ierror, ddeg=ddeg)   
            else:
                (LONG, LATG), MN_FIELDS, GD_FIELDS, GCN, CNT, DT = calc_errors_date_analysis_ensemble(date, expt, variable, insuffix=insuffix, ierror=ierror, ddeg=ddeg, ens_axis=ens_axis)
            MN_LIST.append(MN_FIELDS)
            GD_LIST.append(GD_FIELDS)  
            GC_LIST.append(GCN)
            CN_LIST.append(CNT)
            DT_LIST.append(DT)
    else:
        nproc=min(num_cpus, ndates)
        pool = multiprocessing.Pool(nproc)
        izip = list(zip( dates, itertools.repeat(expt), itertools.repeat(variable)))
        if ( ens_axis == 0 ):
            FULL_RESULT_LIST = pool.starmap(partial(calc_errors_date_analysis, insuffix=insuffix, ierror=ierror, ddeg=ddeg), izip)
        else:
            FULL_RESULT_LIST = pool.starmap(partial(calc_errors_date_analysis_ensemble, insuffix=insuffix, ierror=ierror, ddeg=ddeg, ens_axis=ens_axis), izip)
        pool.close()
        pool.join()
        for idate in range(ndates):
            (LONG, LATG), MN_FIELDS, GD_FIELDS, GCN, CNT, DT = FULL_RESULT_LIST[idate]
            #print(type(MN_FIELDS), type(GD_FIELDS), GCN, type(CNT))
            MN_LIST.append(MN_FIELDS)
            GD_LIST.append(GD_FIELDS)  
            GC_LIST.append(GCN)
            CN_LIST.append(CNT)
            DT_LIST.append(DT)
    if ( dates == DT_LIST ):
        print("DATES ARE CORRECTLY SORTED")
    else:
        print("WARNING/WARNING: DATES ARE NOT CORRECTLY SORTED")
    file=['ANAL_data/U_'+variable+'_'+expt+'.dat', 'ANAL_data/V_'+variable+'_'+expt+'.dat']
    tile=['ANAL_data/U_'+variable+'_'+expt+'.tmp', 'ANAL_data/V_'+variable+'_'+expt+'.tmp']
    for idate, date in enumerate(DT_LIST):
        MN_FIELDS = MN_LIST[idate]
        GCN = GC_LIST[idate]
        dateint=check_date(date, outtype=int)
        for ifld in range(2):
           mn_fields = np.array([FIELD[ifld,0] for FIELD in MN_FIELDS]+[GCN])
           datadatefile.add_to_file(dateint, mn_fields, file=file[ifld], tmpfile=tile[ifld])               
    MNS_FIELDS = add_numpy_list(MN_LIST)
    GDS_FIELDS = add_numpy_list(GD_LIST)
    GCS = sum(GC_LIST)
    CNS = sum(CN_LIST)
    print(GCS)
    timef = time.time() - time0
    print('TIME ELAPSED', timef)
    
    return (LONG, LATG), MNS_FIELDS, GDS_FIELDS, GCS, CNS, MN_LIST, GD_LIST, GC_LIST, CN_LIST              

def calc_mean_field(sum_field, cnt):
    if ( isinstance(sum_field, list) ):
        mean_field = []
        for sum_field_element in sum_field:
            mean_field_element = calc_mean_field(sum_field_element, cnt)
            mean_field.append(mean_field_element)
    elif ( isinstance(cnt, int) or isinstance(cnt, float) ): 
        mean_field=0.0*sum_field
        if ( cnt != 0 ): mean_field=sum_field/cnt
    elif ( isinstance(cnt, np.ndarray) ):
        isfinite = np.where(cnt > 0 )
        mean_field = 0.0 * sum_field 
        mean_field = np.NaN * np.ones(sum_field.shape)
        nd = 0
        #if ( sum_field.ndim == 3 ): nx, ny, nv = sum_field.shape
        if ( sum_field.ndim == 4 ): 
          nx, ny, nv, nd = sum_field.shape
          for iiv in range(nv):
            for iid in range(nd):
              if ( nd > 0 ): 
                iii_sum_field = sum_field[:, :, iiv, iid]
                iii_mean_field = mean_field[:,:,iiv, iid]
                iii_mean_field[isfinite] = iii_sum_field[isfinite] / cnt[isfinite]
                mean_field[:,:,iiv,iid] = iii_mean_field
        if ( sum_field.ndim == 2 ): 
            nv, nt = sum_field.shape
            for iiv in range(nv):
                mean_field[iiv, isfinite] = sum_field[iiv, isfinite] / cnt[isfinite]
    return mean_field

#FIELDS=(obser_padded, field, varobs, varfld, covari, merror, sqrerr, abserr)
def get_correlation_fields(field_list):
    sum_obs = field_list[0]
    sum_fld = field_list[1]
    var_obs = field_list[2]
    var_fld = field_list[3]
    cov_obs_fld = field_list[4]
    return sum_obs, sum_fld, var_obs, var_fld, cov_obs_fld

def calc_variance(ncnt, sum_x, sum_x2):
    if ( isinstance(ncnt, int) or isinstance(ncnt, float) ):
        vari = 0.0
        if ( ncnt > 1 ):
            vari = sum_x2 / ncnt - np.square(sum_x/ncnt)
    elif ( isinstance(ncnt, np.ndarray) ):
        nd=0
        if ( sum_x.ndim == 3 ): 
            nx, ny, nv = sum_x.shape
            pass
        if ( sum_x.ndim == 4 ): 
            nx, ny, nv, nd = sum_x.shape
            vari = 0.0 * sum_x
            isfinite = np.where(ncnt > 1 )
            for iiv in rangve(nv):
                for iid in range(nd):
                    iii_sum_x = sum_x[:,:,iiv,iid]
                    iii_sum_x2 = sum_x2[:,:, iiv, iid]
                    iii_vari = vari[:,:,iiv,iid]
                    var1 = iii_sum_x2[isfinite] / ncnt[isfinite]
                    var2 = np.square(iii_sum_x[isfinite] / ncnt[isfinite])
                    iii_vari[isfinite] = var1 - var2
                    vari[:,:,iiv,iid] = iii_vari
        if ( sum_x.ndim == 2 ): 
            nv, nt = sum_x.shape
            for iiv in range(nv):
                var1 = sum_x2[iiv, isfinite] / ncnt[isfinite]
                var2 = np.square(sum_x[isfinite] / ncnt[isfinite])
                vari[iiv, isfinite] = var1 - var2
    return vari
        
def calc_correlation(ncnt, sum_x, sum_y, sum_xy, sum_x2, sum_y2):
    if ( isinstance(ncnt, int) or isinstance(ncnt, float) ):
        rcorr = 0
        if ( ncnt > 1 ):
            rnum = ncnt * sum_xy - sum_x * sum_y
            rden = np.sqrt((ncnt*sum_x2-(sum_x)**2)*(ncnt*sum_y2-(sum_y)**2))
            rcorr = rnum / rden
    elif ( isinstance(ncnt, np.ndarray) ):
        nd = 0
        if ( sum_x.ndim == 3 ): 
            nx, ny, nv = sum_x.shape
            pass  # not configured yet
        if ( sum_x.ndim == 4 ): 
            nx, ny, nv, nd = sum_x.shape
            rcorr = np.NaN * np.ones(sum_x.shape)
            isfinite = np.where(ncnt > 1 )  # NEED AT LEAST TWO POINTS FOR CORRELATOIN
            for iiv in range(nv):
                for iid in range(nd):
                    iii_sum_xy = sum_xy[:,:,iiv,iid]
                    iii_sum_x = sum_x[:,:,iiv,iid]
                    iii_sum_y = sum_y[:,:,iiv,iid]
                    iii_sum_x2 = sum_x2[:,:,iiv,iid]
                    iii_sum_y2 = sum_y2[:,:,iiv,iid]
                    iii_rcorr = rcorr[:,:,iiv,iid]
                    rnum1 = np.multiply(ncnt[isfinite], iii_sum_xy[isfinite])
                    rnum2 = np.multiply(iii_sum_x[isfinite], iii_sum_y[isfinite])
                    rnum = rnum1 - rnum2
                    rden1 = np.multiply( ncnt[isfinite], iii_sum_x2[isfinite]) - np.square(iii_sum_x[isfinite])
                    rden2 = np.multiply( ncnt[isfinite], iii_sum_y2[isfinite]) - np.square(iii_sum_y[isfinite])
                    rden = np.multiply( np.sqrt(rden1) , np.sqrt(rden2) )
                    nozero = np.where(rden > 0 )
                    noz_rcorr = iii_rcorr[isfinite]
                    noz_rcorr[nozero] = rnum[nozero] / rden[nozero]
                    iii_rcorr[isfinite] = noz_rcorr
                    #print(np.all(iii_rcorr[isfinite][nozero] == noz_rcorr[nozero]))
                    #iii_rcorr[isfinite][nozero] = rnum[nozero] / rden[nozero]
                    #print(np.all(iii_rcorr[isfinite][nozero] == noz_rcorr[nozero]))
                    rcorr[:,:,iiv,iid] = iii_rcorr
        if ( sum_x.ndim == 2 ): 
            nv, nt = sum_x.shape
            rcorr = 0.0 * sum_x
            isfinite = np.where(ncnt > 1 )  # NEED AT LEAST TWO POINTS FOR CORRELATION
            for iiv in range(nv):
                rnum1 = np.multiply(ncnt[isfinite], sum_xy[iiv, isfinite])
                rnum2 = np.multiply(sum_x[iiv, isfinite], sum_y[iiv, isfinite])
                rnum = rnum1 - rnum2
                rden1 = np.multiply( ncnt[isfinite], sum_x2[iiv, isfinite]) - np.square(sum_x[iiv, isfinite])
                rden2 = np.multiply( ncnt[isfinite], sum_y2[iiv, isfinite]) - np.square(sum_y[iiv, isfinite])
                rden = np.multiply( np.sqrt(rden1) , np.sqrt(rden2) )
                nozero = np.where(rden > 0 )
                noz_rcorr = rcorr[iiv, isfinite]
                noz_rcorr[nozero] = rnum[nozero] / rden[nozero]
                rcorr[iiv, isfinite] = noz_rcorr
        if ( sum_x.ndim == 1 ): 
            nt = sum_x.shape[0]
            rcorr = 0.0 * sum_x
            isfinite = np.where(ncnt > 1 )  # NEED AT LEAST TWO POINTS FOR CORRELATION
            rnum1 = np.multiply(ncnt[isfinite], sum_xy[isfinite])
            rnum2 = np.multiply(sum_x[isfinite], sum_y[isfinite])
            rnum = rnum1 - rnum2
            rden1 = np.multiply( ncnt[isfinite], sum_x2[isfinite]) - np.square(sum_x[isfinite])
            rden2 = np.multiply( ncnt[isfinite], sum_y2[isfinite]) - np.square(sum_y[ isfinite])
            rden = np.multiply( np.sqrt(rden1) , np.sqrt(rden2) )
            nozero = np.where(rden > 0 )
            noz_rcorr = rcorr[isfinite]
            noz_rcorr[nozero] = rnum[nozero] / rden[nozero]
            rcorr[isfinite] = noz_rcorr
        if ( sum_x.ndim > 4 ):
            print('NOT CODED')
    return rcorr

#FIELDS=(obser_padded, field, varobs, varfld, covari, merror, sqrerr, abserr)
def find_analysis_errors(date_range, expt, variable, insuffix='GIOPS_orca025_currents.f2', ierror=1, ddeg=4.0, mp=True):

    date_start = date_range[0]
    date_final = date_range[1]
    date_inc = 1
    if ( len(date_range) == 3 ): date_inc = date_range[3]
    if ( len(date_range) > 3 ): 
        date_list = [check_date(idate,outtype=datetime.datetime) for idate in date_range]
    else:
        dates_list = create_dates(date_start, date_final, date_inc)
        
    (LONG, LATG), MNS_FIELDS, GDS_FIELDS, GCS, CNS, MN_LIST, GD_LIST, GC_LIST, CN_LIST = loop_analysis_dates(dates_list, expt, variable, insuffix=insuffix, ierror=ierror, ddeg=ddeg, mp=mp)

    MNA_FIELDS = calc_mean_field(MNS_FIELDS, GCS)
    GDA_FIELDS = calc_mean_field(GDS_FIELDS, CNS)

    sum_obs, sum_fld,  var_obs, var_fld, cov_obs_fld = get_correlation_fields(MNS_FIELDS)
    grd_sum_obs, grd_sum_fld, grd_var_obs, grd_var_fld, grd_cov_obs_fld = get_correlation_fields(GDA_FIELDS)

    RCORR = calc_correlation(GCS, sum_obs, sum_fld, cov_obs_fld, var_obs, var_fld)
    tVAROB = calc_variance(GCS, sum_obs, var_obs)
    tVARFD = calc_variance(GCS, sum_fld, var_fld)
    VAROB = MNA_FIELDS[2] - np.square(MNA_FIELDS[0])
    VARFD = MNA_FIELDS[3] - np.square(MNA_FIELDS[1])
    print('CHECK VARIANCE', VAROB-tVAROB, VARFD-tVARFD)
    SCORR = ( MNA_FIELDS[5] - VAROB - VARFD ) / (2*np.sqrt(VAROB)*np.sqrt(VARFD)) 
    print('CHECK CORRELATION', RCORR - SCORR )
    GRD_RCORR = calc_correlation(CNS, grd_sum_obs, grd_sum_fld, grd_cov_obs_fld, grd_var_obs, grd_var_fld)

    FIELDN=('observations', 'model equiv', 'obs var', 'mod var', 'cov obs var', 'mean error', 'rmse', 'abs error')
    print('GLOBAL MEANS')
    for ierr in range(len(MNA_FIELDS)):
       if ( FIELDN[ierr] == 'rmse' ):
           print(FIELDN[ierr], MNA_FIELDS[ierr], np.sqrt(MNA_FIELDS[ierr]))
       else:
           print( FIELDN[ierr], MNA_FIELDS[ierr])  

    return  MNA_FIELDS, GDA_FIELDS, RCORR, GRD_RCORR

def calculate_mean_rcorr(mn_list, gc_list):
        mna_timeii = calc_mean_field(mn_list, gc_list)
        sum_obs, sum_fld, var_obs, var_fld, cov_obs_fld = get_correlation_fields(mn_list)
        rcorr = calc_correlation(gc_list, sum_obs, sum_fld, cov_obs_fld, var_obs, var_fld)
        mna_timeii.append(rcorr)
        return mna_timeii
              
def compare_analysis_errors(date_range, expts, labels=None, insuffixs=['GIOPS_orca025_currents.f2','GIOPS_orca025_currents.f2'], ens_axes=[0, 0], ierror=1, ddeg=10.0, mp=True, maxtaylor=1.0):

    if ( not labels ): labels=expts
    expt0 = expts[0]
    expt1 = expts[1]

    insuffix0 = insuffixs[0]
    insuffix1 = insuffixs[1]
    
    if ( isinstance(date_range, type(None)) ):
        datelist1=find_dates_in_expt(expts[0], insuffixs[0])
        datelist2=find_dates_in_expt(expts[1], insuffixs[1])
        date_range = find_common_dates(datelist1, datelist2)
        
    dayofweek = check_date(date_range[0], outtype=datetime.datetime).weekday() 
    if ( dayofweek == 1 ):
        pass
    else:
        print('BEST if START DAY IS A TUESDAY')
        print('BEST_ESTIMATE ONLY ON TUESDAY')
    andiff = ( 1-dayofweek ) % 7
    date_atart = check_date(date_range[0], outtype=datetime.datetime)+datetime.timedelta(days=andiff)
    date_start = date_range[0]
    date_final = date_range[1]
    date_inc = 1
    if ( len(date_range) == 3 ): date_inc = date_range[3]
    if ( len(date_range) > 3 ): 
        dates_list = [check_date(idate,outtype=datetime.datetime) for idate in date_range]
    else:
        dates_list = create_dates(date_start, date_final, date_inc)
    dates_last = create_dates(date_atart, dates_list[-1], 7)
    ndates=len(dates_list)

    MNA_LIST = []
    MNA_LITT = []
    GDA_LIST = []
    CNS_LIST = []
    TAY_LIST = []  # obs_var, fld_var, corr, sqr_err
    
    FIELDN=('observations', 'model equiv', 'obs var', 'mod var', 'cov obs var', 'mean error', 'sqr err', 'abs error', 'corr')
    if ( ens_axes[0] != 0 ): FIELDN=('observations', 'model equiv', 'obs var', 'mod var', 'cov obs var', 'mean error', 'sqr err', 'abs error', 'crps', 'corr')
    TAYLOR=('obs var', 'mod var', 'corr', 'sqr err')
    for iexpt, expt in enumerate(expts):
      for variable in ['init_estimate', 'best_estimate']:
        ens_axis = ens_axes[iexpt]
        dates_this = dates_list
        if ( variable == 'best_estimate'): dates_this=dates_last
        if ( variable == 'best_estimate'): 
            if ( ens_axis != 0 ): ens_axis=-1
        
        (LONG, LATG), MNS_FIELDS, GDS_FIELDS, GCS, CNS, MN_LIST, GD_LIST, GC_LIST, CN_LIST = loop_analysis_dates(dates_this, expt, variable, insuffix=insuffixs[iexpt], ierror=ierror, ddeg=ddeg, mp=mp, ens_axis=ens_axis)

        #MNA_FIELDS = calc_mean_field(MNS_FIELDS, GCS)
        MNA_FIELDS = calculate_mean_rcorr(MNS_FIELDS, GCS)
        #GDA_FIELDS = calc_mean_field(GDS_FIELDS, CNS)
        GDA_FIELDS = calculate_mean_rcorr(GDS_FIELDS, CNS)

        #sum_obs, sum_fld,  var_obs, var_fld, cov_obs_fld = get_correlation_fields(MNS_FIELDS)
        #grd_sum_obs, grd_sum_fld, grd_var_obs, grd_var_fld, grd_cov_obs_fld = get_correlation_fields(GDA_FIELDS)

        #RCORR = calc_correlation(GCS, sum_obs, sum_fld, cov_obs_fld, var_obs, var_fld)
        #GRD_RCORR = calc_correlation(CNS, grd_sum_obs, grd_sum_fld, grd_cov_obs_fld, grd_var_obs, grd_var_fld)
        #MNA_FIELDS.append(RCORR)
        #GDA_FIELDS.append(GRD_RCORR)
        
        if ( not mp ):
            time0 = time.time()
            MNA_TIMESS = []
            for idate in range(len(MN_LIST)):
                mna_timeii = calculate_mean_rcorr(MN_LIST[idate], GC_LIST[idate])
                MNA_TIMESS.append(mna_timeii)
            print("SEQ TIME", time.time()-time0)
        else:
            time0 = time.time()
            nproc=min(num_cpus, len(MN_LIST))
            pool = multiprocessing.Pool(nproc)
            izip = list(zip( MN_LIST, GC_LIST ) )
            MNA_TIMESS = pool.starmap(calculate_mean_rcorr, izip)
            pool.close()
            pool.join()
            print("MP TIME", time.time()-time0)
        
        print('LENGTH', len(MNA_FIELDS), len(GDA_FIELDS), len(FIELDN))
        MNA_LIST.append(MNA_FIELDS)
        MNA_LITT.append(MNA_TIMESS)
        GDA_LIST.append(GDA_FIELDS)
        CNS_LIST.append(CNS)
        TAY_LIST.append([ MNA_FIELDS[FIELDN.index(taylor_element)] for taylor_element in TAYLOR ])

        print(expt+' '+variable+' GLOBAL MEANS')
        for ierr in range(len(MNA_FIELDS)):
            if ( FIELDN[ierr] == 'sqr err' ):
                print('rmse', MNA_FIELDS[ierr], np.sqrt(MNA_FIELDS[ierr]))
            else:
                print( FIELDN[ierr], MNA_FIELDS[ierr])  

    TAYLOR_LIST = []
    LABEL_LIST = []
    COLOR_LIST = []

    VARN = ['U', 'V', 'S']
    vari = ['init','best']
    clrs_expt = ['red', 'magenta', 'blue', 'cyan']
    for iexpt in range(len(expts)):
      for ivari in range(len(vari)):
        ielem=2*iexpt+ivari
        TAY_ELE = TAY_LIST[ielem]
        for ifld in range(2):
            taylor_points = [TAY_PTS[ifld,0] for TAY_PTS in TAY_ELE] 
            print(expts[iexpt], vari[ivari], VARN[ifld], 'taylor points', taylor_points)
            TAYLOR_LIST.append(taylor_points)
            variable = vari[ivari] 
            if ( variable == 'best' ):
               variable = 'Analysis'
            elif ( variable == 'init' ):
               variable = 'Trial'   
            LABEL_LIST.append(labels[iexpt]+' '+variable+' '+VARN[ifld])
            COLOR_LIST.append(clrs_expt[ielem])

    dates_string = check_date(dates_list[0], outtype=str)+'_'+check_date(dates_list[-1], outtype=str)
    expts_string = expts[0]+'_'+expts[1]
    outfile = 'ANAL_plots/'+ 'taylor'+'_'+expts_string+'_'+dates_string
    taylor_diagram.make_taylor_figure( TAYLOR_LIST, LABEL_LIST, 
         list_of_colors=COLOR_LIST, pltfile=outfile, maxlim=maxtaylor )    
    outfile = 'ANAL_plots/'+ 'taylorzoom'+'_'+expts_string+'_'+dates_string
    taylor_diagram.make_taylor_figure( TAYLOR_LIST, LABEL_LIST, 
         list_of_colors=COLOR_LIST, pltfile=outfile, maxlim='auto', addvalues=True )    

    #FIELDN=('observations', 'model equiv', 'obs var', 'mod var', 'cov obs var', 'mean error', 'sqr err', 'abs error', 'corr')
    # want 5, 6 & 8 (mean err / sqr err / corr)
    # want 5, 6 & 9 & 8 for ensemble ( mean err / sqr err / ccorr / crps )
    
    GDA_PLOT = []
    GLB_SUMP = []
    TIM_PLOT = []
    if ( ens_axes[0] == 0 ):
        nplts=3
        print('LENGTH GDA_LIST', len(GDA_LIST))
        for iplot in range(len(GDA_LIST)):
            GDA_PLOT.append([GDA_LIST[iplot][5], GDA_LIST[iplot][6], GDA_LIST[iplot][8]])
            GLB_SUMP.append([MNA_LIST[iplot][5], MNA_LIST[iplot][6], MNA_LIST[iplot][8]])
            tim_plot=[]
            print("LENGTH", len(MNA_LITT), len(MNA_LITT[0]), len(MNA_LITT[0][0]), len(MNA_LITT[0][0][6]))
            ndtemp=len(MNA_LITT[iplot])
            for idate in range(ndtemp):
                print('idate', idate)
                print('NUM ERRORS', idate, len(MNA_LITT[iplot][idate]))
                tim_plot.append([MNA_LITT[iplot][idate][5], MNA_LITT[iplot][idate][6], MNA_LITT[iplot][idate][8]])
            print('LENGTH tim_plot', len(tim_plot))
            TIM_PLOT.append(tim_plot)
    else:
        nplts=4
        for iplot in range(len(GDA_LIST)):
            GDA_PLOT.append([GDA_LIST[iplot][5], GDA_LIST[iplot][6], GDA_LIST[iplot][9], GDA_LIST[iplot][8]])
            GLB_SUMP.append([MNA_LIST[iplot][5], MNA_LIST[iplot][6], MNA_LIST[iplot][9], MNA_LIST[iplot][8]])
            tim_plot=[]
            ndtemp=len(MNA_LITT[iplot])
            for idate in range(ndtemp):
                print('idate', idate)
                print('NUM ERRORS', idate, len(MNA_LITT[iplot][idate]))
                tim_plot.append([MNA_LITT[iplot][idate][5], MNA_LITT[iplot][idate][6], MNA_LITT[iplot][idate][9], MNA_LITT[iplot][idate][8]])
            print('LENGTH tim_plot', len(tim_plot))
            TIM_PLOT.append(tim_plot)
        
    GDA_MERR = []
    GDA_RMSE = []
    GDA_CORR = []
    GDA_CRPS = []
    GLB_MERR = []
    GLB_RMSE = []
    GLB_CORR = []
    GLB_CRPS = []
    TIM_MERR = []
    TIM_RMSE = []
    TIM_CORR = []
    TIM_CRPS = []
    # GDA_PLOT index 0, 1 (init/best, expt 1); 2, 3 (init/best, expt 2)
    for ielem in range(nplts):
        print( 'bin count_check', np.all(CNS_LIST[0] == CNS_LIST[2]), np.all(CNS_LIST[1] == CNS_LIST[3]) )
        if ( ielem == 0 ):
            GDA_MERR.append( GDA_PLOT[0][ielem] - GDA_PLOT[2][ielem] )
            GDA_MERR.append( GDA_PLOT[1][ielem] - GDA_PLOT[3][ielem] )
            GLB_MERR.append( [GLB_SUMP[0][ielem], GLB_SUMP[2][ielem]] )
            GLB_MERR.append( [GLB_SUMP[1][ielem], GLB_SUMP[3][ielem]] )
            tim_merr = []
            ndtemp=len(TIM_PLOT[0]) 
            print('LENGTH TIM_PLOT', len(TIM_PLOT), len(TIM_PLOT[0]), len(TIM_PLOT[1]), len(TIM_PLOT[2]), len(TIM_PLOT[3]))
            print('LENGTH TIM_PLOT', ielem, len(TIM_PLOT[0][0]))
            ndtemp=len(TIM_PLOT[0]) 
            for idate in range(ndtemp):
                tim_merr.append( [ TIM_PLOT[0][idate][ielem], TIM_PLOT[2][idate][ielem] ] )
            TIM_MERR.append(tim_merr)
            tim_merr = []
            ndtemp=len(TIM_PLOT[1])
            for idate in range(ndtemp):
                tim_merr.append( [ TIM_PLOT[1][idate][ielem], TIM_PLOT[3][idate][ielem] ] )
            TIM_MERR.append(tim_merr)
        if ( ielem == 1 ):
            GDA_RMSE.append( np.sqrt(GDA_PLOT[0][ielem]) - np.sqrt(GDA_PLOT[2][ielem]) )
            GDA_RMSE.append( np.sqrt(GDA_PLOT[1][ielem]) - np.sqrt(GDA_PLOT[3][ielem]) )
            GLB_RMSE.append( [np.sqrt(GLB_SUMP[0][ielem]), np.sqrt(GLB_SUMP[2][ielem])] )
            GLB_RMSE.append( [np.sqrt(GLB_SUMP[1][ielem]), np.sqrt(GLB_SUMP[3][ielem])] )
            tim_rmse = []
            ndtemp=len(TIM_PLOT[0]) 
            print('LENGTH TIM_PLOT', len(TIM_PLOT), len(TIM_PLOT[0]), len(TIM_PLOT[1]), len(TIM_PLOT[2]), len(TIM_PLOT[3]))
            print('LENGTH TIM_PLOT', ielem, len(TIM_PLOT[0][0]))
            ndtemp=len(TIM_PLOT[0]) 
            for idate in range(ndtemp):
                tim_rmse.append( [ np.sqrt(TIM_PLOT[0][idate][ielem]), np.sqrt(TIM_PLOT[2][idate][ielem]) ] )
            TIM_RMSE.append(tim_rmse)
            tim_rmse = []
            ndtemp=len(TIM_PLOT[1])
            for idate in range(ndtemp):
                tim_rmse.append( [ np.sqrt(TIM_PLOT[1][idate][ielem]), np.sqrt(TIM_PLOT[3][idate][ielem]) ] )
            TIM_RMSE.append(tim_rmse)
        elif (ielem == 2 ):
            GDA_CORR.append( GDA_PLOT[0][ielem] - GDA_PLOT[2][ielem] )
            GDA_CORR.append( GDA_PLOT[1][ielem] - GDA_PLOT[3][ielem] )
            GLB_CORR.append( [ GLB_SUMP[0][ielem], GLB_SUMP[2][ielem]] )
            GLB_CORR.append( [ GLB_SUMP[1][ielem], GLB_SUMP[3][ielem]] )
            tim_corr = []
            ndtemp=len(TIM_PLOT[0])
            for idate in range(ndtemp):
                tim_corr.append( [ TIM_PLOT[0][idate][ielem], TIM_PLOT[2][idate][ielem] ] )
            TIM_CORR.append(tim_corr)
            tim_corr = []
            ndtemp=len(TIM_PLOT[1])
            for idate in range(ndtemp):
                tim_corr.append( [ TIM_PLOT[1][idate][ielem], TIM_PLOT[3][idate][ielem] ] )
            TIM_CORR.append(tim_corr)
        elif (ielem == 3 ):
            GDA_CRPS.append( GDA_PLOT[0][ielem] - GDA_PLOT[2][ielem] )
            GDA_CRPS.append( GDA_PLOT[1][ielem] - GDA_PLOT[3][ielem] )
            GLB_CRPS.append( [ GLB_SUMP[0][ielem], GLB_SUMP[2][ielem]] )
            GLB_CRPS.append( [ GLB_SUMP[1][ielem], GLB_SUMP[3][ielem]] )
            tim_crps = []
            ndtemp=len(TIM_PLOT[0])
            for idate in range(ndtemp):
                tim_crps.append( [ TIM_PLOT[0][idate][ielem], TIM_PLOT[2][idate][ielem] ] )
            TIM_CRPS.append(tim_crps)
            tim_crps = []
            ndtemp=len(TIM_PLOT[1])
            for idate in range(ndtemp):
                tim_crps.append( [ TIM_PLOT[1][idate][ielem], TIM_PLOT[3][idate][ielem] ] )
            TIM_CRPS.append(tim_crps)

    CMAP = cmap_full_field
    AMAP = cmap_anom_field
    BMAP = cmap_anom_rield
    for ivari in range(len(vari)):
        for ifld in range(len(VARN)):
            grd_merr = GDA_MERR[ivari][:,:, ifld, 0]
            grd_rmse = GDA_RMSE[ivari][:,:, ifld, 0]    
            grd_corr = GDA_CORR[ivari][:,:, ifld, 0] 
            glb_merr = [GLB_MERR[ivari][iexpt][ifld] for iexpt in range(len(expts))]
            glb_rmse = [GLB_RMSE[ivari][iexpt][ifld] for iexpt in range(len(expts))]
            glb_corr = [GLB_CORR[ivari][iexpt][ifld] for iexpt in range(len(expts))]
            glb_crps = [0, 0]
            if ( len(GLB_CRPS)  > 0 ):
                grd_crps = GDA_CRPS[ivari][:,:, ifld, 0]    
                glb_crps = [GLB_CRPS[ivari][iexpt][ifld] for iexpt in range(len(expts))]
            
            #print('glb_corr', glb_corr)
            #print('GLB_CORR', GLB_CORR)
            #print('glb_rmse', glb_rmse)
            #print('GLB_RMSE', GLB_RMSE)
            #print('GLB_CORR', GLB_CORR[ivari][ifld][:], [GLB_CORR[ivari][iexpt][ifld] for iexpt in range(2)])
            if ( len(GLB_CRPS) > 0 ):
                print('glb_crps', len(glb_crps), glb_crps[0].shape, glb_crps)
                print('GLB_CRPS', len(GLB_CRPS), GLB_CRPS[0][0].shape, GLB_CRPS)
            rcor_str0 = " (%.4f)" % glb_corr[0]
            rmse_str0 = " (%.4f)" % glb_rmse[0]
            rmne_str0 = " (%.4f)" % glb_merr[0]
            rcor_str1 = " (%.4f)" % glb_corr[1]
            rmse_str1 = " (%.4f)" % glb_rmse[1]
            rmne_str1 = " (%.4f)" % glb_merr[1]
            rcps_str0 = " (%.4f)" % glb_crps[0]
            rcps_str1 = " (%.4f)" % glb_crps[1]
            maxrmse = np.nanmax(np.abs(grd_rmse))
            maxcorr = np.nanmax(np.abs(grd_corr))
            maxmerr = np.nanmax(np.abs(grd_merr))
            if ( len(GLB_CRPS) > 0 ): 
                maxcrps = np.max(np.abs(grd_crps))
            else:
                maxcrps = 0
            LEV_UNIT = np.array([-0.9, -0.7, -0.5, -0.3, -0.1, 0.1, 0.3, 0.5, 0.7, 0.9])
            TIC_UNIT = np.array([-1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
            LEV_MERR = maxmerr*LEV_UNIT
            LEV_RMSE = maxrmse*LEV_UNIT
            LEV_CORR = maxcorr*LEV_UNIT
            LEV_CRPS = maxcrps*LEV_UNIT
            LEV_CORR = 0.1*LEV_UNIT
            LEV_MERR = 0.05*LEV_UNIT
            LEV_RMSE = 0.05*LEV_UNIT
            LEV_CRPS = 0.05*LEV_UNIT
            TIC_CORR = 0.1*TIC_UNIT
            TIC_MERR = 0.05*TIC_UNIT
            TIC_RMSE = 0.05*TIC_UNIT
            TIC_CRPS = 0.05*TIC_UNIT
            print(ivari, ifld, type(grd_rmse), type(grd_corr))
            print(ivari, ifld, grd_rmse.shape, grd_corr.shape)
            print(LONG.shape, LATG.shape) 
            print('MAX', maxrmse, maxcorr,maxcrps)  

            title='MERR difference '+labels[0]+rmne_str0+' - '+labels[1]+rmne_str1+' '+vari[ivari] +' '+VARN[ifld]
            ofile='ANAL_plots/'+ 'merr'+'_'+expts_string+'_'+vari[ivari]+'_'+VARN[ifld]+'_'+dates_string+'.png'
            cplot.pcolormesh(LONG, LATG, grd_merr, levels=LEV_MERR, ticks=TIC_MERR, obar='horizontal', title=title, outfile=ofile, make_global=True, project='PlateCarree', cmap=BMAP, add_gridlines=True)
            title='RMSE difference '+labels[0]+rmse_str0+' - '+labels[1]+rmse_str1+' '+vari[ivari] +' '+VARN[ifld]
            title='RMSE difference '+labels[0]+rmse_str0+' - '+labels[1]+rmse_str1+' '+VARN[ifld]
            ofile='ANAL_plots/'+ 'rmse'+'_'+expts_string+'_'+vari[ivari]+'_'+VARN[ifld]+'_'+dates_string+'.png'
            cplot.pcolormesh(LONG, LATG, grd_rmse, levels=LEV_RMSE, ticks=TIC_RMSE, obar='horizontal', title=title, outfile=ofile, make_global=True, project='PlateCarree', cmap=BMAP, add_gridlines=True)
            title='rcorr difference '+labels[0]+rcor_str0+' - '+labels[1]+rcor_str1+' '+vari[ivari] +' '+VARN[ifld]
            ofile='ANAL_plots/'+ 'corr'+'_'+expts_string+'_'+vari[ivari]+'_'+VARN[ifld]+'_'+dates_string+'.png'
            cplot.pcolormesh(LONG, LATG, grd_corr, levels=LEV_CORR,ticks=TIC_CORR,obar='horizontal', title=title, outfile=ofile, make_global=True, project='PlateCarree', cmap=AMAP, add_gridlines=True)
            if ( len(GLB_CRPS) > 0 ):
                title='crps difference '+labels[0]+rcps_str0+' - '+labels[1]+rcps_str1+' '+vari[ivari] +' '+VARN[ifld]
                ofile='ANAL_plots/'+ 'crps'+'_'+expts_string+'_'+vari[ivari]+'_'+VARN[ifld]+'_'+dates_string+'.png'
                cplot.pcolormesh(LONG, LATG, grd_crps, levels=LEV_CRPS,ticks=TIC_CRPS,obar='horizontal', title=title, outfile=ofile, make_global=True, project='PlateCarree', cmap=BMAP, add_gridlines=True)
            
            tim_merr = []
            tim_rmse = []
            tim_corr = []
            tim_crps = []
            for iexpt in range(2):
                print('Length', ivari, len(TIM_RMSE[ivari]), len(TIM_CORR[ivari]))
                tim_merr.append([ TIM_MERR[ivari][idate][iexpt][ifld] for idate in range(len(TIM_MERR[ivari]))])
                tim_rmse.append([ TIM_RMSE[ivari][idate][iexpt][ifld] for idate in range(len(TIM_RMSE[ivari]))])
                tim_corr.append([ TIM_CORR[ivari][idate][iexpt][ifld] for idate in range(len(TIM_CORR[ivari]))])
                if ( nplts == 4 ): tim_crps.append([ TIM_CRPS[ivari][idate][iexpt][ifld] for idate in range(len(TIM_CRPS[ivari]))])
                
            title='merr timeseries '+labels[0]+rmne_str0+' - '+labels[1]+rmne_str1+' '+vari[ivari] +' '+VARN[ifld]
            ofile='ANAL_plots/'+ 'T_merr'+'_'+expts_string+'_'+vari[ivari]+'_'+VARN[ifld]+'_'+dates_string
            fig, axe = plt.subplots()
            dates_this = dates_list
            if ( vari[ivari] == 'best'): dates_this=dates_last
            print(ivari, vari[ivari], 'length to plot', len(dates_this), len(tim_merr[0]), len(tim_merr[1]))
            axe.plot(dates_this, tim_merr[0], color='r', label=labels[0])
            axe.plot(dates_this, tim_merr[1], color='b', label=labels[1])
            fig.savefig(ofile+'.png')
            fig.savefig(ofile+'.pdf')
            plt.close(fig)
           
            title='rmse timeseries '+labels[0]+rmse_str0+' - '+labels[1]+rmse_str1+' '+vari[ivari] +' '+VARN[ifld]
            ofile='ANAL_plots/'+ 'T_rmse'+'_'+expts_string+'_'+vari[ivari]+'_'+VARN[ifld]+'_'+dates_string
            fig, axe = plt.subplots()
            dates_this = dates_list
            nwin=31
            if ( vari[ivari] == 'best'): 
                dates_this=dates_last
                nwin=5
            print(ivari, ifld, vari[ivari], VARN[ifld], 'length to plot', len(dates_this), len(tim_rmse[0]), len(tim_rmse[1]))
            axe.plot(dates_this, glb_rmse[0]*np.ones(len(tim_rmse[0])), color='r', linestyle='--', linewidth=1.0)
            axe.plot(dates_this, glb_rmse[1]*np.ones(len(tim_rmse[1])), color='b', linestyle='--', linewidth=1.0)
            axe.plot(dates_this, tim_rmse[0], color='r', linewidth=0.1)
            axe.plot(dates_this, tim_rmse[1], color='b', linewidth=0.1)
            print(type(tim_rmse[0]), type(tim_rmse[1]), np.array(tim_rmse[0]).shape, np.array(tim_rmse[1]).shape)
            try:
                axe.plot(dates_this, calc_running_mean(np.array(tim_rmse[0]), nwin=nwin), color='r', linewidth=2.0, label=labels[0])
                axe.plot(dates_this, calc_running_mean(np.array(tim_rmse[1]), nwin=nwin), color='b', linewidth=2.0, label=labels[1])
                axe.legend()
            except:
                print(traceback.print_exc())
            axe.set_title(title)
            fig.savefig(ofile+'.png')
            fig.savefig(ofile+'.pdf')
            plt.close(fig)

            title='rcorr timeseries '+labels[0]+rmse_str0+' - '+labels[1]+rmse_str1+' '+vari[ivari] +' '+VARN[ifld]
            ofile='ANAL_plots/'+ 'T_corr'+'_'+expts_string+'_'+vari[ivari]+'_'+VARN[ifld]+'_'+dates_string
            fig, axe = plt.subplots()
            dates_this = dates_list
            if ( vari[ivari] == 'best'): dates_this=dates_last
            axe.plot(dates_this, glb_corr[0]*np.ones(len(tim_corr[0])), color='r', linestyle='--', linewidth=1.0)
            axe.plot(dates_this, glb_corr[1]*np.ones(len(tim_corr[1])), color='b', linestyle='--', linewidth=1.0)
            axe.plot(dates_this, tim_corr[0], color='r', linewidth=0.1)
            axe.plot(dates_this, tim_corr[1], color='b', linewidth=0.1)
            try:
                axe.plot(dates_this, calc_running_mean(np.array(tim_corr[0]), nwin=nwin), color='r', linewidth=2.0, label=labels[0])
                axe.plot(dates_this, calc_running_mean(np.array(tim_corr[1]), nwin=nwin), color='b', linewidth=2.0, label=labels[1])
                axe.legend()
            except:
                print(traceback.print_exc())
            axe.set_title(title)
            fig.savefig(ofile+'.png')
            fig.savefig(ofile+'.pdf')
            plt.close(fig)

            if ( nplts == 4 ):
                title='crps timeseries '+labels[0]+rmse_str0+' - '+labels[1]+rmse_str1+' '+vari[ivari] +' '+VARN[ifld]
                ofile='ANAL_plots/'+ 'T_corr'+'_'+expts_string+'_'+vari[ivari]+'_'+VARN[ifld]+'_'+dates_string
                fig, axe = plt.subplots()
                dates_this = dates_list
                if ( vari[ivari] == 'best'): dates_this=dates_last
                axe.plot(dates_this, glb_crps[0]*np.ones(len(tim_crps[0])), color='r', linestyle='--', linewidth=1.0)
                axe.plot(dates_this, glb_crps[1]*np.ones(len(tim_crps[1])), color='b', linestyle='--', linewidth=1.0)
                axe.plot(dates_this, tim_crps[0], color='r',linewidth=0.1)
                axe.plot(dates_this, tim_crps[1], color='b',linewidth=0.1)
                try:
                    axe.plot(dates_this, calc_running_mean(np.array(tim_crps[0]),nwin=nwin), color='r', linewidth=2.0, label=labels[0])
                    axe.plot(dates_this, calc_running_mean(np.array(tim_crps[1]),nwin=nwin), color='b', linewidth=2.0, label=labels[1])
                    axe.legend()
                except:
                    print(traceback.print_exc())
                axe.set_title(title)
                fig.savefig(ofile+'.png')
                fig.savefig(ofile+'.pdf')
            plt.close(fig)

    return  MNA_LIST, GDA_LIST, MNA_LITT

#FIELDS=(obser_padded, field, varobs, varfld, covari, merror, sqrerr, abserr, cnt)
def plot_anal_timeseries(expts, labels=None, date_range=None, filepre='ANAL_plots/timeseries', colors=['r','b','m','c']):
    date_min=None
    date_max=None
    if ( date_range ):
        (date_min, date_max) = [check_date(date, outtype=datetime.datetime) for date in date_range]
        
    if ( not labels ):
       labels=expts[:]

    plots=[]
    for variable in ['init', 'best']:
        for fld in ['U', 'V']:
            fig, axe = plt.subplots(3)
            addtitle=' '+variable+' '+fld
            axe[0].set_title('Mean Error'+addtitle)
            axe[1].set_title('RMSE Error'+addtitle)
            axe[2].set_title('Correlation'+addtitle)
            plots.append([fig, axe])
        
        
    for iexpt, expt in enumerate(expts):
        iplot=-1
        label=labels[iexpt]
        color=colors[iexpt%len(colors)]
        for variable in ['init_estimate', 'best_estimate']:
            for fld in ['U', 'V']: 
                iplot=iplot+1
                fig, axe = plots[iplot]
                errfile='ANAL_data/'+fld+'_'+variable+'_'+expt+'.dat'
                print(errfile)
                intdate, errors = datadatefile.read_file(errfile)
                dates = datadatefile.convert_strint_datelist(intdate)

                if ( date_range ):
                    new_dates=[]
                    new_errors=[]
                    (date_min, date_max) = date_range  
                    for idate,date in enumerate(dates):
                        if ( ( date >= date_min ) and (date <= date_max ) ):
                            new_dates.append(date)
                            new_errors.append(errors[:,idate])
                    dates = new_dates
                    errors = np.transpose( np.array(new_errors) )  # np.array will flip the time, variable indices.
                    print(len(dates))
                    
                cnt=errors[-1]
                mne_error, sqr_error, one_error = calc_mean_field(errors, cnt)[[5,6,8],:]
                print(1, '==', np.mean(one_error), '?')
                rms_error = np.sqrt(sqr_error)
                sum_obs, sum_fld,  var_obs, var_fld, cov_obs_fld =  get_correlation_fields(errors)
                rcorr = calc_correlation(cnt, sum_obs, sum_fld, cov_obs_fld, var_obs, var_fld) 
                
                axe[0].plot(dates, mne_error, color=color)
                axe[1].plot(dates, rms_error, color=color)
                axe[2].plot(dates, rcorr, label=label, color=color)

    iplot = -1
    for variable in ['init', 'best']:
        for fld in ['U', 'V']:
            filename=filepre+'_'+variable+'_'+fld
            iplot = iplot+1
            fig, axe = plots[iplot]
            axe[2].legend()
            fig.savefig(filename+'.pdf')
            fig.savefig(filename+'.png')
            plt.close(fig)
       
    return
        
def make_error_grid( NSHAPE, ddeg=4 ):
   grid_lon, grid_lat, lon_bin, lat_bin, grid_sum, grid_cnt = cplot.make_bin_grid(ddeg=ddeg, central_longitude=0)
   nlon, nlat = grid_sum.shape
   grid = np.zeros( (nlon, nlat)+NSHAPE[1:] )
   return (grid_lon, grid_lat), grid, grid_cnt

def make_error_grids( NSHAPE, ddeg=4 ):
   numobs, numvars, numfcsts, numdeps = NSHAPE
   grid_lon, grid_lat, lon_bin, lat_bin, grid_sum, grid_cnt = cplot.make_bin_grid(ddeg=ddeg, central_longitude=0)
   nlon, nlat = grid_sum.shape
   grd_beste = np.zeros(( nlon, nlat, numvars, numdeps ))
   grd_inite = np.zeros(( nlon, nlat, numvars, numdeps ))
   grd_fcstv = np.zeros(( nlon, nlat, numvars, numfcsts, numdeps ))
   grd_persi = np.zeros(( nlon, nlat, numvars, numfcsts, numdeps )) 
   return (grid_lon, grid_lat), (grd_beste, grd_inite, grd_fcstv, grd_persi), grid_cnt

def add_error_grids_system(indir, insuffix, date, ierror=1, ddeg=4, etype='mean', magnitude=False):
    (LON, LAT, depth, plf_code), (obser, beste, inite, fcstv, persi) = load_3errors_from_obsfile(date, indir, insuffix, ierror=ierror, etype=etype, magnitude=magnitude)
    (lon_grid, lat_grid), (Gbeste, Ginite, Gfcstv, Gpersi), CNT = make_error_grids(fcstv.shape, ddeg=ddeg)
    (LONG, LATG), (Gbeste, Ginite, Gfcstv, Gpersi), CNT = add_error_grids( (LON,LAT), (beste, inite, fcstv, persi), (Gbeste, Ginite, Gfcstv, Gpersi), CNT, ddeg=ddeg)     
    return (LONG, LATG), (Gbeste, Ginite, Gfcstv, Gpersi), CNT
       
def add_error_grids( LONLAT, ERRORS, GRDS, CNT, ddeg=4):
    LON, LAT = LONLAT
    
    (Ebeste, Einite, Efcstv, Epersi) = ERRORS
    (Gbeste, Ginite, Gfcstv, Gpersi) = GRDS
    nobs, nvars, nfcst, nlevs = Efcstv.shape

    for iobs in range(nobs):
        LONG, LATG, __, WGT = cplot.binfldsum(np.array([LON[iobs]]),np.array([LAT[iobs]]), np.ones(1), ddeg=ddeg)
        CNT = CNT + WGT
        for ivars in range(nvars):
            for ilevs in range(nlevs):
                Gbeste[:,:,ivars, ilevs] = Gbeste[:,:,ivars, ilevs]+WGT*Ebeste[iobs, ivars, ilevs]
                Ginite[:,:,ivars, ilevs] = Ginite[:,:,ivars, ilevs]+WGT*Einite[iobs, ivars, ilevs]
                for ifcst in range(nfcst):
                    Gfcstv[:,:,ivars,ifcst,ilevs] = Gfcstv[:,:,ivars,ifcst,ilevs]+WGT*Efcstv[iobs, ivars, ifcst, ilevs]
                    Gpersi[:,:,ivars,ifcst,ilevs] = Gpersi[:,:,ivars,ifcst,ilevs]+WGT*Epersi[iobs, ivars, ifcst, ilevs]
    
    return (LONG, LATG), (Gbeste, Ginite, Gfcstv, Gpersi), CNT

def add_error_grid_list( LONLAT, error_list, egrid_list, cnt, ddeg=4):
    LON, LAT = LONLAT
    ngrid_list = []
    for ierror, error in enumerate(error_list):
        egrid = egrid_list[ierror]
        (LONG, LATG), iegrid, icnt = add_error_grid( LONLAT, error, egrid, cnt, ddeg=ddeg)
        ngrid_list.append(iegrid)
        if ( ierror == 0 ): ncnt=icnt
    return (LONG, LATG), ngrid_list, ncnt       
        
def add_error_grid( LONLAT, error, egrid, cnt, ddeg=4):
    LON, LAT = LONLAT
    
    ndims=error.ndim
    if ( ndims == 3 ):
        nobs, nvars, nlevs = error.shape
    if ( ndims == 4 ):
        nobs, nvars, nfcst, nlevs = error.shape
    if ( ndims == 5 ):
        nobs, nvars, nfcst, nensm, nlevs = error.shape

    for iobs in range(nobs):
        LONG, LATG, __, wgt = cplot.binfldsum(np.array([LON[iobs]]),np.array([LAT[iobs]]), np.ones(1), ddeg=ddeg)
        cnt = cnt + wgt
        for ivars in range(nvars):
            for ilevs in range(nlevs):
                if ( ndims == 3 ):
                    egrid[:, :, ivars, ilevs] = egrid[:, :, ivars, ilevs] + wgt*error[iobs, ivars, ilevs]
                else:
                    for ifcst in range(nfcst):
                        if ( ndims == 4 ):
                            egrid[:,:,ivars,ifcst,ilevs] = egrid[:,:,ivars,ifcst,ilevs] + wgt*error[iobs, ivars, ifcst, ilevs]
                        else:
                            for iensm in range(nensm):
                                egrid[:,:,ivars,ifcst,iensm,ilevs] = egrid[:,:,ivars,ifcst,iensm, ilevs] + wgt*error[iobs, ivars, ifcst, iensm, ilevs]
    
    return (LONG, LATG), egrid, cnt

def fin_error_grids( GRDS, CNT, ddeg=4):
    
    (Gbeste, Ginite, Gfcstv, Gpersi) = GRDS
    nlon, nlat, nvars, nfcst, nlevs = Gfcstv.shape

    izer = np.where(CNT == 0)
    NCNT = CNT.copy()
    NCNT[izer] = 1.0   
    for ivars in range(nvars):
        for ilevs in range(nlevs):
            #print(Gbeste[:,:,ivars,ilevs].shape, NCNT.shape)
            Gbeste[:,:,ivars,ilevs] = Gbeste[:,:,ivars,ilevs]/NCNT
            Ginite[:,:,ivars,ilevs] = Ginite[:,:,ivars,ilevs]/NCNT
            for ifcst in range(nfcst):
                Gfcstv[:,:,ivars,ifcst,ilevs] = Gfcstv[:,:,ivars,ifcst,ilevs]/NCNT
                Gpersi[:,:,ivars,ifcst,ilevs] = Gpersi[:,:,ivars,ifcst,ilevs]/NCNT
    
    return (Gbeste, Ginite, Gfcstv, Gpersi)

def fin_error_grid( grid, cnt ):
    
    if ( grid.ndim == 4 ):
        nlon, nlat, nvars, nlevs = grid.shape
    if ( grid.ndim == 5 ):
        nlon, nlat, nvars, nfcst, nlevs = grid.shape
    if ( grid.ndim == 6 ):
        nlon, nlat, nvars, nfcst, nlevs = grid.shape

    izer = np.where(cnt == 0)
    ncnt = cnt.copy()
    ncnt[izer] = 1.0   
    for ivars in range(nvars):
        for ilevs in range(nlevs):
            if ( grid.ndim == 4 ):
                grid[:,:,ivars,ilevs] = grid[:,:,ivars,ilevs]/ncnt
            else:
                for ifcst in range(nfcst):
                    if ( grid.ndim == 5 ):
                        grid[:,:,ivars,ifcst,ilevs] = grid[:,:,ivars,ifcst,ilevs]/ncnt
                    else:
                        for iensm in range(nensm):
                            grid[:,:,ivars,ifcst,iensm,ilevs] = grid[:,:,ivars,ifcst,iensm,ilevs]/NCNT
    
    return grid

indir='CLASS4_currents_CCMEP_FILT2'
insuffix='GIOPS_orca025_currents-filter'

def bin_over_dates_3errors(indir, insuffix, dates, ierror=1, ddeg=4, mp=False, etype='mean', magnitude=False):
    time0=time.time()
    GGLIST = []
    CNLIST = []
    ndates=len(dates)
    if ( not mp ): 
        for date in dates:
            time1=time.time()
            FULL_RESULT = add_error_grids_system(indir, insuffix, date, ierror=ierror, etype=etype,magnitude=magnitude,ddeg=ddeg)
            (LONG, LATG), GLIST, CNT = FULL_RESULT
            GGLIST.append(GLIST)
            CNLIST.append(CNT)
            print(time.time() - time1)
    else:
        nproc=min(num_cpus, ndates)
        pool = multiprocessing.Pool(nproc)
        izip = list(zip( itertools.repeat(indir), itertools.repeat(insuffix), dates))
        FULL_RESULT_LIST = pool.starmap(partial(add_error_grids_system, ierror=ierror, ddeg=ddeg, etype=etype), izip)
        pool.close()
        pool.join()
        for idate in range(ndates):
            (LONG, LATG), GLIST, CNT = FULL_RESULT_LIST[idate]
            GGLIST.append(GLIST)
            CNLIST.append(CNT)
        
    for idate in range(ndates):
        if ( idate == 0 ):
            lon_grid, lat_grid = LONG, LATG
            CNT = CNLIST[idate]
            Gbeste, Ginite, Gfcstv, Gpersi = GGLIST[idate]
        else:
            CNA = CNLIST[idate]
            CNT = CNT + CNA
            Abeste, Ainite, Afcstv, Apersi = GGLIST[idate]
            Gbeste = Gbeste+Abeste
            Ginite = Ginite+Ainite
            Gfcstv = Gfcstv+Afcstv
            Gpersi = Gpersi+Apersi
       
    GLIST = fin_error_grids( (Gbeste, Ginite, Gfcstv, Gpersi), CNT )
    timet = (time.time() - time0) / len(dates)
    print(timet)
    return  (lon_grid, lat_grid),  GLIST     

def make_spatial_plots_3errors(indir, insuffix, date_range, outdir='EPLOTS/', 
                       magnitude=False, levels=np.arange(0, 0.11, 0.01), alevels=np.arange(-.475, 0.5, 0.05), ierror=0, ilev=0, etype='mean', mp=True):

    etypp = etype
    if ( etype == 'square'): etypp='rmse'
    date_start = date_range[0]
    date_final = date_range[1]
    date_inc = 1
    if ( len(date_range) == 3 ): date_inc = date_range[2]
    dates_list = create_dates(date_start, date_final, date_inc)

    (lon_grid, lat_grid),  (Ebeste, Einite, Efcstv, Epersi) = bin_over_dates_3errors(indir, insuffix, dates_list, ierror=ierror, etype=etype, mp=mp, magnitude=magnitude)
    datestr_start = check_date(date_start, outtype=str)
    datestr_final = check_date(date_final, outtype=str)
    datestr_spans = datestr_start+'_'+datestr_final
    datestr_spass = datestr_start+' to '+datestr_final
    
    nlons, nlats, nvars, nfcst, nlevs = Efcstv.shape

    US = ['U', 'V', 'S', 'A']
    BEST=['best', 'init']
    FEST=['fcst', 'pers']
    FFLD=[Efcstv, Epersi]
    
    CMAP = cmap_full_field
    AMAP = cmap_anom_field
    CLEV = levels
    ALEV = alevels
    if ( etype == 'mean' ): 
        CLEV=alevels
        CMAP=cmap_anom_field
    for ifield, field in enumerate([Ebeste, Einite]):
       for ivar in [0, 1]:
           USS = US[ivar+2*magnitude]
           outfile=outdir+USS+'_'+etypp+'_'+BEST[ifield]+'_'+datestr_spans+'.png'
           title = etypp+' error for '+USS+' '+BEST[ifield]+' '+datestr_spass
           plt_fld = field[:, :, ivar, ilev]
           if ( etype == 'square' ): plt_fld = np.sqrt(plt_fld)
           cplot.pcolormesh(lon_grid, lat_grid, plt_fld, levels=CLEV,obar='horizontal', title=title, outfile=outfile, make_global=True, project='PlateCarree', cmap=CMAP, add_gridlines=True)    
           if ( ifield == 1 ):
               outfile=outdir+'D'+USS+'_'+etypp+'_'+BEST[ifield]+'_'+datestr_spans+'.png'
               title = etype+' error difference for '+USS+' '+BEST[ifield]+' '+datestr_spass
               plt_fld = field[:, :, ivar, ilev]-Ebeste[:, :, ivar, ilev]
               #if ( etype == 'square' ): plt_fld = np.sqrt(plt_fld)
               cplot.pcolormesh(lon_grid, lat_grid, plt_fld, levels=ALEV,obar='horizontal', title=title, outfile=outfile, make_global=True, project='PlateCarree', cmap=AMAP, add_gridlines=True)    

           for ifcst in range(nfcst):
               dstr=str(ifcst+1).zfill(2)
               outfile=outdir+USS+'_'+etypp+'_'+FEST[ifield]+'_'+dstr+'_'+datestr_spans+'.png'
               title = etypp+' error for '+USS+' '+FEST[ifield]+' day '+dstr+' '+datestr_spass
               plt_fld = FFLD[ifield][:, :, ivar, ifcst, ilev]
               if ( etype == 'square' ): plt_fld = np.sqrt(plt_fld)
               cplot.pcolormesh(lon_grid, lat_grid, plt_fld, levels=CLEV, obar='horizontal', title=title, outfile=outfile, make_global=True, project='PlateCarree', cmap=CMAP, add_gridlines=True)    

           for ifcst in range(nfcst):
               dstr=str(ifcst+1).zfill(2)
               outfile=outdir+'D'+USS+'_'+etypp+'_'+FEST[ifield]+'_'+dstr+'_'+datestr_spans+'.png'
               title = etype+' error difference for '+USS+' '+FEST[ifield]+' day '+dstr+' '+datestr_spass
               plt_fld = FFLD[ifield][:, :, ivar, ifcst, ilev]-field[:, :, ivar, ilev]
               #if ( etype == 'square' ): plt_fld = np.sqrt(plt_fld)
               cplot.pcolormesh(lon_grid, lat_grid, plt_fld, levels=ALEV, obar='horizontal', title=title, outfile=outfile, make_global=True, project='PlateCarree', cmap=AMAP)    
    return           
    
def make_Dspatial_plots_3errors(indir, insuffix, date_range, outdir='EPLOTS/', 
                       magnitude=False, levels=np.arange(-.095, 0.1, 0.01), ierror=[0,1], ilev=0, etype='mean', mp=True):

    etypp = etype
    if ( etype == 'square'): etypp='sqre'
    date_start = date_range[0]
    date_final = date_range[1]
    date_inc = 1
    if ( len(date_range) == 3 ): date_inc = date_range[2]
    dates_list = create_dates(date_start, date_final, date_inc)

    (lon_grid, lat_grid),  (Ebeste, Einite, Efcstv, Epersi) = bin_over_dates_3errors(indir, insuffix, dates_list, ierror=ierror[0], magnitude=magnitude, etype=etype, mp=True)
    (lon_grid, lat_grid),  (Fbeste, Finite, Ffcstv, Fpersi) = bin_over_dates_3errors(indir, insuffix, dates_list, ierror=ierror[1], magnitude=magnitude, etype=etype, mp=True)

    datestr_start = check_date(date_start, outtype=str)
    datestr_final = check_date(date_final, outtype=str)
    datestr_spans = datestr_start+'_'+datestr_final
    datestr_spass = datestr_start+' to '+datestr_final
    
    nlons, nlats, nvars, nfcst, nlevs = Efcstv.shape

    US = ['U', 'V', 'S', 'A']
    BEST=['best', 'init']
    FEST=['fcst', 'pers']
    FFLD=[Efcstv, Epersi]
    BFLD=[Fbeste, Finite]
    CFLD=[Ffcstv, Fpersi]
    for ifield, field in enumerate([Ebeste, Einite]):
       for ivar in [0, 1]:
           USS = US[ivar+2*magnitude]
           outfile=outdir+USS+'_'+etypp+'_'+BEST[ifield]+'_'+datestr_spans+'.png'
           title = etypp+' error for '+USS+' '+BEST[ifield]+' '+datestr_spass
           plt_fld = field[:, :, ivar, ilev]-BFLD[ifield][:, :, ivar, ilev]
           #if ( etype == 'square' ): plt_fld = np.sqrt(plt_fld)
           CLEV=levels
           CMAP=cmap_anom_field
           cplot.pcolormesh(lon_grid, lat_grid, plt_fld, levels=CLEV, obar='horizontal', title=title, outfile=outfile, make_global=True, project='PlateCarree', cmap=CMAP, add_gridlines=True)    

           for ifcst in range(nfcst):
               dstr=str(ifcst+1).zfill(2)
               outfile=outdir+USS+'_'+etypp+'_'+FEST[ifield]+'_'+dstr+'_'+datestr_spans+'.png'
               title = etypp+' error for '+USS+' '+FEST[ifield]+' day '+dstr+' '+datestr_spass
               plt_fld = FFLD[ifield][:, :, ivar, ifcst, ilev]-CFLD[ifield][:, :, ivar, ifcst, ilev]
               #if ( etype == 'square' ): plt_fld = np.sqrt(plt_fld)
               cplot.pcolormesh(lon_grid, lat_grid, plt_fld, levels=CLEV, obar='horizontal', title=title, outfile=outfile, make_global=True, project='PlateCarree', cmap=CMAP, add_gridlines=True)    

    return           

def make_spatial_plots(indir, insuffix, date_range, outdir='EPLOTS/', 
                       calculate_error=True, magnitude=False, plot_daily=True, plot_time_mean=True, 
                       levels=np.arange(-1.0, 1.1, 0.1), ilev=0, cmap=cmap_anom_fmask):
    date_range=[ check_date(idate,outtype=datetime.datetime) for idate in date_range]
    date_loop = date_range[0]
    isangle=-1
    if ( magnitude ): isangle=1
    error='error'
    if ( not calculate_error ): error='value'
    err_list = ( 'avg', 'abs', 'rms' )

    datestr=date_loop.strftime('%Y%m%d')  
    obsfile=indir+'/'+'class4_'+datestr+'_'+insuffix+'.nc'
    (LONO, LATO, depth), (obser, beste, inite, fcstv, persi, nersi) = read_obsfile_plus(obsfile)
    nobss, nvars, nfcst, nlevs = fcstv.shape
    
    grid_lon, grid_lat, lon_bin, lat_bin, grid_sum, grid_cnt = cplot.make_bin_grid(ddeg=4)
    ERR_LIST = []
    for err_item in err_list:
       ERR_LIST.append( [ grid_sum.copy(), grid_cnt.copy() ])
    ERR_VARS = [ ERR_LIST[:], ERR_LIST[:] ]
    ERR_BEST = [ ERR_LIST[:], ERR_LIST[:] ]
    ERR_INIT = [ ERR_LIST[:], ERR_LIST[:] ]
    ERR_FCST = []
    ERR_PERS = []
    for ifcst in range(nfcst):
        ERR_FCST.append([ ERR_LIST[:], ERR_LIST[:] ]) 
        ERR_PERS.append([ ERR_LIST[:], ERR_LIST[:] ])

    totobs=0
    while ( date_loop <= date_range[1] ):
        datestr=date_loop.strftime('%Y%m%d')  
        obsfile=indir+'/'+'class4_'+datestr+'_'+insuffix+'.nc'
        if ( os.path.isfile(obsfile) ):
            (LONO, LATO, depth), (obser, beste, inite, fcstv, persi, nersi) = read_obsfile_plus(obsfile)
            nobss, __, __, __ = fcstv.shape
            print(datestr, nobss)
            totobs=totobs+nobss
            if ( magnitude ):  (obser, beste, inite, fcstv, persi, nersi) = speed_and_angle_list((obser, beste, inite, fcstv, persi, nersi))
            isangle=-1
            if ( calculate_error ):
                (beste, inite, fcstv, persi, nersi) = calc_error(obser, (beste, inite, fcstv, persi, nersi), isangle=isangle)
        
            if ( plot_daily ):
                for ivar in range(2):
                    if ( ivar == 0 ): var='u'
                    if ( ivar == 1 ): var='v'
                    if ( ( ivar == 0 ) and ( magnitude ) ): var='s'
                    
                    plevels=levels
                    if ( ( magnitude ) and ( ivar == 2 ) ): plevels = np.pi * np.arange(-1.0,1.1, 0.1)
                    FLD=beste[:,ivar,ilev]
                    ofile=outdir+datestr+'_'+'best'+'_'+var+error+'.png'
                    title='Init Estimate for '+datestr
                    cplot.scatter(LONO, LATO, FLD, levels=plevels, obar='horizontal', outfile=ofile, title=title, make_global=True, cmap=cmap, add_gridlines=True)
                    FLD=inite[:,ivar,ilev]
                    ofile=outdir+datestr+'_'+'init'+'_'+var+error+'.png'
                    cplot.scatter(LONO, LATO, FLD, levels=plevels, obar='horizontal', outfile=ofile, title=title, make_global=True, cmap=cmap, add_gridlines=True)
                    if ( not calculate_error ):
                        FLD=obser[:,ivar,ilev] = obser[:,ivar,ilev]
                        ofile=outdir+datestr+'_'+'obsv'+'_'+var+error+'.png'
                        title='Observations for '+datestr
                        cplot.scatter(LONO, LATO, FLD, levels=plevels, obar='horizontal', outfile=ofile, title=title, make_global=True, cmap=cmap, add_gridlines=True)
                    
                    for ifcst in range(nfcst):
                        FLD=fcstv[:,ivar,ifcst, ilev]
                        ofile=outdir+datestr+'_'+'fc'+str(ifcst).zfill(2)+'_'+var+error+'.png'
                        title=str(ifcst).zfill(2)+' Forecast for '+datestr
                        cplot.scatter(LONO, LATO, FLD, levels=plevels, obar='horizontal', outfile=ofile, title=title, make_global=True, cmap=cmap, add_gridlines=True)
                        FLD=persi[:,ivar,ifcst, ilev]
                        ofile=outdir+datestr+'_'+'ps'+str(ifcst).zfill(2)+'_'+var+error+'.png'
                        title=str(ifcst).zfill(2)+' Persist for '+datestr
                        cplot.scatter(LONO, LATO, FLD, levels=plevels, obar='horizontal', outfile=ofile, title=title, make_global=True, cmap=cmap, add_gridlines=True)

            if ( plot_time_mean ):
                for ivar in range(2):
                    if ( ivar == 0 ): var='u'
                    if ( ivar == 1 ): var='v'
                    if ( ( ivar == 0 ) and ( magnitude ) ): var='s'
                    if ( ( ivar == 1 ) and ( magnitude ) ): var='a'
                    iERR_BEST = ERR_BEST[ivar]
                    iERR_INIT = ERR_INIT[ivar]
                    for ierr in range(len(err_list)):
                        grid_sum, grid_cnt = iERR_BEST[ierr]
                        FLD=beste[:, ivar, ilev]
                        if ( err_list[ierr] == 'abs' ):  FLD=np.absolute(FLD)
                        if ( err_list[ierr] == 'rms' ):  FLD=np.square(FLD)
                        grid_sum_new, grid_cnt_new = cplot.binfldsumcum(LONO, LATO, FLD, lon_bin, lat_bin, grid_sum, grid_cnt)
                        print(np.sum(grid_cnt_new), np.sum(grid_cnt), len(FLD))
                        iERR_BEST[ierr] = [grid_sum_new, grid_cnt_new]
                        grid_sum, grid_cnt = iERR_INIT[ierr]
                        FLD=inite[:, ivar, ilev]
                        grid_sum_new, grid_cnt_new = cplot.binfldsumcum(LONO, LATO, FLD, lon_bin, lat_bin, grid_sum, grid_cnt)
                        iERR_INIT[ierr] = [grid_sum_new, grid_cnt_new]
                    ERR_BEST[ivar] = iERR_BEST
                    ERR_INIT[ivar] = iERR_INIT
                for ifcst in range(nfcst):
                    for ivar in range(2):
                        if ( ivar == 0 ): var='u'
                        if ( ivar == 1 ): var='v'
                        if ( ( ivar == 0 ) and ( magnitude ) ): var='s'
                        if ( ( ivar == 1 ) and ( magnitude ) ): var='a'
                        iERR_FCST = ERR_FCST[ifcst][ivar]
                        iERR_PERS = ERR_PERS[ifcst][ivar]
                        for ierr in range(len(err_list)):
                            grid_sum, grid_cnt = iERR_FCST[ierr]
                            FLD=fcstv[:, ivar, ifcst, ilev]
                            if ( err_list[ierr] == 'abs' ):  FLD=np.absolute(FLD)
                            if ( err_list[ierr] == 'rms' ):  FLD=np.square(FLD)
                            grid_sum_new, grid_cnt_new = cplot.binfldsumcum(LONO, LATO, FLD, lon_bin, lat_bin, grid_sum, grid_cnt)
                            iERR_FCST[ierr] = grid_sum_new, grid_cnt_new
                            grid_sum, grid_cnt = iERR_PERS[ierr]
                            FLD=persi[:, ivar, ifcst, ilev]
                            if ( err_list[ierr] == 'abs' ):  FLD=np.absolute(FLD)
                            if ( err_list[ierr] == 'rms' ):  FLD=np.square(FLD)
                            grid_sum_new, grid_cnt_new = cplot.binfldsumcum(LONO, LATO, FLD, lon_bin, lat_bin, grid_sum, grid_cnt)
                            iERR_PERS[ierr] = grid_sum_new, grid_cnt_new
                        ERR_FCST[ifcst][ivar] = iERR_FCST
                        ERR_PERS[ifcst][ivar] = iERR_PERS

        date_loop = date_loop + datetime.timedelta(days=1)
        
    if ( plot_time_mean ):          
        datestr=date_range[0].strftime('%Y%m%d')+'_'+date_range[1].strftime('%Y%m%d')
        for ivar in range(2):
            if ( ivar == 0 ): var='u'
            if ( ivar == 1 ): var='v'
            if ( ( ivar == 0 ) and ( magnitude ) ): var='s'
            if ( ( ivar == 1 ) and ( magnitude ) ): var='a'
            iERR_BEST = ERR_BEST[ivar]
            iERR_INIT = ERR_INIT[ivar]
            for ierr in range(len(err_list)):
                if ( err_list[ierr] == 'avg' ):  error='merror'
                if ( err_list[ierr] == 'abs' ):  error='abserr'
                if ( err_list[ierr] == 'rms' ):  error='rmserr'
                grid_sum, grid_cnt = iERR_BEST[ierr]
                if ( np.sum(grid_cnt) != totobs ):
                    print('ERROR -- BAD OBS #:', totobs, np.sum(grid_cnt))
                else:
                    print('TOTAL OBSERVATIONS #:', totobs, np.sum(grid_cnt))
                grid_plt = cplot.binfldsumFIN(grid_sum, grid_cnt)
                #print(err_list[ierr], np.min(grid_plt))
                if ( err_list[ierr] == 'rms' ): grid_plt = np.sqrt(grid_plt)
                ofile=outdir+datestr+'_'+'best'+'_'+var+error+'.png'
                cplot.pcolormesh(grid_lon, grid_lat, grid_plt, outfile=ofile, levels=levels, cmap=cmap, add_gridlines=True)
                grid_sum, grid_cnt = iERR_INIT[ierr]
                grid_plt = cplot.binfldsumFIN(grid_sum, grid_cnt)
                if ( err_list[ierr] == 'rms' ): grid_plt = np.sqrt(grid_plt)
                ofile=outdir+datestr+'_'+'init'+'_'+var+error+'.png'
                cplot.pcolormesh(grid_lon, grid_lat, grid_plt, outfile=ofile, levels=levels, cmap=cmap, add_gridlines=True)
        for ifcst in range(nfcst):
            for ivar in range(2): 
                if ( ivar == 0 ): var='u'
                if ( ivar == 1 ): var='v'
                if ( ( ivar == 0 ) and ( magnitude ) ): var='s'
                if ( ( ivar == 1 ) and ( magnitude ) ): var='a'
                iERR_FCST = ERR_FCST[ifcst][ivar]
                iERR_PERS = ERR_PERS[ifcst][ivar]
                for ierr in range(len(err_list)):
                    if ( err_list[ierr] == 'avg' ):  error='merror'
                    if ( err_list[ierr] == 'abs' ):  error='abserr'
                    if ( err_list[ierr] == 'rms' ):  error='rmserr'
                    grid_sum, grid_cnt = iERR_FCST[ierr]
                    #if ( np.sum(grid_cnt) != totobs ):
                    #    print('ERROR -- BAD OBS #:', totobs, np.sum(grid_cnt))
                    #else:
                    #    print('TOTAL OBSERVATIONS #:', totobs, np.sum(grid_cnt))
                    grid_plt = cplot.binfldsumFIN(grid_sum, grid_cnt)
                    #print(ierr, err_list[ierr], np.min(grid_plt))
                    if ( err_list[ierr] == 'rms' ): grid_plt = np.sqrt(grid_plt)
                    ofile=outdir+datestr+'_'+'fc'+str(ifcst).zfill(2)+'_'+var+error+'.png'
                    cplot.pcolormesh(grid_lon, grid_lat, grid_plt, outfile=ofile, levels=levels, cmap=cmap, add_gridlines=True)
                    grid_sum, grid_cnt = iERR_PERS[ierr]
                    grid_plt = cplot.binfldsumFIN(grid_sum, grid_cnt)
                    if ( err_list[ierr] == 'rms' ): grid_plt = np.sqrt(grid_plt)
                    ofile=outdir+datestr+'_'+'ps'+str(ifcst).zfill(2)+'_'+var+error+'.png'
                    cplot.pcolormesh(grid_lon, grid_lat, grid_plt, outfile=ofile, levels=levels, cmap=cmap, add_gridlines=True)
    return

def find_date_in_filename(filename):
    idate=filename.find('class4_')+7
    datestr=filename[idate:idate+8]
    date = check_date(datestr, outtype=datetime.datetime)
    return date
    
def find_dates_in_expt(indir, insuffix='GIOPS_orca025_currents.f2'):
    obsfiles = sorted(glob.glob(indir+'/class4_'+'????????'+'_'+insuffix+'.nc'))
    dates = []
    for filename in obsfiles:
        date=find_date_in_filename(filename)
        dates.append(date)
    return dates
    
def find_common_dates(datelist1, datelist2):
    datelistc=[]
    for date in datelist1:
        if ( date in datelist2 ): datelistc.append(date)
    return datelistc
    
def calc_running_mean(field, nwin=31):
    istime=False
    if ( isinstance(field[0], datetime.datetime) ):
       field=time_to_float(field)
       istime=True
    nt=len(field)
    if ( nwin % 2 == 0 ):
        print('window should be odd length to centre')
        nwin=nwin+1
    noff = int( ( nwin - 1 ) / 2 )
    nrunm = np.convolve(field, np.ones(nwin)/nwin, mode='valid')
    new_field = np.nan*np.ones(nt)
    new_field[noff:-1*noff] = nrunm
    if ( istime ): new_field = time_to_float(new_field, inverse=True)
    return(new_field)

def time_to_float(times, inverse=False):
   if ( not inverse ):
       new_times = [ ttime.timestamp() for ttime in times ] 
   else:
      new_times = []
      for ttime in times:
         if ( np.isfinite(ttime) ):
             new_times.append(datetime.datetime.fromtimestamp(ttime))
         else:
             new_times.append(None)
   return(new_times)    
    
def interp_direct_obs( lonlat, field, lonlat_obs, method='2sweep', fill_value=missing):
    if ( isinstance(field, list) ):
        fld_at_obs = []
        for a_field in field:
            a_fld_at_obs = interp_direct_obs( lonlat, a_field, lonlat_obs, method=method, fill_value=fill_value)
            fld_at_obs.append(a_fld_at_obs)
        return fld_at_obs

    lon, lat = lonlat
    lono, lato = lonlat_obs
    this_method = method
    if ( this_method == '2sweeplinear' ): this_method='linear'
    if ( this_method == '2sweepcubic' ): this_method='cubic'
    fld_at_obs = si.griddata( (lon.flatten(), lat.flatten()), field.flatten(), (lono, lato), method=this_method, fill_value=fill_value)
    if ( method[:6] == '2sweep' ):
        fld_nr_obs = si.griddata( (lon.flatten(), lat.flatten()), field.flatten(), (lono, lato), method='nearest', fill_value=fill_value)
        IMISS = np.where( fld_at_obs == fill_value )
        fld_at_obs[IMISS] = fld_nr_obs[IMISS]
    return fld_at_obs
