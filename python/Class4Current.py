import os
import sys
sys.path.insert(0, '/home/dpe000/Class4_Currents/python')

import netCDF4
import shutil 
import numpy as np
import datetime
import time
import subprocess

import find_fcst_file
import stfd
import read_grid
import isoheatcontent
import bering_diurnal
import cplot
import find_hall
import datadatefile
import glob

KCONV=273.16

cmap_full_field='gist_stern_r'
cmap_anom_field='seismic'

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

file_best = '/fs/site4/eccc/mrd/rpnenv/dpe000//tmpdir/gu/ORCA025-CMC-ANAL_1d_grid_T_2021110100.nc'
file_besu = '/fs/site4/eccc/mrd/rpnenv/dpe000//tmpdir/gu/ORCA025-CMC-ANAL_1d_grid_U_2021110100.nc'
file_besv = '/fs/site4/eccc/mrd/rpnenv/dpe000//tmpdir/gu/ORCA025-CMC-ANAL_1d_grid_V_2021110100.nc'
file_fcts = ' tmpdir value'

missing=-999

def UU_ORCA_to_NE(UV):
    hall=find_hall.find_hall()
    fila='/space/'+hall+'/sitestore/eccc/mrd/rpnenv/socn000/env_ubuntu-18.04-skylake-64/datafiles/constants/oce/repository/master/CONCEPTS/orca025/grids/orca025grid_new.std'
    LONA, LATA, TH = stfd.read_fstd_var(fila, 'LAAN')
    UO, VO = UV
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

def speed_and_angle(u, v):
    speed =  np.sqrt(u**2, v**2)
    angle = np.arctan2(v, u)
    return speed, angle
       
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

    ipt, jpt = bering_diurnal.find_nearest_point(LONP, LATP, LONN, LATN)
    
    FDpt = []
    FDpi = []
    for ifld, FDN in enumerate(FD_LIST):
        #ivalid = np.where(FDN.mask == False)
	#print('PRINT', ifld, FDN.shape, LONN.shape, LATN.shape, LONPT, LATPT)
        FDpt.append(FDN[ipt, jpt])
        FDpi.append(bering_diurnal.interpolate_to_point(FDN, LONN, LATN, LONPT, LATPT))
    print('POINT VALUES', FDpt, FDpi)
    return FDpt, FDpi

def find_nearest_points((LONPT, LATPT), LONLAT_LIST):
    IJPTS=[]
    for LONLATin in LONLAT_LIST:
        LONin, LATin = LONLATin
	ipto, jpto = bering_diurnal.find_nearest_point(LONPT, LATPT, LONin, LATin)
	IJPTS.append((ipto,jpto))
    return IJPTS

def msk_flds(FLDSin, (LONin, LATin), maskin=mask0):
    FLDSout=[]
    for FLD in FLDSin:
        LONM, LATM, FLDm = cplot.mask_field(LONin, LATin, FLD, maskin)
	FLDSout.append(FLDm)
    return FLDSout, (LONM, LATM)
        
def put_flds_latlon(FLDSin, (LONin, LATin), ddeg=0.1, method='2sweeplinear'):
    FLDSout=[]
    for FLD in FLDSin:
        LONf, LATf, FLDf = cplot.grdfld(LONin, LATin, FLD, ddeg=ddeg, method=method)
	FLDSout.append(FLDf)
    return FLDSout, (LONf, LATf)
	
def read_anal_file(file, grid='T'):
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
    FLD = np.transpose(np.squeeze(FLD), [0, 2, 1])   # PROBLEM FOR 2D FIELD # NEED TO CONVERT TO KELVIN
    if ( fvar == 'thetao' ): # Convert from Celsius to Kelvin
        FLD = FLD + KCONV 
    FLD = np.ma.masked_array(FLD, mask=1-maskf)
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

grid_dir='/fs/site4/eccc/mrd/rpnenv/dpe000/GDPS/gdps_Class4/'
script_dir='/fs/homeu1/eccc/mrd/ords/rpnenv/dpe000/Class4_Currents/jobscripts/'
template_dir='/fs/homeu1/eccc/mrd/ords/rpnenv/dpe000/Class4_Currents/templates/'
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
    rc=subprocess.call(['bash' , script, '-f='+filter, '-i='+source_file, '-o='+destination_file])
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

def plot_fields(FLDS, (LONFLD, LATFLD), suptitle=None, grid=False, outfile_prefix='PLOT/'):
    [TFLD, UFLD, VFLD] = FLDS
    for ifld, FLD in enumerate(FLDS):
        CLEV = np.arange(-5, 5.5, 0.5)
        if ( ifld == 0 ): CLEV = np.arange(-2, 34, 2)
	if ( ifld == 0 ): FLD = FLD - KCONV
	if ( ifld == 0 ): title='T'
	if ( ifld == 1 ): title='U'
	if ( ifld == 2 ): title='V'
	cmap_use=cmap_anom_field
	if ( ifld == 0 ): cmap_use=cmap_full_field

	outfile=outfile_prefix+title+'.png'
        if ( not grid ):
	    cplot.pcolormesh(LONFLD, LATFLD, FLD, levels=CLEV, cmap=cmap_use, project='PlateCarree', 
	       outfile=outfile, make_global=True, title=title, suptitle=suptitle, 
	       cbar=True, obar='vertical', fontsizes=None)
        else:
	    cplot.grd_pcolormesh(LONFLD, LATFLD, FLD, levels=CLEV, cmap=cmap_use, project='PlateCarree', 
	       outfile=outfile, ddeg=0.2, make_global=True, title=title, suptitle=suptitle, 
	       cbar=True, obar='vertical', fontsizes=None)
	outfile=outfile_prefix+title+'.NP.png'
	cplot.pcolormesh(LONFLD, LATFLD, FLD, levels=CLEV, cmap=cmap_use, project='NorthPolarStereo', 
	       outfile=outfile, make_global=True, title=title, suptitle=suptitle, 
	       cbar=True, obar='vertical', fontsizes=None, box=[-180, 180, 65, 90])
	      
    return

def calc_error(obser, model):

    if ( ( isinstance(model, list) ) or ( isinstance(model, tuple) ) ):
        ERROR = []
	for imodel in model:
	    iERROR = calc_error(obser, imodel)
	    ERROR.append(imodel)
	return ERROR
	
    nobss, nvars, ndeps = obser.shape
    mobss, mvars, mfcst, mdeps = (0, 0, 0, 0)
    if ( model.ndim == 3 ):
        mobss, mvars, mdeps = model.shape
    elif ( model.ndim == 4 ): 
        mobss, mvars, mfcst, mdeps = model.shape

    if ( (nobss != mobss) or (nvars != mvars) or (ndeps != mdeps) ):
        print('SHAPE MISMATCH')
	return None
	
    if ( mfcst == 0 ):
        ERROR = model - obser
    if ( mfcst > 0 ):
        ERROR = 0.0*model.copy()
        for ifcst in range(mfcst):
	    ERROR[:,:,ifcst,:] = model[:,:,ifcst,:] - obser[:,:,:]
	    
    return ERROR

def calc_mean_error(obser, model):
    if ( ( isinstance(model, list) ) or ( isinstance(model, tuple) ) ):
        mean_error = []
	for imodel in model:
	    iERROR = calc_mean_error(obser, imodel)
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

    error = calc_error(obser, model)
    mean_error = np.mean(error, axis=0)
    rmse_error = np.std (error, axis=0)
   
    return mean_error, rmse_error
	
tate=datetime.datetime(2021, 10, 31)
def process_obsfile(date=tate, TEST_SINGLE=False, Plot=False):
    datestr=date.strftime('%Y%m%d')   
    time0=time.time()
    psyfile='CLASS4_currents/class4_'+datestr+'_PSY4V3R1_orca12_currents.nc'
    obsfile1='CLASS4_currents_CCMEP/class4_'+datestr+'_GIOPS_orca025_currents.1.nc'
    obsfile2='CLASS4_currents_CCMEP/class4_'+datestr+'_GIOPS_orca025_currents.2.nc'
    obsfile3='CLASS4_currents_CCMEP/class4_'+datestr+'_GIOPS_orca025_currents.3.nc'
    print('obsfile', obsfile1, obsfile2, obsfile3)
    # READ IN OBS FILE
    (LONO, LATO, depth), (obser, beste, persi, fcstv) = read_obsfile(psyfile)
    obsset = netCDF4.Dataset(psyfile,mode='r')

    neste = beste.copy()    
    Mobser = 0.0*obser.copy()
    Mneste = 0.0*neste.copy()
    Mbeste = 0.0*beste.copy()
    Mfcstv = 0.0*fcstv.copy()
    Mpersi = 0.0*persi.copy()
    
    Mobser[:,0,:], Mobser[:,1, :] = speed_and_angle(obser[:,0,:], obser[:,1,:])
    Mneste[:,0,:], Mneste[:,1, :] = speed_and_angle(neste[:,0,:], neste[:,1,:])
    Mbeste[:,0,:], Mbeste[:,1, :] = speed_and_angle(beste[:,0,:], beste[:,1,:])
    Mfcstv[:,0,:,:], Mfcstv[:,1, :] = speed_and_angle(fcstv[:,0,:,:], fcstv[:,1,:,:])
    Mpersi[:,0,:,:], Mpersi[:,1, :] = speed_and_angle(persi[:,0,:,:], persi[:,1,:,:])

    nobss, nvars, nfcst, ndeps = fcstv.shape
    if ( TEST_SINGLE ):
        nobss_loop = 1
    else:
        nobss_loop = nobss
	
    ## COPY OF OBS FILE VARIABLES FOR EACH OUTPUT
    tmp_neste = 0.0*beste.copy()
    tmp_beste = 0.0*beste.copy()
    tmp_fcstv = 0.0*fcstv.copy()
    tmp_persi = 0.0*persi.copy()
    tmp_neste.mask = False 
    tmp_beste.mask = False
    tmp_fcstv.mask = False
    tmp_persi.mask = False
    NESTE_LIST = [tmp_neste.copy(), tmp_neste.copy(), tmp_neste.copy()]
    BESTE_LIST = [tmp_beste.copy(), tmp_beste.copy(), tmp_beste.copy()]
    FCSTV_LIST = [tmp_fcstv.copy(), tmp_fcstv.copy(), tmp_fcstv.copy()]
    PERSI_LIST = [tmp_persi.copy(), tmp_persi.copy(), tmp_persi.copy()]
    
    tempn=np.zeros((nobss, ndeps))
    tempb=np.zeros((nobss, ndeps))
    tempf=np.zeros((nobss, nfcst, ndeps))
    tempp=np.zeros((nobss, nfcst, ndeps))

    TEMPN_LIST = [tempn.copy(), tempn.copy(), tempn.copy()]
    TEMPB_LIST = [tempb.copy(), tempb.copy(), tempb.copy()]
    TEMPF_LIST = [tempf.copy(), tempf.copy(), tempf.copy()]
    TEMPP_LIST = [tempp.copy(), tempp.copy(), tempp.copy()]

    MERCATOR_MEAN_ERRORS =  calc_mean_error(obser, (neste, beste, fcstv, persi))   
    mr_ERRMn, mr_ERRMb, mr_ERRMf, mr_ERRMp = MERCATOR_MEAN_ERRORS
    mean_ERRMn, rmse_ERRMn = mr_ERRMn
    mean_ERRMb, rmse_ERRMb = mr_ERRMb
    mean_ERRMf, rmse_ERRMf = mr_ERRMf 
    mean_ERRMp, rmse_ERRMp = mr_ERRMp
    MERCATOR_MMEAN_ERRORS = calc_mean_error(Mobser, (Mneste, Mbeste, Mfcstv, Mpersi))
    mr_ERMMn, mr_ERMMb, mr_ERMMf, mr_ERMMp = MERCATOR_MMEAN_ERRORS
    mean_ERMMn, rmse_ERMMn = mr_ERMMn
    mean_ERMMb, rmse_ERMMb = mr_ERMMb
    mean_ERMMf, rmse_ERMMf = mr_ERMMf 
    mean_ERMMp, rmse_ERMMp = mr_ERMMp

    if ( not TEST_SINGLE ):  ## DON"T DO FOR TEST  
        
        write_mean_errors('ERRORS/PSY4a', int(datestr), MERCATOR_MEAN_ERRORS, vec=['u','v'], tmp_prefix='site4/TMP/PSY4a')
        write_mean_errors('ERRORS/PSY4a', int(datestr), MERCATOR_MMEAN_ERRORS, vec=['s','a'], tmp_prefix='site4/TMP/PSY4a')
	#nobss, nvars, nfcst, ndeps = fcstv.shape
        ## Keep 15m 
        MERCA_MEANU_ERRORS=np.concatenate((np.array([mean_ERRMn[0,1], mean_ERRMb[0,1]]), mean_ERRMp[0,:,1].flatten(), mean_ERRMf[0,:,1].flatten()))
        MERCA_MEANV_ERRORS=np.concatenate((np.array([mean_ERRMn[1,1], mean_ERRMb[1,1]]), mean_ERRMp[1,:,1].flatten(), mean_ERRMf[1,:,1].flatten()))
        MERCA_MEANS_ERRORS=np.concatenate((np.array([mean_ERMMn[0,1], mean_ERMMb[0,1]]), mean_ERMMp[0,:,1].flatten(), mean_ERMMf[0,:,1].flatten()))
        MERCA_MEANA_ERRORS=np.concatenate((np.array([mean_ERMMn[1,1], mean_ERMMb[1,1]]), mean_ERMMp[1,:,1].flatten(), mean_ERMMf[1,:,1].flatten()))
        MERCA_RMSEU_ERRORS=np.concatenate((np.array([rmse_ERRMn[0,1], rmse_ERRMb[0,1]]), rmse_ERRMp[0,:,1].flatten(), rmse_ERRMf[0,:,1].flatten()))
        MERCA_RMSEV_ERRORS=np.concatenate((np.array([rmse_ERRMn[1,1], rmse_ERRMb[1,1]]), rmse_ERRMp[1,:,1].flatten(), rmse_ERRMf[1,:,1].flatten()))
        MERCA_RMSES_ERRORS=np.concatenate((np.array([rmse_ERMMn[0,1], rmse_ERMMb[0,1]]), rmse_ERMMp[0,:,1].flatten(), rmse_ERMMf[0,:,1].flatten()))
        MERCA_RMSEA_ERRORS=np.concatenate((np.array([rmse_ERMMn[1,1], rmse_ERMMb[1,1]]), rmse_ERMMp[1,:,1].flatten(), rmse_ERMMf[1,:,1].flatten()))
        datadatefile.add_to_file(int(datestr), MERCA_MEANU_ERRORS, file='ERRORS/PSY4_meanu.dat', tmpfile='site4/TMP/PSY4_meanu.dat')
        datadatefile.add_to_file(int(datestr), MERCA_MEANV_ERRORS, file='ERRORS/PSY4_meanv.dat', tmpfile='site4/TMP/PSY4_meanv.dat')
        datadatefile.add_to_file(int(datestr), MERCA_MEANS_ERRORS, file='ERRORS/PSY4_means.dat', tmpfile='site4/TMP/PSY4_means.dat')
        datadatefile.add_to_file(int(datestr), MERCA_MEANA_ERRORS, file='ERRORS/PSY4_meana.dat', tmpfile='site4/TMP/PSY4_meana.dat')
        datadatefile.add_to_file(int(datestr), MERCA_RMSEU_ERRORS, file='ERRORS/PSY4_rmseu.dat', tmpfile='site4/TMP/PSY4_rmseu.dat')
        datadatefile.add_to_file(int(datestr), MERCA_RMSEV_ERRORS, file='ERRORS/PSY4_rmsev.dat', tmpfile='site4/TMP/PSY4_rmsev.dat')
        datadatefile.add_to_file(int(datestr), MERCA_RMSES_ERRORS, file='ERRORS/PSY4_rmses.dat', tmpfile='site4/TMP/PSY4_rmses.dat')
        datadatefile.add_to_file(int(datestr), MERCA_RMSEA_ERRORS, file='ERRORS/PSY4_rmsea.dat', tmpfile='site4/TMP/PSY4_rmsea.dat')

    timen=time.time()
    timep=timen-time0
    print("TIMING :: Total Observation Preparation Time: "+str(timep))

    time0=timen
    # READ IN MODEL standard FILE
    (LONN, LATN) = (nav_lon, nav_lat)
    bedate=date
    nedate=date + datetime.timedelta(days=1)

    file_nest, __ =find_fcst_file.find_anal_file(nedate, system='gu', var='T',execute=True)
    file_nesu, __ =find_fcst_file.find_anal_file(nedate, system='gu', var='U',execute=False)
    file_nesv, __ =find_fcst_file.find_anal_file(nedate, system='gu', var='V',execute=False)
    print('file_nest', file_nest)
    file_best, __ =find_fcst_file.find_anal_file(bedate, system='gu', var='T',execute=True)
    file_besu, __ =find_fcst_file.find_anal_file(bedate, system='gu', var='U',execute=False)
    file_besv, __ =find_fcst_file.find_anal_file(bedate, system='gu', var='V',execute=False)
    print('file_best', file_best)

    timen=time.time()
    timep=timen-time0
    print("TIMING :: Best Estimate Retieval Time: "+str(timep))
    time0=timen

    ## NEXT DAY ESTIMATE    
    TM, LONT, LATT, LEVT, DAT = read_anal_file(file_nest, grid='T')
    UW, LONU, LATU, LEVU, DAT = read_anal_file(file_nesu, grid='U')
    VW, LONV, LATV, LEVV, DAT = read_anal_file(file_nesv, grid='V')
    UU, LONUU, LATUU = map_to_A_grid(UW, LONU, LATU, grid='U')
    VV, LONVV, LATVV = map_to_A_grid(VW, LONV, LATV, grid='V')
    # Important for U15 calculations that arrays not be masked.
    print(type(TM), type(UU), type(VV))
    if ( isinstance(TM, np.ma.core.MaskedArray) ): TM = TM.data
    if ( isinstance(UU, np.ma.core.MaskedArray) ): UU = UU.data
    if ( isinstance(VV, np.ma.core.MaskedArray) ): VV = VV.data
    print(type(TM), type(UU), type(VV))
    [T00, U00, V00] = [TM[0,:,:], UU[0,:,:], VV[0,:,:]]
    [T15, U15, V15] = calc_m15([TM, UU, VV], e3t, mask)
    if ( Plot ):
        suptitle='0m'+datestr+'Native Grid Nest Estimate'
        prefix='PLOTS/00_NENG_'+datestr+'_'
        plot_fields((T00, U00, V00), (LONN, LATN), suptitle=suptitle, grid=True, outfile_prefix=prefix)
        suptitle='15m'+datestr+'Native Grid Nest Estimate'
        prefix='PLOTS/15_NENG_'+datestr+'_'
        plot_fields((T15, U15, V15), (LONN, LATN), suptitle=suptitle, grid=True, outfile_prefix=prefix)
    [U00, V00] = UU_ORCA_to_NE([U00, V00])
    [U15, V15] = UU_ORCA_to_NE([U15, V15])

    # this actually removes non-ocean points
    (T00m, U00m, V00m, T15m, U15m, V15m), (LONM, LATM) = msk_flds([T00, U00, V00, T15, U15, V15], [LONN, LATN], mask0)
    # now grid to lat long
    (T00f, U00f, V00f, T15f, U15f, V15f), (LONF, LATF) = put_flds_latlon([T00m, U00m, V00m, T15m, U15m, V15m], [LONM, LATM], ddeg=0.2, method='2sweeplinear')
    if ( Plot == True ):
        suptitle='0m'+datestr+'Nest Estimate'
        prefix='PLOTS/00_NEST_'+datestr+'_'
        plot_fields((T00f, U00f, V00f), (LONF, LATF), suptitle=suptitle, grid=False, outfile_prefix=prefix)
        suptitle='15m'+datestr+'Nest Estimate'
        prefix='PLOTS/15_NEST_'+datestr+'_'
        plot_fields((T15f, U15f, V15f), (LONF, LATF), suptitle=suptitle, grid=False, outfile_prefix=prefix)
 
    NEST=[
          [[T00,  U00,  V00],  [T15,  U15,  V15]],
	  [[T00f, U00f, V00f], [T15f, U15f, V15f]]
	 ]
    
    ## BEST ESTIMATE    
    TM, LONT, LATT, LEVT, DAT = read_anal_file(file_best, grid='T')
    UW, LONU, LATU, LEVU, DAT = read_anal_file(file_besu, grid='U')
    VW, LONV, LATV, LEVV, DAT = read_anal_file(file_besv, grid='V')
    UU, LONUU, LATUU = map_to_A_grid(UW, LONU, LATU, grid='U')
    VV, LONVV, LATVV = map_to_A_grid(VW, LONV, LATV, grid='V')
    # Important for U15 calculations that arrays not be masked.
    print(type(TM), type(UU), type(VV))
    if ( isinstance(TM, np.ma.core.MaskedArray) ): TM = TM.data
    if ( isinstance(UU, np.ma.core.MaskedArray) ): UU = UU.data
    if ( isinstance(VV, np.ma.core.MaskedArray) ): VV = VV.data
    print(type(TM), type(UU), type(VV))
    [T00, U00, V00] = [TM[0,:,:], UU[0,:,:], VV[0,:,:]]
    [T15, U15, V15] = calc_m15([TM, UU, VV], e3t, mask)
    if ( Plot ):
        suptitle='0m'+datestr+'Native Grid Best Estimate'
        prefix='PLOTS/00_BENG_'+datestr+'_'
        plot_fields((T00, U00, V00), (LONN, LATN), suptitle=suptitle, grid=True, outfile_prefix=prefix)
        suptitle='15m'+datestr+'Native Grid Best Estimate'
        prefix='PLOTS/15_BENG_'+datestr+'_'
        plot_fields((T15, U15, V15), (LONN, LATN), suptitle=suptitle, grid=True, outfile_prefix=prefix)
    [U00, V00] = UU_ORCA_to_NE([U00, V00])
    [U15, V15] = UU_ORCA_to_NE([U15, V15])

    # this actually removes non-ocean points
    (T00m, U00m, V00m, T15m, U15m, V15m), (LONM, LATM) = msk_flds([T00, U00, V00, T15, U15, V15], [LONN, LATN], mask0)
    # now grid to lat long
    (T00f, U00f, V00f, T15f, U15f, V15f), (LONF, LATF) = put_flds_latlon([T00m, U00m, V00m, T15m, U15m, V15m], [LONM, LATM], ddeg=0.2, method='2sweeplinear')
    if ( Plot == True ):
        suptitle='0m'+datestr+'Best Estimate'
        prefix='PLOTS/00_BEST_'+datestr+'_'
        plot_fields((T00f, U00f, V00f), (LONF, LATF), suptitle=suptitle, grid=False, outfile_prefix=prefix)
        suptitle='15m'+datestr+'Best Estimate'
        prefix='PLOTS/15_BEST_'+datestr+'_'
        plot_fields((T15f, U15f, V15f), (LONF, LATF), suptitle=suptitle, grid=False, outfile_prefix=prefix)
 
    BEST=[
          [[T00,  U00,  V00],  [T15,  U15,  V15]],
	  [[T00f, U00f, V00f], [T15f, U15f, V15f]]
	 ]
    timen=time.time()
    timep=timen-time0
    print("TIMING :: Best Estimate Processing Time: "+str(timep))
    time0=timen
    timel=time0
    
     
    FCST=[]
    PERS=[]
    for ifcst in range(nfcst):
        jfcst=ifcst+1
        fchour=jfcst*24
        fcdate = date - datetime.timedelta(days=jfcst)
        file_best, __ = find_fcst_file.find_anal_file(fcdate, system='gu', var='T',execute=True)
        file_besu, __ = find_fcst_file.find_anal_file(fcdate, system='gu', var='U',execute=True)
        file_besv, __ = find_fcst_file.find_anal_file(fcdate, system='gu', var='V',execute=True)
        file_fcst     = find_fcst_file.find_fcst_file('OPD', fcdate, fchour, 0, src='ocn', execute=True)
	filter_standard_file(template_dir+'/flt.ksh', file_fcst, file_fcst)
	# Produce a cst interpolated grid
	## Can't afford memory for 0.1deg grids -- try 0.2deg grids (should be 4 times less).
	file_intr = file_fcst+'.Z02'
	if ( not os.path.isfile(file_intr) ):
	    rc = cst_interpolation(file_fcst, file_intr, ref_grid=grid_dir+'dest_grid.std.2')
        print('file_best', file_best, os.path.isfile(file_best))
        print('file_fcst', file_fcst, os.path.isfile(file_fcst))
        print('file_intr', file_intr, os.path.isfile(file_intr))
    
        TM, LONT, LATT, LEVT, DAT = read_anal_file(file_best, grid='T')
        UW, LONU, LATU, LEVU, DAT = read_anal_file(file_besu, grid='U')
        VW, LONV, LATV, LEVV, DAT = read_anal_file(file_besv, grid='V')
        UU, LONUU, LATUU = map_to_A_grid(UW, LONU, LATU, grid='U')
        VV, LONVV, LATVV = map_to_A_grid(VW, LONV, LATV, grid='V')
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

        PERS.append([
	             [[T00,  U00,  V00],[T15,  U15,  V15]],
	             [[T00f, U00f, V00f],[T15f, U15f, V15f]]
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

        [T00i, U00i, V00i] = [TI[0,:,:], UI[0,:,:], VI[0,:,:]]
        [T15i, U15i, V15i] = calc_m15([TI, UI, VI], e3ti, MI)
	print("Interpolated Shapes", T00i.shape, T15i.shape, LONG.shape, LATG.shape)
	
        if ( Plot == True ):
            suptitle='0m'+datestr+'Native Grid Lead '+str(ifcst+1)
            prefix='PLOTS/00_NG'+str(ifcst+1)+'_'+datestr+'_'
            plot_fields((T00f, U00f, V00f), (LONF, LATF), suptitle=suptitle, grid=False, outfile_prefix=prefix)
            suptitle='15m'+datestr+'Native Grid Lead '+str(ifcst+1)
            prefix='PLOTS/15_NG'+str(ifcst+1)+'_'+datestr+'_'
            plot_fields((T15f, U15f, V15f), (LONF, LATF), suptitle=suptitle, grid=False, outfile_prefix=prefix)

	## ROTATE NATIVE GRID VECTORS TO EAST and NORTH
        [U00, V00] = UU_ORCA_to_NE([U00, V00])
        [U15, V15] = UU_ORCA_to_NE([U15, V15])

        # this actually removes non-ocean points
        (T00m, U00m, V00m, T15m, U15m, V15m), (LONM, LATM) = msk_flds([T00, U00, V00, T15, U15, V15], [LONN, LATN], mask0)
        # now grid to lat long
        (T00f, U00f, V00f, T15f, U15f, V15f), (LONF, LATF) = put_flds_latlon([T00m, U00m, V00m, T15m, U15m, V15m], [LONM, LATM], ddeg=0.2, method='2sweeplinear')

        if ( Plot == True ):
            suptitle='0m'+datestr+'Forecast Lead '+str(ifcst+1)
            prefix='PLOTS/00_FC'+str(ifcst+1)+'_'+datestr+'_'
            plot_fields((T00f, U00f, V00f), (LONF, LATF), suptitle=suptitle, grid=False, outfile_prefix=prefix)
            suptitle='15m'+datestr+'Forecast Lead '+str(ifcst+1)
            prefix='PLOTS/15_FC'+str(ifcst+1)+'_'+datestr+'_'
            plot_fields((T15f, U15f, V15f), (LONF, LATF), suptitle=suptitle, grid=False, outfile_prefix=prefix)
            suptitle='0m'+datestr+'cstintr Forecast Lead '+str(ifcst+1)
            prefix='PLOTS/00_FI'+str(ifcst+1)+'_'+datestr+'_'
            plot_fields((T00i, U00i, V00i), (LONG, LATG), suptitle=suptitle, grid=False, outfile_prefix=prefix)
            suptitle='15m'+datestr+'cstintr Forecast Lead '+str(ifcst+1)
            prefix='PLOTS/15_FI'+str(ifcst+1)+'_'+datestr+'_'
            plot_fields((T15i, U15i, V15i), (LONG, LATG), suptitle=suptitle, grid=False, outfile_prefix=prefix)

        FCST.append([
	             [[T00,  U00,  V00],[T15,  U15,  V15]],
	             [[T00f, U00f, V00f],[T15f, U15f, V15f]],
		     [[T00i, U00i, V00i],[T15i, U15i, V15i]]
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
	IJPTS = find_nearest_points((LONP, LATP), [(LONN, LATN), (LONF, LATF), (LONG, LATG) ])
	IJPT, IJPF, IJPI = IJPTS
	print('IJPTS', IJPTS)

	for igrid, GRID in enumerate(NEST):
	    
	    neste = NESTE_LIST[igrid]
	    tempn = TEMPN_LIST[igrid]
	    IJPT = IJPTS[igrid]
	    
	    for kk, KLEV in enumerate(GRID):
	        nestl = []
	        for ifld, FLD in enumerate(KLEV):
		    if ( IJPT[1] != 1 ):
		        print(FLD.shape, IJPT)
		        FLP=FLD[IJPT]
	            else: 
			if ( FLD.ndim == 1 ):
		          FLP=FLD[IJPT[0]]
			else:
			  print ("Error:  Mismatch in points and field shapes")
			  print( "NEST", igrid, kk, ifld, FLD.shape, IJPT ) 
			  FLP=missing
		    nestl.append(FLP)
		  
		Tpi, Upi, Vpi = nestl   ## These are the best values for level kk
		## Repeatedly fill in the nest values.  Initially Orca grid, then orca grid masked, then lat lon grid.
	        neste[iobs, 0, kk] = Upi
	        neste[iobs, 1, kk] = Vpi
	        tempn[iobs, kk] = Tpi
		print (iobs, igrid, kk, 'Nest Values', Upi, Vpi, Tpi, obser[iobs, :, kk] )
            NESTE_LIST[igrid] = neste
	    TEMPN_LIST[igrid] = tempn
        # no 3rd IGRID.  COPY LAST
	NESTE_LIST[2] = neste
	TEMPN_LIST[2] = tempn
	
	for igrid, GRID in enumerate(BEST):
	    beste = BESTE_LIST[igrid]
	    tempb = TEMPB_LIST[igrid]
	    IJPT = IJPTS[igrid]
	    
	    for kk, KLEV in enumerate(GRID):
	        bestl = []
	        for ifld, FLD in enumerate(KLEV):
		    if ( IJPT[1] != 1 ):
		        print(FLD.shape, IJPT)
		        FLP=FLD[IJPT]
	            else: 
			if ( FLD.ndim == 1 ):
		          FLP=FLD[IJPT[0]]
			else:
			  print ("Error:  Mismatch in points and field shapes")
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
	BESTE_LIST[2] = beste
	TEMPB_LIST[2] = tempb 
	
        for ip, PRED in enumerate([FCST, PERS]): 
          for ld, LEAD in enumerate(PRED):
	    for ig, GRID in enumerate(LEAD): 
	      IJPT = IJPTS[ig]
	      if ( ip == 0 ):
	        fcstv = FCSTV_LIST[ig]
	        tempf = TEMPF_LIST[ig]
	      if ( ip == 1 ):
	        persi = PERSI_LIST[ig]
	        tempp = TEMPP_LIST[ig]
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
              if ( ip == 0 ):
	        FCSTV_LIST[ig] = fcstv
	        TEMPF_LIST[ig] = tempf
              if ( ip == 1 ):
	        PERSI_LIST[ig] = persi
	        TEMPP_LIST[ig] = tempp
            # no 3rd IGRID for persi.  COPY LAST
	    if ( ip == 1 ):
	      PERSI_LIST[2] = persi
	      TEMPP_LIST[2] = tempp

        print('IOBS', iobs, obser[iobs,:, :].data, beste[iobs,:, :].data)

    timen = time.time()
    timet = timen-time0
    print('TIMING :: Total Observation Processing Time: '+str(timet)+' seconds')

    for ig in range(3):
        neste = NESTE_LIST[ig]
	beste = BESTE_LIST[ig]
	fcstv = FCSTV_LIST[ig]
	persi = PERSI_LIST[ig]
	tempn = TEMPN_LIST[ig]
	tempb = TEMPB_LIST[ig]
	tempf = TEMPF_LIST[ig]
	tempp = TEMPP_LIST[ig]
		
        Mobser = 0.0*obser.copy()
        Mneste = 0.0*neste.copy()
        Mbeste = 0.0*beste.copy()
        Mfcstv = 0.0*fcstv.copy()
        Mpersi = 0.0*persi.copy()

        Mobser[:,0,:], Mobser[:,1, :] = speed_and_angle(obser[:,0,:], obser[:,1,:])
        Mneste[:,0, ], Mneste[:,1, :] = speed_and_angle(neste[:,0,:], neste[:,1,:])
        Mbeste[:,0, ], Mbeste[:,1, :] = speed_and_angle(beste[:,0,:], beste[:,1,:])
        Mfcstv[:,0,:,:], Mfcstv[:,1, :] = speed_and_angle(fcstv[:,0,:,:], fcstv[:,1,:,:])
        Mpersi[:,0,:,:], Mpersi[:,1, :] = speed_and_angle(persi[:,0,:,:], persi[:,1,:,:])

        CCMEP_MEAN_ERRORS = calc_mean_error(obser, (neste, beste, fcstv, persi))
        mr_ERRCn, mr_ERRCb, mr_ERRCf, mr_ERRCp = CCMEP_MEAN_ERRORS
        mean_ERRCn, rmse_ERRCn = mr_ERRCn
	mean_ERRCb, rmse_ERRCb = mr_ERRCb
	mean_ERRCf, rmse_ERRCf = mr_ERRCf
	mean_ERRCp, rmse_ERRCp = mr_ERRCp
	CCMEP_MMEAN_ERRORS = calc_mean_error(Mobser, (Mneste, Mbeste, Mfcstv, Mpersi))
        mr_ERMCn, mr_ERMCb, mr_ERMCf, mr_ERMCp = CCMEP_MMEAN_ERRORS
        mean_ERMCn, rmse_ERMCn = mr_ERMCn
	mean_ERMCb, rmse_ERMCb = mr_ERMCb
	mean_ERMCf, rmse_ERMCf = mr_ERMCf
	mean_ERMCp, rmse_ERMCp = mr_ERMCp

        print(datestr+' RESULTS')
    
        print('Best Error Mercator', mean_ERRMb, rmse_ERRMb, mean_ERMMb, rmse_ERMMb)
        print(ig, 'Best Error CCMEP', mean_ERRCb, rmse_ERRCb, mean_ERMCb, rmse_ERMCb)
        print(ig, 'Nest Error CCMEP', mean_ERRCn, rmse_ERRCn, mean_ERMCn, rmse_ERMCn)
    
        print('Persistence Error Mercator', mean_ERRMp, rmse_ERRMp, mean_ERMMp, rmse_ERMMp)
        print(ig, 'Persistence Error CCMEP', mean_ERRCp, rmse_ERRCp, mean_ERMCp, rmse_ERMCp)
    
        print('Forecast Error Mercator', mean_ERRMf, rmse_ERRMf, mean_ERMMf, rmse_ERMMf)
        print(ig, 'Forecast Error CCMEP', mean_ERRCf, rmse_ERRCf, mean_ERMCf, rmse_ERMCf)

        if ( not TEST_SINGLE ):  ## DON"T DO FOR TEST   
	    if ( ig == 0 ) : obsfile_ig=obsfile1
	    if ( ig == 1 ) : obsfile_ig=obsfile2
	    if ( ig == 2 ) : obsfile_ig=obsfile3
	    write_model_obsfile(obsfile_ig, psyfile, (beste, fcstv, persi))
    
        if ( not TEST_SINGLE ):  ## DON"T DO FOR TEST  

            if ( ig == 0 ) : add='.1'
            if ( ig == 1 ) : add='.2'
            if ( ig == 2 ) : add='.3'
            write_mean_errors('ERRORS/GIOPSa'+add, int(datestr), CCMEP_MEAN_ERRORS, vec=['u','v'], tmp_prefix='site4/TMP/GIOPSa'+add)
            write_mean_errors('ERRORS/GIOPSa'+add, int(datestr), CCMEP_MMEAN_ERRORS, vec=['s','a'], tmp_prefix='site4/TMP/GIOPSa'+add)
            CCMEP_MEANU_ERRORS=np.concatenate((np.array([mean_ERRCn[0,1], mean_ERRCb[0,1]]), mean_ERRCp[0,:,1].flatten(), mean_ERRCf[0,:,1].flatten()))
            CCMEP_MEANV_ERRORS=np.concatenate((np.array([mean_ERRCn[1,1], mean_ERRCb[1,1]]), mean_ERRCp[1,:,1].flatten(), mean_ERRCf[1,:,1].flatten()))
            CCMEP_MEANS_ERRORS=np.concatenate((np.array([mean_ERMCn[0,1], mean_ERMCb[0,1]]), mean_ERMCp[0,:,1].flatten(), mean_ERMCf[0,:,1].flatten()))
            CCMEP_MEANA_ERRORS=np.concatenate((np.array([mean_ERMCn[1,1], mean_ERMCb[1,1]]), mean_ERMCp[1,:,1].flatten(), mean_ERMCf[1,:,1].flatten()))
            CCMEP_RMSEU_ERRORS=np.concatenate((np.array([rmse_ERRCn[0,1], rmse_ERRCb[0,1]]), rmse_ERRCp[0,:,1].flatten(), rmse_ERRCf[0,:,1].flatten()))
            CCMEP_RMSEV_ERRORS=np.concatenate((np.array([rmse_ERRCn[1,1], rmse_ERRCb[1,1]]), rmse_ERRCp[1,:,1].flatten(), rmse_ERRCf[1,:,1].flatten()))
            CCMEP_RMSES_ERRORS=np.concatenate((np.array([rmse_ERMCn[0,1], rmse_ERMCb[0,1]]), rmse_ERMCp[0,:,1].flatten(), rmse_ERMCf[0,:,1].flatten()))
            CCMEP_RMSEA_ERRORS=np.concatenate((np.array([rmse_ERMCn[1,1], rmse_ERMCb[1,1]]), rmse_ERMCp[1,:,1].flatten(), rmse_ERMCf[1,:,1].flatten()))
            datadatefile.add_to_file(int(datestr), CCMEP_MEANU_ERRORS, file='ERRORS/GIOPS'+add+'_meanu.dat', tmpfile='site4/TMP/GIOPS'+add+'_meanu.dat')
            datadatefile.add_to_file(int(datestr), CCMEP_MEANV_ERRORS, file='ERRORS/GIOPS'+add+'_meanv.dat', tmpfile='site4/TMP/GIOPS'+add+'_meanv.dat')
            datadatefile.add_to_file(int(datestr), CCMEP_MEANS_ERRORS, file='ERRORS/GIOPS'+add+'_means.dat', tmpfile='site4/TMP/GIOPS'+add+'_means.dat')
            datadatefile.add_to_file(int(datestr), CCMEP_MEANS_ERRORS, file='ERRORS/GIOPS'+add+'_meana.dat', tmpfile='site4/TMP/GIOPS'+add+'_meana.dat')
            datadatefile.add_to_file(int(datestr), CCMEP_RMSEU_ERRORS, file='ERRORS/GIOPS'+add+'_rmseu.dat', tmpfile='site4/TMP/GIOPS'+add+'_rmseu.dat')
            datadatefile.add_to_file(int(datestr), CCMEP_RMSEV_ERRORS, file='ERRORS/GIOPS'+add+'_rmsev.dat', tmpfile='site4/TMP/GIOPS'+add+'_rmsev.dat')
            datadatefile.add_to_file(int(datestr), CCMEP_RMSES_ERRORS, file='ERRORS/GIOPS'+add+'_rmses.dat', tmpfile='site4/TMP/GIOPS'+add+'_rmses.dat')
            datadatefile.add_to_file(int(datestr), CCMEP_RMSEA_ERRORS, file='ERRORS/GIOPS'+add+'_rmsea.dat', tmpfile='site4/TMP/GIOPS'+add+'_rmsea.dat')
    
    return	
    	
date_start = datetime.datetime(2019, 07, 01)
date_final = datetime.datetime(2019, 11, 01)

def plot_errors(EXPTS, date_range=(date_start, date_final)):
    (date_min, date_max) = date_range
    plt.rc('font', family='serif')
    plt.rc('text', usetex=True)
    file='ERRORS/GIOPS_rmseu.dat'
    intdate, errors = datadatefile(file)
    dates = datadatefile.convert_strint_datelist(intdate)
    
    new_dates=[]
    new_errors=[]
    
    for idate,date in dates:
        if ( ( date >= date_min ) and (date <= date_max ) ):
	    new_date.append(date)
	    new_error.append(errors[:,idate])
	  
    dates = new_dates
    errors = np.transpose( np.array(new_errors) )
    
    errors_list = errors.tolist()
    tfig, taxe = plt.subplots()
    fcst = []
    pers = []
    for ierr, error in enumerate(errors_list):
        if ( ierr == 0 ):  # BEST ESTIMATE
           clr='k'
           best=np.mean(error)
        elif ( ierr == 1 ):  # CLIMATATOLOGY:  IGNORE
           clr='w'	
           clim=np.mean(error)
        elif ( ierr <= 11 ): # FORECAST
           clr='b'
           fcst.append(np.mean(error))	
        elif ( ierr <= 21 ): # PERSIST
           clr='g'
           pers.append(np.mean(error))
        if ( ierr != 1 ):	
           taxe.plot(dates, error, color=clr)
    afig.save('PLOTS/u_time.png')

    lfig, laxe = plt.subplots()
    fcst.insert(0, best)
    pers.insert(0, best)
    laxe.plot(range(len(fcst)), fcst, color='b')
    laxe.plot(range(len(pers)), pers, color='g')
    lfig.save('PLOTS/u_lead.png')
    
    return

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

    return (LONO, LATO, depth), (obser, beste, persi, fcstv)    

def write_model_obsfile(obsfile, tplfile, (beste, fcstv, persi)):
    # cp template file to new obsfile.  Then enter model data.
    shutil.copy(tplfile, obsfile)
    obsset = netCDF4.Dataset(obsfile,mode='r+')
    obsset['best_estimate'][:]=beste
    obsset['forecast'][:]=fcstv
    obsset['persistence'][:]=persi
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
    


def write_mean_errors(out_prefix, dateint, MEAN_ERRORS, vec=['u','v'], tmp_prefix='NULL'):

    if ( isinstance(dateint, datetime.datetime) or isinstance(dateint, datetime.date) ): dateint=date.strftime("%Y%m%d")
    if ( isinstance(dateint, str) ): dateint=int(dateint)    
    (mean_neste, rmse_neste), (mean_beste, rmse_beste), (mean_fcstv, rmse_fcstv), (mean_persi, rmse_persi) = MEAN_ERRORS
    #KEEP 15m errors
    MEANU_ERRORS=np.concatenate((np.array([mean_neste[0,1], mean_beste[0,1]]), mean_persi[0,:,1].flatten(), mean_fcstv[0,:,1].flatten()))
    MEANV_ERRORS=np.concatenate((np.array([mean_neste[1,1], mean_beste[1,1]]), mean_persi[1,:,1].flatten(), mean_fcstv[1,:,1].flatten()))
    RMSEU_ERRORS=np.concatenate((np.array([rmse_neste[0,1], rmse_beste[0,1]]), rmse_persi[0,:,1].flatten(), rmse_fcstv[0,:,1].flatten()))
    RMSEV_ERRORS=np.concatenate((np.array([rmse_neste[1,1], rmse_beste[1,1]]), rmse_persi[1,:,1].flatten(), rmse_fcstv[1,:,1].flatten()))

    file1=out_prefix+'_mean'+vec[0]+'.dat'
    file2=out_prefix+'_mean'+vec[1]+'.dat'   
    file3=out_prefix+'_rmse'+vec[0]+'.dat'
    file4=out_prefix+'_rmse'+vec[1]+'.dat'   
    
    if ( tmp_prefix == 'NULL' ): tmp_prefix=out_prefix
    tile1=tmp_prefix+'_mean'+vec[0]+'.tmp'
    tile2=tmp_prefix+'_mean'+vec[1]+'.tmp'
    tile3=tmp_prefix+'_rmse'+vec[0]+'.tmp'
    tile4=tmp_prefix+'_rmse'+vec[1]+'.tmp'
     
    datadatefile.add_to_file(dateint, MEANU_ERRORS, file=file1, tmpfile=tile1)
    datadatefile.add_to_file(dateint, MEANV_ERRORS, file=file2, tmpfile=tile2)
    datadatefile.add_to_file(dateint, RMSEU_ERRORS, file=file3, tmpfile=tile3)
    datadatefile.add_to_file(dateint, RMSEV_ERRORS, file=file4, tmpfile=tile4)
    
    return
   
psyfile='CLASS4_currents/class4_20211031_PSY4V3R1_orca12_currents.nc'
obsfile='CLASS4_currents_CCMEP/class4_20211031_GIOPS_orca025_currents.nc'
def write_mean_errors_from_obsfile(date, indir, insuffix, out_prefix, tmp_prefix='NULL'):
    if ( isinstance(date, datetime.datetime) or isinstance(date, datetime.date) ):  datestr=date.strftime("%Y%m%d")
    if ( isinstance(date, int) ): datestr=str(date)
    if ( isinstance(date, str) ): datestr=date
    obsfile=indir+'/class4_'+datestr+'_'+insuffix+'.nc'
    (LONO, LATO, depth), (obser, beste, persi, fcstv) = read_obsfile(obsfile)
    neste = beste.copy()  # NESTE is actually the best estimate from the day AFTER the observation.  

    Mobser = 0.0*obser.copy()
    Mneste = 0.0*neste.copy()
    Mbeste = 0.0*beste.copy()
    Mfcstv = 0.0*fcstv.copy()
    Mpersi = 0.0*persi.copy()

    Mobser[:,0,:], Mobser[:,1, :] = speed_and_angle(obser[:,0,:], obser[:,1,:])
    Mneste[:,0, ], Mneste[:,1, :] = speed_and_angle(neste[:,0,:], neste[:,1,:])
    Mbeste[:,0, ], Mbeste[:,1, :] = speed_and_angle(beste[:,0,:], beste[:,1,:])
    Mfcstv[:,0,:,:], Mfcstv[:,1, :] = speed_and_angle(fcstv[:,0,:,:], fcstv[:,1,:,:])
    Mpersi[:,0,:,:], Mpersi[:,1, :] = speed_and_angle(persi[:,0,:,:], persi[:,1,:,:])

    MEAN_ERRORS = calc_mean_error(obser, (neste, beste, persi, fcstv))
    MMEAN_ERRORS = calc_mean_error(Mobser, (Mneste, Mbeste, Mpersi, Mfcstv))
    write_mean_errors(out_prefix, datestr, MEAN_ERRORS, vec=['u','v'], tmp_prefix=tmp_prefix)
    write_mean_errors(out_prefix, datestr, MMEAN_ERRORS, vec=['s','a'], tmp_prefix=tmp_prefix)
    return
    
def write_mean_errors_from_obsfiles(dates, indir, insuffix, out_prefix, tmp_prefix='NULL'):
    if ( isinstance(dates, str) ):
        if ( ( dates == 'all') or (dates == 'ALL') or ( dates == 'ls' ) ):
	    datestr='????????'
	    obsfiles = sorted(glob.glob(indir+'/class4_'+datestr+'_'+insuffix+'.nc'))
	    date_list = []
	    for obsfile in obsfiles:
	        date_start = obsfile.find('class4_')+7
		date_final = date_start+8
		date_list.append(obsfile[date_start:date_final])
    elif ( isinstance(dates, list) ):
        date_list=dates
        obsfiles=[]
	for datestr in dates:  
	   obsfiles.append(indir+'/'+indir+'/class4_'+datestr+'_'+insuffix+'.nc')
	  
    for date in date_list: 
        write_mean_errors_from_obsfile(date, indir, insuffix, out_prefix, tmp_prefix='NULL')

    return

