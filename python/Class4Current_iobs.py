import os
import sys
sys.path.insert(0, '/home/dpe000/Class4_Currents/python')
#from importlib import reload

import numpy as np
import scipy.interpolate as si
import time
import datetime
import shutil 
import netCDF4

import Class4Current
import find_fcst_file
import get_archive
import stfd
import read_grid
import find_hall
import check_date
import cplot

mask = read_grid.read_mask(var='tmask')
maskt = read_grid.read_mask(var='tmask')
masku = read_grid.read_mask(var='umask')
maskv = read_grid.read_mask(var='vmask')
mask0 = np.squeeze(mask[0,:,:])
e3t = read_grid.read_e3t_mesh(var='e3t_0')
e1t = read_grid.read_mesh_var('e1t')
e2t = read_grid.read_mesh_var('e2t')
nav_lon, nav_lat, grid_area = read_grid.read_coord(grid='T')

site=find_hall.get_site()
datadir='/fs/'+site+'/eccc/mrd/rpnenv/dpe000/'
tempdir=datadir+'/'+'tmpdir'

missing=-999

date_ic2 = Class4Current.date_ic2
date_ic3 = Class4Current.date_ic3
date_mg1 = Class4Current.date_mg1

def interp_direct_obs( lonlat, field, lonlat_obs, method='2sweep', fill_value=missing):
    if ( isinstance(field, list) ):
        fld_at_obs = []
        for a_field in field:
            a_fld_at_obs, IMISS = interp_direct_obs( lonlat, a_field, lonlat_obs, method=method, fill_value=fill_value)
            fld_at_obs.append(a_fld_at_obs)
        return fld_at_obs, IMISS
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
    return fld_at_obs, IMISS

def process_uvanalfiles(files, thedate, LONLAT_OBS):
    LONO, LATO = LONLAT_OBS
    fileu, filev = files
    UW, LONU, LATU, LEVU, DAT = Class4Current.read_anal_file(fileu, grid='U', dates=[thedate]); UW=np.squeeze(UW)
    VW, LONV, LATV, LEVV, DAT = Class4Current.read_anal_file(filev, grid='V', dates=[thedate]); VW=np.squeeze(VW)
    UU, LONUU, LATUU = Class4Current.map_to_A_grid(np.squeeze(UW), LONU, LATU, grid='U')
    VV, LONVV, LATVV = Class4Current.map_to_A_grid(np.squeeze(VW), LONV, LATV, grid='V')
    # Important for U15 calculations that arrays not be masked.
    if ( isinstance(UU, np.ma.core.MaskedArray) ): UU = UU.data
    if ( isinstance(VV, np.ma.core.MaskedArray) ): VV = VV.data
    [U00, V00] = [UU[0,:,:], VV[0,:,:]]
    [U15, V15] = Class4Current.calc_m15([UU, VV], e3t, mask)
    [U00, V00] = Class4Current.UU_ORCA_to_NE([U00, V00])
    [U15, V15] = Class4Current.UU_ORCA_to_NE([U15, V15])
    # this actually removes non-ocean points
    (U00m, V00m, U15m, V15m), (LONM, LATM) = Class4Current.msk_flds([U00, V00, U15, V15], [nav_lon, nav_lat], mask0)
    # now grid to observation locations
    (UOOo, V00o, U15o, V15o), IMISS = interp_direct_obs( (LONM, LATM), [U00m, V00m, U15m, V15m], [LONO, LATO], method='2sweeplinear')
    return (UOOo, V00o), (U15o, V15o), IMISS 

def process_uvfcstfiles(file_here,  LONLAT_OBS, file_altn=None):
    LONO, LATO = LONLAT_OBS
    file_3hrs = None
    file_fcst = None
    get_from_A = True
    if ( isinstance(file_here, list) ): 
        file_3hrs = file_here
        file_fcst = None
        get_from_A = False
    else:
        file_3hrs = None
        file_fcst = file_here
        get_from_A = True
    #    
    if ( isinstance(file_3hrs, type(None)) and isinstance(file_altn, list) ):
        file_3hrs = file_altn
    if ( isinstance(file_fcst, type(None)) and isinstance(file_altn, str) ):
        file_fcst = file_altn
    #        
    # IF 3D files
    U00A, V00A, U15A, V15A = None, None, None, None
    if ( file_fcst ):
        try:
            leu, UF = stfd.read_fstd_multi_lev(file_fcst, 'UUW',vfreq=24, typvar='P@')
            lev, VF = stfd.read_fstd_multi_lev(file_fcst, 'VVW',vfreq=24, typvar='P@')
            # Important for U15 calculations that arrays not be masked.
            ## BUT NEED THE INTERPOLATED MASK -- which isn't actually masked.
            if ( isinstance(UF, np.ma.core.MaskedArray) ): UF = UF.data
            if ( isinstance(VF, np.ma.core.MaskedArray) ): VF = VF.data
            [U00A, V00A] = [UF[0,:,:], VF[0,:,:]]
            [U15A, V15A] = Class4Current.calc_m15([UF, VF], e3t, mask)
            print("SUCCESSFUL retrieval U15 from U3D")
        except:
            U00A, V00A, U15A, V15A = None, None, None, None
            print("UNSUCCESSFUL retrieval U15 from U3D")
    # IF 3H 2D files
    U00B, V00B, U15B, V15B = None, None, None, None
    if (file_3hrs):
        try:
            U15B, V15B = Class4Current.calc_u15_from_3hr(file_3hrs)
            print("SUCCESSFUL retrieval U15 from U10/20")
        except:
            try:
                U15B, V15B, T15B = Class4Current.calc_u00_from_3hr(file_3hrs, 15.0)
                print("SUCCESSFUL retrieval U15 from U15")
            except:
                U15B, V15B = None, None
                print("UNSUCCESSFULL retrieval of U15 from 3hr files")
        try:
            U00B, V00B, T00B = Class4Current.calc_u00_from_3hr(file_3hrs, 0.0)
        except:
            U00B, V00B, T00B = None, None, None
    #
    U00, V00, U15, V15 = None, None, None, None
    ## THIS SHOULD JUST ABOUT TAKE ALL INSTANCES INTO CONSIDERATION
    if ( get_from_A and (not isinstance(U15A, type(None) ) ) ): 
        print('SETTING VELOCITY FROM A')
        U00, V00, U15, V15 = U00A, V00A, U15A, V15A
    if ( ( not get_from_A ) and ( not isinstance(U00B, type(None)) ) ):
        print('SETTING 00 VELOCITY FROM B')
        U00, V00 = U00B, V00B
    if ( ( not get_from_A ) and ( not isinstance(U15B, type(None)) ) ):
        print('SETTING 15 VELOCITY FROM B')
        U15, V15 = U15B, V15B
    if ( isinstance(U00, type(None)) and ( not isinstance(U00B, type(None) ) ) ):
        print('USING ALTERNATIVE VELOCITY FROM B')
        U00, V00 = U00B, V00B
    if ( isinstance(U00, type(None)) and ( not isinstance(U00A, type(None) ) ) ):
        print('USING ALTERNATIVE VELOCITY FROM A')
        U00, V00 = U00A, V00A
    if ( isinstance(U15, type(None)) and ( not isinstance(U15B, type(None) ) ) ):
        print('USING ALTERNATIVE VELOCITY FROM B')
        U15, V15 = U15B, V15B
    if ( isinstance(U15, type(None)) and ( not isinstance(U15A, type(None) ) ) ):
        print('USING ALTERNATIVE VELOCITY FROM A')
        U15, V15 = U15A, V15A
    U00o, V00o, U15o, V15o = None, None, None, None
    try:
        # this actually removes non-ocean points
        (U00m, V00m, U15m, V15m), (LONM, LATM) = Class4Current.msk_flds([U00, V00, U15, V15], [nav_lon, nav_lat], mask0)
        # now grid to observation locations
        (U00o, V00o, U15o, V15o), IMISS = interp_direct_obs( (LONM, LATM), [U00m, V00m, U15m, V15m], [LONO, LATO], method='2sweeplinear')
        print('OBSERVATION GRID SIZES', U15o.shape, V15o.shape)
    except:
        print('UNSUCCESSFUL INTERPOLATION TO OBSERVATIONS')
    
    return (U00o, V00o), (U15o, V15o), IMISS

def preretrieve_forecast(thisdate, bdays, ffreq=3):
    # GET A WHOLE BUNCH OF FILES AHEAD OF REQUIRED TIME:  LESS CALLS TO RARC?  From thisdate -- back bdays (length of forecast)
    SYS='NUL'    
    if ( thisdate-datetime.timedelta(days=bdays) > date_ic2 ): 
        SYS='OPD'
        tmpdir=tempdir+'/'+SYS+'/'
        branch='operation.forecasts.giops.prog.glboce'
    if ( thisdate < date_ic2 ): 
        SYS='PSD'
        tmpdir=tempdir+'/'+SYS+'/'
        branch='parallel.forecasts.giops.prog.glboce'
    if ( SYS != 'NUL'):
        rc = get_archive.get_archive_leads(tmpdir, branch, [thisdate-datetime.timedelta(days=bdays), thisdate], (np.arange(0,24*bdays,ffreq)+ffreq).tolist(), ensnum=None, execute=True)
    return

def process_giops(datein, CHARLY=True, filter=True):
    date=check_date.check_date(datein, outtype=datetime.datetime)
    datestr=date.strftime('%Y%m%d')   
    time0=time.time()

    if ( not CHARLY ):
        psyfile='CLASS4_currents/class4_'+datestr+'_PSY4V3R1_orca12_currents.nc'
        obsfile='CLASS4_currents_CCMEP/class4_'+datestr+'_GIOPS_orca025_currents.n.nc'
    elif ( CHARLY ):
      if ( not filter ): 
        psyfile='CLASS4_currents_CHARLY/class4_'+datestr+'_PSY4V3R1_orca12_currents.nc'
        obsfile='CLASS4_currents_CCMEP_UFIL/class4_'+datestr+'_GIOPS_orca025_currents.n.nc'
      elif ( filter ): 
        psyfile='CLASS4_currents_CHARLY/class4_'+datestr+'_PSY4V3R1_orca12_currents-filtr.nc'
        obsfile='CLASS4_currents_CCMEP_FILT/class4_'+datestr+'_GIOPS_orca025_currents.fn.nc'

    (LONO, LATO, depth), (obser, beste, fcstv, persi) = Class4Current.read_obsfile(psyfile)
    tidev, stokv = Class4Current.read_stokes_variables(psyfile)
    inite = beste.copy()
    MERCATOR_MEAN_ERRORS = Class4Current.calc_mean_error(obser, [Class4Current.adjust_velocities(field, (tidev, stokv), ierror=1) for field in (beste, inite, fcstv, persi)])
    nobss, nvars, nfcst, ndeps = fcstv.shape

    ## COPY OF OBS FILE VARIABLES FOR EACH OUTPUT
    beste = 0.0*beste.copy()
    inite = 0.0*beste.copy()
    fcstv = 0.0*fcstv.copy()
    persi = 0.0*persi.copy()
    beste.mask = False
    inite.mask = False
    fcstv.mask = False
    persi.mask = False

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
    file_besu, __ =find_fcst_file.find_anal_file(andate, system=GD, anal=anal, var='U',execute=True)
    file_besv, __ =find_fcst_file.find_anal_file(andate, system=GD, anal=anal, var='V',execute=False)
    file_iniu, __ =find_fcst_file.find_anal_file(bedate, system=GU, anal=True, var='U',execute=True)
    file_iniv, __ =find_fcst_file.find_anal_file(bedate, system=GU, anal=True, var='V',execute=False)
    print('file_best', file_besu, file_besv)

    timen=time.time()
    timep=timen-time0
    print("TIMING :: Best Estimate Retieval Time: "+str(timep))
    time0=timen

    ## BEST ESTIMATE   
    (UOOb, V00b), (U15b, V15b), IMISS = process_uvanalfiles((file_besu, file_besv), bedate, [LONO, LATO])
    if ( ndeps == 2 ):
        BEST=[[U00b, V00b], [U15b, V15b]]
    elif ( ndeps == 1 ):
        BEST=[[U15b, V15b]]

    ## INIT ESTIMATE    
    (UOOi, V00i), (U15i, V15i), __ = process_uvanalfiles((file_iniu, file_iniv), bedate, [LONO, LATO])
    if ( ndeps == 2 ):
        INIT=[[U00i, V00i], [U15i, V15i]]
    elif ( ndeps == 1 ):
        INIT=[[U15i, V15i]]

    timen=time.time()
    timep=timen-time0
    print("TIMING :: Best Estimate Processing Time: "+str(timep))
    time0=timen
    timel=time0
    
    FCST=[]
    PERS=[]

    # GET A WHOLE BUNCH OF FILES AHEAD OF REQUIRED TIME:  LESS CALLS TO RARC?  From thisdate -- back bdays (length of forecast)
    preretrieve_forecast(bedate, nfcst, ffreq=3)

    for ifcst in range(nfcst):
        jfcst=ifcst+1
        fchour=jfcst*24
        fchours=(np.arange(fchour-24, fchour, 3)+3).tolist()
        fcdate = bedate - datetime.timedelta(days=jfcst)
        andiff = (2 - fcdate.weekday() ) % 7
        andate = fcdate + datetime.timedelta(days=andiff) 
        anal=False
        if ( andiff == 0 ): anal=True
        print( 'Persistence files', fcdate, andate)
        CHOOSE='best'
        if ( CHOOSE == 'init' ):
            GU='gu'
            if ( fcdate <= date_ic2 ):  GU='pu'
            file_besu, __ = find_fcst_file.find_anal_file(fcdate, system=GU, var='U',execute=True)
            file_besv, __ = find_fcst_file.find_anal_file(fcdate, system=GU, var='V',execute=False)
        elif ( CHOOSE == 'best' ):
            GD='gd'
            if ( andate <= date_ic2 - datetime.timedelta(days=7) ):  GD='pd'
            file_besu, __ = find_fcst_file.find_anal_file(andate, system=GD, anal=anal, var='U',execute=True)
            file_besv, __ = find_fcst_file.find_anal_file(andate, system=GD, anal=anal, var='V',execute=False)
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
        # And SO SHOULD THIS NOW BE FINDING A FILE ALREADY ON THE TMPDIR.
        file_3hrs=[]
        for fc3hour in fchours:
            file_3hr=find_fcst_file.find_fcst_file(SYS, fcdate, fc3hour, 0, src='ocn', execute=True)
            file_3hrs.append(file_3hr)
        print('file_best', file_besu, os.path.isfile(file_besu), file_besv, os.path.isfile(file_besv))
        print('file_fcst', file_fcst, os.path.isfile(file_fcst))
        for file_3hr in file_3hrs:
           print('file_3hr', file_3hr, os.path.isfile(file_3hr))  
    
        (UOOp, V00p), (U15p, V15p), __ = process_uvanalfiles((file_besu, file_besv), fcdate, [LONO, LATO])
        if ( ndeps == 2 ):
            PERS.append([[U00p, V00p],[U15p, V15p]])
        elif ( ndeps == 1 ):
            PERS.append([[U15p, V15p]])


        (U00F, V00F), (U15F, V15F), __ = process_uvfcstfiles(file_3hrs, [LONO, LATO], file_altn=file_fcst) 
        
        if ( ndeps == 2 ):
            FCST.append([[UOOF, VOOF],[U15F, V15F]])
        elif ( ndeps == 1 ):
            FCST.append([[U15F, V15F]])
        timen=time.time()
        timep=timen-timel
        timet=timen-time0
        print("TIMING :: LEAD "+str(ifcst+1)+" Processing Time: "+str(timep))
        print("TIMING :: Cumalative Forecast Processing and Retrieval Time "+str(timet))
        timel=timen
        
    #nobss, nvars, nfcst, ndeps = fcstv.shape
    for kk in range(ndeps):
        UUK, VVK = BEST[kk]
        beste[:, 0, kk] = UUK
        beste[:, 1, kk] = VVK
    for kk in range(ndeps):
        UUK, VVK = INIT[kk]
        inite[:, 0, kk] = UUK
        inite[:, 1, kk] = VVK
    for ld in range(nfcst):
        for kk in range(ndeps):
            UUK, VVK = FCST[ld][kk]
            fcstv[:, 0, ld, kk] = UUK
            fcstv[:, 1, ld, kk] = VVK
    for ld in range(nfcst):
        for kk in range(ndeps):
            UUK, VVK = PERS[ld][kk]
            persi[:, 0, ld, kk] = UUK
            persi[:, 1, ld, kk] = VVK

    write_model_obsfile(obsfile, psyfile, (beste, inite, fcstv, persi))

    CCMEP_MEAN_ERRORS = Class4Current.calc_mean_error(obser, [ Class4Current.adjust_velocities(field, (tidev, stokv), ierror=1) for field in (beste, inite, fcstv, persi)])
    
    for iierror in range(4):
        if ( iierror == 0 ): print('best errors')
        if ( iierror == 1 ): print('init errors')
        if ( iierror == 2 ): print('fcst errors')
        if ( iierror == 3 ): print('pers errors')
        
        print('Mercator', MERCATOR_MEAN_ERRORS[iierror])
        print('CCMER', CCMEP_MEAN_ERRORS[iierror])
    
    iidepth=0
    if ( ndeps == 2 ): iidepth1
    best_rmse = np.squeeze(np.sqrt( np.sum( np.square( obser - Class4Current.adjust_velocities(beste, (tidev, stokv), ierror=1)), axis=1 ) )[:, iidepth])
    
    IGOOD = [ iobs for iobs in range(nobss) if iobs not in IMISS[0] ]
    
    suptitle='0m'+datestr+'Native Grid Init Estimate'
    outfile='PLOTS/scatter_'+datestr+'.png'
    cplot.scatterdots([(LONO[IGOOD], LATO[IGOOD], best_rmse[IGOOD], 'b'), (LONO[IMISS], LATO[IMISS], best_rmse[IMISS],'r')], labels=['GOOD', 'EXTRA'], make_global=True, suptitle=suptitle, outfile=outfile, legend_title=['', 'velocity (m/s)'])
    outfile='PLOTS/badpoints_'+datestr+'.png'
    cplot.scatterdots([(LONO[IMISS], LATO[IMISS], best_rmse[IMISS],'r')], labels=['GOOD', 'EXTRA'], make_global=True, suptitle=suptitle, outfile=outfile, legend_title=['', 'velocity (m/s)'])
    print(IMISS)
    for imiss in IMISS:
       print(LON0(imiss), LAT0(imis), best_rmse(imiss))
    return
    
def write_model_obsfile(obsfile, tplfile, fields):
    (beste, inite, fcstv, persi)=fields
    # cp template file to new obsfile.  Then enter model data.
    shutil.copy(tplfile, obsfile)
    obsset = netCDF4.Dataset(obsfile,mode='r+')

    inite_var = obsset.createVariable('init_estimate', np.float32,
                                      ('numobs', 'numvars', 'numdeps'),
                                      fill_value=obsset['best_estimate']._FillValue)
    inite_var.long_name = "Model Initial Estimate"
    inite_var.units = "m s-1"

    obsset['best_estimate'][:]=beste
    obsset['forecast'][:]=fcstv
    obsset['persistence'][:]=persi
    inite_var[:] = inite
    obsset.close()
    return
