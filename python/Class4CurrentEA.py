import os
import sys
sys.path.insert(0, '/home/dpe000/Class4_Currents/python')
#from importlib import reload
       
import netCDF4
import shutil 
import numpy as np
import datetime
import time
import subprocess
import glob

import multiprocessing
import itertools
from functools import partial
num_cpus = len(os.sched_getaffinity(0))

import stfd
import read_grid
import isoheatcontent
import find_value_at_point
import Class4Current
import find_fcst_file
import shapiro

KCONV=273.16

mask = read_grid.read_mask(var='tmask')
maskt = read_grid.read_mask(var='tmask')
masku = read_grid.read_mask(var='umask')
maskv = read_grid.read_mask(var='vmask')
mask0 = np.squeeze(mask[0,:,:])
e3t = read_grid.read_e3t_mesh(var='e3t_0')
e1t = read_grid.read_mesh_var('e1t')
e2t = read_grid.read_mesh_var('e2t')
nav_lon, nav_lat, grid_area = read_grid.read_coord(grid='T')

mdir5='/fs/site5/eccc/mrd/rpnenv/dpe000/maestro_archives'
mdir6='/fs/site5/eccc/mrd/rpnenv/dpe000/maestro_archives'
hdir5='/fs/site5/eccc/mrd/rpnenv/dpe000/maestro_hpcarchives'
hdir6='/fs/site5/eccc/mrd/rpnenv/dpe000/maestro_hpcarchives'
eg_anal_file=mdir5+'/GIOPS_T0/SAM2/20210602/DIA/ORCA025-CMC-ANAL_1d_grid_T_2021060200.nc'
tate=datetime.datetime(2022, 1, 5)

date_ic2 = datetime.datetime(2019,7,3)
date_ic3 = datetime.datetime(2021,12,1,12)
date_mg1 = datetime.datetime(2020,1,21,12)

def get_anal_dates():
    dates = Class4Current.create_dates(20210602, 20220601, date_inc=7)
    return dates

dates_AVAIL=get_anal_dates()

def decide_anal(anal):
    if ( not isinstance(anal, bool) ):
        if ( isinstance(anal, str) ):
            if ( ( anal == 'ANAL' ) or ( anal == 'anal' ) ): anal=True
            if ( ( anal == 'TRIAL' ) or ( anal == 'trial' ) ): anal=False
    if ( not isinstance(anal, bool) ):
        print('Do Not Understand anal = ', anal)
        return 'Null'
    if ( anal ): ANAL='ANAL'
    if ( not anal ): ANAL='TRIAL'
    return ANAL

def find_manal_file(andate, ens, var='T', anal=True, expt='GIOPS_T', ddir=mdir5+'/'):
    andate=Class4Current.check_date(andate, outtype=datetime.datetime)
    file='Null'
    date_str=Class4Current.check_date(andate)
    date_s10=date_str+'00'
    estr=str(ens)
    ANAL=decide_anal(anal)
    SDIR=expt+estr+'/SAM2/'+date_str+'/DIA'
    GRID='grid_'+var
    if ( var == 'S' ): GRID='grid_T'
    bile='ORCA025-CMC-'+ANAL+'_1d_'+GRID+'_'+date_s10+'.nc'
    file=ddir+'/'+SDIR+'/'+bile
    # NOW CHECK IF FILE EXISTS
    if ( os.path.isfile(file) ):
        pass
    else:
        #CHECK IF FILE SHOULD EXIST
        if ( andate in dates_AVAIL ):
            print('File SHOULD BE AVAILABLE', file)
            file='MISS'
        else:
            print('File NOT AVAILABLE')
            file='FNaN'
    return file

eg_fcst_file='/fs/site5/eccc/mrd/rpnenv/dpe000/maestro_archives/GEPS_STO22/gridpt/prog/ens.glboce/2022060100_360_020'
def find_mfcst_file(fcdate, fchour, ens, ddir=mdir5, expt='GEPS_STO2X'):
    fcst_file='Null'
    midpath='gridpt/prog/ens.glboce'
    date_str=Class4Current.check_date(fcdate)
    date_s10=date_str+'00'    
    lead_str=str(fchour).zfill(3)
    ensm_str=str(ens).zfill(3)
    fcst_file=ddir+'/'+expt+'/'+midpath+'/'+date_s10+'_'+lead_str+'_'+ensm_str
    # NOW CHECK IF FILE EXISTS
    if ( os.path.isfile(fcst_file) ):
        pass
    else:
        #CHECK IF FILE SHOULD EXIST
        if ( fcdate in dates_AVAIL ):
            print('File SHOULD BE AVAILABLE', fcst_file)
            fcst_file='MISS'
        else:
            print('File NOT AVAILABLE')
            fcst_file='FNaN'
    return fcst_file

def read_manal_file(file, grid='T', timname='time_instant', dates=[]):
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
    FLD, LON, LAT, lev, DAT = read_grid.read_netcdf(file, fvar, depname=depname, timname=timname, tunits='seconds', DAY0=datetime.datetime(1950, 1, 1, 0) )
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
            index = DAT.index(date)
            new_DAT.append(DAT[index])
            indices.append(index)
        DAT=new_DAT
        if ( FLD.ndim == 4 ):
            FLD=FLD[indices,:,:,:]
    LON=np.transpose(LON)
    LAT=np.transpose(LAT)
    return FLD, LON,  LAT, lev, DAT
    
def process_enan_obs(date=tate, ens_list=[0], filter=True, expt='GEPS_STO2X', nshapiro=0):

    nfcst=15
    nenss=len(ens_list)
    ens_str=str(ens_list[0]).zfill(3)
    
    datestr=date.strftime('%Y%m%d')   

    EDIR=expt
    if ( expt == 'ENAN' ): expt='GEPS_STO2X'
    if ( 'GEPS_' in expt):  EDIR=expt.replace('GEPS_','')
    if ( expt == 'GEPS_STO2X' or expt == 'GEPS_STO21' or expt == 'GEPS_STO22' ): EDIR="ENAN" 
    if ( nshapiro > 0 ): EDIR=EDIR+'_'+str(nshapiro)
    if ( not filter ): 
        obsfile1='CLASS4_currents_CCMEP_UFIL/class4_'+datestr+'_GIOPS_orca025_currents.1.nc'
        obsfile2='CLASS4_currents_CCMEP_UFIL/class4_'+datestr+'_GIOPS_orca025_currents.2.nc'
        oesfile1='CLASS4_currents_'+EDIR+'_UFIL/class4_'+datestr+'_'+EDIR+'_orca025_currents.1.'+ens_str+'.nc'
        oesfile2='CLASS4_currents_'+EDIR+'_UFIL/class4_'+datestr+'_'+EDIR+'_orca025_currents.2.'+ens_str+'.nc'
    elif ( filter ): 
        obsfile1='CLASS4_currents_CCMEP_FILT_OLD/class4_'+datestr+'_GIOPS_orca025_currents.f1.nc'
        obsfile2='CLASS4_currents_CCMEP_FILT_OLD/class4_'+datestr+'_GIOPS_orca025_currents.f2.nc'
        oesfile1='CLASS4_currents_'+EDIR+'_FILT/class4_'+datestr+'_'+EDIR+'_orca025_currents.f1.'+ens_str+'.nc'  # Nearest neighbour from ORCA025 grid
        oesfile2='CLASS4_currents_'+EDIR+'_FILT/class4_'+datestr+'_'+EDIR+'_orca025_currents.f2.'+ens_str+'.nc'  # Nearest neighbour from 0.2 lat/lon grid.
    OBSFILE_LIST = [obsfile1, obsfile2]
    OESFILE_LIST = [oesfile1, oesfile2] 

    FCSTV_LIST = []
    BESTE_LIST = []
    INITE_LIST = []
    PERSI_LIST = []
    print(obsfile1)
    (LONO, LATO, depth), (obser, beste, inite, fcstv, persi, nersi) = Class4Current.read_obsfile_plus(obsfile1)
    nobss, nvars, ofcst, ndeps = fcstv.shape
    fcstv = np.zeros((nobss, nvars, nfcst, nenss, ndeps))
    fcstv = np.ma.array(fcstv)
    inite = np.zeros((nobss, nvars, nenss, ndeps))
    inite = np.ma.array(inite)
    beste = np.zeros((nobss, nvars, ndeps))
    beste = np.ma.array(beste)
    persi = np.zeros((nobss, nvars, nfcst, nenss, ndeps))
    persi = np.ma.array(persi)
    
    for obsfile in [obsfile1, obsfile2]:
        FCSTV_LIST.append(fcstv.copy())
        BESTE_LIST.append(beste.copy())
        INITE_LIST.append(inite.copy())
        PERSI_LIST.append(persi.copy())

    (LONN, LATN) = (nav_lon, nav_lat)
    bedate=date + datetime.timedelta(days=1)
    #find days to next gd analysis
    andiff = ( 2-bedate.weekday() ) % 7
    andate = bedate + datetime.timedelta(days=andiff)
    anal=False   # retrieve TRIAL
    if ( andiff == 0 ): anal=True   # retrieve anal

    GD='gd'
    if ( andate <= date_ic2-datetime.timedelta(days=7) ): GD='pd'
    file_best, __ =find_fcst_file.find_anal_file(andate, system=GD, anal=anal, var='T',execute=True)
    file_besu, __ =find_fcst_file.find_anal_file(andate, system=GD, anal=anal, var='U',execute=False)
    file_besv, __ =find_fcst_file.find_anal_file(andate, system=GD, anal=anal, var='V',execute=False)
    print('file_best', file_best)

    ## BEST ESTIMATE    
    #TM, LONT, LATT, LEVT, DAT = Class4Current.read_anal_file(file_best, grid='T', dates=[bedate]); TM=np.squeeze(TM)
    UW, LONU, LATU, LEVU, DAT = Class4Current.read_anal_file(file_besu, grid='U', dates=[bedate]); UW=np.squeeze(UW)
    VW, LONV, LATV, LEVV, DAT = Class4Current.read_anal_file(file_besv, grid='V', dates=[bedate]); VW=np.squeeze(VW)
    UU, LONUU, LATUU = Class4Current.map_to_A_grid(np.squeeze(UW), LONU, LATU, grid='U')
    VV, LONVV, LATVV = Class4Current.map_to_A_grid(np.squeeze(VW), LONV, LATV, grid='V')
    # Important for U15 calculations that arrays not be masked.
    print(type(UU), type(VV))
    if ( isinstance(UU, np.ma.core.MaskedArray) ): UU = UU.data
    if ( isinstance(VV, np.ma.core.MaskedArray) ): VV = VV.data
    print(type(UU), type(VV))
    [U00,  V00] = [UU[0,:,:], VV[0,:,:]]
    [U15,  V15] = Class4Current.calc_m15([UU, VV], e3t, mask)
    [U00F, V00F] = Class4Current.UU_ORCA_to_NE([U00, V00])
    [U15F, V15F] = Class4Current.UU_ORCA_to_NE([U15, V15])

    if ( nshapiro > 0 ):
        time0s=time.time()
        [U00F, V00F, U15F, V15F] = do_shapiro(nshapiro, [U00F, V00F, U15F, V15F], mask0)
        timefs=time.time()-time0s
        print("TIME FOR SHAPIRO N=", nshapiro, timefs)
    # this actually removes non-ocean points
    (U00m, V00m, U15m, V15m), (LONM, LATM) = Class4Current.msk_flds([U00F, V00F, U15F, V15F], [LONN, LATN], mask0)
    # now grid to lat long
    (U00f, V00f, U15f, V15f), (LONF, LATF) = Class4Current.put_flds_latlon([U00m, V00m, U15m, V15m], [LONM, LATM], ddeg=0.2, method='2sweeplinear')
    BEST=[
          [[U15F, V15F]],
          [[U15f, V15f]]
         ]

    ## INIT ESTIMATE
    INIT=[]
    for iensm in ens_list:
        file_init = find_manal_file(andate, iensm, anal=anal, var='T')
        file_iniu = find_manal_file(andate, iensm, anal=anal, var='U')
        file_iniv = find_manal_file(andate, iensm, anal=anal, var='V')
        if ( expt == 'PDAT' ):
            file_init = find_manal_file(andate, 0, anal=anal, var='T')
            file_iniu = find_manal_file(andate, 0, anal=anal, var='U')
            file_iniv = find_manal_file(andate, 0, anal=anal, var='V')
        print('file_init', file_init, os.path.isfile(file_init), os.path.isfile(file_iniu), os.path.isfile(file_iniv))
        #TM, LONT, LATT, LEVT, DAT = read_manal_file(file_init, grid='T', dates=[bedate]); TM=np.squeeze(TM)
        UW, LONU, LATU, LEVU, DAT = read_manal_file(file_iniu, grid='U', dates=[bedate]); UW=np.squeeze(UW)
        VW, LONV, LATV, LEVV, DAT = read_manal_file(file_iniv, grid='V', dates=[bedate]); VW=np.squeeze(VW)
        UU, LONUU, LATUU = Class4Current.map_to_A_grid(np.squeeze(UW), LONU, LATU, grid='U')
        VV, LONVV, LATVV = Class4Current.map_to_A_grid(np.squeeze(VW), LONV, LATV, grid='V')
        # Important for U15 calculations that arrays not be masked.
        print(type(UU), type(VV))
        if ( isinstance(UU, np.ma.core.MaskedArray) ): UU = UU.data
        if ( isinstance(VV, np.ma.core.MaskedArray) ): VV = VV.data
        print(type(UU), type(VV))
        [U00,  V00] = [UU[0,:,:], VV[0,:,:]]
        [U15,  V15] = Class4Current.calc_m15([UU, VV], e3t, mask)
        [U00F, V00F] = Class4Current.UU_ORCA_to_NE([U00, V00])
        [U15F, V15F] = Class4Current.UU_ORCA_to_NE([U15, V15])
        if ( nshapiro > 0 ):
            time0s=time.time()
            [U00F, V00F, U15F, V15F] = do_shapiro(nshapiro, [U00F, V00F, U15F, V15F], mask0)
            timefs=time.time()-time0s
            print("TIME FOR SHAPIRO N=", nshapiro, timefs)
        # this actually removes non-ocean points
        (U00m, V00m, U15m, V15m), (LONM, LATM) = Class4Current.msk_flds([U00F, V00F, U15F, V15F], [LONN, LATN], mask0)
        # now grid to lat long
        (U00f, V00f, U15f, V15f), (LONF, LATF) = Class4Current.put_flds_latlon([U00m, V00m, U15m, V15m], [LONM, LATM], ddeg=0.2, method='2sweeplinear')
        INIT.append(
          [
            [[U15F, V15F]],
            [[U15f, V15f]]
          ]
                   )

    FCST=[]
    PERS=[]
    for ifcst in range(nfcst):
      LEAD=[]
      LEAP=[]
      jfcst=ifcst+1
      fchour=jfcst*24
      fcdate = bedate - datetime.timedelta(days=jfcst)
      andiff = (2 - fcdate.weekday() ) % 7
      andate = fcdate + datetime.timedelta(days=andiff) 
      anal=False
      if ( andiff == 0 ): anal=True
      for iensm in ens_list:
        file_fcst  = find_mfcst_file(fcdate, fchour, iensm, expt=expt)  # had optional arguements if needed.
        print('file_fcst', file_fcst, os.path.isfile(file_fcst))
        if ( file_fcst == 'FNaN' ):   # NO FORECASTS FOR THESE DAYS  ## NOTE:  REPEATED FOR ALL ENSEMBLE MEMBERS (I HOPE).
            ENSM = [ [[None, None]], [[None, None]] ]  ##  NOT SURE IF None type will work here.  Try.
        else:
          # Read in Forecast File
          LONS, LATS, TF = stfd.read_fstd_var(file_fcst, 'TM', typvar='P@')
          #leu, UF = stfd.read_fstd_multi_lev(file_fcst, 'UU2W',vfreq=24, typvar='P@')
          #lev, VF = stfd.read_fstd_multi_lev(file_fcst, 'VV2W',vfreq=24, typvar='P@')
          if ( fcdate < date_mg1 ): 
            print('FIELDS NOT AVAILABLE BEFORE', date_mg1)
            return 99
          newnorm=True
          if ( fcdate < date_ic3 ): newnorm=False
          (U15, V15) = Class4Current.calc_u15_from_3hr([file_fcst], vfreq=24, newnorm=newnorm)
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
            (U15, V15) = Class4Current.calc_u15_from_3hr([file_fcst], vfreq=24, newnorm=True)
          [U15F, V15F] =  Class4Current.UU_ORCA_to_NE([U15, V15])
        
          # NOT ENOUGH WALLTIME
          #if ( nshapiro > 0 ):
          #    [U15F, V15F] = do_shapiro(nshapiro, [U15F, V15F], mask0)
          # this actually removes non-ocean points
          (U15m, V15m), (LONM, LATM) =  Class4Current.msk_flds([U15F, V15F], [LONN, LATN], mask0)
          # now grid to lat long
          (U15f, V15f), (LONF, LATF) =  Class4Current.put_flds_latlon([U15m, V15m], [LONM, LATM], ddeg=0.2, method='2sweeplinear')
          ENSM = [ [[U15F, V15F]], [[U15f, V15f]] ]
        LEAD.append(ENSM)   ## EXTRA BRACKETS ARE HISTORICAL BAGGAGE.
        
        # NOW DO SAME FOR PERSISTENCE
        file_init = find_manal_file(andate, iensm, anal=anal, var='T')
        file_iniu = find_manal_file(andate, iensm, anal=anal, var='U')
        file_iniv = find_manal_file(andate, iensm, anal=anal, var='V')
        print('file_init', file_init, os.path.isfile(file_init), os.path.isfile(file_iniu), os.path.isfile(file_iniv), fcdate, andate)
        #TM, LONT, LATT, LEVT, DAT = read_manal_file(file_init, grid='T', dates=[fcdate]); TM=np.squeeze(TM)
        UW, LONU, LATU, LEVU, DAT = read_manal_file(file_iniu, grid='U', dates=[fcdate]); UW=np.squeeze(UW)
        VW, LONV, LATV, LEVV, DAT = read_manal_file(file_iniv, grid='V', dates=[fcdate]); VW=np.squeeze(VW)
        UU, LONUU, LATUU = Class4Current.map_to_A_grid(np.squeeze(UW), LONU, LATU, grid='U')
        VV, LONVV, LATVV = Class4Current.map_to_A_grid(np.squeeze(VW), LONV, LATV, grid='V')
        # Important for U15 calculations that arrays not be masked.
        print(type(UU), type(VV))
        if ( isinstance(UU, np.ma.core.MaskedArray) ): UU = UU.data
        if ( isinstance(VV, np.ma.core.MaskedArray) ): VV = VV.data
        print(type(UU), type(VV))
        [U00,  V00] = [UU[0,:,:], VV[0,:,:]]
        [U15,  V15] = Class4Current.calc_m15([UU, VV], e3t, mask)
        [U00F, V00F] = Class4Current.UU_ORCA_to_NE([U00, V00])
        [U15F, V15F] = Class4Current.UU_ORCA_to_NE([U15, V15])
        # NOT ENOUGH WALLTIME
        #if ( nshapiro > 0 ):
        #    [U00F, V00F, U15F, V15F] = do_shapiro(nshapiro, [U00F, V00F, U15F, V15F], mask0)
        # this actually removes non-ocean points
        (U00m, V00m, U15m, V15m), (LONM, LATM) = Class4Current.msk_flds([U00F, V00F, U15F, V15F], [LONN, LATN], mask0)
        # now grid to lat long
        (U00f, V00f, U15f, V15f), (LONF, LATF) = Class4Current.put_flds_latlon([U00m, V00m, U15m, V15m], [LONM, LATM], ddeg=0.2, method='2sweeplinear')
        ## WHY DO I USE MASK FIELDS FOR GDPS AND ORCA FIELDS FOR GEPS???
        ENSP = [ [[U15F, V15F]], [[U15f, V15f]] ]
        LEAP.append(ENSP)
      FCST.append(LEAD)   ## INDICES:  nfcst=16, nens=21, ngrid=2, ideps=1, nflds=2
      PERS.append(LEAP)   ## INDICES:  nfcst=16, nens=21, ngrid=2, ideps=1, nflds=2
        
    for iobs in range(nobss):
        time0o = time.time()
        [LONP, LATP] = LONO[iobs], LATO[iobs]
        # NOTE:  GRIDS ARE 0: ORCA025, 1: 0.2 LAT/LON
        IJPTS = Class4Current.find_nearest_points((LONP, LATP), [(LONN, LATN), (LONF, LATF)])
        IJPT, IJPF = IJPTS
        print('IJPTS', IJPTS)
        
        for ig, GRIB in enumerate(BEST):
          IJPT = IJPTS[ig]
          beste = BESTE_LIST[ig].copy()
          for kk, KLEB in enumerate(GRIB):
            bestl = []
            for ifld, FLB in enumerate(KLEB):
              bestl.append(FLB[IJPT])
            Upb, Vpb = bestl
            print( 'best', iobs, Upb, Vpb)
            beste[iobs, 0, kk] = Upb
            beste[iobs, 1, kk] = Vpb
          BESTE_LIST[ig] = beste.copy()

        for ie, ENSI in enumerate(INIT):
          for ig, GRII in enumerate(ENSI):
            IJPT = IJPTS[ig]
            inite = INITE_LIST[ig].copy()
            for kk, KLEI in enumerate(GRII):
              initl = []
              for ifld, FLI in enumerate(KLEI):
                initl.append(FLI[IJPT])
              Upi, Vpi = initl
              print( 'init', iobs, ie, Upi, Vpi)
              inite[iobs, 0, ie, kk] = Upi
              inite[iobs, 1, ie, kk] = Vpi
            INITE_LIST[ig] = inite.copy()

        for ld, LEAD in enumerate(FCST): 
          LEAP = PERS[ld]
          for ie, ENSM in enumerate(LEAD):
            ENSP = LEAP[ie]
            for ig, GRID in enumerate(ENSM):
              GRIP = ENSP[ig]
              IJPT = IJPTS[ig]
              fcstv = FCSTV_LIST[ig].copy()
              persi = PERSI_LIST[ig].copy()
              for kk, KLEV in enumerate(GRID):
                KLEP = GRIP[kk]
                fcstl = []
                persl = []
                for ifld, FLF in enumerate(KLEV):
                    FLP = KLEP[ifld]
                    # NOT ALL FORECAST FIELDS ARE AVAILABLE
                    if ( isinstance(FLF, type(None)) ):
                      FLFpt = np.NaN
                    else:
                      FLFpt=FLF[IJPT]
                    FLPpt=FLP[IJPT]
                    fcstl.append(FLFpt)
                    persl.append(FLPpt)
                Upf, Vpf = fcstl   ## These are the best values for level kk, grid ig, ensemble ie, lead ld
                Upp, Vpp = persl
                fcstv[iobs, 0, ld, ie, kk] = Upf
                fcstv[iobs, 1, ld, ie, kk] = Vpf
                persi[iobs, 0, ld, ie, kk] = Upp
                persi[iobs, 1, ld, ie, kk] = Vpp
              FCSTV_LIST[ig] = fcstv.copy()
              PERSI_LIST[ig] = persi.copy()

    for ig in range(len(FCSTV_LIST)):
        fcstv = FCSTV_LIST[ig].copy()
        persi = PERSI_LIST[ig].copy()
        beste = BESTE_LIST[ig].copy()
        inite = INITE_LIST[ig].copy()
        obsfile = OBSFILE_LIST[ig]
        oesfile = OESFILE_LIST[ig]
        print(obsfile, oesfile, fcstv.shape)
        rc = write_model_obsfile_ensemble_analysis(oesfile, obsfile, beste, inite, fcstv, persi, clobber=True)

    return 0

#/fs/site5/eccc/cmd/e/kch001/maestro_archives/drz_tests/drz_toe_DS_smooth_23/SAM2/20220309/DIA/ORCA025-CMC-TRIAL_1d_gridT-RUN-crs_20220302-20220308.nc
#/fs/site5/eccc/cmd/e/kch001/maestro_archives/drz_tests/drz_toe_DS_smooth_23/SAM2/20220309/DIA/ORCA025-CMC-TRIAL_1d_gridU-RUN-crs_20220302-20220308.nc
#/fs/site5/eccc/cmd/e/kch001/maestro_archives/drz_tests/drz_toe_DS_smooth_23/SAM2/20220309/DIA/ORCA025-CMC-TRIAL_1d_gridV-RUN-crs_20220302-20220308.nc
#/fs/site5/eccc/cmd/e/kch001/maestro_archives/drz_tests/drz_toe_DS_smooth_23/SAM2/20220309/DIA/ORCA025-CMC-TRIAL_1d_grid_T_2022030900.nc
#/fs/site5/eccc/cmd/e/kch001/maestro_archives/drz_tests/drz_toe_DS_smooth_23/SAM2/20220309/DIA/ORCA025-CMC-TRIAL_1d_grid_U_2022030900.nc
#/fs/site5/eccc/cmd/e/kch001/maestro_archives/drz_tests/drz_toe_DS_smooth_23/SAM2/20220309/DIA/ORCA025-CMC-TRIAL_1d_grid_V_2022030900.nc
eg_file='/fs/site5/eccc/mrd/rpnenv/dpe000/maestro_archives/GIOPS_T0/SAM2/20220601/DIA/ORCA025-CMC-TRIAL_1d_grid_U_2022060100.nc'
# BUT THERE MIGHT BE daily 2d 16m depth data here
e2_file='/fs/site5/eccc/mrd/rpnenv/dpe000/maestro_archives/GIOPS_T0/SAM2/20220601/DIA/ORCA025-CMC-TRIAL_1d_grid_U_2D_2022060100.nc'
er_file='/fs/site5/eccc/mrd/rpnenv/roy000/maestro/giops/gx_dev-3.4.0_DA_y2_MDT2/hub/gridpt/trial_netcdf/2019031300/ORCA025-CMC-TRIAL_1d_grid_U.nc'
def construct_anal_file(date, expt, var, is3D=True, anal=True, freq='1d', ddir=mdir5):
    trial='TRIAL'
    if ( anal ): trial='ANAL'
    grid='grid_'+var
    if ( not is3D ):
        grid=grid+'_2D'
    datestr=Class4Current.check_date(date)
    if ( not 'roy000/maestro/' in ddir ):
        file=ddir+'/'+expt+'/SAM2/'+datestr+'/DIA/'+'ORCA025-CMC-'+trial+'_'+freq+'_'+grid+'_'+datestr+'00.nc'
    else:
        file=ddir+'/'+expt+'/hub/gridpt/trial_netcdf/'+datestr+'00/ORCA025-CMC-'+trial+'_'+freq+'_'+grid+'.nc'
    return file
 
def construct_anal_stfd(date, expt, ddir=mdir5):
    datestr=Class4Current.check_date(date)
    dpth=ddir+'/'+expt+'/hub/sorties_rpn/NEMO/'
    file=dpth+datestr+'00_000'
    return file
   
def process_anal_obs(date=tate, filter=True, expt='GIOPS_T0', ddir=mdir5, nshapiro=0):

    nfcst=0
    datestr=date.strftime('%Y%m%d')   

    EDIR=expt
    if ( not filter ): 
        psyfile='CLASS4_currents_CHARLY/class4_'+datestr+'_PSY4V3R1_orca12_currents.nc'
        obsfile1='CLASS4_currents_CCMEP_UFIL_JUNE/class4_'+datestr+'_GIOPS_orca025_currents.1.nc'
        obsfile2='CLASS4_currents_CCMEP_UFIL_JUNE/class4_'+datestr+'_GIOPS_orca025_currents.2.nc'
    elif ( filter ): 
        psyfile='CLASS4_currents_CHARLY/class4_'+datestr+'_PSY4V3R1_orca12_currents-filtr.nc'
        obsfile1=EDIR+'/class4_'+datestr+'_GIOPS_orca025_currents.f1.nc'
        obsfile2=EDIR+'/class4_'+datestr+'_GIOPS_orca025_currents.f2.nc'
    OBSFILE_LIST = [obsfile1, obsfile2]

    BESTE_LIST = []   # USE ANAL FIELDS (every 7 days)
    INITE_LIST = []   # USE TRIAL FIELDS (every day)
    print(psyfile, obsfile1, obsfile2)
    (LONO, LATO, depth), (obser, beste, fcstv, persi) = Class4Current.read_obsfile(psyfile)
    nobss, nvars, ndeps = beste.shape
    inite = np.zeros((nobss, nvars, ndeps))
    inite = np.ma.array(inite)
    beste = np.zeros((nobss, nvars, ndeps))
    beste = np.ma.array(beste)
    
    for obsfile in [obsfile1, obsfile2]:
        BESTE_LIST.append(beste.copy())
        INITE_LIST.append(inite.copy())

    (LONN, LATN) = (nav_lon, nav_lat)
    bedate=date + datetime.timedelta(days=1)
    #find days to next gd analysis
    andiff = ( 2-bedate.weekday() ) % 7
    andate = bedate + datetime.timedelta(days=andiff)
    anal=False   # retrieve TRIAL
    if ( andiff == 0 ): anal=True   # retrieve anal

    # THERE IS daily 3D data here
    eg_file='/fs/site5/eccc/mrd/rpnenv/dpe000/maestro_archives/GIOPS_T0/SAM2/20220601/DIA/ORCA025-CMC-TRIAL_1d_grid_U_2022060100.nc'
    # BUT THERE MIGHT BE daily 2d 15m depth data here
    e2_file='/fs/site5/eccc/mrd/rpnenv/dpe000/maestro_archives/GIOPS_T0/SAM2/20220601/DIA/ORCA025-CMC-TRIAL_1d_grid_U_2D_2022060100.nc'
    # ASSUME ALL FILES ARE ALREADY ON DISK.
    file_best = construct_anal_file(andate, expt, var='T', anal=True, is3D=True, ddir=ddir)
    file_besu = construct_anal_file(andate, expt, var='U', anal=True, is3D=True, ddir=ddir)
    file_besv = construct_anal_file(andate, expt, var='V', anal=True, is3D=True, ddir=ddir)
    file_init = construct_anal_file(andate, expt, var='T', anal=False, is3D=True, ddir=ddir)
    file_iniu = construct_anal_file(andate, expt, var='U', anal=False, is3D=True, ddir=ddir)
    file_iniv = construct_anal_file(andate, expt, var='V', anal=False, is3D=True, ddir=ddir)
    isstfd = False
    if ( not os.path.isfile(file_iniu) ):
       file_iniu = construct_anal_stfd(bedate, expt, ddir=ddir) 
       print(file_iniu, os.path.isfile(file_iniu))
       if ( os.path.isfile(file_iniu) ): isstfd = True
    if ( not os.path.isfile(file_iniv) ):
       file_iniv = construct_anal_stfd(bedate, expt, ddir=ddir) 
       print(file_iniv, os.path.isfile(file_iniv))
       if ( os.path.isfile(file_iniv) ): isstfd = True
    print('best/init', file_best, file_init)

    ## BEST ESTIMATE    
    #TM, LONT, LATT, LEVT, DAT = Class4Current.read_anal_file(file_best, grid='T', dates=[bedate]); TM=np.squeeze(TM)
    if ( os.path.isfile(file_besu) ):
        UW, LONU, LATU, LEVU, DAU = Class4Current.read_anal_file(file_besu, grid='U', dates=[bedate]); UW=np.squeeze(UW)
    else:
        DAU = None
    if ( os.path.isfile(file_besu) ):
        VW, LONV, LATV, LEVV, DAV = Class4Current.read_anal_file(file_besv, grid='V', dates=[bedate]); VW=np.squeeze(VW)
    else:
        DAV = None
    print(bedate, DAU, DAV)
    if ( DAU and DAV ):
        UU, LONUU, LATUU = Class4Current.map_to_A_grid(np.squeeze(UW), LONU, LATU, grid='U')
        VV, LONVV, LATVV = Class4Current.map_to_A_grid(np.squeeze(VW), LONV, LATV, grid='V')
        # Important for U15 calculations that arrays not be masked.
        print(type(UU), type(VV))
        if ( isinstance(UU, np.ma.core.MaskedArray) ): UU = UU.data
        if ( isinstance(VV, np.ma.core.MaskedArray) ): VV = VV.data
        print(type(UU), type(VV))
        [U00,  V00] = [UU[0,:,:], VV[0,:,:]]
        [U15,  V15] = Class4Current.calc_m15([UU, VV], e3t, mask)
        [U00F, V00F] = Class4Current.UU_ORCA_to_NE([U00, V00])
        [U15F, V15F] = Class4Current.UU_ORCA_to_NE([U15, V15])
        if ( nshapiro > 0 ):
            [U00F, V00F, U15F, V15F] = do_shapiro(nshapiro, [U00F, V00F, U15F, V15F], mask0)
        # this actually removes non-ocean points
        (U00m, V00m, U15m, V15m), (LONM, LATM) = Class4Current.msk_flds([U00F, V00F, U15F, V15F], [LONN, LATN], mask0)
        # now grid to lat long
        (U00f, V00f, U15f, V15f), (LONF, LATF) = Class4Current.put_flds_latlon([U00m, V00m, U15m, V15m], [LONM, LATM], ddeg=0.2, method='2sweeplinear')
        BEST=[
              [[U15F, V15F]],
              [[U15f, V15f]]
             ]
    else:
        BEST=[
              [[None, None]],
              [[None, None]]
             ]

    ## INIT ESTIMATE    
    #TM, LONT, LATT, LEVT, DAT = Class4Current.read_anal_file(file_init, grid='T', dates=[bedate]); TM=np.squeeze(TM)
    if ( not isstfd ):
        print('READING NETCDF FILE')
        UW, LONU, LATU, LEVU, DAU = Class4Current.read_anal_file(file_iniu, grid='U', dates=[bedate]); UW=np.squeeze(UW)
        VW, LONV, LATV, LEVV, DAV = Class4Current.read_anal_file(file_iniv, grid='V', dates=[bedate]); VW=np.squeeze(VW)
        if ( DAU and DAV ):
            UU, LONUU, LATUU = Class4Current.map_to_A_grid(np.squeeze(UW), LONU, LATU, grid='U')
            VV, LONVV, LATVV = Class4Current.map_to_A_grid(np.squeeze(VW), LONV, LATV, grid='V')
            # Important for U15 calculations that arrays not be masked.
            print(type(UU), type(VV))
            print(UU.shape, VV.shape)
            if ( isinstance(UU, np.ma.core.MaskedArray) ): UU = UU.data
            if ( isinstance(VV, np.ma.core.MaskedArray) ): VV = VV.data
            print(type(UU), type(VV))
            [U00,  V00] = [UU[0,:,:], VV[0,:,:]]
            [U15,  V15] = Class4Current.calc_m15([UU, VV], e3t, mask)
            [U00F, V00F] = Class4Current.UU_ORCA_to_NE([U00, V00])
            [U15F, V15F] = Class4Current.UU_ORCA_to_NE([U15, V15])
            if ( nshapiro > 0 ):
                [U00F, V00F, U15F, V15F] = do_shapiro(nshapiro, [U00F, V00F, U15F, V15F], mask0)
            # this actually removes non-ocean points
            (U00m, V00m, U15m, V15m), (LONM, LATM) = Class4Current.msk_flds([U00F, V00F, U15F, V15F], [LONN, LATN], mask0)
            # now grid to lat long
            (U00f, V00f, U15f, V15f), (LONF, LATF) = Class4Current.put_flds_latlon([U00m, V00m, U15m, V15m], [LONM, LATM], ddeg=0.2, method='2sweeplinear')
            INIT=[
              [[U15F, V15F]],
              [[U15f, V15f]]
             ]
        else:
            INIT=[
              [[None, None]],
              [[None, None]]
             ]
    else:
        print('READING RPNSTD FILE')
        leu, UF = stfd.read_fstd_multi_lev(file_iniu, 'UUW', typvar='P@')
        lev, VF = stfd.read_fstd_multi_lev(file_iniv, 'VVW', typvar='P@')
        if ( isinstance(UF, np.ma.core.MaskedArray) ): UF = UF.data
        if ( isinstance(VF, np.ma.core.MaskedArray) ): VF = VF.data
        print(UF.shape, VF.shape)
        [U00, V00] = [UF[0,:,:], VF[0,:,:]]
        [U15, V15] = Class4Current.calc_m15([UF, VF], e3t, mask)
        [U00F, V00F] = Class4Current.UU_ORCA_to_NE([U00, V00])
        [U15F, V15F] = Class4Current.UU_ORCA_to_NE([U15, V15])
        if ( nshapiro > 0 ): 
            [U00F, V00F, U15F, V15F] = do_shapiro(nshapiro, [U00F, V00F, U15F, V15F], mask0)
        # this actually removes non-ocean points
        (U00m, V00m, U15m, V15m), (LONM, LATM) = Class4Current.msk_flds([U00F, V00F, U15F, V15F], [LONN, LATN], mask0)
        # now grid to lat long
        (U00f, V00f, U15f, V15f), (LONF, LATF) = Class4Current.put_flds_latlon([U00m, V00m, U15m, V15m], [LONM, LATM], ddeg=0.2, method='2sweeplinear')

        INIT=[
              [[U15, V15]],
              [[U15f, V15f]]
             ]
        
    for iobs in range(nobss):
        time0o = time.time()
        [LONP, LATP] = LONO[iobs], LATO[iobs]
        # NOTE:  GRIDS ARE 0: ORCA025, 1: 0.2 LAT/LON
        IJPTS = Class4Current.find_nearest_points((LONP, LATP), [(LONN, LATN), (LONF, LATF)])
        IJPT, IJPF = IJPTS
        #print('IJPTS', IJPTS)
        
        for ig, GRIB in enumerate(BEST):
          IJPT = IJPTS[ig]
          beste = BESTE_LIST[ig].copy()
          for kk, KLEB in enumerate(GRIB):
            bestl = []
            for ifld, FLB in enumerate(KLEB):
                if ( FLB is None ):
                    FLBpt = np.NaN
                else:
                    FLBpt = FLB[IJPT]
                bestl.append(FLBpt)
            Upb, Vpb = bestl
            #print( 'best', iobs, Upb, Vpb)
            beste[iobs, 0, kk] = Upb
            beste[iobs, 1, kk] = Vpb
          BESTE_LIST[ig] = beste.copy()

        for ig, GRII in enumerate(INIT):
          IJPT = IJPTS[ig]
          inite = INITE_LIST[ig].copy()
          for kk, KLEI in enumerate(GRII):
            initl = []
            for ifld, FLI in enumerate(KLEI):
                if ( FLI is None ):
                    FLIpt = np.NaN
                else:
                    FLIpt = FLI[IJPT]
                initl.append(FLIpt)
            Upi, Vpi = initl
            #print( 'init', iobs, Upi, Vpi)
            inite[iobs, 0, kk] = Upi
            inite[iobs, 1, kk] = Vpi
          INITE_LIST[ig] = inite.copy()

    for ig in range(len(INITE_LIST)):
        beste = BESTE_LIST[ig].copy()
        inite = INITE_LIST[ig].copy()
        obsfile = OBSFILE_LIST[ig]
        print(obsfile, psyfile, beste.shape, inite.shape)
        rc = write_model_obsfile_analysis_only(obsfile, psyfile, beste, inite, clobber=True)

    return 0

def do_shapiro(nshapiro, FLDS, mask):
    if ( nshapiro == 0 ):
        return FLDS
    FLDT = []
    return_list = True
    if ( not isinstance(FLDS, list) ): return_list = False
    if ( not isinstance(FLDS, list) ): FLDS=[FLDS]
    wtc=None
    for FLD in FLDS:
        # SEEM TO BE PASSING UNMASKED DATA
        FLDM = np.ma.array(FLD, mask=(1-mask).astype(bool))
        FLDN, wtc = shapiro.shapiro2D(FLDM,npass=nshapiro, wtc=wtc)
        FLDT.append(FLDN.data)
    if ( not return_list ): FLDT=FLDT[0]
    return FLDT
    
def write_model_obsfile_ensemble_analysis(obsfile, tplfile, beste, inite, fcstv, persi, clobber=True, fullpath='/fs/site6/eccc/mrd/rpnenv/dpe000/Class4_Currents'):

    #print('fcstv SHAPE', np.shape(fcstv))
    if ( fcstv.ndim == 5 ):
      nobs, nvars, nfcsts, nenss, ndeps = np.shape(fcstv)
    elif (fcstv.ndim == 4 ):
      nobs, nvars, nfcsts, ndeps = np.shape(fcstv)
    if ( persi.ndim == 5 ):
      __, __, nfcstp, nensp, __ = np.shape(persi)
    elif (fcstv.ndim == 4 ):
      __, __, nfcstp, __ = np.shape(persi)
    if ( inite.ndim == 4 ):
      __, __, nensi, __ = np.shape(inite)

    # ncks template file to new obsfile -- but remove forecast -- which we will need to re-create.
    # removing all three numfcsts length variables SHOULD remove numfcst dimension as well.  [THIS MAY FAIL?]
    # Using same module, we should be able to recreate persistence with a longer time line too.
    # Fullpath added for use with ssh.  Not needed for local host execution.
    if ( clobber ):
      #rc=subprocess.call(['ncks','-O','-x','-v','init_estimate,forecast,persistence,negative_persistence',fullpath+'/'+tplfile, fullpath+'/'+obsfile]) 
      #print('Job RC = ', rc)
      rc=subprocess.call(['ssh', 'ppp6', '/usr/bin/ncks','-O','-x','-v','init_estimate,forecast,persistence,negative_persistence',fullpath+'/'+tplfile, fullpath+'/'+obsfile]) 
      print('Job RC = ', rc)
      #rc=subprocess.POpen(['ssh', 'ppp6', '/usr/bin/ncks','-O','-x','-v','init_estimate,forecast,persistence,negative_persistence',fullpath+'/'+tplfile, fullpath+'/'+obsfile], stdout=subprocess.PIPE, stderr=subprocess.PIPE) 
      #print('Job RC = ', rc.returncode)
      #print(rc.stdout.read(), rc.stderr.read())
      #print(' '.join(['ssh', 'ppp6', '/usr/bin/ncks','-O','-x','-v','forecast,persistence,negative_persistence',fullpath+'/'+tplfile, fullpath+'/'+obsfile]))

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

    if ( persi.ndim == 5 ):
        persi_var = obsset.createVariable('persistence', np.float32,
                                      ('numobs', 'numvars', 'numfcsts', 'numens', 'numdeps'),
                                      fill_value=obsset['best_estimate']._FillValue)
    if ( persi.ndim == 4 ):
        persi_var = obsset.createVariable('persistence', np.float32,
                                      ('numobs', 'numvars', 'numfcsts', 'numdeps'),
                                      fill_value=obsset['best_estimate']._FillValue)
   
    persi_var.units = "m s-1"
    persi_var.long_name = "Model persistence counterpart of obs. value" ;

    if ( inite.ndim == 4 ):
        inite_var = obsset.createVariable('init_estimate', np.float32,
                                      ('numobs', 'numvars', 'numens', 'numdeps'),
                                      fill_value=obsset['best_estimate']._FillValue)
    if ( inite.ndim == 3 ):
        inite_var = obsset.createVariable('init_estimate', np.float32,
                                      ('numobs', 'numvars', 'numdeps'),
                                      fill_value=obsset['best_estimate']._FillValue)
   
    inite_var.units = "m s-1"
    inite_var.long_name = "Model initial estimate counterpart of obs. value" ;

    fcstv_var[:] = fcstv
    persi_var[:] = persi
    inite_var[:] = inite

    # ADD EXISTING (UNCHANGED DIMENSION) VARIABLES    
    obsset['best_estimate'][:]=beste
    obsset.close()
    return

def write_model_obsfile_analysis_only(obsfile, tplfile, beste, inite, clobber=True, fullpath='/fs/site6/eccc/mrd/rpnenv/dpe000/Class4_Currents'):

    nobs, nvars, ndeps = np.shape(beste)
    if ( inite.ndim == 4 ):
        __, __, nens, __ = np.shape(inite)

    # ncks template file to new obsfile -- but remove forecast -- which we will need to re-create.
    # removing all three numfcsts length variables SHOULD remove numfcst dimension as well.  [THIS MAY FAIL?]
    # Using same module, we should be able to recreate persistence with a longer time line too.
    # Fullpath added for use with ssh.  Not needed for local host execution.
    if ( clobber ):
      rc=subprocess.call(['ncks','-O','-x','-v','forecast,persistence',fullpath+'/'+tplfile, fullpath+'/'+obsfile], stdout=subprocess.PIPE, stderr=subprocess.PIPE) 
      print('Job RC = ', rc)

    obsset = netCDF4.Dataset(obsfile,mode='r+')


    if ( inite.ndim == 4 ):
        inite_var = obsset.createVariable('init_estimate', np.float32,
                                      ('numobs', 'numvars', 'numens', 'numdeps'),
                                      fill_value=obsset['best_estimate']._FillValue)
    if ( inite.ndim == 3 ):
        inite_var = obsset.createVariable('init_estimate', np.float32,
                                      ('numobs', 'numvars', 'numdeps'),
                                      fill_value=obsset['best_estimate']._FillValue)
   
    inite_var.units = "m s-1"
    inite_var.long_name = "Model initial estimate counterpart of obs. value" ;

    inite_var[:] = inite

    # ADD EXISTING (UNCHANGED DIMENSION) VARIABLES    
    obsset['best_estimate'][:]=beste
    obsset.close()
    return

def create_obsset_variable( obsset, varname, dimension_names, var_fill_value, var_units, var_long_name):
    this_var = obsset.createVariable(varname, np.float32,
                    dimension_names,
                    fill_value=var_fill_value)
    this_var.units = units
    this_var.long_name = var_long_name
    return this_var

def cat_date_ensembleEA(date, indir='CLASS4_currents_ENAN_FILT', insuffix='ENAN_orca025_currents.f2', clobber=True, nens=21):
    if ( isinstance(date, list) ):
        for idate in date:
            cat_date_ensemble(idate, indir=indir, insuffix=insuffix, clobber=clobber, nens=nens)
    datestr=check_date(date)
    obsfile=indir+'/class4_'+datestr+'_'+insuffix+'.nc'
    tplfile=indir+'/class4_'+datestr+'_'+insuffix+'.000.nc'
    cat_obsfile_ensemble(obsfile, tplfile, nens=nens, clobber=clobber)
    return

def cat_obsfile_ensembleEA(obsfile, tplfile, nens=21, VARS=['best_estimate', 'init_estimate', 'forecast', 'persistence'], ENSS=[False, True, True, True], clobber=True):
    FLDS_FULL=[]
    FLDS_MEAN=[]
    for ivar, var in enumerate(VARS):
      ENS=ENSS[ivar]
      if ( ENS ):
        for iens in range(nens):
          ensstr=str(iens).zfill(3)
          ensfile=obsfile.replace('.nc', '.'+ensstr+'.nc')
          (LONO, LATO, depth), field_sngl = Class4Current.read_obsfile_variable(ensfile, var)
          if ( iens == 0 ):
            if ( field_sngl.ndim == 5 ):
              nobs, nvars, nfcsts, nenss, ndeps = np.shape(field_sngl)
              field_full=np.zeros((nobs, nvars, nfcsts, nens, ndeps))
            elif ( field_sngl.ndim == 4 ):
              nobs, nvars, nenss, ndeps = np.shape(field_sngl)
              field_full=np.zeros((nobs, nvars, nens, ndeps))
          #endif ( iens == 0)
          if ( field_sngl.ndim == 5 ):
            field_full[:,:,:,iens,:] = field_sngl[:,:,:,0,:]
          elif ( field_sngl.ndim == 4 ):
            field_full[:,:,iens,:] = field_sngl[:,:,0,:]
        # end for iens in range(nens):
        if ( field_sngl.ndim == 5 ):
          field_mean=np.mean(field_full, axis=3)
        elif ( field_sngl.ndim == 4 ):
          field_mean=np.mean(field_full, axis=2)
      else:
        ensstr=str(0).zfill(3)
        ensfile=obsfile.replace('.nc', '.'+ensstr+'.nc')
        (LONO, LATO, depth), field_sngl = Class4Current.read_obsfile_variable(ensfile, var)
        field_full = field_sngl
        field_mean = field_sngl
      #endif ( ENS ):
      FLDS_FULL.append(field_full)
      FLDS_MEAN.append(field_mean)
    # end for ivar, var in enumerate(VARS):

    beste_full, inite_full, fcstv_full, persi_full = FLDS_FULL  ##  REDO IF ORDER EVER CHANGES
    beste_mean, inite_mean, fcstv_mean, persi_mean = FLDS_MEAN  ##  REDO IF ORDER EVER CHANGES   
    obsmile=obsfile.replace('.nc', '.enm.nc')
    write_model_obsfile_ensemble_analysis(obsfile, tplfile, beste_full, inite_full, fcstv_full, persi_full, clobber=clobber)
    write_model_obsfile_ensemble_analysis(obsmile, tplfile, beste_mean, inite_mean, fcstv_mean, persi_mean, clobber=clobber)
    return

def assemble_ensembleEA_date(date, obspre='CLASS4_currents_ENAN_FILT/class4', obssuf='ENAN_orca025_currents', iters=['f1','f2'], nens=21, clobber=True):
    datestr=Class4Current.check_date(date)
    for itstr in iters:
        obsfile=obspre+'_'+datestr+'_'+obssuf+'.'+itstr+'.nc'
        print(obsfile)
        tplfile=obsfile.replace('nc','000.nc')
        cat_obsfile_ensembleEA(obsfile, tplfile, nens=nens, clobber=clobber)

def write_model_obsfile_ensembleEA(obsfile, tplfile, FLDS, VARS, clobber=True, fullpath='/fs/site6/eccc/mrd/rpnenv/dpe000/Class4_Currents'):
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

    
