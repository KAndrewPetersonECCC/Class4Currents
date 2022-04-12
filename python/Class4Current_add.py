
date=datetime.datetime(2022,1,1)

def process_geps_obs(date=tate, Plot=False, filter=True):

    nfcst=16
    nenss=21
    
    datestr=date.strftime('%Y%m%d')   

    if ( not filter ): 
        obsfile1='CLASS4_currents_CCMEP_UFIL/class4_'+datestr+'_GIOPS_orca025_currents.1.nc'
        obsfile2='CLASS4_currents_CCMEP_UFIL/class4_'+datestr+'_GIOPS_orca025_currents.2.nc'
        oesfile1='CLASS4_currents_CCMEP_UFIL/class4_'+datestr+'_GEPS_orca025_currents.1.nc'
        oesfile2='CLASS4_currents_GEPS_UFIL/class4_'+datestr+'_GEPS_orca025_currents.2.nc'
    elif ( filter ): 
        obsfile1='CLASS4_currents_CCMEP_FILT/class4_'+datestr+'_GIOPS_orca025_currents.f1.nc'
        obsfile2='CLASS4_currents_CCMEP_FILT/class4_'+datestr+'_GIOPS_orca025_currents.f2.nc'
        oesfile1='CLASS4_currents_GEPS_FILT/class4_'+datestr+'_GEPS_orca025_currents.f1.nc'
        oesfile2='CLASS4_currents_GEPS_FILT/class4_'+datestr+'_GEPS_orca025_currents.f2.nc'
    OBSFILE_LIST = [obsfile1, obsfile2]
    OESFILE_LIST = [obsfile1, obsfile2] 

    FCSTV_LIST = []
    (LONO, LATO, depth), (obser, beste, inite, fcstv, persi, nersi) = read_obsfile_plus(obsfile1)
    nobss, nvars, ofcst, ndeps = fcstv.shape
    fcstv = np.zeros(nobss, nvars, nfcst, nenss, ndeps)
    fcstv = np.ma.array(fcstv)
    
    for obsfile in [obsfile1, obsfile2]:
        FCSTV_LIST.append(fcstv.copy())

    (LONN, LATN) = (nav_lon, nav_lat)
    bedate=date + datetime.timedelta(days=1)
    FCST=[]

    nenss=21
    for ifcst in range(nfcst):
      jfcst=ifcst+1
      fchour=jfcst*24
      fcdate = bedate - datetime.timedelta(days=jfcst)
      SYS='OP'
      tmpdir=tempdir+'/'+SYS+'/'+'O'+'/'
      branch='operation.ensemble.prog.ens.glboce'
      rc = get_archive.get_archive(tmpdir, branch, fcdate, fchour, ensnum=range(nenss), execute=True) 
      ENSM=[]
      for iensm in range(nenss):
        # THIS SHOULD NOW BE FINDING A FILE ALREADY ON THE TMPDIR.
        file_fcst     = find_fcst_file.find_fcst_file(SYS, fcdate, fchour, iensm, src='ocn', execute=True)
        print('file_fcst', file_fcst, os.path.isfile(file_fcst))        let, TF = stfd.read_fstd_multi_lev(file_fcst, 'TM',vfreq=24, typvar='P@')
        # Read in Forecast File
        LONS, LATS, TF = stfd.read_fstd_var(file_fcst, 'TM', typvar='P@')
        leu, UF = stfd.read_fstd_multi_lev(file_fcst, 'UU2W',vfreq=24, typvar='P@')
        lev, VF = stfd.read_fstd_multi_lev(file_fcst, 'VV2W',vfreq=24, typvar='P@')
	newnorm=True
	if ( fcdate < datetime.datetime(2021,12,1,12) ): newnorm=False
        (U15, V15) = calc_u15_from_3hr([file_fcst], vfreq=24, newnorm=newnorm)
	if ( np.max(U15) > 20.0 ):  
	    print("LOOKS LIKE NORMALIZATION WRONG")
	    print("DATE", fcdate)
	    (U15, V15) = calc_u15_from_3hr([file_fcst], vfreq=24, newnorm=False)
	if ( np.max(U15) < 0.05 ):  
	    print("LOOKS LIKE NORMALIZATION WRONG")
	    print("DATE", fcdate)
	    (U15, V15) = calc_u15_from_3hr([file_fcst], vfreq=24, newnorm=True)
	[U15F, V15F] = UU_ORCA_to_NE([U15, V15])
	
        # this actually removes non-ocean points
        (U15m, V15m), (LONM, LATM) = msk_flds([U15F, V15F], [LONN, LATN], mask0)
        # now grid to lat long
        (U15f, V15f), (LONF, LATF) = put_flds_latlon([U15m, V15m], [LONM, LATM], ddeg=0.2, method='2sweeplinear')

	ENSM.append=[ [[U15F, V15F]], [[U15f, V15f]] ]   ## EXTRA BRACKETS ARE HISTORICAL BAGGAGE.
	
      FCST.APPEND[ENSM]   ## INDICES:  nens=21, nfcst=16, ngrid=2, ideps=1, nflds=2
	
    for iobs in range(nobss_loop):

        [LONP, LATP] = LONO[iobs], LATO[iobs]
	IJPTS = find_nearest_points((LONP, LATP), [(LONN, LATN), (LONF, LATF)])
	IJPT, IJPF, IJPI, IJPG = IJPTS
	print('IJPTS', IJPTS)

        for ie, ENSM in enumerate(FCST): 
          for ld, LEAD in enumerate(ENSM):
	    for ig, GRID in enumerate(LEAD):
	      IJPT = IJPTS[ig]
	      fcstv = FCSTV_LIST[ig].copy()
	      for kk, klev in enumerate(GRID):
	        fcstv[iobs, 0, ld, kk] = Upi
	        fcstv[iobs, 1, ld, kk] = Vpi
	
	      for kk, KLEV in enumerate(GRID):
	        fcstl = []
	        for ifld, FLD in enumerate(KLEV):
                    FLP=FLD[IJPT]
		    initl.append(FLP)
		Upi, Vpi = fcstl   ## These are the best values for level kk
                fcstv[iobs, 0, ld, ie, kk] = Upi
                fcstv[iobs, 1, ld, ie, kk] = Vpi
              FCSTV_LIST[ig] = fcstv.copy()

    for	ig in range(FCSTV_LIST):
        fcstv = FCSTV_LIST[ig].copy()
	obsfile = OBSFILE_LIST[ig]
	oesfile = OESFILE_LIST[ig]
	rc = write_model_obsfile_ensemble(oesfile, obsfile, fcstv, clobber=True)

     return
	
def write_model_obsfile_ensembles(obsfile, tplfile, fcstv, clobber=True):
    if ( fcstv.ndim == 5 ):
      nobs, nvars, nfcsts, nenss, ndeps = np.shape(fcstv)
    # ncks template file to new obsfile -- but remove forecast -- which we will need to re-create.
    # removing all three numfcsts length variables SHOULD remove numfcst dimension as well.  [THIS MAY FAIL?]
    # Using same module, we should be able to recreate persistence with a longer time line too.
    if ( clobber ):
      rc=subprocess.call(['ncks','-O','-x','-v','forecast,persistence,negative_persistence',tplfile, obsfile]) 

    obsset = netCDF4.Dataset(obsfile,mode='r+')
    nfcsts_d = obsset.createDimension('numfcsts', nfcsts)
    nens_d = obsset.createDimension('numens', nenss)
    fcstv_var = obsset.createVariable('forecast', np.float32,
                                      ('numobs', 'numvars', 'numfcsts', 'numens', 'numdeps'),
				      fill_value=obsset['best_estimate']._FillValue)
    
    fcstv_var.units = "m s-1"
    fcstv_var.long_name = "Model forecast counterpart of obs. value" ;
    fcstv_var[:] = fcstv
    obsset.close()
    return

