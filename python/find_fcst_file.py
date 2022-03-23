import sys
import traceback
import os

import find_hall
import get_archive
import tmpcp


main_host=find_hall.get_main_host()
site=find_hall.get_site()
datadir='/fs/'+site+'/eccc/mrd/rpnenv/dpe000/'
scratch=os.environ['SCRATCH']
workdir=os.environ['TMPDIR']
archive=datadir+'/tmpdir/'

print 'WORKDIR = '+workdir
print 'ARCHIVE = '+archive

d_geps_20190128='/home/smco500/.suites/geps_20190128/e1/hub/'+main_host+'/'
d_geps_20190703='/home/smco500/.suites/geps_20190703/e1/hub/'+main_host+'/'
d_geps_20191231='/home/smco500/.suites/geps_20191231/e1/hub/'+main_host+'/'
d_geps_20210312='/home/smco501/.suites/geps_20210312/e1/hub/'+main_host+'/'

a_geps=archive+'/PS'
t_geps=workdir+'/PS'
d_geps_op='/home/smco500/.suites/geps/forecast/hub/'+main_host
a_geps_op=archive+'/OP'
t_geps_op=workdir+'/OP'
d_gdps='/home/smco501/.suites/gdps/g1/hub/'+main_host+'/gridpt/prog/oce'
a_gdps=archive+'/PSD'
t_gdps=workdir+'/PSD'
d_gdps_op='/home/smco500/.suites/gdps/g1/hub/'+main_host+'/gridpt/prog/oce'
a_gdps_op=archive+'/OPD'
t_gdps_op=workdir+'/OPD'
m_geps_20191231='/home/smco500/.suites/geps_20191231/e1/hub/eccc-ppp4'
a_geps_mg=archive+'/MG'
aa_geps_mg=archive+'/MGA'
a_geps_om=archive+'/OM'
aa_geps_om=archive+'/OMA'


def systems_ensembles(system):
    nens = 21
    if ( system == 'PS' ):
        nens = 21
    if ( system == 'MG' ):
        nens = 21
    if ( system == 'PSD' ):
        nens = 1
    if ( system == 'OPD' ): 
        nens = 1
    if ( system == 'OP' ):
        nens = 21
    return nens

def systems_leads(system):
    leads = 32 
    if ( system == 'PS' ):
	leads = 32
    if ( system == 'MG' ):
	leads = 32
    if ( system == 'PSD' ):
	leads = 10
    if ( system == 'OPD' ): 
	leads = 10
    if ( system == 'OP' ):
	leads = 32
    return leads


def find_fcst_file(system, date, leadhrs, ie, src='ocn', execute=False):
 
    this_date_str = date.strftime('%Y%m%d%H')
    lead_hour_str = str(leadhrs).zfill(3)
    ens_str = str(ie)
    en3_str = str(ie).zfill(3)
    weekday = date.weekday()

    file_fcst='Null'

    if ( system == 'MG' ):
        if ( src == 'ocn' ):
	    for e_geps in [ 'oce', 'ens.glboce' ]:
	        file_fcst=d_geps_20191231+'/gridpt/prog/'+e_geps+'/'+this_date_str+'_'+lead_hour_str+'_'+en3_str
		if ( os.path.isfile(file_fcst) ): 
		    break
        if ( src == 'atm' ):
	    file_fcst=m_geps_20191231+'/gridpt/prog/ens.glbmodel/'+this_date_str+'_'+lead_hour_str+'_'+en3_str
	    arc_fcst=a_geps_mg+'/A'
	    file_fcst=tmpcp.tmpcp(file_fcst, tmpdir=arc_fcst)
            branch = 'parallel.ensemble.ens.glbmodel'
	    if ( not os.path.isfile(file_fcst) ):
                file_fcst=arc_fcst+'/'+branch+'/'+this_date_str+'_'+lead_hour_str+'_'+en3_str
		# WHEN OCN DATA ARCHIVED.  REORDER.
	        if ( not os.path.isfile(file_fcst) ):
		    if ( ( weekday == 3 ) or ( leadhrs < 385 ) ): # ONLY LONGER LEADS ON THURSDAY
		        if ( ie == 0 ): # Now gets full ensemble if ens==0.
		            enslst = range(systems_ensembles(system))
			    rc=get_archive.get_archive(arc_fcst, branch=branch, date=date, fcst_hour=leadhrs, ensnum=enslst, execute=execute)
		        else:
		            rc=get_archive.get_archive(arc_fcst, branch=branch, date=date, fcst_hour=leadhrs, ensnum=ie, execute=execute)

    if ( system == 'PS_old' ):
        if ( src == 'ocn' ): 
	    for d_geps in [d_geps_20190703, d_geps_20190128]:
	        for e_geps in [ 'oce', 'ens.glboce' ]:
	            file_fcst=d_geps+'/gridpt/prog/'+e_geps+'/'+this_date_str+'_'+lead_hour_str+'_'+en3_str
		    if ( os.path.isfile(file_fcst) ): 
		        break
        if ( src == 'atm' ):
	    for d_geps in [d_geps_20190703, d_geps_20190128]:
	        file_fcst=d_geps+'/gridpt/prog/ens.glbmodel/'+this_date_str+'_'+lead_hour_str+'_'+en3_str
		if ( os.path.isfile(file_fcst) ): 
		    break
	    branch='parallel.ensemble.ens.glbmodel'
	    arc_fcst=a_geps+'/A'
	    if ( not os.path.isfile(file_fcst) ):
                file_fcst=arc_fcst+'/'+branch+'/'+this_date_str+'_'+lead_hour_str+'_'+en3_str
		# WHEN OCN DATA ARCHIVED.  REORDER.
	        if ( not os.path.isfile(file_fcst) ):
		    if ( ( weekday == 3 ) or ( leadhrs < 385 ) ): # ONLY LONGER LEADS ON THURSDAY
		        if ( ie == 0 ): # Now gets full ensemble if ens==0.
		            enslst = range(systems_ensembles(system))
			    rc=get_archive.get_archive(arc_fcst, branch=branch, date=date, fcst_hour=leadhrs, ensnum=enslst, execute=execute)
		        else:
		            rc=get_archive.get_archive(arc_fcst, branch=branch, date=date, fcst_hour=leadhrs, ensnum=ie, execute=execute)
    if ( system == 'PS' ):
        if ( src == 'ocn' ): 
	    for d_geps in [d_geps_20210312]:
	        file_fcst=d_geps+'/gridpt/prog/ens.glboce/'+this_date_str+'_'+lead_hour_str+'_'+en3_str
		if ( os.path.isfile(file_fcst) ): 
		    break
            branch = 'parallel.ensemble.prog.ens.glboce'
	    arc_fcst=a_geps+'/O'
        if ( src == 'atm' ):
	    for d_geps in [d_geps_20210312]:
	        file_fcst=d_geps+'/gridpt/prog/ens.glbmodel/'+this_date_str+'_'+lead_hour_str+'_'+en3_str
		if ( os.path.isfile(file_fcst) ): 
		    break
	    branch='parallel.ensemble.ens.glbmodel'
	    arc_fcst=a_geps+'/A'
        if ( not os.path.isfile(file_fcst) ):
            file_fcst=arc_fcst+'/'+branch+'/'+this_date_str+'_'+lead_hour_str+'_'+en3_str
            # WHEN OCN DATA ARCHIVED.  REORDER.
	    if ( not os.path.isfile(file_fcst) ):
		if ( ( weekday == 3 ) or ( leadhrs < 385 ) ): # ONLY LONGER LEADS ON THURSDAY
		    if ( ie == 0 ): # Now gets full ensemble if ens==0.
		        enslst = range(systems_ensembles(system))
			rc=get_archive.get_archive(arc_fcst, branch=branch, date=date, fcst_hour=leadhrs, ensnum=enslst, execute=execute)
		    else:
		        rc=get_archive.get_archive(arc_fcst, branch=branch, date=date, fcst_hour=leadhrs, ensnum=ie, execute=execute)

    if ( system == 'OP' ):
        if ( src == 'ocn' ):
	    file_fcst=d_geps_op+'/gridpt/prog/ens.glboce/'+this_date_str+'_'+lead_hour_str+'_'+en3_str
            branch = 'operation.ensemble.prog.ens.glboce'
	    arc_fcst=a_geps_op+'/O'
        if ( src == 'atm' ):
	    file_fcst=d_geps_op+'/gridpt/prog/ens.glbmodel/'+this_date_str+'_'+lead_hour_str+'_'+en3_str
            branch = 'operation.ensemble.prog.ens.glbmodel'
	    arc_fcst=a_geps_op+'/A'
        if ( not os.path.isfile(file_fcst) ):
            file_fcst=arc_fcst+'/'+branch+'/'+this_date_str+'_'+lead_hour_str+'_'+en3_str
            # WHEN OCN DATA ARCHIVED.  REORDER.
	    if ( not os.path.isfile(file_fcst) ):
		if ( ( weekday == 3 ) or ( leadhrs < 385 ) ): # ONLY LONGER LEADS ON THURSDAY
		    if ( ie == 0 ): # Now gets full ensemble if ens==0.
		        enslst = range(systems_ensembles(system))
			rc=get_archive.get_archive(arc_fcst, branch=branch, date=date, fcst_hour=leadhrs, ensnum=enslst, execute=execute)
		    else:
		        rc=get_archive.get_archive(arc_fcst, branch=branch, date=date, fcst_hour=leadhrs, ensnum=ie, execute=execute)
	    
    if ( ( system == 'PSD' ) or ( system == 'OPD' ) ):
        if ( system == 'PSD' ): 
	    suite='smco501'
	    branch='parallel.forecasts.giops.prog.glboce'
	    dir_fcst=d_gdps
	    arc_fcst=a_gdps
        if ( system == 'OPD' ): 
	    suite='smco500'
	    branch='operation.forecasts.giops.prog.glboce'
	    dir_fcst=d_gdps_op
	    arc_fcst=a_gdps_op

        #dir_fcst='/home/'+suite+'/.suites/gdps/forecast/hub/eccc-ppp1/gridpt/prog/oce'
        file_fcst=dir_fcst+'/'+this_date_str+'_'+lead_hour_str
	if ( not os.path.isfile(file_fcst) ): 
            file_fcst=arc_fcst+'/'+branch+'/'+this_date_str+'_'+lead_hour_str
	    if ( not os.path.isfile(file_fcst) ):
	        # GET FROM ARCHIVE WHEN MODULE WRITTEN.
	        rc=get_archive.get_archive(arc_fcst, branch=branch, date=date, fcst_hour=leadhrs, execute=execute)
	        if ( ( not os.path.isfile(file_fcst) ) and execute ):
	            raise Exception('File '+file_fcst+' does not exist' )

    if ( not os.path.isfile(file_fcst) and execute ):
        print "find_fcst_file", system, date, leadhrs, ie, src, execute
	raise Exception('File '+file_fcst+' does not exist' )

    return file_fcst

def find_anal_file(date, system='gd', anal=True, var='T',execute=False):
    if ( not isinstance(anal, bool) ):
        if ( isinstance(anal, str) ):
	    if ( ( anal == 'ANAL' ) or ( anal == 'anal' ) ): anal=True
	    if ( ( anal == 'TRIAL' ) or ( anal == 'trial' ) ): anal=False
    if ( not isinstance(anal, bool) ):
        print('Do Not Understand anal = ', anal)
	return 'Null', 'Null'
         
    #if ( system == 'pd' ): 
    #    print('Not configured ', system)
    #	return 'Null', 'Null'
    this_date_str1 = date.strftime('%Y%m%d')
    this_date_str2 = this_date_str1+'00'
    weekday = date.weekday()

    sample_file='eccc-ppp4:/home/smco501/.suites/giops/gu/hub/eccc-ppp4/SAM2/20190515/DIA/ORCA025-CMC-ANAL_1d_grid_U_2019051500.nc.gz'
    dir_psan='eccc-ppp4:/home/smco501/.suites/giops/'+'gu'+'/hub/eccc-ppp4/SAM2/'+this_date_str1+'/DIA'
    dir_opan='eccc-ppp4:/home/smco500/.suites/giops/'+system+'/hub/eccc-ppp4/SAM2/'+this_date_str1+'/DIA'
    #dir_mgan='eccc-ppp4:/home/smco500/.suites/giops_20191231/'+system+'/hub/eccc-ppp4/SAM2/'+this_date_str1+'/DIA'
    #dir_mgbn='eccc-ppp4:/home/smco500/.suites/giops_20191231/'+system+'/hub/eccc-ppp4/SAM2/'+this_date_str1+'/DIA'

    if ( anal ): ANAL='ANAL'
    if ( not anal ): ANAL='TRIAL'
    ps_file_anal= dir_psan+'/ORCA025-CMC-'+ANAL+'_1d_grid_'+var+'_'+this_date_str2+'.nc.gz'
    ps_file_ssh= dir_psan+'/ORCA025-CMC-'+ANAL+'_1h_grid_T_2D_'+this_date_str2+'.nc.gz'
    op_file_anal= dir_opan+'/ORCA025-CMC-'+ANAL+'_1d_grid_'+var+'_'+this_date_str2+'.nc.gz'
    op_file_ssh= dir_opan+'/ORCA025-CMC-'+ANAL+'_1h_grid_T_2D_'+this_date_str2+'.nc.gz'
    #mg_file_anal= dir_mgan+'/ORCA025-CMC-'+ANAL+'_1d_grid_'+var+'_'+this_date_str2+'.nc.gz'
    #mg_file_ssh= dir_mgan+'/ORCA025-CMC-'+ANAL+'_1h_grid_T_2D_'+this_date_str2+'.nc.gz'

    if ( ( system == 'gu' ) or ( system == 'gd' ) ):
        file_anal = tmpcp.tmpcpgz(op_file_anal, tmpdir=datadir+'/tmpdir/'+system)
        file_ssh  = tmpcp.tmpcpgz(op_file_ssh,  tmpdir=datadir+'/tmpdir/'+system)

    if ( ( system == 'pu' ) or ( system == 'pd' ) ):
        file_anal = tmpcp.tmpcpgz(ps_file_anal, tmpdir=datadir+'/tmpdir/'+system)
        file_ssh  = tmpcp.tmpcpgz(ps_file_ssh,  tmpdir=datadir+'/tmpdir/'+system)
 
    if ( ( not os.path.isfile(file_anal) ) or ( not os.path.isfile(file_ssh) ) ):
        get_archived_analysis(date, system=system, ANAL=ANAL, execute=execute)
    
    if ( not os.path.isfile(file_anal) ):
	raise Exception('File '+file_anal+' does not exist' )
    if ( not os.path.isfile(file_ssh) ):
	raise Exception('File '+file_ssh+' does not exist' )

    return file_anal, file_ssh

def get_cmcarc_analysis(date, system='gu'):
    dir_opan='eccc-ppp4:/home/smco500/.suites/giops/'+system+'/hub/eccc-ppp4/SAM2/'
    dir_mgan='eccc-ppp4:/home/smco500/.suites/giops_20191231/'+system+'/hub/eccc-ppp4/SAM2/'
    this_date_str1 = date.strftime('%Y%m%d')
    this_date_str2 = this_date_str1+'00'
    file_cmcarc = dir_mgan+this_date_str2+'_.ca'
    return
    
def get_archived_analysis(date, system='gd', ANAL='ANAL', execute=False):
    date_string = date.strftime('%Y%m%d%H')
    arc_anal=archive+'/'+system
    branch='operation.diagnostics.giops.'+system
    if ( system == 'pu' ): branch='parallel.diagnostics.giops.gu'
    if ( system == 'pd' ): branch='parallel.diagnostics.giops.gd'
    rc=get_archive.get_archive(arc_anal, branch=branch, date=date, fcst_hour=0, ext_string='.ca', execute=execute)
    arc_file=arc_anal+'/'+branch+'/'+date_string+'_.ca'
    arguements=['-v', '-p', '-x']
    files=[]
    for grid in ['T', 'U', 'V']:
        file='ORCA025-CMC-'+ANAL+'_1d_grid_'+grid+'_'+date_string+'.nc.gz'
        add='DIA/'+file
	arguements.append(add)
	files.append(arc_anal+'/'+file)
    file='ORCA025-CMC-'+ANAL+'_1h_grid_T_2D_'+date_string+'.nc.gz'
    add='DIA/'+file
    arguements.append(add)
    files.append(arc_anal+'/'+file)	
    tmpcp.cdcmcarc(arc_file, arc_anal, arguements=arguements,rm=True)
    for file in files:
        nile=tmpcp.gunzip(file)
    return

def find_cmcsst_file(date, cmc_suite, execute=False):
    if ( cmc_suite == 'smco501' ): 
        branch = 'parallel.analyses.glb.netcdf.anal.glbsfc6'
        arc_fcst=a_gdps
    if ( cmc_suite == 'smco500' ): 
        branch = 'operation.analyses.glb.netcdf.anal.glbsfc6'
        arc_fcst=a_gdps_op
	
    this_date_str = date.strftime('%Y%m%d')+'00'  
    #file='/home/smco501/.suites/gdps/surface/final/hub/eccc-ppp1/netcdf/anal/sfc/2019061400_sst_latlon0.1x0.1.nc'
    basefile=this_date_str+'_sst_latlon0.1x0.1.nc'
    dir_anal='/home/'+cmc_suite+'/.suites/gdps/surface/final/hub/'+main_host+'/netcdf/anal/sfc/'
    file_anal=dir_anal+basefile
    if ( not os.path.isfile(file_anal) ):
        file_anal=arc_fcst+'/'+branch+'/'+basefile
	if ( not os.path.isfile(file_anal) ):
	    # GET FROM ARCHIVE WHEN MODULE WRITTEN.
	    rc=get_archive.get_archive(arc_fcst, branch=branch, date=date, fcst_hour=None, ext_string='sst_latlon0.1x0.1.nc', execute=execute)
	    if ( not os.path.isfile(file_anal) ):
	        raise Exception('File '+file_anal+' does not exist' )

    if ( not os.path.isfile(file_anal) ):
	raise Exception('File '+file_anal+' does not exist' )
    return file_anal
