import sys
import subprocess
import os
import numpy as np
import datetime
        
def get_archive(destdir, branch, date, fcst_hour, ensnum=None, ext_string=None, execute=False, filter=[], filter_file=None):
    # DO NOT EXECUTE RARC UNTIL BUGS WORKED OUT.
    #  NOTE:  Date can be a single date, or range of dates
    add_filter=False
    use_filter=False
    if ( len(filter) > 0 ): add_filter=True
    if ( isinstance(filter_file, str) ): use_filter=True

    rc=subprocess.call(['mkdir','-p', destdir])
    if ( isinstance(date, list) ):
        date_string = date[0].strftime('%Y,%m,%d')
        date_string = date_string+','+date[-1].strftime('%Y,%m,%d')
        date_hour_string = date[0].strftime('%H')
    else:
        date_string = date.strftime('%Y,%m,%d')
        date_hour_string = date.strftime('%H')
    if ( isinstance(ext_string, type(None) ) ):
        if ( isinstance(ensnum, type(None) ) ): 
            ext_string = str(fcst_hour).zfill(3)
            fil_string = ext_string
        elif ( isinstance(ensnum, list) or isinstance(ensnum, np.ndarray) ):
            fil_string = str(fcst_hour).zfill(3)+'_'+'lst'
            for ie,ensmem in enumerate(ensnum):
                if ( ie == 0 ):
                    ext_string = str(fcst_hour).zfill(3)+'_'+str(ensmem).zfill(3)
                else:
                    ext_string = ext_string+','+str(fcst_hour).zfill(3)+'_'+str(ensmem).zfill(3)
        else:
            ext_string = str(fcst_hour).zfill(3)+'_'+str(ensnum).zfill(3)
            fil_string = ext_string
    else:
        fil_string = ext_string
    if ( isinstance(date,list) ):
        date_alll_string = date[0].strftime('%Y%m%d%H')+'_'+fil_string
    else:
        date_alll_string = date.strftime('%Y%m%d%H')+'_'+fil_string
    fetch_file=destdir+'/'+date_alll_string+'.rarc'
    log_file=destdir+'/'+date_alll_string+'.log'
    file_fetch = open(fetch_file,"w+")
    if ( not isinstance(filter_file, str) ):
        filter_file=destdir+'/'+date_alll_string+'.flt.ksh'
    file_fetch.write('target = '+destdir+'\n')
    if ( add_filter or use_filter ):
        file_fetch.write('filter = '+filter_file+'\n')
    else:
        file_fetch.write('filter = copy'+'\n')
    file_fetch.write('postprocess = nopost'+'\n')
    file_fetch.write('date = '+date_string+'\n')
    file_fetch.write('branche = '+branch+'\n')
    file_fetch.write('ext = '+ext_string+'\n')
    file_fetch.write('heure = '+date_hour_string+'\n')
    file_fetch.write('priority = online'+'\n')
    file_fetch.write('inc = 1'+'\n')
   
    file_fetch.close()

    if ( add_filter ):
        file_filter = open(filter_file, "w+")
        file_filter.write("#!/bin/ksh"+"\n")
        file_filter.write("#-*-ksh-*-"+"\n")
        file_filter.write("editfst -s $1 -d $2 -l << EOF"+"\n")
        file_filter.write("desire(-1,['^^'])"+"\n")
        file_filter.write("desire(-1,['>>'])"+"\n")
        file_filter.write("desire(-1,['^>'])"+"\n")
        for field in filter:
            file_filter.write("desire(-1,['"+field+"'])"+"\n")
        file_filter.write("EOF"+"\n")
        file_filter.close()
        rc=subprocess.call(['chmod','755',filter_file])

    if ( execute ):
        #rc=subprocess.call(['rarc','-i',fetch_file, '-m', log_file])
        command='/home/dpe000/GEOPS/jobscripts/rarc_retrieval.sh'
        total_command=[command, '-i='+fetch_file,'-l='+log_file]
        print(total_command)
        rc=subprocess.call(total_command)
    else:
        rc=0
    return rc

def get_archive_leads(destdir, branch, date, fcst_list, ensnum=None, execute=False, filter=[], filter_file=None):
    # DO NOT EXECUTE RARC UNTIL BUGS WORKED OUT.
    
    add_filter=False
    use_filter=False
    if ( len(filter) > 0 ): add_filter=True
    if ( isinstance(filter_file, str) ): use_filter=True

    rc=subprocess.call(['mkdir','-p', destdir])
    if ( isinstance(date, datetime.datetime) ):
        date_string = date.strftime('%Y,%m,%d')
        date_hour_string = date.strftime('%H')
    if ( isinstance(date, list) ):
        date_string = date[0].strftime('%Y,%m,%d')+','+date[1].strftime('%Y,%m,%d')
        date_hour_strin0 = date[0].strftime('%H')
        date_hour_strin1 = date[1].strftime('%H')
        if ( date_hour_strin1 == date_hour_strin0 ): 
            date_hour_string = date_hour_strin0
        else:
            date_hour_string = date_hour_strin0+','+date_hour_strin1
    ie=0
    for fcst_hour in fcst_list:    
        if ( isinstance(ensnum, type(None) ) ): 
            if ( ie == 0 ):
                ext_string = str(fcst_hour).zfill(3)
            else:
                ext_string = ext_string+','+str(fcst_hour).zfill(3)
            ie=ie+1
        elif ( isinstance(ensnum, list) or isinstance(ensnum, np.ndarray) ):
            for ensmem in ensnum:
                if ( ie == 0 ):
                    ext_string = str(fcst_hour).zfill(3)+'_'+str(ensmem).zfill(3)
                else:
                    ext_string = ext_string+','+str(fcst_hour).zfill(3)+'_'+str(ensmem).zfill(3)
                ie=ie+1
        else:
            if ( ie == 0 ):
                ext_string = str(fcst_hour).zfill(3)+'_'+str(ensnum).zfill(3)
            else:
                ext_string = ext_string+','+str(fcst_hour).zfill(3)+'_'+str(ensnum).zfill(3)
            ie=ie+1

    if ( isinstance(date, datetime.datetime) ):
        date_alll_string = date.strftime('%Y%m%d%H')
    if ( isinstance(date, list) ):
        date_alll_string = date[0].strftime('%Y%m%d%H')+"_"+date[1].strftime('%Y%m%d%H')
    fetch_file=destdir+'/'+date_alll_string+'.rarc'
    log_file=destdir+'/'+date_alll_string+'.log'
    if ( not isinstance(filter_file, str) ):
        filter_file=destdir+'/'+date_alll_string+'.flt.ksh'
    file_fetch = open(fetch_file,"w+")

    file_fetch.write('target = '+destdir+'\n')
    if ( add_filter or use_filter ):
        file_fetch.write('filter = '+filter_file+'\n')
    else:
        file_fetch.write('filter = copy'+'\n')
    file_fetch.write('postprocess = nopost'+'\n')
    file_fetch.write('date = '+date_string+'\n')
    file_fetch.write('branche = '+branch+'\n')
    file_fetch.write('ext = '+ext_string+'\n')
    file_fetch.write('heure = '+date_hour_string+'\n')
    file_fetch.write('priority = online'+'\n')
    file_fetch.write('inc = 1'+'\n')
   
    file_fetch.close()

    if ( add_filter ):
        file_filter = open(filter_file, "w+")
        file_filter.write("#!/bin/ksh"+"\n")
        file_filter.write("#-*-ksh-*-"+"\n")
        file_filter.write("editfst -s $1 -d $2 -l << EOF"+"\n")
        file_filter.write("desire(-1,['^^'])"+"\n")
        file_filter.write("desire(-1,['>>'])"+"\n")
        file_filter.write("desire(-1,['^>'])"+"\n")
        for field in filter:
            file_filter.write("desire(-1,['"+field+"'])"+"\n")
        file_filter.write("EOF"+"\n")
        file_filter.close()
        rc=subprocess.call(['chmod','755',filter_file])

    if ( execute ):
        #rc=subprocess.call(['rarc','-i',fetch_file, '-m', log_file])
        rc=subprocess.call(['cat', filter_file])
        rc=subprocess.call(['cat', fetch_file])
        command='/home/dpe000/GEOPS/jobscripts/rarc_retrieval.sh'
        total_command=[command, '-i='+fetch_file,'-l='+log_file]
        print(total_command)
        rc=subprocess.call(total_command)
    else:
        rc=0
    return rc

def make_filter_file(filter_file, filter):
    file_filter = open(filter_file, "w+")
    file_filter.write("#!/bin/ksh"+"\n")
    file_filter.write("#-*-ksh-*-"+"\n")
    file_filter.write("editfst -s $1 -d $2 -l << EOF"+"\n")
    file_filter.write("desire(-1,['^^'])"+"\n")
    file_filter.write("desire(-1,['>>'])"+"\n")
    for field in filter:
        file_filter.write("desire(-1,['"+field+"'])"+"\n")
    file_filter.write("EOF"+"\n")
    file_filter.close()
    rc=subprocess.call(['chmod','755',filter_file])
    return rc
    
def call_filter(src_file, filter_file, filter, dest_file='null'): 
    devnull = open(os.devnull, 'w')
    if ( dest_file == 'null' ):  dest_file=src_file+'.filtered' 
    rc=make_filter_file(filter_file, filter)
    if ( rc == 0 ):
        command=[filter_file, src_file, dest_file]
        rc=subprocess.call(command, stdout=devnull, stderr=devnull)
    return rc    
    
  
