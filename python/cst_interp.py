import subprocess


grid_dir='/fs/site6/eccc/mrd/rpnenv/dpe000/GDPS/gdps_Class4/'
script_dir='/fs/homeu2/eccc/mrd/ords/rpnenv/dpe000/Class4_Currents/jobscripts/'
template_dir='/fs/homeu2/eccc/mrd/ords/rpnenv/dpe000/Class4_Currents/templates/'
grid_ref='dest_grid.std.25'

def cst_interpolation(source_file, destination_file, ref_grid=grid_dir+grid_ref):
    script=script_dir+'perform_interpolation_25.sh'
    rc = subprocess.call(['bash', script, '-fs='+source_file, '-fd='+destination_file, '-fr='+ref_grid])
    return rc

def filter_standard_file(filter, source_file, destination_file):
    script=script_dir+'perform_filter.sh'
    clobber=''
    if ( destination_file == source_file ): 
        clobber="--clobber"
        destination_file = source_file+'.tmp'
    rc=subprocess.call(['bash' , script, '-f='+filter, '-i='+source_file, '-o='+destination_file, clobber])
    return rc
    
