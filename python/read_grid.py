# NOTE netCDF4 needs to be imported before any other modules at top of python script
import sys
#sys.path.insert(0,'/fs/ssm/eccc/cmd/cmdm/satellite/master/sat-pyenv_1.2.2_ubuntu-14.04-amd64-64/envs/sat-pyenv_1.2.2/lib/python2.7/site-packages')
import netCDF4
#import netcdftime
#del sys.path[0]
import numpy as np
import datetime
import find_hall

hall=find_hall.find_hall()

def ubuntu(hall=hall):
    ubu='env_ubuntu-14.04-amd64-64'
    if ( (hall == 'hall3') or (hall == 'hall4') ): ubu='env_ubuntu-18.04-skylake-64'
    if ( (hall == 'hall6') or (hall == 'hall5') ): ubu='env_rhel-8-icelake-64'
    return ubu
      
def read_coord(hall=hall, grid='T'):
    site='/space/'+hall+'/sitestore'
    ubu=ubuntu(hall=hall)
    dire=site+'/eccc/mrd/rpnenv/socn000/'+ubu+'/datafiles/constants/oce/repository/master/CONCEPTS/orca025/grids'
    file=dire+'/coordinates_ORCA025_LIM.nc'
    #print('Coord file = '+file)
    dataset = netCDF4.Dataset(file) 
    if ( grid == 'T' ):
        e1t=np.transpose( np.squeeze(dataset.variables['e1t'][:]) )
        e2t=np.transpose( np.squeeze(dataset.variables['e2t'][:]) )
        nav_lon=np.transpose( np.squeeze(dataset.variables['nav_lon'][:]) )
        nav_lat=np.transpose( np.squeeze(dataset.variables['nav_lat'][:]) )
    if ( grid == 'U' ):
        e1t=np.transpose( np.squeeze(dataset.variables['e1u'][:]) )
        e2t=np.transpose( np.squeeze(dataset.variables['e2u'][:]) )
        nav_lon=np.transpose( np.squeeze(dataset.variables['glamu'][:]) )
        nav_lat=np.transpose( np.squeeze(dataset.variables['gphiu'][:]) )
    if ( grid == 'V' ):
        e1t=np.transpose( np.squeeze(dataset.variables['e1v'][:]) )
        e2t=np.transpose( np.squeeze(dataset.variables['e2v'][:]) )
        nav_lon=np.transpose( np.squeeze(dataset.variables['glamv'][:]) )
        nav_lat=np.transpose( np.squeeze(dataset.variables['gphiv'][:]) )
    if ( grid == 'F' ):
        e1t=np.transpose( np.squeeze(dataset.variables['e1f'][:]) )
        e2t=np.transpose( np.squeeze(dataset.variables['e2f'][:]) )
        nav_lon=np.transpose( np.squeeze(dataset.variables['glamf'][:]) )
        nav_lat=np.transpose( np.squeeze(dataset.variables['gphif'][:]) )
    return nav_lon, nav_lat, e1t*e2t
 
def read_mask(var='tmask', hall=hall):
    site='/space/'+hall+'/sitestore'
    ubu=ubuntu(hall=hall)
    dire=site+'/eccc/mrd/rpnenv/socn000/'+ubu+'/datafiles/constants/oce/repository/master/CONCEPTS/orca025/grids'
    file=dire+'/mask_float.nc'
    #print('Mask file = '+ file)
    dataset = netCDF4.Dataset(file) 
    mask=np.squeeze(dataset.variables[var][:])
    mask=np.moveaxis(mask, [0, 1, 2], [0, 2, 1])
    return mask

def vars_in_mesh(hall=hall):
    site='/space/'+hall+'/sitestore'
    dire=site+'/eccc/mrd/rpnenv/dpe000/NEMO_MESH_MASK'
    file=dire+'/mesh_mask.nc'
    #print('Mesh Mask file = ', file)
    dataset = netCDF4.Dataset(file) 
    variables = dataset.variables.keys()
    return variables
    
def read_mesh_var(var, hall=hall):
    site='/space/'+hall+'/sitestore'
    dire=site+'/eccc/mrd/rpnenv/dpe000/NEMO_MESH_MASK'
    file=dire+'/mesh_mask.nc'
    #print('Mesh Mask file = '+file)
    dataset = netCDF4.Dataset(file) 
    var=np.squeeze( dataset.variables[var][:] )
    # transpose last x and y positions (for adherence with standard file)
    var=np.transpose(var, list(range(var.ndim-2))+[-1, -2])
    return var
 
def read_coord_mesh(hall=hall,grid='T'):
    site='/space/'+hall+'/sitestore'
    dire=site+'/eccc/mrd/rpnenv/dpe000/NEMO_MESH_MASK'
    file=dire+'/mesh_mask.nc'
    #print('Mesh Mask file = ', file)
    dataset = netCDF4.Dataset(file) 
    if ( grid == 'T' ):
        e1t=np.transpose( np.squeeze(dataset.variables['e1t'][:]) )
        e2t=np.transpose( np.squeeze(dataset.variables['e2t'][:]) )
        nav_lon=np.transpose( np.squeeze(dataset.variables['nav_lon'][:]) )
        nav_lat=np.transpose( np.squeeze(dataset.variables['nav_lat'][:]) )
    if ( grid == 'U' ):
        e1t=np.transpose( np.squeeze(dataset.variables['e1u'][:]) )
        e2t=np.transpose( np.squeeze(dataset.variables['e2u'][:]) )
        nav_lon=np.transpose( np.squeeze(dataset.variables['glamu'][:]) )
        nav_lat=np.transpose( np.squeeze(dataset.variables['gphiu'][:]) )
    if ( grid == 'V' ):
        e1t=np.transpose( np.squeeze(dataset.variables['e1v'][:]) )
        e2t=np.transpose( np.squeeze(dataset.variables['e2v'][:]) )
        nav_lon=np.transpose( np.squeeze(dataset.variables['glamv'][:]) )
        nav_lat=np.transpose( np.squeeze(dataset.variables['gphiv'][:]) )
    if ( grid == 'F' ):
        e1t=np.transpose( np.squeeze(dataset.variables['e1f'][:]) )
        e2t=np.transpose( np.squeeze(dataset.variables['e2f'][:]) )
        nav_lon=np.transpose( np.squeeze(dataset.variables['glamf'][:]) )
        nav_lat=np.transpose( np.squeeze(dataset.variables['gphif'][:]) )
    return nav_lon, nav_lat, e1t*e2t

def read_e3t_mesh(hall=hall, var='e3t_0'):
    e3t=read_mesh_var(var, hall)
    return e3t

def read_tmask_mesh(hall=hall, var='tmask'):
    tmask=read_mesh_var(var, hall)
    return tmask

def bottom_depth_from_e3t(e3t, mask):
    Depth=np.sum(e3t*mask, axis=0)
    return Depth
    
def bottom_depth_mesh(hall=hall):
    tmask=read_mesh_var('tmask', hall)
    e3t=read_mesh_var('e3t_0', hall)
    Depth=np.sum(e3t*tmask, axis=0)
    return Depth

def read_netcdf(file, varname, latname='nav_lat', lonname='nav_lon', depname='deptht', timname='time_counter',tunits='hours', DAY0=datetime.datetime(1950, 1, 1, 0) ):
    dataset = netCDF4.Dataset(file) 
    # HANDY TO BE ABLE TO READ IN ONLY THE COORDS
    if ( isinstance(varname, type(None) ) ):
        FLD=[]
    else:
        FLD=dataset.variables[varname][:]
    if ( isinstance(lonname, type(None) ) ):
       LON=[]
    else:
       LON=dataset.variables[lonname][:]
    if ( isinstance(latname, type(None) ) ):
        LAT=[]
    else:
        LAT=dataset.variables[latname][:]
    if ( isinstance(depname, type(None) ) ):
        DEP=[]
    else:
        DEP=dataset.variables[depname][:]
    TIM = []
    if ( timname != None ):
        TIM=dataset.variables[timname][:]
        t_unit = dataset.variables[timname].units
        print('t_unit', t_unit)  ## FIGURE OUT HOW TO USE THIS
        try :
            t_cal = dataset.variables[timname].calendar
        except AttributeError : # Attribute doesn't exist
            t_cal = u"gregorian" # or standard

    DATE = []
    #DAY0 = datetime.datetime(1950, 1, 1, 0)
    if ( tunits == 'hours' ):
        tmult=3600
    elif ( tunits == 'seconds' ):
        tmult=1
    elif ( tunits == 'days' ):
        tmult=86400
    else:
        tmult=1
    for tim in TIM: DATE.append( DAY0+ datetime.timedelta(seconds=tim*tmult))
    dataset.close()
    return FLD, LON, LAT, DEP, DATE

def read_angle(hall=hall, grid='T'):
    site='/space/'+hall+'/sitestore'
    ubun=ubuntu(hall=hall)
    dire=site+'/eccc/mrd/rpnenv/socn000/'+ubun+'/datafiles/constants/oce/repository/master/CONCEPTS/orca025/grids'
    file=dire+'/coordinates_ORCA025_LIM.nc'
    print(file)
    dataset = netCDF4.Dataset(file) 
    angle = []
    return angle

def read_bathymetry(hall=hall):
    site='/space/'+hall+'/sitestore'
    ubu=ubuntu(hall=hall)
    dire=site+'/eccc/mrd/rpnenv/socn000/'+ubu+'/datafiles/constants/oce/repository/master/CONCEPTS/orca025/grids'
    file=dire+'/bathy_ORCA025_LIM.nc'
    #print('Bathy file = ', file)
    dataset = netCDF4.Dataset(file) 
    bathy=dataset['Bathymetry'][:]
    bathy=np.transpose( bathy )
    return bathy

