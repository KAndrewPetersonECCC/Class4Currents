'Standard File Routines'
'Global ice-ocean prediction system (GIOPS)'
'Requires . r.load.dot eccc/mrd/rpn/MIG/ENV/py/2.7/rpnpy/2.0.4 '
'Requires ? . ssmuse-sh -d /fs/ssm/eccc/mrd/rpn/OCEAN/cstint-3.0.4 '
        
import datetime
import fnmatch
import logging
import os.path
import sys

import matplotlib.pyplot as plt
import scipy.constants as constants
import numpy as np

import rpnpy.librmn.all as rmn
import rpnpy.rpndate as rpndate

## REMOVE IF YOU WANT MORE verbose librmn diagnostic output.
## OTHER OPTIONS
## FSTOPI_MSG_DEBUG,   FSTOPI_MSG_INFO,  FSTOPI_MSG_WARNING,
## FSTOPI_MSG_ERROR,   FSTOPI_MSG_FATAL, FSTOPI_MSG_SYSTEM,
## FSTOPI_MSG_CATAST
rmn.fstopt(rmn.FSTOP_MSGLVL, rmn.FSTOPI_MSG_ERROR)

SHALLOWEST_DEPTH = .494025
GRID_SHAPE = 1442, 1021
file='/home/rmu001/data/ppp2/maestro/rm_hind_month_025/Sorties_3/nemo_0204/CICE/1996020500_000'
        
def file_query(file, nomvar=None):

    try:
        funit = rmn.fstopenall(file)
    except:
        sys.stderr.write("Problem opening the file: %s\n" % file)

    try:    
        if ( nomvar ): 
            keylist = rmn.fstinl(funit, nomvar=nomvar)
            print("Found %d records matching %s" % (len(keylist), nomvar) )
        else:
            keylist = rmn.fstinl(funit)
            print("Found %d records" % (len(keylist)) )

        # Get every record meta data
        VARS=[]
        VARF=[]
        for k in keylist:
            m = rmn.fstprm(k)
            print("%s (%d, %d, %d, %s)" % (m['nomvar'], m['ip1'], m['ip2'], m['ip3'], m['datev']))
            VARS.append(m['nomvar'])
            VARF.append(m)
    except:
        pass
    finally:
        # Close file even if an error occured above
        rmn.fstcloseall(funit)
    return VARS, VARF

def read_fstd_grid(fstd, **kwargs):
    '''Read lat/lon fields in GIOPS FSTD.

    Returns
    -------
    lon : ndarray
        longitudes
    lat : ndarray
        latitudes
    area : grid area = e1t*e2t

    '''

    logger = logging.getLogger(__name__)

    funit = rmn.fstopenall(fstd)

    try:
        lon = rmn.fstinl(funit, nomvar='>>')
        lon = rmn.fstluk(lon[0], dtype=np.float32)
        lon = lon['d']
    except:
        lon = None

    try:
        lat = rmn.fstinl(funit, nomvar='^^')
        lat = rmn.fstluk(lat[0], dtype=np.float32)
        lat = lat['d']
    except:
        lat = None

    e1t_rec = rmn.fstlir(funit, nomvar='e1t')
    e1t = e1t_rec['d']

    e2t_rec = rmn.fstlir(funit, nomvar='e2t')
    e2t = e2t_rec['d']

    rmn.fstcloseall(funit)

    return lon, lat, e1t*e2t

def read_fstd_gid(fstd, nomvar, **kwargs):
    '''Read GRID-ID from FSTD.

    Returns
    -------
    gid : grid_id

    '''

    logger = logging.getLogger(__name__)

    funit = rmn.fstopenall(fstd)
    var_recs = rmn.fstinl(funit, nomvar=nomvar)
    ## ONLY NEED ONE FIELD TO FIND GRID
    fldstr = rmn.fstluk(var_recs[0])

    var_gridU = rmn.readGrid(funit, fldstr)
    var_grid0 = var_gridU['subgrid'][0]
    var_grid1 = var_gridU['subgrid'][1]

    ## U GRID ID
    gid = var_gridU['id']

    gid0 = []
    gid1 = []
    if ( var_gridU['grtyp'] == 'U' ):
        ## Z SUBGRID ID's
        gid0 = var_grid0['id']
        gid1 = var_grid1['id']

    gid_tuple = (gid, gid0, gid1) 
    flds=[]
    for var_rec in var_recs:
        fldstr = rmn.fstluk(var_rec)
        fld = fldstr['d'].astype(np.float32)
        flds.append(fld)

    rmn.fstcloseall(funit)

    # NOTE:  If field occurs multiple times:  This may not be the field you want.  
    return gid_tuple, flds

def read_fstd_var(fstd, nomvar, **kwargs):
    '''Read field in GIOPS FSTD.

    Returns
    -------
    lon : ndarray
        longitudes
    lat : ndarray
        latitudes
    var : MaskedArray
        data

    '''

    logger = logging.getLogger(__name__)

    funit = rmn.fstopenall(fstd)

    try:
        lon = rmn.fstinl(funit, nomvar='>>')
        lon = rmn.fstluk(lon[0], dtype=np.float32)
        lon = lon['d']
    except:
        lon = None

    try:
        lat = rmn.fstinl(funit, nomvar='^^')
        lat = rmn.fstluk(lat[0], dtype=np.float32)
        lat = lat['d']
    except:
        lat = None

    var = rmn.fstlir(funit, nomvar=nomvar, **kwargs)
    datev = rpndate.RPNDate(var['datev']).toDateTime()
    logger.debug('datev: %s', datev.strftime('%c'))
    var = var['d']

    kwargs.pop('typvar', None)
    try:
        var_mask = rmn.fstlir(funit, nomvar=nomvar, typvar='@@', **kwargs)
        var_mask = var_mask['d']
    except:
        var_mask = np.logical_not(np.zeros(var.shape))

    rmn.fstcloseall(funit)

    var = np.ma.array(var, mask=np.logical_not(var_mask))

    return lon, lat, var

def read_fstd_multi(fstd, nomvar, **kwargs):
    '''Read field in GIOPS FSTD.

    Returns
    -------
    date : ndarray
        list of datetimes
    var : MaskedArray
        data

    '''

    logger = logging.getLogger(__name__)

    funit = rmn.fstopenall(fstd)

    record = rmn.fstinl(funit, nomvar=nomvar, **kwargs)
    var=[]
    date=[]
    for krec in record:
        varstr=rmn.fstluk(krec)
        datev = rpndate.RPNDate(varstr['datev']).toDateTime()
        date.append(datev)
        #print datev
        logger.debug('datev: %s', datev.strftime('%c'))
        var.append(varstr['d'])

    rmn.fstcloseall(funit)

    var = np.array(var)
    var.shape

    return date, var

def read_fstd_multi_lev(fstd, nomvar, vfreq=0, **kwargs):
    '''Read field in GIOPS FSTD.

    Returns
    -------
    level : ndarray
        list of levels
    var : MaskedArray
        data

    '''

    logger = logging.getLogger(__name__)

    funit = rmn.fstopenall(fstd)

    record = rmn.fstinl(funit, nomvar=nomvar, **kwargs)
    var=[]
    lev=[]
    for krec in record:
        varstr=rmn.fstluk(krec)
        ## exclude variables that don't have the required output frequency (default 24)
        if ( ( varstr['ip2'] - varstr['ip3' ] ) == vfreq or ( vfreq == 0 ) ):
            (rp1, rp2, rp3) = rmn.convertIPtoPK(varstr['ip1'], varstr['ip2'], varstr['ip3'])
            if ( ( rmn.FLOATIPtoList(rp1)[2] == 0) or (rmn.FLOATIPtoList(rp1)[2] == 4 ) ):  # MAKE SURE ON DEPTH NOT SIGMA LEVEL
               var.append(varstr['d'])
               lev.append(rmn.FLOATIPtoList(rp1)[0])

    rmn.fstcloseall(funit)

    var = np.array(var)
    #var.shape

    return lev, var

def find_date_in_dates(dates, this_date):
    for ii, list_date in enumerate(dates):
        if ( list_date == this_date ):
            jj=ii
            date=list_date
    return jj, date

def read_ip1(forecast, nomvar):
    'Read IP1 vertical levels in RPN Standard File'

    funit = rmn.fstopenall(forecast)
    rec = rmn.fstinl(funit, nomvar=nomvar, typvar='P@')
    ip1 = [None]*len(rec)
    for i, rec_i in enumerate(rec):
        ip1[i] = rmn.fstprm(rec_i)['ip1']
    rmn.fstcloseall(funit)
    return ip1

def shift_orca(lon, lat, var):
    'Shift ORCA grid for Basemap'

    # Code provided by Christopher Subich

    minlons = np.amin(lon, 0)
    shifts = np.zeros((np.shape(lon)[1]), dtype=int)
    for col in range(np.shape(lon)[1]):
        shifts[col] = plt.mlab.find(lon[:, col] == minlons[col])[0]

    def myroll(inarray):
        'Edit ORCA grid dataset for compatibility with Basemap.'
        outarray = inarray.copy()
        for col in range(np.shape(lon)[1]):
            outarray[:, col] = np.roll(inarray[:, col], -shifts[col], 0)
        return outarray

    rlon = myroll(lon)[:, :-1]
    rlat = myroll(lat)[:, :-1]
    rvar = myroll(var)[:, :-1]

    return rlon, rlat, rvar

def area_weights(scale_fstd_):
    """Calculate grid point area weights.

    See Also
    --------

    iris.analysis.cartography.area_weights : Similar function in
        SciTools Iris.

    Examples
    --------
    >>> from mpl_toolkits.basemap import Basemap
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> import os.path
    >>> import datetime
    >>> from scipy.constants import kilo
    >>> glboce = os.path.join(seq_exp_home_gridpt(), 'prog', 'glboce')
    >>> dateo = datetime.datetime(2014, 7, 3)
    >>> lon, lat, area_weights_ = area_weights(scale_fstd(glboce, dateo))
    >>> m = Basemap(projection='nplaea',boundinglat=45,lon_0=270,resolution='l')
    >>> m.drawcoastlines()
    >>> m.drawparallels(np.arange(-90.,120.,30.), labels=[1, 0, 0, 0])
    >>> m.drawmeridians(np.arange(0.,360.,60.), labels=[0, 0, 0, 1])
    >>> lon_s, lat_s, area_weights_s = shift_orca(lon, lat, area_weights_)
    >>> area_weights_s /= kilo**2
    >>> m.pcolormesh(lon_s, lat_s, area_weights_s, latlon=True, vmin=100, vmax=500)
    >>> cbar = plt.colorbar()
    >>> cbar.set_label(r'grid-point area (km$^2$)')

    """

    lon, lat, e1t = read_fstd_var(scale_fstd_, "e1t")
    e2t = read_fstd_var(scale_fstd_, "e2t")[2]

    return lon, lat, e1t*e2t

def read_latlon(fstd, nomlon='>>', nomlat='^^', **kwargs):
    '''Read field in lat-lon variables.

    Returns
    -------
    lon : ndarray
        longitudes
    lat : ndarray
        latitudes

    '''

    logger = logging.getLogger(__name__)

    funit = rmn.fstopenall(fstd)

    lon = rmn.fstinl(funit, nomvar=nomlon)
    lon = rmn.fstluk(lon[0], dtype=np.float32)
    lon = lon['d']

    lat = rmn.fstinl(funit, nomvar=nomlat)
    lat = rmn.fstluk(lat[0], dtype=np.float32)
    lat = lat['d']

    rmn.fstcloseall(funit)

    return lon, lat

def read_yy_fstd_var(fstd, nomvar, skip_date=True, **kwargs):
    '''Read field in GIOPS FSTD.

    Returns
    -------
    lon : ndarray
        longitudes
    lat : ndarray
        latitudes
    var : MaskedArray
        data

    '''

    logger = logging.getLogger(__name__)

    funit = rmn.fstopenall(fstd, rmn.FST_RO)

    var_rec = rmn.fstlir(funit, nomvar=nomvar)
    # NOTE:  Date information is not returned!
    if ( not skip_date ):
        datev = rpndate.RPNDate(var_rec['datev']).toDateTime()
        logger.debug('datev: %s', datev.strftime('%c'))
    var = var_rec['d']

    # Prefered method to get grid lat, lon. Works on any RPNSTD grid
    # type (except 'X')
    #try:
    var_gridU = rmn.readGrid(funit, var_rec)
    var_grid0 = var_gridU['subgrid'][0]
    var_grid1 = var_gridU['subgrid'][1]
    gridLatLon0 = rmn.gdll(var_grid0)
    gridLatLon1 = rmn.gdll(var_grid1)
    lat0 = gridLatLon0['lat']
    lon0 = gridLatLon0['lon']
    lat1 = gridLatLon1['lat']
    lon1 = gridLatLon1['lon']
    
    nx0, ny0 = lat0.shape
    nx1, ny1 = lat1.shape
    nx, ny = var.shape
    #print 'NY',ny0, ny1, ny
    
    var0 = var[:, 0:ny0]
    var1 = var[:, ny0:ny]
    #print lat0.shape, lat1.shape
    lat = np.append(lat0, lat1, axis=1)
    lon = np.append(lon0, lon1, axis=1)
    #print 'shape', lat.shape, lon.shape

    rmn.fstcloseall(funit)

    G0 = [var0, lon0, lat0]
    G1 = [var1, lon1, lat1]
    GG = [var , lon , lat] 
    return  GG, G0, G1
    
def read_yy_fstd_multi(fstd, nomvar, **kwargs):
    '''Read field in GIOPS FSTD.

    Returns
    -------
    lon : ndarray
        longitudes
    lat : ndarray
        latitudes
    var : MaskedArray
        data

    '''
    logger = logging.getLogger(__name__)

    funit = rmn.fstopenall(fstd, rmn.FST_RO)

    var_recs = rmn.fstinl(funit, nomvar=nomvar, **kwargs)
    #print var_recs

    var = []
    date = []
    for krec in var_recs:
        varstr=rmn.fstluk(krec)
        datev = rpndate.RPNDate(varstr['datev']).toDateTime()
        date.append(datev)
        #print datev
        #logger.debug('datev: %s', datev.strftime('%c'))
        var.append(varstr['d'])

    # Prefered method to get grid lat, lon. Works on any RPNSTD grid
    # type (except 'X')
    #try:
    varstr=rmn.fstluk(var_recs[0])
    var_gridU = rmn.readGrid(funit, varstr)
    var_grid0 = var_gridU['subgrid'][0]
    var_grid1 = var_gridU['subgrid'][1]
    gridLatLon0 = rmn.gdll(var_grid0)
    gridLatLon1 = rmn.gdll(var_grid1)
    lat0 = gridLatLon0['lat']
    lon0 = gridLatLon0['lon']
    lat1 = gridLatLon1['lat']
    lon1 = gridLatLon1['lon']
    lat = np.append(lat0, lat1, axis=1)
    lon = np.append(lon0, lon1, axis=1)
    
    nx0, ny0 = lat0.shape
    nx1, ny1 = lat1.shape
    nx, ny = var[0].shape
    #print 'NY',ny0, ny1, ny

    var0 = []
    var1 = []
    for irec, ivar in enumerate(var):    
        var0.append(ivar[:, 0:ny0])
        var1.append(ivar[:, ny0:ny])


    rmn.fstcloseall(funit)

    G0 = [var0, date, lon0, lat0]
    G1 = [var1, date, lon1, lat1]
    GG = [var , date, lon , lat] 
    return  GG, G0, G1

def read_yy_fstd_multi_lev(fstd, nomvar, **kwargs):
    '''Read field in GIOPS FSTD.

    Returns
    -------
    lon : ndarray
        longitudes
    lat : ndarray
        latitudes
    var : MaskedArray
        data

    '''
    logger = logging.getLogger(__name__)

    funit = rmn.fstopenall(fstd, rmn.FST_RO)

    var_recs = rmn.fstinl(funit, nomvar=nomvar, **kwargs)
    #print var_recs

    var = []
    lev = []
    date = []
    for krec in var_recs:
        varstr=rmn.fstluk(krec)
        datev = rpndate.RPNDate(varstr['datev']).toDateTime()
        date.append(datev)
        #print datev
        #logger.debug('datev: %s', datev.strftime('%c'))
        var.append(varstr['d'])
        (rp1, rp2, rp3) = rmn.convertIPtoPK(varstr['ip1'], varstr['ip2'], varstr['ip3'])
        lev.append(rmn.FLOATIPtoList(rp1)[0])

    # Prefered method to get grid lat, lon. Works on any RPNSTD grid
    # type (except 'X')
    #try:
    varstr=rmn.fstluk(var_recs[0])
    var_gridU = rmn.readGrid(funit, varstr)
    var_grid0 = var_gridU['subgrid'][0]
    var_grid1 = var_gridU['subgrid'][1]
    gridLatLon0 = rmn.gdll(var_grid0)
    gridLatLon1 = rmn.gdll(var_grid1)
    lat0 = gridLatLon0['lat']
    lon0 = gridLatLon0['lon']
    lat1 = gridLatLon1['lat']
    lon1 = gridLatLon1['lon']
    lat = np.append(lat0, lat1, axis=1)
    lon = np.append(lon0, lon1, axis=1)
    
    nx0, ny0 = lat0.shape
    nx1, ny1 = lat1.shape
    nx, ny = var[0].shape
    #print 'NY',ny0, ny1, ny

    var0 = []
    var1 = []
    for irec, ivar in enumerate(var):    
        var0.append(ivar[:, 0:ny0])
        var1.append(ivar[:, ny0:ny])


    rmn.fstcloseall(funit)

    G0 = [var0, date, lev, lon0, lat0]
    G1 = [var1, date, lev, lon1, lat1]
    GG = [var , date, lev, lon , lat] 
    return  GG, G0, G1

def make_subgrid_fields(FLD, gid, nomvar='TT'):

      if ( isinstance(gid, str) ):
          gid_tuple, __ = read_fstd_gid(gid, nomvar)
      elif ( isinstance(gid, int) ):
          sub_tuple = rmn.decodeGrid(gid)['subgridid']
          gid_tuple = [gid] + sub_tuple
      elif ( isinstance(gid, tuple) ):
          gid_tuple = gid

      (gid, gid0, gid1) = gid_tuple 
      grid_U = rmn.decodeGrid(gid)
      grid_0 = rmn.decodeGrid(gid0)
      grid_1 = rmn.decodeGrid(gid1)
      
      ny = grid_U['nj']
      ny0 = grid_0['nj']
      ny1 = grid_1['nj']
      
      FLD0 = None
      FLD1 = None
      if ( isinstance(FLD, list) ):
          FLD0 = []
          FLD1 = []
          for fld in FLD:
              fld0 = fld[:, 0:ny0]
              fld1 = fld[:, ny0:ny]
              FLD0.append(fld0)
              FLD1.append(fld1)
      elif ( isinstance(FLD, np.ndarray) ):
          FLD0 = FLD[:, 0:ny0]
          FLD1 = FLD[:, ny0:ny]

      return FLD0, FLD1
