import sys
import os
import glob
sys.path.insert(0, '/home/dpe000/GEOPS/python')
import subprocess

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import scipy.interpolate
import netCDF4
import geopy.distance
import rpnpy.librmn.all as rmn

import datetime
import pytz
import time

import stfd
import cplot
import read_grid
import datafiles
import datadatefile

base_dir='/fs/homeu1/eccc/mrd/ords/rpnenv/dpe000/GEOPS/GEOPS'
data_dir='/fs/site4/eccc/mrd/rpnenv/dpe000//GEPS_OPER/geps/000.glboce.ATM'
dato_dir='/fs/site4/eccc/mrd/rpnenv/dpe000//GEPS_OPER/geps/000.glboce.TM'
dats_dir='/fs/site4/eccc/mrd/rpnenv/dpe000//GEPS_OPER/geps/000.glboce.SFC'
datw_dir='/fs/site4/eccc/mrd/rpnenv/dpe000//GEPS_OPER/geps/000.glboce.WND'

plot_dir='/fs/site4/eccc/mrd/rpnenv/dpe000/GEOPS/DIURNAL'
NPbox=[-180, 180, 65, 90]
HBbox=[-97, -47, 44, 69]
ECbox=[-75, -45, 45, 65]
PNbox=[130, 220, 45, 75]

missing = -999.

NPStereo='NorthPolarStereo'
PLCarree='PlateCarree'
PCCarree='PacificCarree'

CLA=np.arange(0.0, 4.0, 0.2)
CLT=np.arange(0.0, 4.0, 0.4)
CLF=np.arange(-2.0, 16.0, 2.0)
CLD=np.arange(-1.9, 2.1, 0.2)

cmap='gist_stern_r'
cmaa='seismic'

def convert_knots_SI(Vel_in_knots):
    Vel_in_ms = Vel_in_knots * 0.514444
    return Vel_in_ms
    
def get_data_dir():
    return data_dir, dato_dir, dats_dir, plot_dir
    
def plt_diurnal(date, lead=24):
    lead_str=str(lead).zfill(3)
    date_str=date.strftime("%Y%m%d%H")
    file=data_dir+'/'+date_str+'_'+lead_str+'_000'
    #print file, os.path.isfile(file)                             
    (TMD, dates, lon, lat) , __, __ =  stfd.read_yy_fstd_multi(file, 'TM')                                                                              

    nx, ny = TMD[0].shape

    TMT = TMD[1] - TMD[5]
    TMS=np.zeros((nx,ny))
    for ix in range(nx):
      for iy in range(ny):
        TMP=[TMD[i][ix, iy] for i in range(len(TMD))]
        MIN=min(TMP)
        MAX=max(TMP)
        TMS[ix, iy] = MAX - MIN

    outfile=plot_dir+'/D.'+date_str+'_'+lead_str+'_'+'pac.png'
    cplot.grd_pcolormesh(lon, lat, TMS, outfile=outfile, levels=CLA, ticks=CLT, cmap=cmap, project='Percator', box=PNbox,ddeg=0.25)
    return

def plt_date(date, lead, pdir=plot_dir):
    mask_nemo = read_grid.read_tmask_mesh(var='tmask')[0]
    lead_str=str(lead).zfill(3)
    date_str=date.strftime("%Y%m%d")
    flat=dato_dir+'/'+'latlon'
    lonn, latn = stfd.read_latlon(flat) 
    for shr in ['00','12']:
        bile=date_str+shr+'_'+lead_str+'_000'
        file=data_dir+'/'+bile
        filo=dato_dir+'/'+bile
	print file, os.path.isfile(file)
	print filo, os.path.isfile(filo)
        (TMD, dates, lon, lat) , __, __ =  stfd.read_yy_fstd_multi(file, 'TM')
        TMA = sum(TMD) / len(TMD)
	TMA = TMA - 273.15
        if ( lead != 0 ):  # No SST in 0h lead file
            LONN, LATN, TMO = stfd.read_fstd_var(filo, 'TM', typvar='P@', ip1=15728640)
	    TMO = TMO - 273.15
            TAO = scipy.interpolate.griddata((lon.flatten(),lat.flatten()),TMA.flatten() , (lonn%360, latn), method='nearest', fill_value=0)
            TAO = np.ma.array(TAO, mask=1-mask_nemo)
        else:
            TMO=None
	
        outfile=pdir+'/TMA.'+date_str+shr+'_'+lead_str+'_'+'pac.png'
        cplot.grd_pcolormesh(lon, lat, TMA, outfile=outfile, levels=CLF, ticks=None, cmap=cmap, project='Percator', box=PNbox,ddeg=0.25,title=shr+'Z SST on Atmospheric Grid')
        if ( shr == '00' ): TMA00=TMA.copy()
        if ( shr == '12' ): TMA12=TMA.copy()
        if ( lead > 0 ):
            TMD = TAO - TMO
            outfile=pdir+'/TMO.'+date_str+shr+'_'+lead_str+'_'+'pac.png'
            cplot.grd_pcolormesh(lonn, latn, TMO, outfile=outfile, levels=CLF, ticks=None, cmap=cmap, project='Percator', box=PNbox,ddeg=0.25,title=shr+'Z SST on Ocean Grid')
            outfile=pdir+'/TMD.'+date_str+shr+'_'+lead_str+'_'+'pac.png'
            cplot.grd_pcolormesh(lonn, latn, TMD, outfile=outfile, levels=CLD, ticks=None, cmap=cmaa, project='Percator', box=PNbox,ddeg=0.25,title=shr+'Z Difference SST A-O')
            if ( shr == '00' ): 
	        TMO00=TMO.copy()
		TMD00=TMD.copy()
            if ( shr == '12' ): 
	        TMO12=TMO.copy()
		TMD12=TMD.copy()
    shr='DD'
    outfile=pdir+'/TMA.'+date_str+shr+'_'+lead_str+'_'+'pac.png'
    cplot.grd_pcolormesh(lon, lat, TMA12-TMA00, outfile=outfile, levels=CLF, ticks=None, cmap=cmap, project='Percator', box=PNbox,ddeg=0.25,title='12Z-00Z SST Diff on Atmospheric Grid')
    if ( lead > 0 ):
        outfile=pdir+'/TMO.'+date_str+shr+'_'+lead_str+'_'+'pac.png'
        cplot.grd_pcolormesh(lonn, latn, TMO12-TMO00, outfile=outfile, levels=CLD, ticks=None, cmap=cmap, project='Percator', box=PNbox,ddeg=0.25,title='12Z-00Z SST diff on Ocean Grid')
        outfile=pdir+'/TMD.'+date_str+shr+'_'+lead_str+'_'+'pac.png'
        cplot.grd_pcolormesh(lonn, latn, TMD12-TMD00, outfile=outfile, levels=CLD, ticks=None, cmap=cmaa, project='Percator', box=PNbox,ddeg=0.25,title='12Z-00Z double diff SST A-O')
    return
    
def find_nearest_point(lon_pt, lat_pt, lon_grid, lat_grid):
    ## lon: 0 -> 360
    ## lat: -90 -> 90
    
    min_distance=[]
    icc = []
    jcc = []
    for cc in [0, -1, 1]:
        lon_cc = lon_pt + cc*360.0 
        distance = np.square(lat_grid-lat_pt)+np.square(lon_grid-lon_cc)
        if ( distance.ndim > 1 ):
	    ipt, jpt = np.unravel_index(np.argmin(distance), distance.shape)
	    min_distance.append(distance[ipt, jpt])
	else:
	    ipt = np.argmin(distance)
	    jpt = 1
	    min_distance.append(distance[ipt])
	icc.append(ipt)
	jcc.append(jpt)
    #print min_distance
    val, idx = min((val, idx) for (idx, val) in enumerate(min_distance))
    #print idx, val
    ipt=icc[idx]
    jpt=jcc[idx]
    return ipt, jpt

def find_nearest_glbpt(lon_pt, lat_pt, lon_grid, lat_grid):
    ## lon: 0 -> 360
    ## lat: -90 -> 90
    
    min_distance=[]
    icc = []
    jcc = []
    for cc in [0, -1, 1]:
        lon_cc = lon_pt + cc*360.0 
        distance = np.square(lat_grid-lat_pt)+np.square(np.cos((np.pi/180.0)*(0.5*lat_grid+0.5*lat_pt)))*np.square(lon_grid-lon_cc)
        ipt, jpt = np.unravel_index(np.argmin(distance), distance.shape)
	min_distance.append(distance[ipt, jpt])
	icc.append(ipt)
	jcc.append(jpt)
    #print min_distance
    val, idx = min((val, idx) for (idx, val) in enumerate(min_distance))
    #print idx, val
    ipt=icc[idx]
    jpt=jcc[idx]
    return ipt, jpt

def find_nearest_glcpt(lon_pt, lat_pt, lon_grid, lat_grid):
    ## lon: 0 -> 360
    ## lat: -90 -> 90
    
    min_distance=[]
    icc = []
    jcc = []
    for cc in [0, -1, 1]:
        lon_cc = lon_pt + cc*360.0 
	converged_lon_cc = lon_cc * np.cos( ( np.pi / 180.0 ) * lat_pt )
	converged_lon_gd = lon_grid * np.cos( ( np.pi / 180.0 ) *  lat_grid )
	distance = np.square(lat_grid-lat_pt)+np.square(converged_lon_gd - converged_lon_cc )
        ipt, jpt = np.unravel_index(np.argmin(distance), distance.shape)
	min_distance.append(distance[ipt, jpt])
	icc.append(ipt)
	jcc.append(jpt)
    #print min_distance
    val, idx = min((val, idx) for (idx, val) in enumerate(min_distance))
    #print idx, val
    ipt=icc[idx]
    jpt=jcc[idx]
    return ipt, jpt

def grid_geopy_distance(lon_pt, lat_pt, lon_grid, lat_grid):
    nx, ny = lon_grid.shape
    distance_grid = np.zeros((nx,ny))
    distance_min = 40e6 # circumference of earth in metres ( minimum distance SHOULD be smaller than this!) 
    for ix in range(nx):
        for iy in range(ny):
	    distance_pt = geopy.distance.distance((lat_grid[ix,iy],lon_grid[ix,iy]),(lat_pt, lon_pt)).m
	    if ( distance_pt < distance_min ):
	        distance_min = distance_pt
		imin = ix
		jmin = iy
	    distance_grid[ix,iy] = distance_pt
    return distance_grid, distance_min, imin, jmin
	    
def find_nearest_geopt(lon_pt, lat_pt, lon_grid, lat_grid):
    ## lon: 0 -> 360
    ## lat: -90 -> 90
    
    min_distance=[]
    icc = []
    jcc = []
    distance_grid, distance_min, imin, jmin = grid_geopy_distance(lon_pt, lat_pt, lon_grid, lat_grid)
    return imin, jmin

def plot_timeseries_point(lon, lat, dates, lead,pdir='DIURNAL'):
    marker='.'
    lon=lon%360  ## Not strickly neccessary
    all_dates = []
    hor_dates = []
    avg_dates = []
    ocn_dates = []
    all_TM = []
    all_TO = []
    hor_TM = []
    avg_TM = []
    fig, ax = plt.subplots()	
    lead_str=str(lead).zfill(3)
    lalo_str_words="%.1fE / %.1fN"%(lon%360,lat)
    lalo_str="%.iE%iN"%(lon%360,lat)
    file=dato_dir+'/'+'latlon'
    lonn, latn = stfd.read_latlon(file) 
    first00=True
    first12=True
    for date in dates:
        hour_str=date.strftime("%H")
	lead_date=date+datetime.timedelta(hours=lead)
	color='b'
	if ( hour_str == '12' ): color='r'
        date_str=date.strftime("%Y%m%d%H")
        file=data_dir+'/'+date_str+'_'+lead_str+'_000'
        #print file, os.path.isfile(file)                             
        (TMD, fcst_dates, lon_grid, lat_grid) , __, __ =  stfd.read_yy_fstd_multi(file, 'TM') 
	if ( date == dates[0] ):
	    ipt, jpt = find_nearest_point(lon, lat, lon_grid, lat_grid) 
	    ipo, jpo = find_nearest_point(lon, lat, lonn, latn)
	    print 'Using ipt/jpt = ', ipt, jpt, lon_grid[ipt,jpt], lat_grid[ipt,jpt]
	    print 'Using ipo/jpo = ', ipo, jpo, lonn[ipo,jpo], latn[ipo,jpo]
	all_dates.extend(fcst_dates)
	hor_dates.append(fcst_dates[-1])
	avg_dates.append( avg_time( fcst_dates ) )
	TM_point = [TMD[idx][ipt, jpt]-273.15 for idx in range(len(TMD))]
	all_TM.extend( TM_point ) 
	hor_TM.append( TM_point[-1] )
	avg_TM.append( sum(TM_point) / len(TM_point) )

        file=dato_dir+'/'+date_str+'_'+lead_str+'_000'
        #print file, os.path.isfile(file)                             
	if ( lead != 0 ):  # No SST in 0h lead file
	    LONN, LATN, TMO = stfd.read_fstd_var(file, 'TM', typvar='P@', ip1=15728640)
	    TO_point = TMO[ipo, jpo]-273.15
            all_TO.append( TO_point )
	    ocn_dates.append( lead_date - datetime.timedelta(hours=12) )

	if ( first00 ):
	   ax.plot(fcst_dates, TM_point, marker=marker, color=color, label='00Z') 
           first00 = False
	elif ( first12 ):
	   ax.plot(fcst_dates, TM_point, marker=marker, color=color, label='12Z') 
	   first12 = False
	else:
	   ax.plot(fcst_dates, TM_point, marker=marker, color=color)	   
	ax.plot(fcst_dates[-1], TM_point[-1], marker='x', color=color)

    myFmt = mdates.DateFormatter('%m/%d')
    ax.xaxis.set_major_formatter(myFmt)
    ax.plot(hor_dates, hor_TM, color='g', label='lead')
    if ( len(all_TO) > 0 ):
        ax.plot(ocn_dates, all_TO, color='m', label='ocn')
    ax.legend()
    ax.set_title('Temperature at '+lalo_str_words)
    fig.savefig(pdir+'/'+'TM'+lalo_str+'_'+lead_str+'.pdf')
    fig.savefig(pdir+'/'+'TM'+lalo_str+'_'+lead_str+'.png')
    plt.close(fig)
    
    fig, ax = plt.subplots()	
    ax.plot( avg_dates, avg_TM, color='b', marker='x', label='avg')
    ax.plot( hor_dates, hor_TM, color='g', marker='x', label='lead')
    if ( len(all_TO) > 0 ):
        ax.plot( ocn_dates, all_TO, color='m', marker='x', label='ocn')
    ax.legend()
    ax.xaxis.set_major_formatter(myFmt)
    ax.set_title('Temperature at '+lalo_str_words)
    fig.savefig(pdir+'/'+'TO'+lalo_str+'_'+lead_str+'.pdf')
    fig.savefig(pdir+'/'+'TO'+lalo_str+'_'+lead_str+'.png')
    plt.close(fig)
    return

def avg_time(dates):
    date_ref = dates[0]
    delta = [ (date - date_ref).total_seconds() for date in dates ]
    avgdt = sum(delta) / len(dates)
    date_avg = date_ref + datetime.timedelta(seconds=avgdt)
    return date_avg
 
#def make_TZaware(dates):
#    new_dates=[]
#    for date in dates:
#        new_date = date.replace(tzinfo=pytz.UTC)
#	new_dates.append(new_date)
#    return new_dates

def make_TZaware(date):
    if ( isinstance(date, list) ):
        new_date=[]
        for idate in date:
	    new_idate = make_TZaware(idate)
	    new_date.append(new_idate)
    elif ( isinstance(date, datetime.datetime) ):
        new_date = date.replace(tzinfo=pytz.UTC)
    else:
        new_date = []
    return new_date

       
def read_drone_latlon(drone=1034):
    drone_str=str(drone).zfill(4)
    if ( drone_str == '0000' ):
        position_file=base_dir+'/'+'NOAA_SailDrones/LatLon_1034.txt'
    else:
        position_file=base_dir+'/'+'NOAA_SailDrones/LatLon_'+drone_str+'.txt'
    dronedata = np.loadtxt(position_file, unpack=True, delimiter=',') 
    M = dronedata[0]
    D = dronedata[1]
    lat = dronedata[2]
    lon = dronedata[3]
    if ( drone_str == '0000' ): 
        lat = np.zeros(lat.shape) +  65.33333333
	lon = np.zeros(lon.shape) +  2.333333333
    dates=[]
    for ii in range(len(M)):
        dates.append(datetime.datetime(2019, int(M[ii]), int(D[ii]), 12))
    return dates, lat, lon

def read_drone_latlon2(drone=1034):
    files=[base_dir+'/'+'NOAA_SailDrones/Lat_sd.txt',base_dir+'/'+'NOAA_SailDrones/Lon_sd.txt']
    fileLat=files[0]
    fileLon=files[1]
    #Attached are two files for daily lat/lon along the track for the period 6-1 to 9-30 (122 days total).
    #Columns are for the boats (1034, 1035, 1036, and 1037) (written in the file).
    #Please let me know if you have any questions.
    #Thank you!
    #Muyin

    drone=int(drone)
    Date0=datetime.datetime(2019,6,1,12)
    dates=[]
    lat=[]
    lon=[]
    for idate in range(122):
        dates.append(Date0+datetime.timedelta(days=idate))
    LatLoad=np.loadtxt(fileLat)
    LonLoad=np.loadtxt(fileLon)
    Drones=LatLoad[0].astype(int).tolist()
    Dronec=LonLoad[0].astype(int).tolist()
    if ( Drones != Dronec ): return dates, lat, long
    LatDrone=LatLoad[1:,:]
    LonDrone=LonLoad[1:,:]
    if ( drone != 0 ):
        idrone=Drones.index(drone)
        lat=LatDrone[:,idrone]
        lon=LonDrone[:,idrone]
        inan = np.where( np.isnan(lat) )
        lat[inan] = missing
        lon[inan] = missing
    else:
        lat=np.zeros(122) + 65.33333333
	lon=np.zeros(122) +  2.333333333
    return dates, lat, lon


def double_dates(dates, lon, lat):
    ddates = [dates[0]]
    dlon = [lon[0]]
    dlat = [lat[0]]
    for idate in range(1,len(dates)):
        middate = dates[idate-1] + ( dates[idate] - dates[idate-1] ) / 2
	ddates.extend( [middate, dates[idate]] )
	if ( ( lon[idate-1] != missing ) and ( lon[idate] != missing ) ):
	    midlon = 0.5 * ( lon[idate-1] + lon[idate] )
	    midlat = 0.5 * ( lat[idate-1] + lat[idate] )
	else:
	    midlon = missing
	    midlat = missing
	dlon.extend( [midlon, lon[idate] ])
	dlat.extend( [midlat, lat[idate] ])
    return ddates, dlon, dlat
    
def double_datee(dates, lon, lat): 
    ddates = [ dates[0]- (dates[1]-dates[0])/2 ] 
    dlon = [lon[0]-0.5*(lon[1]-lon[0])] 
    dlat = [lat[0]-0.5*(lat[1]-lat[0])] 
    for idate in range(len(dates)): 
        if ( idate == len(dates)-1 ): 
	    middate = dates[idate] + ( dates[idate] - dates[idate-1] ) / 2 
	    if ( ( lon[idate] != missing ) and ( lon[idate-1] != missing ) ):
	        midlon = lon[idate] + 0.5 * ( lon[idate]-lon[idate-1] ) 
	        midlat = lat[idate] + 0.5 * ( lat[idate]-lat[idate-1] ) 
	    else:
	        midlon = missing
		madlat = missing
	else:
            middate = dates[idate] + ( dates[idate+1] - dates[idate] ) / 2 
	    if ( ( lon[idate] != missing ) and ( lon[idate+1] != missing ) ):
	        midlon = 0.5 * ( lon[idate] + lon[idate+1] ) 
	        midlat = 0.5 * ( lat[idate] + lat[idate+1] ) 
	    else:
	        midlon = missing
		midlat = missing
	ddates.extend( [dates[idate], middate] ) 
	dlon.extend( [lon[idate], midlon]) 
	dlat.extend( [lat[idate], midlat]) 
    return ddates, dlon, dlat

def interpolate_to_point(TM_grid, lon_grid, lat_grid, lon_pt, lat_pt, method='linear', convlon=False):
    if ( isinstance(TM_grid, np.ma.core.MaskedArray) ):
        ivalid = np.where(TM_grid.mask == False)
	TM_flat = TM_grid[ivalid]
	lon_flat = lon_grid[ivalid]%360
	lat_flat = lat_grid[ivalid]
    else:
        TM_flat = TM_grid.flatten()
	lon_flat = lon_grid.flatten()%360
	lat_flat = lat_grid.flatten()

    lon_in = lon_pt%360
    if ( convlon ):
        lon_flat = lon_flat * np.cos(lat_flat * np.pi / 180.0 )
        lon_in = lon_in * np.cos(lat_pt * np.pi / 180.0 )
    TM_pt = float(scipy.interpolate.griddata( (lon_flat,lat_flat), TM_flat, (lon_in, lat_pt), method=method, fill_value=missing))
    return TM_pt
    
def ezinterpolate_to_point(lon_pt, lat_pt, TM_grid, file, field, src='A', grid_lat=None, grid_lon=None):
    if ( src == '0' ):
        TM_pt = np.NaN
    if ( ( src == 'A' ) or ( src == 'S' ) or ( src == 'W' ) ):  # Yin Yan grids
        gid_tuple, flds = stfd.read_fstd_gid(file, field)
        (gid, gid0, gid1) = gid_tuple


    if ( ( src == 'A' ) or ( src == 'S' ) ): # Scalar fields.  Direct call to gdllsval should work
        TM_pt  = float( rmn.gdllsval(gid,  lat_pt,  lon_pt,  TM_grid)  )  

    if ( src == 'W' ):  # Vector points.  Need gdllvval -- but also need to work on subgrids.
        if ( isinstance(grid_lon, type(None) ) ):  # Can either pass grid_lon or read from file.
            (__, __, grid_lon, grid_lat) , __, __ =  stfd.read_yy_fstd_multi(file, field)
        UU_grid = TM_grid[0]
        VV_grid = TM_grid[1]

        # Need to break these up into the 2 sub-grids 
	
        ( UU_grid0, UU_grid1 ) = stfd.make_subgrid_fields(UU_grid, gid_tuple)
	( VV_grid0, VV_grid1 ) = stfd.make_subgrid_fields(VV_grid, gid_tuple)
	( lon0, lon1 ) = stfd.make_subgrid_fields(grid_lon, gid_tuple)
	( lat0, lat1 ) = stfd.make_subgrid_fields(grid_lat, gid_tuple)

        (uupt0,vvpt0)  = rmn.gdllvval(gid0,  lat_pt,  lon_pt,  UU_grid0, VV_grid0)
        (uupt1,vvpt1)  = rmn.gdllvval(gid1,  lat_pt,  lon_pt,  UU_grid1, VV_grid1)

        spdpt0 = np.zeros(1).astype(np.float32)
        dirpt0 = np.zeros(1).astype(np.float32)
        spdpt1 = np.zeros(1).astype(np.float32)
        dirpt1 = np.zeros(1).astype(np.float32)
        # This should now give the meteorological direction relative to true north (plus speed).
        ier = rmn.c_gdwdfuv(gid0, spdpt0, dirpt0, uupt0, vvpt0, np.array(lat_pt).astype(np.float32), np.array(lon_pt).astype(np.float32), 1)
        ier = rmn.c_gdwdfuv(gid1, spdpt1, dirpt1, uupt1, vvpt1, np.array(lat_pt).astype(np.float32), np.array(lon_pt).astype(np.float32), 1)
        
	uept0 = -1.0 * spdpt0 * np.sin(dirpt0 * np.pi / 180.0) 
	vnpt0 = -1.0 * spdpt0 * np.cos(dirpt0 * np.pi / 180.0) 
	uept1 = -1.0 * spdpt1 * np.sin(dirpt1 * np.pi / 180.0) 
	vnpt1 = -1.0 * spdpt1 * np.cos(dirpt1 * np.pi / 180.0) 
	print 'COMPONENTS', uept0, vnpt0, uept1, vnpt1
        uept0, vnpt0 = wind_components_from_direction( spdpt0, dirpt0)
        uept1, vnpt1 = wind_components_from_direction( spdpt1, dirpt1)
	print 'COMPONENTS', uept0, vnpt0, uept1, vnpt1

        # Only one of these is correct
        ipt, jpt = find_nearest_point(lon_pt, lat_pt, grid_lon, grid_lat)
	
	in0 = len( np.where( ( lat0 == grid_lat[ipt, jpt] ) & ( lon0 == grid_lon[ipt, jpt] ) )[0]) > 0
	in1 = len( np.where( ( lat1 == grid_lat[ipt, jpt] ) & ( lon1 == grid_lon[ipt, jpt] ) )[0]) > 0
	
	if ( in0 ): 
	    uept = uept0
	    vnpt = vnpt0
	    spdpt = spdpt0
	    dirpt = dirpt0
	if ( in1): 
	    uept = uept1
	    vnpt = vnpt1
	    spdpt = spdpt1
	    dirpt = dirpt1
	TM_pt = ( float(uept), float(vnpt), float(spdpt), float(dirpt))
    return TM_pt

def wind_components_from_direction( spd, dir):
    uu = -1.0 * spd * np.sin(dir * np.pi / 180.0)   
    vv = -1.0 * spd * np.cos(dir * np.pi / 180.0)
    return uu, vv
    
def generate_NS_EW_velocities( UU_tuple, grid_lat, grid_lon, file):

    gid_tuple, flds = stfd.read_fstd_gid(file, 'UU')
    (gid, gid0, gid1) = gid_tuple
    
    (uu, vv)  = UU_tuple
    
    uu0, uu1 = stfd.make_subgrid_fields(uu, gid_tuple)
    vv0, vv1 = stfd.make_subgrid_fields(vv, gid_tuple)
    
    (nx0, ny0) = uu0.shape
    (nx1, ny1) = uu1.shape
    
    lat0, lat1 = stfd.make_subgrid_fields(grid_lat, gid_tuple)
    lon0, lon1 = stfd.make_subgrid_fields(grid_lon, gid_tuple)
    
    spdllout0 = np.zeros(uu0.shape).astype(np.float32)
    dirllout0 = np.zeros(uu0.shape).astype(np.float32)
    spdllout1 = np.zeros(uu1.shape).astype(np.float32)
    dirllout1 = np.zeros(uu1.shape).astype(np.float32)
    
    ier = rmn.c_gdwdfuv(gid0, spdllout0, dirllout0, uu0.flatten(), vv0.flatten(), lat0.flatten(), lon0.flatten(), nx0*ny0)
    ier = rmn.c_gdwdfuv(gid1, spdllout1, dirllout1, uu1.flatten(), vv1.flatten(), lat1.flatten(), lon1.flatten(), nx1*ny1)
    spdllout = np.append(spdllout0, spdllout1, axis=1)
    dirllout = np.append(dirllout0, dirllout1, axis=1)
    
    ## WHAT I REALLY WANT IS east/west north/south velocities
    uue0 = -1.0 * spdllout0 * np.sin(dirllout0 * np.pi / 180.0)
    vvn0 = -1.0 * spdllout0 * np.cos(dirllout0 * np.pi / 180.0)
    uue1 = -1.0 * spdllout1 * np.sin(dirllout1 * np.pi / 180.0)
    vvn1 = -1.0 * spdllout1 * np.cos(dirllout1 * np.pi / 180.0)

    uue0, vvn0 = wind_components_from_direction( spdllout0, dirllout0)
    uue1, vvn1 = wind_components_from_direction( spdllout1, dirllout1)
    
    uue = np.append(uue0, uue1, axis=1)
    vvn = np.append(vvn0, vvn1, axis=1)
    
    UEN = (uue, vvn, spdllout, dirllout)

    return UEN
    
def get_drone_SST(lead, drone=1034, src='O',field='TM',method=0):
    fill_value = missing
    ip1=-1
    if ( src == 'O' ):
        if ( field == 'TM' ): ip1=15728640
    if ( src == 'A' ):
        if ( field == 'TM' ): ip1=0
    if ( src == 'S' ):
        if ( field == 'TM' ): ip1=0
	if ( field == 'TT' ): ip1=76696048
	if ( field == 'HU' ): ip1=76696048
	if ( field == 'PN' ): ip1=0
    if ( src == 'W' ):
	ip1=75597472
	
    # OCEAN TEMPERATURE IS A DAILY AVERAGE VALID 12H BEFORE LEAD TIME.
    print src
    if ( ( src != 'O' ) or ( src != 'A' )  or ( src != 'S' ) ):
       print 'src is O or A or S'
    if ( src == 'O' ):  
        # I'd rather handle this manually
	#real_lead = lead - 12
        real_lead = lead
    else:
        real_lead = lead

    if ( ( method == 0 ) or ( method == 10 ) or ( method == 100 )  or ( method == 110 ) or ( method == 200 ) or ( method == 210 ) ):
        dates, lat, lon = read_drone_latlon(drone)
    elif ( ( method == 1 ) or ( method == 11 ) or ( method == 101 ) or ( method == 111 ) or  (method == 201) or ( method == 211 ) ):
        dates, lat, lon = read_drone_latlon2(drone)

    ddates, dlat, dlon = double_dates(dates, lat, lon)
    if ( ( method == 10 ) or ( method == 11 ) or ( method == 110 ) or ( method == 111 ) or ( method == 210 ) or ( method == 211 ) ):
        ddates, dlat, dlon = double_datee(dates, lat, lon)

    T=[]
    M=[]
    sdates=[]
    for idate, date in enumerate(ddates):
        time0 = time.time()
        start_date = date - datetime.timedelta(hours=real_lead)
	date_str=start_date.strftime("%Y%m%d%H")
	lead_str=str(lead).zfill(3)
	bile=date_str+'_'+lead_str+'_000'	
	print 'Using '+bile+' for '+date.strftime("%Y%m%d%H")+' at '+lead_str+' ( '+start_date.strftime("%Y%m%d%H")+' )'
	exist=False
	drone_exist=True
	if ( dlat[idate] == missing ): drone_exist=False
	if ( src == 'O' ): 
            file=dato_dir+'/'+'latlon'
            grid_lon, grid_lat = stfd.read_latlon(file) 
	    file=dato_dir+'/'+bile
	    if ( os.path.isfile(file) ):
		exist=True
                __, __, TM = stfd.read_fstd_var(file, field, typvar='P@', ip1=ip1)
	if ( src == 'A' ): 
	    file=data_dir+'/'+bile
 	    if ( os.path.isfile(file) ):
	        exist=True
                (TMD, rates, grid_lon, grid_lat) , __, __ =  stfd.read_yy_fstd_multi(file, field, ip1=ip1)	
	        TM = TMD[-1]
	        rate=rates[-1]
	        print 'Drone date '+date.strftime("%Y%m%d%H")+' / '+'Valid date '+rate.strftime("%Y%m%d%H")
	if ( src == 'S' ): 
	    file=dats_dir+'/'+bile
 	    if ( os.path.isfile(file) ):
	        exist=True
                (TMD, rates, grid_lon, grid_lat) , __, __ =  stfd.read_yy_fstd_multi(file, field, ip1=ip1)	
	        TM = TMD[-1]
	        rate=rates[-1]
	        print 'Drone date '+date.strftime("%Y%m%d%H")+' / '+'Valid date '+rate.strftime("%Y%m%d%H")
	if ( src == 'W' ): 
	    file=datw_dir+'/'+bile
 	    if ( os.path.isfile(file) ):
	        exist=True
		print(file)
		GU, GU0, GU1 = stfd.read_yy_fstd_multi(file, 'UU', ip1=ip1)
		(UUD, rates, grid_lon, grid_lat) = GU
		## Standard File velocities are in KNOTS!!!!
	        UU = convert_knots_SI(UUD[-1])
		GV, GV0, GV1 = stfd.read_yy_fstd_multi(file, 'VV', ip1=ip1)
		(VVD, rates, grid_lon, grid_lat) = GV
		## Standard File velocities are in KNOTS!!!!
	        VV = convert_knots_SI(VVD[-1])
		# GENERATE N/S and E/W velocities
		(UUE, VVN, SPD, DIR) = generate_NS_EW_velocities( (UU, VV), grid_lat, grid_lon, file)
		
	        rate=rates[-1]
	        print 'Drone date '+date.strftime("%Y%m%d%H")+' / '+'Valid date '+rate.strftime("%Y%m%d%H")
		
        if ( exist and drone_exist ):
	    ipt, jpt = find_nearest_point(dlon[idate], dlat[idate], grid_lon, grid_lat)
            if ( src == 'W' ):
	        UEpt = UUE[ipt, jpt]
		VNpt = VVN[ipt, jpt]
		SPpt = SPD[ipt, jpt]
		DIpt = DIR[ipt, jpt]
	        TMpt = ( UEpt, VNpt, SPpt, DIpt )
	    else:
	        TMpt = TM[ipt, jpt]
	    TMpi = np.nan
	    TMpe = np.nan
	    if ( src == 'W' ):
	      TMpi = (np.nan, np.nan, np.nan, np.nan)
	      TMpe = (np.nan, np.nan, np.nan, np.nan)
	    if ( method > 99 ):
	      if ( src == 'W' ):
	        UEpi = interpolate_to_point(UUE, grid_lon, grid_lat, dlon[idate], dlat[idate])
		VNpi = interpolate_to_point(VVN, grid_lon, grid_lat, dlon[idate], dlat[idate])
		SPpi = interpolate_to_point(SPD, grid_lon, grid_lat, dlon[idate], dlat[idate])
		DIpi = interpolate_to_point(DIR, grid_lon, grid_lat, dlon[idate], dlat[idate])
	        TMpi = ( UEpi, VNpi, SPpi, DIpi )
	      else:
	          TMpi = interpolate_to_point(TM, grid_lon, grid_lat, dlon[idate], dlat[idate])
	      if ( src == 'O' ):  # DON'T KNOW HOW TO USE EZinterp on ocean grid yet.
	          TMpe = TMpi
              elif ( src == 'W' ):
	          TMpe = ezinterpolate_to_point(dlon[idate], dlat[idate], (UU, VV), file, field, src=src,  grid_lat=grid_lat, grid_lon=grid_lon)
	      else:
	          TMpe = ezinterpolate_to_point(dlon[idate], dlat[idate], TM, file, field, src=src)
	    print 'Using '+str(grid_lon[ipt,jpt])+'/'+str(grid_lat[ipt,jpt])+' for '+str(dlon[idate])+'/'+str(dlat[idate])
	    if ( field == 'TM' ): # Convert to Celsius from Kelvin
	        TMpt = TMpt - 273.15
		TMpi = TMpi - 273.15
		TMpe = TMpe - 273.15
	        print 'SST is '+str(TMpt)+' or interpolated '+str(TMpi)+' or ez interpolated '+str(TMpe)
	    else:
	        print field+' is '+str(TMpt)+' or interpolated '+str(TMpi)+' or ez interpolated '+str(TMpe)
            if ( method < 100 ):   # Use Nearest Point
	        T.append( TMpt )
	    elif ( method < 200 ):  # Use Interpolated Value
	        T.append( TMpi )
	    else:  # Use E-Z Interpolated Value
	        T.append( TMpe )
	    if ( src == 'W' ):
	        M.append( (False, False, False, False) )
	    else:
	        M.append( False )
	    sdates.append(start_date)
	else:
	    print 'NO DATA'
	    if ( src == 'W' ):
	      T.append( (fill_value, fill_value, fill_value, fill_value) )
	      M.append( (True, True, True, True) )
	    else:
	      T.append( fill_value )
	      M.append(True)
	    sdates.append(start_date)
	timeE = time.time() - time0
	print 'Elapsed Time', timeE
        sys.stdout.flush()
	
    T_ma = np.ma.array(T, mask=M, fill_value=fill_value)
    return ddates, sdates, dlat, dlon, T_ma
 
 
def read_drone_Temp(file='NOAA_SailDrones/TEMP_SBE37_MEAN.1034.nc'):
    dataset = netCDF4.Dataset(file) 
    tsec=dataset.variables['TAX1'][:]
    TEMP=dataset.variables['TEMP_SBE37_MEAN'][:]
    tref=datetime.datetime(1970,1,1,0,tzinfo=pytz.UTC)
    time=[]
    for isec in tsec:
        time.append( tref + datetime.timedelta(seconds=isec) )
    return time, TEMP

def bin_time(time, time_bins):
    nt=len(time_bins)
    bins = []
    for it in range(nt-1):
        this_bin = []
        for iit, this_time in enumerate(time):
            in_bin=( this_time > time_bins[it] ) and ( this_time <= time_bins[it+1] )
            if ( in_bin ): this_bin.append(iit)
	bins.append(this_bin)
    return bins
    
def bin_temp(TEMP, time, time_bins):
    bins = bin_time(time, time_bins)
    print len(bins)
    binTEMP = []
    binTIME = []
    for bin in bins:
        if ( len(bin) != 0 ):
            thisTEMP = np.mean(TEMP[bin])
	    thisTIME = avg_time([time[index] for index in bin])
	else:
	    thisTEMP = missing
	    thisTIME = missing
	binTEMP.append(thisTEMP)
	binTIME.append(thisTIME)
    npTEMP=np.ma.array(binTEMP)
    return binTIME, npTEMP
    
def create_binned_time(date_start, date_final, dt_in_hours):
     date_start = make_TZaware(date_start)
     date_final = make_TZaware(date_final)
     loop_time = date_start
     binned_time=[]
     while (loop_time <= date_final ):
         binned_time.append(loop_time)
	 loop_time = loop_time + datetime.timedelta(hours=dt_in_hours)   
     return binned_time

drones = [1034, 1035, 1036, 1037, 0000]
leads = [0, 24, 120, 240, 360]

def get_out_dir(method=0):
    mstr=str(method)
    out_dir=base_dir+'/'+'DRONES_'+mstr
    subprocess.call([ 'mkdir', '-p', out_dir])
    return out_dir
    
def analyse_drone_timeseries(drone, lead, method=0):
    drone_str=str(drone).zfill(4)
    print drone_str
    lead_str=str(lead).zfill(3)

    out_dir = get_out_dir(method=method)
        
    temperature_file = 'NOAA_SailDrones/TEMP_SBE37_MEAN.'+drone_str+'.nc'

    x = sys.stdout
    f = open(out_dir+'/Check_DRONE_'+drone_str+'_'+lead_str+'.out', 'w')
    sys.stdout = f
    OCEAN_RESULTS =  get_drone_SST(lead, drone=drone, src='O', method=method)
    ATMOS_RESULTS =  get_drone_SST(lead, drone=drone, src='A', method=method)
    AIRTT_RESULTS =  get_drone_SST(lead, drone=drone, src='S', field='TT', method=method)
    AMSLP_RESULTS =  get_drone_SST(lead, drone=drone, src='S', field='PN', method=method)
    AHUMI_RESULTS =  get_drone_SST(lead, drone=drone, src='S', field='HU', method=method)
    AWIND_RESULTS =  get_drone_SST(lead, drone=drone, src='W', field='UU', method=method)
    sys.stdout = x
    
    Odates, OSdates, dlat, dlon, TO =  OCEAN_RESULTS
    Adates, ASdates, dlat, dlon, TA =  ATMOS_RESULTS
    Tdates, TSdates, __, __, TT =  AIRTT_RESULTS
    Pdates, PSdates, __, __, PN =  AMSLP_RESULTS
    Hdates, HSdates, __, __, HU =  AHUMI_RESULTS
    Wdates, WSdates, __, __, WD =  AWIND_RESULTS

    Rdates = [ date-datetime.timedelta(hours=12) for date in Odates ]
    # create ascii time series to send to Chidong
    Odata=np.array([TO, dlat, dlon])
    Adata=np.array([TA, dlat, dlon])
    Tdata=np.array([TT, dlat, dlon])
    Pdata=np.array([PN, dlat, dlon])
    Hdata=np.array([HU, dlat, dlon])
    print 'SHAPE', len(HU), Hdata.shape
    print 'SHAPE', len(WD), len(WD[0]), type(WD[0][0]), len(dlat), len(dlon)
    #Wlist = WD[:] ; Wlist.append(dlat) ; Wlist.append(dlon)
    print 'SHAPE', np.array(WD).shape, np.array([dlat, dlon]).shape
    Wdata=np.append( np.transpose(np.array(WD)), np.array([dlat, dlon]), axis=0 )
    print 'SHAPE', Wdata.shape

    Ofile= out_dir+'/TMO_'+drone_str+'_'+lead_str+'.dat'
    Afile= out_dir+'/TMA_'+drone_str+'_'+lead_str+'.dat'
    Tfile= out_dir+'/T2M_'+drone_str+'_'+lead_str+'.dat'
    Pfile= out_dir+'/SLP_'+drone_str+'_'+lead_str+'.dat'
    Hfile= out_dir+'/H2M_'+drone_str+'_'+lead_str+'.dat'
    Wfile= out_dir+'/W10_'+drone_str+'_'+lead_str+'.dat'

    Odates_str, Odates_int = datadatefile.convert_datelist_strint(Odates)
    OSdates_str, OSdates_int = datadatefile.convert_datelist_strint(OSdates)
    Adates_str, Adates_int = datadatefile.convert_datelist_strint(Adates)
    ASdates_str, ASdates_int = datadatefile.convert_datelist_strint(ASdates)
    Tdates_str, Tdates_int = datadatefile.convert_datelist_strint(Tdates)
    TSdates_str, TSdates_int = datadatefile.convert_datelist_strint(TSdates)
    Pdates_str, Pdates_int = datadatefile.convert_datelist_strint(Pdates)
    PSdates_str, PSdates_int = datadatefile.convert_datelist_strint(PSdates)
    Hdates_str, Hdates_int = datadatefile.convert_datelist_strint(Hdates)
    HSdates_str, HSdates_int = datadatefile.convert_datelist_strint(HSdates)
    Wdates_str, Wdates_int = datadatefile.convert_datelist_strint(Wdates)
    WSdates_str, WSdates_int = datadatefile.convert_datelist_strint(WSdates)

    datadatefile.write_2dates_file(Odates_int, OSdates_int, Odata, file=Ofile)
    datadatefile.write_2dates_file(Adates_int, ASdates_int, Adata, file=Afile)
    datadatefile.write_2dates_file(Tdates_int, TSdates_int, Tdata, file=Tfile)
    datadatefile.write_2dates_file(Pdates_int, PSdates_int, Pdata, file=Pfile)
    datadatefile.write_2dates_file(Hdates_int, HSdates_int, Hdata, file=Hfile)
    datadatefile.write_2dates_file(Wdates_int, WSdates_int, Wdata, file=Wfile)

    #Ocean time stamp is actually 12h earlier (daily mean of previous 12h).    
    O12dates=[]
    for idate in Odates:
        O12dates.append(idate - datetime.timedelta(hours=12))

    dronetime, TEMP = read_drone_Temp(file=temperature_file)
    ## THIS SHOULD(?) be a time-series with the same length as the model time-series
    ##    Which starts on 6/1 and ends 8/31
    if ( ( method == 0 ) or ( method == 100 ) ):
        binned_time=create_binned_time(datetime.datetime(2019,6,1,6), datetime.datetime(2019,9,1,0), 12)
    if ( ( method == 1 ) or ( method == 101 ) ):
        binned_time=create_binned_time(datetime.datetime(2019,6,1,6), datetime.datetime(2019,10,1,0), 12)
    if ( ( method == 10 ) or ( method == 110 ) ):
        binned_time=create_binned_time(datetime.datetime(2019,5,31,18), datetime.datetime(2019,9,1,6), 12)
    if ( ( method == 11 ) or ( method == 111 ) ):
        binned_time=create_binned_time(datetime.datetime(2019,6,1,6), datetime.datetime(2019,10,1,6), 12)
    bin_TIME, bin_TEMP = bin_temp(TEMP, dronetime, binned_time)

    ivalio = np.where(TO.mask == False)
    ivalia = np.where(TA.mask == False)
    AVGO = np.mean(bin_TEMP[ivalio] - TO[ivalio])
    AVGA = np.mean(bin_TEMP[ivalia] - TA[ivalia]) 
    AVGD = np.mean(TO[ivalio] - TA[ivalio])  ## A will always be valid where O is valid.  Not vice-versa

    RMSO = np.sqrt( np.mean(np.square(bin_TEMP[ivalio] - TO[ivalio])))
    RMSA = np.sqrt( np.mean(np.square(bin_TEMP[ivalia] - TA[ivalia])))
    RMSD = np.sqrt( np.mean(np.square(TO[ivalio] - TA[ivalio])))

    if ( len(ivalio[0]) > 0 ):
        MO, BO = np.polyfit(TT[ivalio], TO[ivalio], 1) 
    else:
        MO = 0.0
	BO = 0.0	

    if ( len(ivalia[0]) > 0 ):
        MA, BA = np.polyfit(TT[ivalia], TA[ivalia], 1) 
    else:
        MA = 0.0
	BA = 0.0	
    
    print 'SAILDRONE: '+drone_str 
    print 'Models at leadtime: '+lead_str
    print 'Mean Error :', AVGO, AVGA, AVGD, 'Ocn/Atm/DOA'
    print 'RMS Error :', RMSO, RMSA, RMSD, 'Ocn/Atm/DOA'
    print 'Regression :', str(MO)+'TT+'+str(BO), str(MA)+'TT+'+str(BA)

    outfile = out_dir+'/TD_'+drone_str+'_'+lead_str
    fig, ax = plt.subplots()
    ax.plot(Odates, TO, color='c',label='Ocean TM at midtime')
    ax.plot(Rdates, TO, color='b',label='Ocean TM')
    ax.plot(Adates, TA, color='r',label='Atmos TM')
    ax.plot(dronetime, TEMP, color='k', linewidth=0.1)
    ax.plot(bin_TIME, bin_TEMP, linewidth=2.0, color='k',label='Drone')
    ax.legend()
    ax.set_xlabel('Date')
    ax.set_ylabel('Sea Surface Temperature')
    fig.savefig(outfile+'.png')
    fig.savefig(outfile+'.pdf')
    plt.close(fig)

    outfile = out_dir+'/TA_'+drone_str+'_'+lead_str
    fig, ax = plt.subplots()
    ax.plot([bin_TIME[iivalid] for iivalid in ivalio[0]], TO[ivalio]-bin_TEMP[ivalio], color='b',label='Ocean TM Error wrt Drone')
    ax.plot([bin_TIME[iivalid] for iivalid in ivalia[0]], TA[ivalia]-bin_TEMP[ivalia], color='r',label='Atmos TM Error wrt Drone')
    ax.legend()
    ax.set_xlabel('Date')
    ax.set_ylabel('Sea Surface Temperature Error')
    fig.savefig(outfile+'.png')
    fig.savefig(outfile+'.pdf')
    plt.close(fig)

    outfile = out_dir+'/TR_'+drone_str+'_'+lead_str
    fig, ax = plt.subplots()
    ax.scatter(TT, TO, color='b',label='Ocean')
    if ( MO > 0.0 ):
        ax.plot(TT , MO*TT +BO, color='b')
    ax.scatter(TT, TA, color='r',label='Atmos')
    ax.plot(TT , MA*TT +BA, color='r')
    ax.set_xlabel('2m Air Temperature')
    ax.set_ylabel('Sea Surface Temperature')
    ax.legend()
    fig.savefig(outfile+'.png')
    fig.savefig(outfile+'.pdf')
    plt.close(fig)

    outfile = out_dir+'/TD_'+drone_str+'_'+lead_str+'_AUG'
    fig, ax = plt.subplots()
    ax.plot(Odates, TO, color='c',label='Ocean TM at midtime')
    ax.plot(Rdates, TO, color='b',label='Ocean TM')
    ax.plot(Adates, TA, color='r',label='Ocean TM')
    ax.plot(dronetime, TEMP, color='k', linewidth=0.1)
    ax.plot(bin_TIME, bin_TEMP, linewidth=2.0, color='k',label='drone')
    ax.set_xlim([datetime.datetime(2019,8,1, tzinfo=pytz.UTC), datetime.datetime(2019,9,1, tzinfo=pytz.UTC)])
    ax.legend()
    ax.set_xlabel('Date')
    ax.set_ylabel('Sea Surface Temperature')
    fig.savefig(outfile+'.png')
    fig.savefig(outfile+'.pdf')
    plt.close(fig)

    outfile = out_dir+'/TA_'+drone_str+'_'+lead_str+'_AUG'
    fig, ax = plt.subplots()
    ax.plot([bin_TIME[iivalid] for iivalid in ivalio[0]], TO[ivalio]-bin_TEMP[ivalio], color='b',label='Ocean TM Error wrt Drone')
    ax.plot([bin_TIME[iivalid] for iivalid in ivalia[0]], TA[ivalia]-bin_TEMP[ivalia], color='r',label='Atmos TM Error wrt Drone')
    ax.set_xlim([datetime.datetime(2019,8,1, tzinfo=pytz.UTC), datetime.datetime(2019,9,1, tzinfo=pytz.UTC)])
    ax.legend()
    ax.set_xlabel('Date')
    ax.set_ylabel('Sea Surface Temperature Error')
    fig.savefig(outfile+'.png')
    fig.savefig(outfile+'.pdf')
    plt.close(fig)

    ERROR_ARRAY = (AVGO, AVGA, AVGD, RMSO, RMSA, RMSD)
    sys.stdout.flush()

    return OCEAN_RESULTS, ATMOS_RESULTS, ERROR_ARRAY

def loop_leads(drone=drones[0], argleads=leads, method=0):
    drone_str=str(drone).zfill(4)
    out_dir=get_out_dir(method=method)
    ALL_TO = []
    ALL_TA = []
    ALL_ER = []
    for lead in argleads:
        OCEAN_RESULTS, ATMOS_RESULTS, ERROR_ARRAY = analyse_drone_timeseries(drone, lead, method=method)
        (AVGO, AVGA, AVGD, RMSO, RMSA, RMSD) = ERROR_ARRAY
        (Odates, OSdates, dlat, dlon, TO) =  OCEAN_RESULTS
        (Adates, ASdates, dlat, dlon, TA) =  ATMOS_RESULTS
	ALL_TO.append(TO)
	ALL_TA.append(TA)
	ALL_ER.append(ERROR_ARRAY)
    ALL_TO = np.ma.array(ALL_TO)
    ALL_TA = np.ma.array(ALL_TA)
    ALL_ER = np.transpose(np.array(ALL_ER))
    print 'SHAPES', ALL_TO.shape, ALL_TA.shape, ALL_ER.shape
    Ofile=out_dir+'/TMO_'+drone_str+'_'+'all_leads'+'.dat'
    Afile=out_dir+'/TMA_'+drone_str+'_'+'all_leads'+'.dat'
    Odates_str, Odates_int = datadatefile.convert_datelist_strint(Odates)
    Adates_str, Adates_int = datadatefile.convert_datelist_strint(Adates)
    datadatefile.write_file(Odates_int, ALL_TO, file=Ofile)
    datadatefile.write_file(Adates_int, ALL_TA, file=Afile)
    Rfile=out_dir+'/ERR_'+drone_str+'.dat'   
    datafiles.write_nvarfile(argleads, ALL_ER, file=Rfile)
    
    return Rfile

def one_lead(drone=drones[0], lead=leads[0], method=0):

    drone_str=str(drone).zfill(4)
    lead_str=str(lead).zfill(3)
    out_dir=get_out_dir(method=method)
    OCEAN_RESULTS, ATMOS_RESULTS, ERROR_ARRAY = analyse_drone_timeseries(drone,lead,method=method)
    (AVGO, AVGA, AVGD, RMSO, RMSA, RMSD) = ERROR_ARRAY
    Rfile=out_dir+'/ERR_'+drone_str+'.dat' 
    nerr=len(ERROR_ARRAY)
    datafiles.test_and_create(nerr, leads, Rfile)
    datafiles.addto_nvarfile(lead, np.reshape(ERROR_ARRAY, (nerr, 1) ) , file=Rfile)
    return
    
def loop_drones(argdrones=drones,method=0):
    out_dir = get_out_dir(method=method)
    for drone in argdrones:
        drone_str=str(drone).zfill(4)
        Rfile = loop_leads(drone=drone, argleads=leads,method=method)
	post_drone(drone, method=method)
    return
    
def post_drone(drone, method=0):
    out_dir = get_out_dir(method=method)
    drone_str=str(drone).zfill(4)
    Rfile = out_dir+'/ERR_'+drone_str+'.dat'
    print 'Rfile', Rfile
    Ofile = out_dir+'/TERR_'+drone_str
    leadtimes, data = datafiles.read_nvarfile(file=Rfile)
    AVGO = data[0,:]
    AVGA = data[1,:]
    AVGD = data[2,:]
    RMSO = data[3,:]
    RMSA = data[4,:]
    RMSD = data[5,:]
    fig, ax = plt.subplots()
    ax.plot(leadtimes, RMSO, color='b', label='Ocean Grid')
    ax.plot(leadtimes, RMSA, color='r', label='Atmos Grid')
    ax.plot(leadtimes, RMSD, color='g', label='A - O Grid')
    ax.plot(leadtimes, AVGO, color='b', linestyle='--')
    ax.plot(leadtimes, AVGA, color='r', linestyle='--')
    ax.plot(leadtimes, AVGD, color='g', linestyle='--')
    ax.legend()
    fig.savefig(Ofile+'.png')
    fig.savefig(Ofile+'.pdf')
    plt.close(fig) 
    return
   
