import os
import sys
sys.path.insert(0, '/home/dpe000/Class4_Currents/python')

import matplotlib.pyplot as plt
import matplotlib.colors as clr
import cartopy.crs as ccrs
import scipy.interpolate as si
import scipy.stats as ss
import cartopy.feature as cfeature

import datetime
import numpy as np

import Class4Current
import cplot

cmap_full_field='gist_stern_r'
cmap_anom_field='seismic'
cmap_anom_field='RdYlBu_r'


cmap=cmap_anom_field
levels=np.arange(-1.0, 1.1, 0.1)

thisdate=datetime.datetime(2020,02,01)
while ( thisdate < datetime.datetime(2020,02,29) ):
  datestr=thisdate.strftime('%Y%m%d')  
  obsfile='CLASS4_currents_CCMEP_FILT/class4_'+datestr+'_GIOPS_orca025_currents.f1.nc'
  (LONO, LATO, depth), (obser, beste, inite, fcstv, persi, nersi) = Class4Current.read_obsfile_plus(obsfile)
  thisdate=thisdate+datetime.timedelta(days=1)
  print datestr, np.min(fcstv)

thisdate=datetime.datetime(2020,2,25)  
datestr=thisdate.strftime('%Y%m%d')  
obsfile='CLASS4_currents_CCMEP_FILT/class4_'+datestr+'_GIOPS_orca025_currents.f1.nc'
(LONO, LATO, depth), (obser, beste, inite, fcstv, persi, nersi) = Class4Current.read_obsfile_plus(obsfile)

Mobser = Class4Current.speed_and_angle_easy(obser)
Mbeste = Class4Current.speed_and_angle_easy(beste)
Minite = Class4Current.speed_and_angle_easy(inite)
Mfcstv = Class4Current.speed_and_angle_easy(fcstv)
Mpersi = Class4Current.speed_and_angle_easy(persi)
Mnersi = Class4Current.speed_and_angle_easy(nersi)

Merror = Mobser - Mbeste

fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(1,1,1, projection=ccrs.PlateCarree())
ax.set_global()
ax.add_feature(cfeature.COASTLINE, edgecolor="green")
ax.gridlines()

Ncolors=plt.get_cmap(cmap).N
norm = clr.BoundaryNorm(levels, ncolors=Ncolors, clip=True)

plt.scatter(x=LONO, y=LATO, c=Merror[:,0,0], s=5, alpha=0.5, transform=ccrs.PlateCarree(), cmap=cmap, norm=norm) ## Important
plt.colorbar(orientation='horizontal')

plt.show()

grid_lon, grid_lat, lon_bin, lat_bin, grid_sum, grid_cnt = cplot.make_bin_grid(ddeg=5)
thisdate=datetime.datetime(2020,03,01)
while ( thisdate < datetime.datetime(2020,03,31) ):
  datestr=thisdate.strftime('%Y%m%d')  
  obsfile='CLASS4_currents_CCMEP_FILT/class4_'+datestr+'_GIOPS_orca025_currents.f1.nc'
  (LONO, LATO, depth), (obser, beste, inite, fcstv, persi, nersi) = Class4Current.read_obsfile_plus(obsfile)
  thisdate=thisdate+datetime.timedelta(days=1)
  print datestr
  Mobser = Class4Current.speed_and_angle_easy(obser)
  Mbeste = Class4Current.speed_and_angle_easy(beste)
  Merror = Mobser - Mbeste
  Serror = Merror[:,0,0] 
  grid_sum, grid_cnt = cplot.binfldsumcum(LONO, LATO, Serror, lon_bin, lat_bin, grid_sum, grid_cnt)
  print(np.sum(grid_cnt), len(obser))

grid_plt = cplot.binfldsumFIN(grid_sum, grid_cnt)
cplot.pcolormesh(grid_lon, grid_lat, grid_plt, outfile='bin.png', levels=levels, cmap=cmap)

grid_lon, grid_lat, grid_plt1 = cplot.binfld(LONO, LATO, Serror, ddeg=5)
cplot.pcolormesh(grid_lon, grid_lat, grid_plt1, outfile='bin1.png', levels=levels, cmap=cmap)

