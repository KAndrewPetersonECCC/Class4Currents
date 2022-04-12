import matplotlib.pyplot as plt
import matplotlib.colors as clr
import cartopy.crs as ccrs
import scipy.interpolate as si
import scipy.stats as ss

import numpy as np
import math

plt.rc('text', usetex=True)

dmap = 'RdBu_r'
missing = -999.

def make_projections(**kwargs):
    projections = {}
    projections['Mollweide'] = ccrs.Mollweide(**kwargs)
    projections['PlateCarree'] = ccrs.PlateCarree(**kwargs)
    projections['PacificCarree'] = ccrs.PlateCarree(central_longitude=180)
    projections['Mercator'] = ccrs.Mercator()
    projections['Percator'] = ccrs.Mercator(central_longitude=180)
    projections['Miller'] = ccrs.Miller()
    projections['NorthPolarStereo'] = ccrs.NorthPolarStereo(**kwargs)
    projections['SouthPolarStereo'] = ccrs.SouthPolarStereo(**kwargs)
    projections['Orthographic'] = ccrs.Orthographic(**kwargs)
    projections['EqualEarth'] = ccrs.EqualEarth(**kwargs)
    pcarree=ccrs.PlateCarree()
    return projections, pcarree

def pcolormesh(lon, lat, field, levels=None, ticks=None, cmap=dmap, project='PlateCarree', outfile='plt.png', 
               box=[-180, 180, -90, 90], make_global=False, title='', suptitle=None, 
	       cbar=True, obar='vertical', fontsizes=None, **kwargs):

    title_fontsize = fontsizes
    cbar_fontsize = fontsizes
    fontsize=None
    if ( isinstance(fontsizes, int) or isinstance(fontsizes, float) ):
        fontsize=fontsizes
    if ( not isinstance(fontsizes, type(None)) ):
        if ( len(fontsizes) == 2 ):
	    title_fontsize = fontsizes[0]
	    cbar_fontsize = fontsizes[1]
	    fontsize=fontsizes[0]
    #print(fontsize, title_fontsize, cbar_fontsize)
    fig = plt.figure()
    projections, pcarree = make_projections(**kwargs)
    ax = plt.subplot(projection=projections[project])

    if ( isinstance(levels, type(None) ) ):
        norm=None
    else:
        Ncolors=plt.get_cmap(cmap).N
        norm = clr.BoundaryNorm(levels, ncolors=Ncolors, clip=True)

    if ( make_global ):
        ax.set_global()
    else:
        ax.set_extent(box, crs=ccrs.PlateCarree())
    if ( isinstance(norm, type(None) ) ):
        mesh = ax.pcolormesh(lon, lat, field, cmap=cmap,transform=pcarree)
    else:
        mesh = ax.pcolormesh(lon, lat, field, norm=norm, cmap=cmap,transform=pcarree)
    mesh.set_cmap(cmap)
    if ( cbar ): 
        #print('cbar', cbar_fontsize)
        cbar_fig=fig.colorbar(mesh, format='%.2f',orientation=obar, ticks=ticks)
	cbar_fig.ax.tick_params(labelsize=cbar_fontsize)
    ax.coastlines()
    #print('title', title_fontsize)
    if ( suptitle != None ): fig.suptitle(suptitle, fontsize=title_fontsize)
    ax.set_title(title, fontsize=title_fontsize)
    if ( isinstance(outfile, str) ):  outfile_list = [ outfile ]
    if ( isinstance(outfile, list) ): outfile_list = outfile
    for ofile in outfile_list:
        fig.savefig(ofile,bbox_inches='tight')
    plt.close(fig)
    return

def cycle_lon(lon, central_longitude=0):
    new_lon = lon.copy()
    icycle = np.where(lon < -180.00 + central_longitude) 
    new_lon[icycle] = lon[icycle]+360
    icycle = np.where(lon > 180.0 + central_longitude)
    new_lon[icycle] = lon[icycle]-360
    return new_lon

def grdfld(lon, lat, field, ddeg=0.5, method='nearest', central_longitude=0):
    regu_lon = np.arange(-180.0+ddeg/2.0, 180.0, ddeg)+central_longitude
    regu_lat = np.arange(-90.0+ddeg/2.0, 90.0, ddeg)
    grid_lat, grid_lon = np.meshgrid(regu_lat, regu_lon)

    lon = cycle_lon(lon, central_longitude=central_longitude)

    if ( method[0:6] == '2sweep' ):
        first=method[6:]
        grid_fld = si.griddata( (lon.flatten(),lat.flatten()), field.flatten(), (grid_lon, grid_lat), method='linear', fill_value=missing)
	grid_fld2 = si.griddata( (lon.flatten(),lat.flatten()), field.flatten(), (grid_lon, grid_lat), method='nearest')
	replace = np.where(grid_fld == missing)
	grid_fld[replace] = grid_fld2[replace]
    else:
        grid_fld = si.griddata( (lon.flatten(),lat.flatten()), field.flatten(), (grid_lon, grid_lat), method=method)
    
    return grid_lon, grid_lat, grid_fld

def binfld(lon, lat, field, ddeg=2.0,central_longitude=0):
    lon_bin=np.arange(-180,180+ddeg,ddeg)+central_longitude
    lat_bin=np.arange(-90, 90+ddeg, ddeg)
    lon_mid = (lon_bin[1:] + lon_bin[:-1]) / 2
    lat_mid = (lat_bin[1:] + lat_bin[:-1]) / 2
    grid_lat, grid_lon = np.meshgrid(lat_mid, lon_mid) 
    lon = cycle_lon(lon, central_longitude=central_longitude)
    grid_fld, xedges, yedges, binnumber = ss.binned_statistic_2d(lon.flatten(), lat.flatten(), values=field.flatten(), 
       statistic='mean', bins=[lon_bin, lat_bin])
       
    return grid_lon, grid_lat, grid_fld

def binfldsum(lon, lat, field, ddeg=2.0, central_longitude=0):
    lon_bin=np.arange(-180,180+ddeg,ddeg)+central_longitude
    lat_bin=np.arange(-90, 90+ddeg, ddeg)
    lon_mid = (lon_bin[1:] + lon_bin[:-1]) / 2
    lat_mid = (lat_bin[1:] + lat_bin[:-1]) / 2
    grid_lat, grid_lon = np.meshgrid(lat_mid, lon_mid) 
    lon = cycle_lon(lon, central_longitude=central_longitude)
    grid_sum, xedges, yedges, binnumber = ss.binned_statistic_2d(lon.flatten(), lat.flatten(), values=field.flatten(), 
       statistic='sum', bins=[lon_bin, lat_bin])
    grid_cnt, xedges, yedges, binnumber = ss.binned_statistic_2d(lon.flatten(), lat.flatten(), values=field.flatten(), 
       statistic='count', bins=[lon_bin, lat_bin])
       
    return grid_lon, grid_lat, grid_sum, grid_cnt

def make_bin_grid(ddeg=2, central_longitude=0):
    lon_bin=np.arange(-180,180+ddeg,ddeg)+central_longitude
    lat_bin=np.arange(-90, 90+ddeg, ddeg)
    lon_mid = (lon_bin[1:] + lon_bin[:-1]) / 2
    lat_mid = (lat_bin[1:] + lat_bin[:-1]) / 2
    grid_lat, grid_lon = np.meshgrid(lat_mid, lon_mid) 
    grid_sum = np.zeros(grid_lat.shape)
    grid_cnt = np.zeros(grid_lat.shape)
    return grid_lon, grid_lat, lon_bin, lat_bin, grid_sum, grid_cnt
    
def binfldsumcum(lon, lat, field, lon_bin, lat_bin, grid_sum, grid_cnt, central_longitude=0):
    lon = cycle_lon(lon, central_longitude=central_longitude)
    grid_sum_add, xedges, yedges, binnumber = ss.binned_statistic_2d(lon.flatten(), lat.flatten(), values=field.flatten(), 
       statistic='sum', bins=[lon_bin, lat_bin])
    grid_cnt_add, xedges, yedges, binnumber = ss.binned_statistic_2d(lon.flatten(), lat.flatten(), values=field.flatten(), 
       statistic='count', bins=[lon_bin, lat_bin])
    grid_sum_new = grid_sum + grid_sum_add
    grid_cnt_new = grid_cnt + grid_cnt_add
    return grid_sum_new, grid_cnt_new

def binfldsumFIN(grid_sum, grid_cnt):
    grid_msk=np.ones(grid_sum.shape).astype(bool)
    grid_plt=np.ma.array(grid_sum.copy(),mask=grid_msk)
    iplot = np.where(grid_cnt > 0)
    grid_plt[iplot].mask = False
    grid_plt[iplot] = grid_sum[iplot] / grid_cnt[iplot]
    return grid_plt
        
def grd_pcolormesh(lon, lat, field, levels=None, ticks=None,
                   cmap=dmap, project='PlateCarree', 
		   outfile='plt.png', box=[-180, 180, -90, 90], 
		   make_global=False, title='', suptitle=None, ddeg=0.5, 
		   cbar=True, interp_method='nearest', fontsizes=None,
		   **kwargs
		   ):

    central_longitude=0.0
    if ( ( project == 'PacificCarree' ) or ( project == 'Percator' ) ): central_longitude=180.0
    grid_lon, grid_lat, grid_fld = grdfld(lon, lat, field, ddeg=ddeg, method=interp_method,central_longitude=central_longitude)
    pcolormesh(grid_lon, grid_lat, grid_fld, levels=levels, ticks=ticks, cmap=cmap, project=project, outfile=outfile, 
               box=box, make_global=make_global, title=title, suptitle=suptitle, cbar=cbar, fontsizes=fontsizes, **kwargs)
    return
    
def bin_pcolormesh(lon, lat, field, levels=None, ticks=None, 
                   cmap=dmap, project='PlateCarree', 
		   outfile='plt.png', box=[-180, 180, -90, 90], 
		   make_global=False, title='',  suptitle=None, ddeg=2.0, cbar=True, obar='vertical', fontsizes=None,
		   **kwargs
		   ):

    grid_lon, grid_lat, grid_fld = binfld(lon, lat, field, ddeg=ddeg)
    pcolormesh(grid_lon, grid_lat, grid_fld, levels=levels, ticks=ticks, cmap=cmap, project=project, outfile=outfile, 
               box=box, make_global=make_global, title=title, suptitle=suptitle, 
	       cbar=cbar, obar=obar,fontsizes=fontsizes
	       **kwargs)
    return
   
def mask_field(lon, lat, fld, mask):
    imask = np.where(mask == 1)
    msk_lon = lon[imask]
    msk_lat = lat[imask]
    msk_fld = fld[imask]
    return msk_lon, msk_lat, msk_fld
    
def msk_grd_pcolormesh(lon, lat, field,  mask, levels=None, ticks=None, 
                   cmap=dmap, project='PlateCarree', 
		   outfile='plt.png', box=[-180, 180, -90, 90], 
		   make_global=False, title='', suptitle=None, ddeg=0.5, cbar=True, obar='vertical',
		   addmask=False, interp_method='nearest', fontsizes=None,
		   **kwargs
		   ):

    mask_lon, mask_lat, mask_fld = mask_field(lon, lat, field, mask)
    central_longitude=0.0
    if ( ( project == 'PacificCarree' ) or ( project == 'Percator' ) ): central_longitude=180.0
    grid_lon, grid_lat, grid_fld = grdfld(mask_lon, mask_lat, mask_fld, ddeg=ddeg,method=interp_method,central_longitude=central_longitude)
    grid_mlon, grid_mlat, grid_msk = grdfld(lon, lat, mask, ddeg=ddeg, method='nearest',central_longitude=central_longitude)  ## ONLY NEAREST NEIGHBOUR MAKES SENSE HERE?
    if ( not addmask ):
        pcolormesh(grid_lon, grid_lat, grid_fld, levels=levels, ticks=ticks, cmap=cmap, project=project, outfile=outfile, 
	           box=box, make_global=make_global, title=title, suptitle=suptitle, cbar=cbar,
		   obar=obar, fontsizes=fontsizes, **kwargs)
    else:
        pcolormesh(grid_lon, grid_lat, np.ma.array(grid_fld, mask=np.logical_not(grid_msk)),
	           levels=levels, ticks=ticks, cmap=cmap, project=project, outfile=outfile, 
		   box=box, make_global=make_global, title=title, suptitle=suptitle, cbar=cbar, 
		   obar=obar, fontsizes=fontsizes,**kwargs)
    return
    
def msk_bin_pcolormesh(lon, lat, field, mask, levels=None, ticks=None, 
                   cmap=dmap, project='PlateCarree', 
		   outfile='plt.png', box=[-180, 180, -90, 90], 
		   make_global=False, title='', suptitle=None, ddeg=2.0, cbar=True, obar='vertial', fontsizes=None, **kwargs):

    mask_lon, mask_lat, mask_fld = mask_field(lon, lat, field, mask)
    grid_lon, grid_lat, grid_fld = binfld(mask_lon, mask_lat, mask_fld, ddeg=ddeg)
    pcolormesh(grid_lon, grid_lat, grid_fld, levels=levels, ticks=ticks, cmap=cmap, project=project, 
               outfile=outfile, box=box, make_global=make_global, title=title, suptitle=suptitle, 
	       cbar=cbar,obar=obar, fontsizes=fontsizes, **kwargs)
    return
   
def quiver(lon, lat, ufield, vfield, levels=None, cmap=dmap, project='PlateCarree', outfile='plt.png', 
           box=[-180, 180, -90, 90], make_global=False, title='', suptitle=None, 
	   cbar=True, obar='vertical', Nskip=1,**kwargs):

    projections,pcarree = make_projections(**kwargs)

    field = np.sqrt(np.square(ufield)+np.square(vfield))
    skip=(slice(None,None,Nskip),slice(None,None,Nskip))    
    fig = plt.figure()
    ax = plt.subplot(projection=projections[project])

    if ( isinstance(levels, type(None) ) ):
        norm=None
    else:
        Ncolors=plt.get_cmap(cmap).N
        norm = clr.BoundaryNorm(levels, ncolors=Ncolors, clip=True)

    if ( make_global ):
        ax.set_global()
    else:
        ax.set_extent(box, crs=ccrs.PlateCarree())
    if ( isinstance(norm, type(None) ) ):
        mesh = ax.pcolormesh(lon, lat, field, cmap=cmap,transform=pcarree)
    else:
        mesh = ax.pcolormesh(lon, lat, field, norm=norm, cmap=cmap,transform=pcarree)
    mesh.set_cmap(cmap)
    quiver = ax.quiver(lon[skip], lat[skip], ufield[skip], vfield[skip])
    if ( cbar ): fig.colorbar(mesh, format='%.2f',orientation=obar)
    ax.coastlines()
    ax.set_title(title)
    if ( suptitle != None ): fig.suptitle(suptitle)
    fig.savefig(outfile,bbox_inches='tight')
    plt.close(fig)
    return
        
def scatter(lon, lat, field, levels=None, ticks=None, cmap=dmap, project='PlateCarree', outfile='plt.png', 
               box=[-180, 180, -90, 90], make_global=False, title='', suptitle=None, 
	       cbar=True, obar='vertical', fontsizes=None, **kwargs):

    title_fontsize = fontsizes
    cbar_fontsize = fontsizes
    fontsize=None
    if ( isinstance(fontsizes, int) or isinstance(fontsizes, float) ):
        fontsize=fontsizes
    if ( not isinstance(fontsizes, type(None)) ):
        if ( len(fontsizes) == 2 ):
	    title_fontsize = fontsizes[0]
	    cbar_fontsize = fontsizes[1]
	    fontsize=fontsizes[0]
    #print(fontsize, title_fontsize, cbar_fontsize)
    fig = plt.figure()
    projections, pcarree = make_projections(**kwargs)
    ax = plt.subplot(projection=projections[project])

    if ( isinstance(levels, type(None) ) ):
        norm=None
    else:
        Ncolors=plt.get_cmap(cmap).N
        norm = clr.BoundaryNorm(levels, ncolors=Ncolors, clip=True)

    if ( make_global ):
        ax.set_global()
    else:
        ax.set_extent(box, crs=ccrs.PlateCarree())
    if ( isinstance(norm, type(None) ) ):
        #mesh = ax.scatter(lon.flatten(), lat.flatten(), field.flatten(), cmap=cmap,transform=pcarree)
        scat = ax.scatter(x=lon, y=lat, c=field, s=5, alpha=0.5, transform=ccrs.PlateCarree(), cmap=cmap,) ## Important
    else:
        #mesh = ax.scatter(lon, lat, field, norm=norm, cmap=dmap,transform=pcarree)
        scat = ax.scatter(x=lon.flatten(), y=lat.flatten(), c=field.flatten(), s=5, alpha=0.5, transform=ccrs.PlateCarree(), cmap=cmap, norm=norm) ## Important
    scat.set_cmap(cmap)
    if ( cbar ): 
        #print('cbar', cbar_fontsize)
        cbar_fig=fig.colorbar(scat, format='%.2f',orientation=obar, ticks=ticks)
	cbar_fig.ax.tick_params(labelsize=cbar_fontsize)
    ax.coastlines()
    #print('title', title_fontsize)
    if ( suptitle != None ): fig.suptitle(suptitle, fontsize=title_fontsize)
    ax.set_title(title, fontsize=title_fontsize)
    if ( isinstance(outfile, str) ):  outfile_list = [ outfile ]
    if ( isinstance(outfile, list) ): outfile_list = outfile
    for ofile in outfile_list:
        fig.savefig(ofile,bbox_inches='tight')
    plt.close(fig)
    return
