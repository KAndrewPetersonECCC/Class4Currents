import numpy as np
import scipy.interpolate as si

missing = -999.

def grdfld(lon, lat, field, ddeg=0.25, method='nearest', central_longitude=0):
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

def cycle_lon(lon, central_longitude=0):
    new_lon = lon.copy()
    icycle = np.where(lon < -180.00 + central_longitude) 
    new_lon[icycle] = lon[icycle]+360
    icycle = np.where(lon > 180.0 + central_longitude)
    new_lon[icycle] = lon[icycle]-360
    return new_lon

