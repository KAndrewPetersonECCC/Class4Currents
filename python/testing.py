import os
import sys
sys.path.insert(0, '/home/dpe000/Class4_Currents/python')

import netCDF4
import shutil 
import numpy as np
import datetime

import find_fcst_file
import stfd
import read_grid
import isoheatcontent
import bering_diurnal

oile='CLASS4_currents/class4_20191130_PSY4V3R1_orca12_currents.nc'
file='CLASS4_currents_CCMEP/class4_20191130_GIOPS_orca025_currents.nc'
oile='CLASS4_currents/class4_20211031_PSY4V3R1_orca12_currents.nc'
file='CLASS4_currents_CCMEP/class4_20211031_GIOPS_orca025_currents.nc'

mask = read_grid.read_mask(var='tmask')
maskt = read_grid.read_mask(var='tmask')
masku = read_grid.read_mask(var='umask')
maskv = read_grid.read_mask(var='vmask')
e3t = read_grid.read_e3t_mesh(var='e3t_0')
e3u = read_grid.read_e3t_mesh(var='e3u_0')
e3v = read_grid.read_e3t_mesh(var='e3v_0')
e1t = read_grid.read_mesh_var('e1t')
e2t = read_grid.read_mesh_var('e2t')
e1u = read_grid.read_mesh_var('e1u')
e2u = read_grid.read_mesh_var('e2u')
e1v = read_grid.read_mesh_var('e1v')
e2v = read_grid.read_mesh_var('e2v')
THE = read_grid.read_angle('T')
#shutil.copy(oile, file)

fila='/space/hall4/sitestore/eccc/mrd/rpnenv/socn000/env_ubuntu-18.04-skylake-64/datafiles/constants/oce/repository/master/CONCEPTS/orca025/grids/orca025grid_new.std'
LONA, LATA, TH = stfd.read_fstd_var(fila, 'LAAN')

obsset = netCDF4.Dataset(file)

def UU_ORCA_to_NE(UV):
    fila='/space/hall4/sitestore/eccc/mrd/rpnenv/socn000/env_ubuntu-18.04-skylake-64/datafiles/constants/oce/repository/master/CONCEPTS/orca025/grids/orca025grid_new.std'
    LONA, LATA, TH = stfd.read_fstd_var(fila, 'LAAN')
    UO, VO = UV
    UE = UO * np.cos(TH) - VO * np.sin(TH)
    VN = UO * np.sin(TH) + VO * np.cos(TH)
    return [UE, VN]

def calc_u15(UVin, e3tin, maskin):
    UUin, VVin = UVin
    u_20, __ = isoheatcontent.depth_integral(UUin, e3tin, maskin, depth=20)
    u_10, __ = isoheatcontent.depth_integral(UUin, e3tin, maskin, depth=10)
    u15 = (u_20 - u_10)/10

    v_20, __ = isoheatcontent.depth_integral(VVin, e3tin, maskin, depth=20)
    v_10, __ = isoheatcontent.depth_integral(VVin, e3tin, maskin, depth=10)
    v15 = (v_20 - v_10)/10
    
    return [u15, v15]
    
file_fcst=find_fcst_file.find_fcst_file('OPD', datetime.datetime(2021,10,30), 24, 0, src='ocn', execute=True)
LONN, LATN, UW = stfd.read_fstd_var(file_fcst, 'UUW', typvar='P@')
__, __, VW = stfd.read_fstd_var(file_fcst, 'VVW', typvar='P@')
__, __, TW = stfd.read_fstd_var(file_fcst, 'TM', typvar='P@')

leu, UU = stfd.read_fstd_multi_lev(file_fcst, 'UUW',vfreq=24, typvar='P@')
lev, VV = stfd.read_fstd_multi_lev(file_fcst, 'VVW',vfreq=24, typvar='P@')
let, TM = stfd.read_fstd_multi_lev(file_fcst, 'TM',vfreq=24, typvar='P@')

[u15, v15] = calc_u15([UU, VV], e3t, mask]
u15e, v15e = UU_ORCA_to_NE([u15, v15])

LON=obsset.variables['longitude']
LAT=obsset.variables['latitude']

file_interp='/fs/site4/eccc/mrd/rpnenv/dpe000//tmpdir/OPD/2021103000_024.dest.1'
lfu, UF = stfd.read_fstd_multi_lev(file_interp, 'UUW',vfreq=24, typvar='P@')
lfv, VF = stfd.read_fstd_multi_lev(file_interp, 'VVW',vfreq=24, typvar='P@')
lft, TF = stfd.read_fstd_multi_lev(file_interp, 'TM',vfreq=24, typvar='P@')
gridf='/fs/site4/eccc/mrd/rpnenv/dpe000/tmpdir/OPD/WGT/dest_grid.std'
LONF, LATF = stfd.read_latlon(gridf, nomlat='lat', nomlon='lon')

def fine_grid_mask(FD, e3tin):
    maskf=np.ones(FD.shape)

    e3t0=np.max(e3tin, axis=(1,2))
    e3tf=np.zeros(FD.shape)

    nzf, nxf, nyf = e3tf.shape
    for ix in range(nxf):
        for iy in range(nyf):
            e3tf[:,ix,iy] = e3t0
	    
    return maskf, e3tf

maskf, e3tf = find_grid_mask(TF, e3t)
    
u15f, v15f = calc_u15([UF, VF]. e3tf, maskf)

LONPT, LATPT = LON[0], LAT[0]
ipt, jpt = bering_diurnal.find_nearest_point(LONPT, LATPT, LONN, LATN)
ipf, jpf = bering_diurnal.find_nearest_point(LONPT, LATPT, LONF, LATF)

TMpt = TM[0,ipt, jpt]
TMpi = bering_diurnal.interpolate_to_point(TM[0,:,:], LONN, LATN, LONPT, LATPT)
TMpf = TF[0,ipf, jpf]
TMpg = bering_diurnal.interpolate_to_point(TF[0,:,:], LONF, LATF, LONPT, LATPT)

print TMpt, TMpi, TMpf, TMpg
