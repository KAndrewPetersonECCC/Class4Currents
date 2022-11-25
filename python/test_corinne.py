#from importlib import reload
import sys
sys.path.insert(0, '/home/dpe000/Class4_Currents/python')
import read_rotate_velocities
import cst_interp
import numpy as np

TMPDIR='/fs/site5/eccc/mrd/rpnenv/dpe000/tmpdir/GEPS_TMP'

# Forecasts usually disappear from operational directory after they can be verified (i.e. by their validity date)
# Date below should be as slose to realtime as possible.
test_file = '/fs/site5/eccc/prod/ops/suites/geps_20220621/e1/gridpt/prog/ens.glboce/2022110700_384_000'
dest_file='/fs/site5/eccc/mrd/rpnenv/dpe000/tmpdir/2022110700_384_000.z25'

date = 2022111000
lhr = 24

# Control Member
ens=0
lon, lat, UWC, VWC = read_rotate_velocities.read_velocities_timestamp(date, lhr, ens, ip1_str='sfc')    

# Full Ensemble
lon, lat, UW_list, VW_list = read_rotate_velocities.read_ensemble_velocities(date, lhr, enslist=list(range(21)), ip1_str='sfc')

# ROTATE VELOCITIES to East and North
UE_list, VN_list = read_rotate_velocities.UU_list_ORCA_to_NE(UW_list, VW_list)

# GRID VELOCITIES ON 0.25 grid
lon_grid, lat_grid, UE_grdlist = read_rotate_velocities.grd_list(lon, lat, UE_list)
lon_grid, lat_grid, VN_grdlist = read_rotate_velocities.grd_list(lon, lat, VN_list)
lon_grir, lat_grir, UW_grdlist = read_rotate_velocities.grd_list(lon, lat, UW_list)   # UNROTATED (for check)
lon_grir, lat_grir, VW_grdlist = read_rotate_velocities.grd_list(lon, lat, VW_list)   # UNROTATED (for check)

# CALCULATE ENSEMBLE MEAN AND VARIANCE
# ON ORIGINAL GRID    
UE_ENM, UE_VAR = read_rotate_velocities.ens_mean_variance(UE_list)
VN_ENM, VN_VAR = read_rotate_velocities.ens_mean_variance(VN_list)
# ON 0.25 GRID    
UG_ENM, UG_VAR = read_rotate_velocities.ens_mean_variance(UE_grdlist)
VG_ENM, VG_VAR = read_rotate_velocities.ens_mean_variance(VN_grdlist)
UW_ENM, UW_VAR = read_rotate_velocities.ens_mean_variance(UW_grdlist)   # UNROTATED (for check)
VW_ENM, VW_VAR = read_rotate_velocities.ens_mean_variance(VW_grdlist)   # UNROTATED (for check)

# OR JUST USE CSTINTRP to rotate velocities and place on 0.25 grid

#Assuming this run correctly (You need a tmpdir to place the file in for one).
read_rotate_velocities.create_interp_files(date, lhr, tmpdir=TMPDIR)  

#WHICH YOU CAN THEN SIMPLY READ IN
lon_cstin, lat_cstin, UI_list, VI_list = read_rotate_velocities.read_ensemble_velocities(date, lhr, dir=TMPDIR, suffix='.z25')
lat_cstgr, lon_cstgr = np.meshgrid(lat_cstin, lon_cstin)


# AND PRODUCE AN ENSMEBLE MEAL
UI_ENM, UI_VAR = read_rotate_velocities.ens_mean_variance(UI_list)
VI_ENM, VI_VAR = read_rotate_velocities.ens_mean_variance(VI_list)

import cplot
import numpy as np

lon_cstgr = cplot.cycle_lons(cstgr)
cplot.pcolormesh(lon_grid, lat_grid, UG_ENM, project='NorthPolarStereo', outfile='UM_GRID.png', obar='horizontal', box=[-180, 180, 50, 90], levels=np.arange(-0.195, 0.2, 0.01) )
cplot.pcolormesh(lon_grid, lat_grid, VG_ENM, project='NorthPolarStereo', outfile='VM_GRID.png', obar='horizontal', box=[-180, 180, 50, 90], levels=np.arange(-0.195, 0.2, 0.01) )

cplot.pcolormesh(lon, lat, UE_ENM, project='NorthPolarStereo', outfile='UM_NEMO.png', obar='horizontal', box=[-180, 180, 50, 90], levels=np.arange(-0.195, 0.2, 0.01) )
cplot.pcolormesh(lon, lat, VN_ENM, project='NorthPolarStereo', outfile='VM_NEMO.png', obar='horizontal', box=[-180, 180, 50, 90], levels=np.arange(-0.195, 0.2, 0.01) )

cplot.pcolormesh(lon_cstgr, lat_cstgr, UI_ENM, project='NorthPolarStereo', outfile='UM_CSTI.png', obar='horizontal', box=[-180, 180, 50, 90], levels=np.arange(-0.195, 0.2, 0.01) )
cplot.pcolormesh(lon_cstgr, lat_cstgr, VI_ENM, project='NorthPolarStereo', outfile='VM_CSTI.png', obar='horizontal', box=[-180, 180, 50, 90], levels=np.arange(-0.195, 0.2, 0.01) )

#FINALLY GRID UNROTATED VELOCITIES (SEE IF YOU CAN SPOT THE DIFFERENCE)
cplot.pcolormesh(lon_grid, lat_grid, UW_ENM, project='NorthPolarStereo', outfile='UM_ORIG.png', obar='horizontal', box=[-180, 180, 50, 90], levels=np.arange(-0.195, 0.2, 0.01) )
cplot.pcolormesh(lon_grid, lat_grid, VW_ENM, project='NorthPolarStereo', outfile='VM_ORIG.png', obar='horizontal', box=[-180, 180, 50, 90], levels=np.arange(-0.195, 0.2, 0.01) )

# THIS MIGHT SEEM ODD TO DO -- BUT UNFORTUNATELY the two grids are unaligned by 0.125 deg -- as well as one begin 0-360 and the other -180-180.
lon_grih, lat_grih, UH_ENM = read_rotate_velocities.grd_fld(lon_cstgr, lat_cstgr, UI_ENM)
lon_grih, lat_grih, VH_ENM = read_rotate_velocities.grd_fld(lon_cstgr, lat_cstgr, VI_ENM)

# BUT NOW THEY CAN BE COMPARED DIRECTLY.
SMALLE = np.sum( np.square(UH_ENM-UG_ENM) + np.square(VH_ENM-VG_ENM) )
ROTATE = np.sum( np.square(UH_ENM-UW_ENM) + np.square(VH_ENM-VW_ENM) )

print(SMALLE, ROTATE)
