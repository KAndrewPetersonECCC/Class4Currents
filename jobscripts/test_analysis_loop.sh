#!/bin/bash
# ord_soumet /fs/homeu2/eccc/mrd/ords/rpnenv/dpe000/Class4_Currents/jobscripts/test_analysis_loop.sh -cpus 80 -mpi -cm 2000M -t 10800 -shell=/bin/bash
#
#bash /fs/homeu2/eccc/mrd/ords/rpnenv/dpe000/Class4_Currents/jobscripts/test_analysis_loop.sh
#

cd /home/dpe000/Class4_Currents
source jobscripts/prepython.sh

cd /home/dpe000/Class4_Currents
source jobscripts/prepython.sh

python << EOD
from importlib import reload
import os
import sys
sys.path.insert(0, '/home/dpe000/Class4_Currents/python')
import matplotlib as mpl
mpl.use('Agg')
import traceback
import datetime
import numpy as np
import datadatefile
import Class4Current


expt0 = 'GX340-U2-CONTROLE-LR'
variable='init_estimate'

dates = Class4Current.create_dates(20210105, 20220204)
expt = 'dev-3.4.0-gx-7dIAU'

dates = Class4Current.create_dates(20190402, 20200331)
expt='gx_dev-3.4.0_DA_y2_MDT2' ; label='MDT'
#expt='gx_dev-3.4.0_DA_y2_core' ; label='core'
expt0 = 'gx_dev-3.4.0_DA_y2'

test1=False
test2=True

if test1:
  (LONG, LATG), MNS_FIELDS, GDS_FIELDS, GCS, CNS = Class4Current.loop_analysis_dates(dates, expt, variable, insuffix='GIOPS_orca025_currents.f2', ierror=1, ddeg=4.0, mp=True)

  MNA_FIELDS = Class4Current.calc_mean_field(MNS_FIELDS, GCS)
  GDA_FIELDS = Class4Current.calc_mean_field(GDS_FIELDS, CNS)

  sum_obs, sum_fld,  var_obs, var_fld, cov_obs_fld = Class4Current.get_correlation_fields(MNS_FIELDS)
  grd_sum_obs, grd_sum_fld, grd_var_obs, grd_var_fld, grd_cov_obs_fld = Class4Current.get_correlation_fields(GDA_FIELDS)

  rcorr = Class4Current.calc_correlation(GCS, sum_obs, sum_fld, cov_obs_fld, var_obs, var_fld)
  grd_rcorr = Class4Current.calc_correlation(CNS, grd_sum_obs, grd_sum_fld, grd_cov_obs_fld, grd_var_obs, grd_var_fld)

  print(rcorr)
  print(np.mean(grd_rcorr))

if test2:
  Class4Current.compare_analysis_errors(dates, [expt, expt0], [label, 'ctrl'], maxtaylor=1.0)
EOD
