#!/bin/bash
# ord_soumet /fs/homeu2/eccc/mrd/ords/rpnenv/dpe000/Class4_Currents/jobscripts/test_Class4_ENAN.sh -cpus 80 -mpi -cm 2000M -t 10800 -shell=/bin/bash

WDIR=/fs/homeu2/eccc/mrd/ords/rpnenv/dpe000/Class4_Currents

cd ${WDIR}
## Adding NCO operations
. ssmuse-sh -d eccc/cmd/cmds/ext/20220331
## Adding RPNPY
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
import Class4CurrentEA

Class4CurrentEA.process_enan_obs(date=datetime.datetime(2022,1,20), ens_list=[0])

EOD
