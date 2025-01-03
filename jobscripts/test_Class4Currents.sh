#!/bin/bash
# ord_soumet /fs/homeu2/eccc/mrd/ords/rpnenv/dpe000/Class4_Currents/jobscripts/test_Class4Currents.sh  -cpus 80 -cm 180000M -t 21600 -shell=/bin/bash

TDIR=/fs/homeu2/eccc/mrd/ords/rpnenv/dpe000/Class4_Currents
cd ${TDIR}
if [[ $? -ne 0 ]] ; then 
    echo "CD TDIR NO WORK"
    exit 99
fi

source jobscripts/prepython.sh

python3 << EOD 

import os
import sys
sys.path.insert(0, '/home/dpe000/Class4_Currents/python')
import matplotlib as mpl
mpl.use('Agg')
import traceback
import datetime
import time

import Class4Current
import Class4CurrentEA

time0=time.time(); Class4CurrentEA.process_enan_obs(date=datetime.datetime(2022,1,1), ens_list=[0], filter=True, expt='GEPS_STO2X', nshapiro=30), print('PROCESS TIME ', time.time()-time0)

time0=time.time(); Class4CurrentEA.process_enan_obs(date=datetime.datetime(2022,2,1), ens_list=[0], filter=True, expt='GEPS_STO2X', nshapiro=30), print('PROCESS TIME ', time.time()-time0)

EOD
