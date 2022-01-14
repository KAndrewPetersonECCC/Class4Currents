#!/bin/bash
# ord_soumet /fs/homeu1/eccc/mrd/ords/rpnenv/dpe000/Class4_Currents/jobscripts/test_Class4Currents.sh -cpus 1 -mpi -cm 64000M -t 7200 -shell=/bin/bash

WDIR=/fs/homeu1/eccc/mrd/ords/rpnenv/dpe000/Class4_Currents
cd ${WDIR}

source jobscripts/preconda.sh
source activate metcarto

python python/Class4Current.py

