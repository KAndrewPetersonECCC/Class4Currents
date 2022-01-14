#!/bin/bash

# ord_soumet /home/dpe000/Class4_Currents/jobscripts/perform_interpolation.sh -cpus 1 -mpi -cm 64000M -t 3600 -shell=/bin/bash

. ssmuse-sh -d /fs/ssm/eccc/mrd/rpn/OCEAN/cstint-3.2.7

GRID_DIR=/fs/site4/eccc/mrd/rpnenv/dpe000/GDPS/gdps_Class4/
REFERENCE=${GRID_DIR}/dest_grid.std.2

## 0.1 deg grid
rm ${REFERNCE} 
cstcreateZll -lomin 0 -lomax 360 -lamin -80 -lamax 90 -dlon 0.2 -dlat 0.2 -rpnfile ${REFERENCE} -nomvar TM 
