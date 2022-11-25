#!/bin/bash

# ord_soumet /home/dpe000/Class4_Currents/jobscripts/perform_interpolation.sh -cpus 1 -mpi -cm 64000M -t 3600 -shell=/bin/bash

r.load /fs/ssm/eccc/mrd/rpn/OCEAN/cstint/3.2.13

GRID_DIR=/fs/site6/eccc/mrd/rpnenv/dpe000/GDPS/gdps_Class4/
REFERENCE=${GRID_DIR}/dest_grid.std.25

## 0.25 deg grid
rm ${REFERNCE} 
cstcreateZll -lomin 0 -lomax 360 -lamin -80 -lamax 90 -dlon 0.25 -dlat 0.25 -rpnfile ${REFERENCE} -nomvar TM 
