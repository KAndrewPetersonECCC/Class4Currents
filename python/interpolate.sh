#!/bin/bash
# ord_soumet /home/dpe000/Class4_Currents/python/interpolate.sh -cpus 1 -mpi -cm 64000M -t 3600 -shell=/bin/bash


. ssmuse-sh -d /fs/ssm/eccc/mrd/rpn/OCEAN/cstint-3.2.7

ODIR=/fs/site4/eccc/mrd/rpnenv/dpe000/GDPS/gdps_Class4
DDIR=/fs/site4/eccc/mrd/rpnenv/dpe000/GDPS/gdps_Class4

## 0.1 deg grid
rm $DDIR/dest_grid.std.1
cstcreatell -lomin 0 -lomax 360 -lamin -80 -lamax 90 -dlon 0.1 -dlat 0.1 -rpnfile $DDIR/dest_grid.std.1 -nomvar TM 

## use -owgts T
rm ${DDIR}/2021103000_024.dest.1
cstintrp -ns TM -fs ${ODIR}/2021103000_024.orca025.prune.1  -fr ${DDIR}/dest_grid.std.1 -nr TM -fstype rpn -fd $DDIR/2021103000_024.dest.1 -intyp bilin -owgts T -status ${DDIR}/file_status.1 > ${DDIR}/output_interp.1.log



## 0.05 deg grid
rm ${DDIR}/dest_grid.std.05
cstcreatell -lomin 0 -lomax 360 -lamin -80 -lamax 90 -dlon 0.05 -dlat 0.05 -rpnfile ${DDIR}/dest_grid.std.05 -nomvar TT 


##  use -owgts T
rm ${DDIR}/2021103000_024.dest.05
cstintrp -fs ${ODIR}/2021103000_024.orca025.prune.1 -fr ${DDIR}/dest_grid.std.05 -nr TT -fstype rpn -fd ${DDIR}/2021103000_024.dest.05 -intyp bilin -owgts T -status ${DDIR}/file_status.05 > ${DDIR}/output_interp.05.log
