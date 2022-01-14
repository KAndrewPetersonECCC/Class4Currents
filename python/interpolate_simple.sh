#!/bin/bash
# ord_soumet /home/dpe000/Class4_Currents/python/interpolate_simple.sh -cpus 1 -mpi -cm 64000M -t 3600 -shell=/bin/bash


. ssmuse-sh -d /fs/ssm/eccc/mrd/rpn/OCEAN/cstint-3.2.7
. ssmuse-sh -p eccc/crd/ccmr/EC-CAS/master/fstd2nc_0.20190903.0

ODIR=/fs/site4/eccc/mrd/rpnenv/dpe000/GDPS/gdps_Class4
DDIR=/fs/site4/eccc/mrd/rpnenv/dpe000/GDPS/gdps_Class4

cd ${DDIR} 

## 0.1 deg grid
rm $DDIR/dest_grid.std.2
cstcreateZll -lomin 0 -lomax 360 -lamin -80 -lamax 90 -dlon 0.2 -dlat 0.2 -rpnfile $DDIR/dest_grid.std.2 -nomvar TM 

rm ${DDIR}/2021103000_024.dest.2
cstintrp -fs ${ODIR}/2021103000_024.orca025.prune  -fr ${DDIR}/dest_grid.std.2 -nr TM -fstype custom -fdtype rpn -fd $DDIR/2021103000_024.dest.2 -intyp bilin -status ${DDIR}/file_status.2 -glbsrc T
