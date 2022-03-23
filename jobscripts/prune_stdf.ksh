#!/bin/bash
# ord_soumet /home/dpe000/Class4_Currents/jobscripts/prune_stdf.ksh -cpus 1 -mpi -cm 4000M -t 900 -shell=/bin/bash

# LOAD editfst
. r.load.dot rpn/utils/19.5

ODIR=/fs/site4/eccc/mrd/rpnenv/dpe000/GDPS/gdps_Class4
file=${ODIR}/2021103000_024.orca025

rm ${file}.prune
editfst -s  ${file} -d ${file}.prune -l << EOF
desire(-1,['^^','>>'])
desire(-1,['TM','UUW', 'VVW'], -1, -1, -1, -1)
EOF

rm ${file}.prune.1
editfst -s  ${file} -d ${file}.prune.1 -l << EOF
desire(-1,['^^','>>'])
desire(-1,['TM'], -1, -1, 10979785, -1)
EOF
