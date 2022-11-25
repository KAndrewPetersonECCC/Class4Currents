#!/bin/bash
# ord_soumet /home/dpe000/Class4_Currents/jobscripts/perform_interpolation.sh -cpus 1 -mpi -cm 64000M -t 3600 -shell=/bin/bash

GRID_DIR=/fs/site6/eccc/mrd/rpnenv/dpe000/GDPS/gdps_Class4
REFERENCE=${GRID_DIR}/dest_grid.std.25

for i in "$@"
do
case $i in
    -fs=*|--source=*)
    SOURCE="${i#*=}"
    shift # past argument=value
    ;;
    -fd=*|--destination=*)
    DESTINATION="${i#*=}"
    shift # past argument=value
    ;;
    -fr=*|--reference=*)
    REFERENCE="${i#*=}"
    shift # past argument=value
    ;;
    -i|--ignore)
    IGNORE=True
    shift # past argument=value
    ;;
    *)
    echo ${USAGE}
          # unknown option
    ;;
esac
done

. r.load.dot /fs/ssm/eccc/mrd/rpn/OCEAN/cstint/3.2.13

rm ${DESTINATION}
cstintrp -fs ${SOURCE} -fr ${REFERENCE} -nr TM -fstype custom -fdtype rpn -fd ${DESTINATION} -intyp bilin -glbsrc T
