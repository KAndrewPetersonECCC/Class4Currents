#!/bin/bash
# ord_soumet /fs/homeu1/eccc/mrd/ords/rpnenv/dpe000/Class4_Currents/jobscripts/final_rename.sh -cpus 1 -mpi -cm 8000M -t 18000 -shell=/bin/bash
#bash /fs/homeu1/eccc/mrd/ords/rpnenv/dpe000/Class4_Currents/jobscripts/final_rename.sh

# THIS WOULD NOT BE NECESSARY IF IT WAS DONE IN final_product.

SDATE=????????
WDIR=/fs/site4/eccc/mrd/rpnenv/dpe000/Class4_Currents
HDIR=/fs/homeu1/eccc/mrd/ords/rpnenv/dpe000/Class4_Currents
SOURCE=CLASS4_currents_CCMEP_FILT2

OPWD=$(pwd)
cd ${WDIR}

for file in ${SOURCE}/class4_${SDATE}_GIOPS_orca025_currents.nc.gz ; do 
    mv ${file} ${file/currents.nc/currents-filtr.nc}
done

cd ${OPWD}
