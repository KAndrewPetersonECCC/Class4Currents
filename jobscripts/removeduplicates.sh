#!/bin/bash
# ord_soumet /fs/homeu2/eccc/mrd/ords/rpnenv/dpe000/Class4_Currents/jobscripts/removeduplicates.sh -cpus 1 -mpi -cm 8000M -t 18000 -shell=/bin/bash
# bash /fs/homeu2/eccc/mrd/ords/rpnenv/dpe000/Class4_Currents/jobscripts/removeduplicates.sh


SCRATCH=/fs/site6/eccc/mrd/rpnenv/dpe000/Class4_Currents/SCRATCH
DIR1=/fs/site6/eccc/mrd/rpnenv/dpe000/Class4_Currents/CLASS4_currents_CCMEP_FILT4.APR/
DIR2=/fs/site6/eccc/mrd/rpnenv/dpe000/Class4_Currents/CLASS4_currents_CCMEP_FILT4.JUNE/

cd ${SCRATCH}

for file in ${DIR1}/*.nc ; do 
  bile=$(basename ${file}) ; 
  echo ${bile}
  rm [12].dump
  if [[ -e  ${DIR1}/${bile} ]]; then 
    ncdump ${DIR1}/${bile} > 1.dump ; 
  else
    ls -l ${DIR1}/${bile}
  fi
  if [[ -e  ${DIR2}/${bile} ]]; then 
    ncdump ${DIR2}/${bile} > 2.dump ; 
  fi
  diff 1.dump 2.dump ; 
  if [[ $? -eq 0 ]] ; then 
    echo "identical ${bile}" ; 
    rm ${file} ;
  fi ; 
done
