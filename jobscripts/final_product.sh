#!/bin/bash
# ord_soumet /fs/homeu2/eccc/mrd/ords/rpnenv/dpe000/Class4_Currents/jobscripts/final_product.sh -cpus 1 -mpi -cm 8000M -t 18000 -shell=/bin/bash
# bash /fs/homeu2/eccc/mrd/ords/rpnenv/dpe000/Class4_Currents/jobscripts/final_product.sh

WDIR=/fs/site6/eccc/mrd/rpnenv/dpe000/Class4_Currents
HDIR=/fs/homeu2/eccc/mrd/ords/rpnenv/dpe000/Class4_Currents

## Adding NCO operations
. ssmuse-sh -d eccc/cmd/cmds/ext/20220331

cd ${WDIR}

SDATE=202304??
FFS=(f4)

SOURCE=CLASS4_currents_CCMEP_FILT
for FF in ${FFS[*]} ; do
  case ${FF} in
    f1) dest=FILT1 ;;
    f2) dest=FILT2 ;;
    f3) dest=FILT3 ;;
    f4) dest=FILT4 ;;
  esac
  case ${FF} in
    f1) suff="-filtr" ;;
    f2) suff="-filtr" ;;
    f3) suff="-filtr" ;;
    f4) suff="-filtr" ;;
  esac
  
  DEST=CLASS4_currents_CCMEP_${dest}
  if [[ ! -d ${DEST} ]] ; then mkdir ${DEST} ; fi
  if [[ ! -s ${HDIR}/${DEST} ]] ; then ln -s ${WDIR}/${DEST} ${HDIR}/${DEST} ; fi
  
  
  for file in ${SOURCE}/class4_${SDATE}_GIOPS_orca025_currents.${FF}.nc ; do
      bile=$(basename ${file})
      dile=${bile/$FF/tmp}
      eile=${dile/.tmp.nc/${suff}.nc}
      date=$(echo ${bile} | cut -c 8-15)
      echo ${date}
      
      system=GIOPSv3.0.0
      if [[ ${date} -gt 20211201 ]]; then system=GIOPSv3.3.0 ; fi 
      ncks -O -x -v negative_persistence ${file} -o ${DEST}/${dile}
      ncatted -O  -a contact,global,o,c,"Andrew.Peterson@ec.gc.ca" ${DEST}/${dile}
      ncatted -O  -a init_estimate_description,global,o,c,"analysis produced in near real time" ${DEST}/${dile}
      ncatted -O  -a system,global,o,c,"${system}"  ${DEST}/${dile}
      ncatted -O  -a institution,global,o,c,"CCMEP" ${DEST}/${dile}
      ncks -O ${DEST}/${dile} -o ${DEST}/${eile}
      ncatted -O -h -a history,global,d,, ${DEST}/${eile}
      if [[ ! -e ${eile}.gz ]] ; then 
          gzip ${eile}
      fi
      rm ${DEST}/${dile}
   done
done

