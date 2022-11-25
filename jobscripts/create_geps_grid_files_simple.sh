#!/bin/bash

script_dir=/fs/homeu2/eccc/mrd/ords/rpnenv/dpe000/Class4_Currents/jobscripts

INDIR=/fs/site5/eccc/prod/ops/suites/geps_20220621/e1/gridpt/prog/ens.glboce/
OUTDIR=/fs/site5/eccc/mrd/rpnenv/dpe000/tmpdir/GEPS_TMP

for i in "$@"
do
case $i in
    -d=*|--date=*)
    DATE="${i#*=}"
    shift # past argument=value
    ;;
    -l=*|-h=*|--lead=*)
    LDHR=("${i#*=}")
    shift # past argument=value
    ;;
    -e=*|--ensm=*)
    ENSM=("${i#*=}")
    shift # past argument=value
    ;;
    -i=*|--indir=*)
    INDIR="${i#*=}"
    shift # past argument=value
    ;;
    -o=*|--outdir=*)
    OUTDIR="${i#*=}"
    shift # past argument=value
    ;;
    *)
    echo ${USAGE}
          # unknown option
    ;;
esac
done

if [[ -z ${DATE} ]] ; then
    echo "NEED DATE CCYYMMDDHH"
    exit 99 
fi
if [[ ${#DATE} -eq 8 ]] ; then
  echo "NEED DATE CCYYMMDDHH -- ADDING 00"
  DATE=${DATE}00
fi
if [[ ${#DATE} -ne 10 ]] ; then
  echo "NEED DATE CCYYMMDDHH"
  exit 99
fi

if [[ -z ${LDHR} ]] ; then 
    echo "NEED LDHR ASSUME ALL"
    LDHR=($(seq -f '%03g' 24 24 768))
fi

if [[ -z ${ENSM} ]] ; then
   echo "NEED ENSM ASSUME ALL"
   ENSM=($(seq -f '%03g' 0 20))
fi
if [[ ${LDHR} == ALL ]] ; then 
    LDHR=($(seq -f '%03g' 24 24 768))
fi

if [[ ${ENSM} == ALL ]] ; then
   ENSM=($(seq -f '%03g' 0 20))
fi

if [[ ! -d ${OUTDIR} ]] ; then 
    mkdir -p ${OUTDIR} 
fi

stime=$(date +%s)

for LDH in ${LDHR[*]} ; do
   LDH=$(printf "%03g" ${LDH})
   for ENS in ${ENSM[*]} ; do
       ENS=$(printf "%03g" ${ENS})
       BASEF=${DATE}_${LDH}_${ENS}
       INPUT=${INDIR}/${BASEF}
       if [[ ! -f ${INPUT} ]] ; then 
           echo "No file \${INPUT}"
       else
           OUTINT=${OUTDIR}/${BASEF}.z25
           OUTNCF=${OUTDIR}/${BASEF}.z25.nc
           bash ${script_dir}/perform_interpolation_25.sh -fs=${INPUT} -fd=${OUTINT}
           bash ${script_dir}/create_netcdf.sh -i=${OUTINT} -o=${OUTNCF}
       fi
    done
done 


