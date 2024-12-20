#!/bin/bash
# ord_soumet /fs/homeu2/eccc/mrd/ords/rpnenv/dpe000/Class4_Currents/jobscripts/test_analysis_loop.sh -cpus 80 -mpi -cm 2000M -t 10800 -shell=/bin/bash
#
#bash /fs/homeu2/eccc/mrd/ords/rpnenv/dpe000/Class4_Currents/jobscripts/test_analysis_loop.sh
#
USAGE="USAGE:  GEPSdate_Class4Currents.sh -d=CCYYMMDD"

SUBMIT=False
SUFF=GIOPS_orca025_currents.f2
SUFR=GIOPS_orca025_currents.f2
ENSA="[0,0]"

for i in "$@"
do
case $i in
    -E=*|--expt=*|--EX=*)
    EXPT="${i#*=}"
    shift # past argument=value
    ;;
    -R=*|--reference=*|--control=*|--ctrl=*)
    REFE="${i#*=}"
    shift # past argument=value
    ;;
    -N=*|--ensemble=*|--ens=*)
    ENSA="${i#*=}"
    shift # past argument=value
    ;;
    -L=*|--label=*)
    LABE="${i#*=}"
    shift # past argument=value
    ;;
    -C=*|--ref_label=*|--reference_label=*)
    RLAB="${i#*=}"
    shift # past argument=value
    ;;
    --CC)
    RLAB="ctrl"
    shift # past argument=value
    ;;
    -s=*|--start=*)
    SATE="${i#*=}"
    shift # past argument=value
    ;;
    -f=*|--final=*)
    FATE="${i#*=}"
    shift # past argument=value
    ;;
    -p=*|--prefix=*)
    SUFF="${i#*=}"
    shift # past argument=value
    ;;
    -q=*|--prefix_ref=*)
    SUFR="${i#*=}"
    shift # past argument=value
    ;;
    -d|--do|--submit)
    SUBMIT=True
    shift # past argument=value
    ;;
    *)
    echo ${USAGE}
          # unknown option
    ;;
esac
done

if [[ -z ${EXPT} ]] ; then
    echo "NEED EXPT"
    echo ${USAGE}
    exit 99 
fi
if [[ -z ${REFE} ]] ; then
    echo "NEED REFERENCE"
    echo ${USAGE}
    exit 99 
fi
if [[ -z ${SATE} ]] ; then
    echo "NEED START DATE"
    echo ${USAGE}
    exit 99 
fi
if [[ -z ${FATE} ]] ; then
    echo "NEED FINAL DATE"
    echo ${USAGE}
    exit 99 
fi
if [[ -z ${LABE} ]] ; then
    LABE=${EXPT}
    echo ${LABE}
fi
if [[ -z ${RLAB} ]] ; then
    RLAB=${REFE}
    echo ${RLAB}
fi

cd /home/dpe000/Class4_Currents
source jobscripts/prepython.sh

ASUF=$(echo ${SUFF} | rev | cut -c-2 | rev)
BASE=JOBS/do_anal_${EXPT}_${REFE}_${SATE}_${FATE}.${ASUF}
BJOB=${BASE}.sh
PJOB=${BASE}.py
SJOB="ord_soumet ${BJOB} -cpus 80 -mpi -cm 2000M -t 10800 -shell=/bin/bash"

cat > ${BJOB} << EOB 
#!/bin/bash 
# ${SJOB}

cd /home/dpe000/Class4_Currents
source jobscripts/prepython.sh
python ${PJOB}
EOB

cat > ${PJOB} << EOP
from importlib import reload
import os
import sys
sys.path.insert(0, '/home/dpe000/Class4_Currents/python')
import matplotlib as mpl
mpl.use('Agg')
import traceback
import datetime
import numpy as np
import datadatefile
import Class4Current


expt = '${EXPT}'
expt0 = '${REFE}'

insuffixs=['GIOPS_orca025_currents.f2','GIOPS_orca025_currents.f2']
insuffixs=['${SUFF}', '${SUFR}']

if ( ${SATE} ):
    dates = Class4Current.create_dates(${SATE}, ${FATE})
else:
    dates = None

MNA_LIST, GDA_LIST, MNA_LITT = Class4Current.compare_analysis_errors(dates, [expt, expt0], ['${LABE}', '${RLAB}'], insuffixs=insuffixs, maxtaylor=1.0, ens_axes=${ENSA})
EOP

if [[ ${SUBMIT} == True ]] ; then 
  ${SJOB}
fi
