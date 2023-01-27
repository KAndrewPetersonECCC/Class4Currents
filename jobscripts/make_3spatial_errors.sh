#!/bin/bash
#
#bash /fs/homeu2/eccc/mrd/ords/rpnenv/dpe000/Class4_Currents/jobscripts/make_3spatial_errors.sh --expt=GDPS -c=1 --t=square

USAGE="USAGE:  date_Class4Currents.sh --begin=CCYYMMDD --end=CCYYMMDD"

IERROR=1
ERTYPE='square'
MAG=False
#INDIR="CLASS4_currents_CCMEP_FILT2"
#INSUF="GIOPS_orca025_currents-filter"
#EXPT="PSY"

START_DATE=20210101
FINAL_DATE=20211231

for i in "$@"
do
case $i in
    -b=*|--begin=*)
    START_DATE="${i#*=}"
    shift # past argument=value
    ;;
    -e=*|--end=*)
    FINAL_DATE="${i#*=}"
    shift # past argument=value
    ;;
    -o=*|--expt=*|--model=*)
    EXPT="${i#*=}"
    shift # past argument=value
    ;;
    -c=*|--correct=*|--ierror=*)
    IERROR="${i#*=}"
    shift # past argument=value
    ;;
    -t=*|--type=*)
    ERTYPE="${i#*=}"
    shift # past argument=value
    ;;
    -m|--magnitude|--angle)
    MAG=True
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
   exit 99
fi

OUTDIR="EPLOTS/${EXPT}_"

WDIR=/fs/homeu2/eccc/mrd/ords/rpnenv/dpe000/Class4_Currents
cd ${WDIR}

BJOB=${WDIR}/JOBS/makemaps.${EXPT}.${START_DATE}_${FINAL_DATE}.${ERTYPE}.${IERROR}${MAG:0:1}.sh
PJOB=${WDIR}/JOBS/makemaps.${EXPT}.${START_DATE}_${FINAL_DATE}.${ERTYPE}.${IERROR}.${MAG:0:1}.py
SJOB="ord_soumet ${BJOB} -cpus 80 -mpi -cm 2000M -t 21600 -shell=/bin/bash"
cat > ${BJOB} << EOJ
#!/bin/bash
####${SJOB}
echo "STARTING JOB for DATE ${DATE}"
cd ${WDIR}
source jobscripts/prepython.sh
python ${PJOB}
echo "FINISHED JOB for ${EXPT}"
EOJ

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

datestart_str="${START_DATE}"
datefinal_str="${FINAL_DATE}"
datestart=Class4Current.check_date(datestart_str, outtype=datetime.datetime)
datefinal=Class4Current.check_date(datefinal_str, outtype=datetime.datetime)
daterange=[datestart, datefinal]

ierror=${IERROR}
etype="${ERTYPE}"
mag=${MAG}

expt='${EXPT}'
if ( expt == 'GIOPS' or expt == 'GDPS' ):
    indir='CLASS4_currents_CCMEP_FILT2'
    suffix="GIOPS_orca025_currents-filter"
    outdir='EPLOTS/${EXPT}${IERROR}_'
      
if ( expt == 'GEPS' ):
    indir='CLASS4_currents_GEPS_FILT2'
    suffix="GEPS_orca025_currents-filter.enm"
    outdir='EPLOTS/${EXPT}${IERROR}_'

if ( expt == 'PSY' ):
    indir='CLASS4_currents_CHARLY'
    suffix='PSY4V3R1_orca12_currents-filtr'
    outdir='EPLOTS/${EXPT}${IERROR}_'
               
Class4Current.make_spatial_plots_3errors(indir, suffix, daterange, outdir=outdir, magnitude=mag, levels=np.arange(0.0, 0.55, 0.05), alevels=np.arange(-.475, 0.5, 0.05), ierror=ierror, etype=etype, mp=True)

EOP

${SJOB}
