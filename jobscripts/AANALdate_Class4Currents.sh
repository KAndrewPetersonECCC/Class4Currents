#!/bin/bash
#bash /fs/homeu2/eccc/mrd/ords/rpnenv/dpe000/Class4_Currents/jobscripts/ANALdate_Class4Currents.sh --date=CCYYMMDD

USAGE="USAGE:  ANALdate_Class4Currents.sh -d=CCYYMMDD --filter"

SUBMIT=False
FILTER=True

KDIR=/fs/site5/eccc/cmd/e/kch001/maestro_archives/
DDIR=/fs/site6/eccc/mrd/rpnenv/dpe000/maestro_archives/
DDRU=/fs/site6/eccc/mrd/rpnenv/dpe000/maestro_archives/
DAT6=/fs/site6/eccc/mrd/rpnenv/dpe000/Class4_Currents

for i in "$@"
do
case $i in
    -E=*|--experiment=*|--expt=*)
    EXPT="${i#*=}"
    shift # past argument=value
    ;;
    -D=*|--data_dir=*|--ddir=*)
    DDIR="${i#*=}"
    shift # past argument=value
    ;;
    -d=*|--date=*)
    DATE="${i#*=}"
    shift # past argument=value
    ;;
    -n=*|--next=*)
    NEXT="${i#*=}"
    shift # past argument=value
    ;;
    -e=*|--ensemble=*)
    ENSM="${i#*=}"
    shift # past argument=value
    ;;
    -s|--submit)
    SUBMIT=True
    shift # past argument=value
    ;;
    -f|--filter)
    FILTER=True
    shift # past argument=value
    ;;
    -u|--not_filter)
    FILTER=False
    shift # past argument=value
    ;;
    *)
    echo ${USAGE}
          # unknown option
    ;;
esac
done

if [[ -z ${DATE} ]] ; then
    echo "NEED DATE"
    echo ${USAGE}
    exit 99 
fi
if [[ ${#DATE} -ne 8 ]] ; then
  echo "NEED DATE CCYYMMDD"
  exit 99
fi
if [[ -z ${EXPT} ]] ; then
    echo "NEED EXPT"
    echo ${USAGE}
    exit 99 
fi
if [[ -z ${DDIR} ]] ; then
    echo "NEED DDIR"
    echo ${USAGE}
    exit 99 
fi
if [[ ${DDIR} == KAMEL ]] ; then
    DDIR=${KDIR}
fi
if [[ ${DDIR} == DREW ]] ; then
    DDIR=${DDRU}
fi
if [[ ! -d ${DAT6}/${EXPT} ]] ; then 
    mkdir ${DAT6}/${EXPT}
    ln -s ${DAT6}/${EXPT} . 
fi

WDIR=/fs/homeu2/eccc/mrd/ords/rpnenv/dpe000/Class4_Currents
cd ${WDIR}

BJOB=${WDIR}/JOBS/ANAL_Class4Currents.${DATE}.${EXPT}.${FILTER:0:1}.sh
NJOB=${WDIR}/JOBS/ANAL_Class4Currents.${NEXT}.${EXPT}.${FILTER:0:1}.sh
PJOB=${WDIR}/JOBS/ANAL_Class4Currents.${DATE}.${EXPT}.${FILTER:0:1}.py
SJOB="ord_soumet ${BJOB} -cpus 1 -mpi -cm 64000M -t 21600 -shell=/bin/bash"
CJOB="ord_soumet ${NJOB} -cpus 1 -mpi -cm 64000M -t 21600 -shell=/bin/bash"

cat > ${BJOB} << EOJ
#!/bin/bash
####${SJOB}
echo "STARTING JOB for DATE ${DATE}"
cd ${WDIR}
## Adding NCO operations
. ssmuse-sh -d eccc/cmd/cmds/ext/20220331
## Adding RPNPY
source jobscripts/prepython.sh
python ${PJOB}
echo "FINISHED JOB for DATE ${DATE}"
EOJ

if [[ ! -z ${NEXT} ]] ; then
if [[ ${#NEXT} -eq 8 ]] ; then
cat >> ${BJOB} << EOJ
echo "Adding NEXT JOB for ${NEXT}"
${CJOB}
EOJ
fi
fi

cat > ${PJOB} << EOP
import os
import sys
sys.path.insert(0, '/home/dpe000/Class4_Currents/python')
import matplotlib as mpl
mpl.use('Agg')
import traceback
import datetime
import time
import datadatefile
import Class4Current
import Class4CurrentEAA
datestr="${DATE}"
date=datadatefile.convert_strint_date(datestr)

try:
    # NEED TO ADD NON FILTER OPTIONS.
    t0 = time.time()
    Class4CurrentEAA.process_anal_obs(date=date, filter=${FILTER}, expt='${EXPT}', ddir='${DDIR}')
    te = time.time() - t0
    print("PROCESS SUCCESS:  TIME ELAPSED ", te)
except:
    print(traceback.format_exc())
    sys.stdout.flush()
EOP

if [[ ${SUBMIT} == True ]] ; then 
    ${SJOB}
fi

