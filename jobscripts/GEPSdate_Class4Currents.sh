#!/bin/bash
#bash /fs/homeu2/eccc/mrd/ords/rpnenv/dpe000/Class4_Currents/jobscripts/GEPSdate_Class4Currents.sh --date=CCYYMMDD

USAGE="USAGE:  GEPSdate_Class4Currents.sh -d=CCYYMMDD"

SUBMIT=False
FILTER=True
for i in "$@"
do
case $i in
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
if [[ -z ${ENSM} ]] ; then
    echo "NEED ENSM"
    echo ${USAGE}
    exit 99 
fi

WDIR=/fs/homeu2/eccc/mrd/ords/rpnenv/dpe000/Class4_Currents
cd ${WDIR}

BJOB=${WDIR}/JOBS/GEPS_Class4Currents.${DATE}.${ENSM}.${FILTER:0:1}.sh
NJOB=${WDIR}/JOBS/GEPS_Class4Currents.${NEXT}.${ENSM}.${FILTER:0:1}.sh
PJOB=${WDIR}/JOBS/GEPS_Class4Currents.${DATE}.${ENSM}.${FILTER:0:1}.py
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
import datadatefile
import Class4Current
datestr="${DATE}"
date=datadatefile.convert_strint_date(datestr)
if not ( '${ENSM}' == 'A' ):
    ens_list=[${ENSM}]
    print("PROCESSING DATE", datestr, date, ens_list)
    try:
        Class4Current.process_geps_obs(date=date, ens_list=ens_list, filter=${FILTER})
    except:
        print(traceback.format_exc())
        sys.stdout.flush()
else:
    print("ASSEMBLING DATE", datestr, date)
    try:
        # NEED TO ADD NON FILTER OPTIONS.
        Class4Current.assemble_ensemble_date(date, obspre='CLASS4_currents_GEPS_FILT/class4', obssuf='GEPS_orca025_currents', iters=['f1','f2'], nens=21, clobber=True)
    except:
        print(traceback.format_exc())
        sys.stdout.flush()
EOP

if [[ ${SUBMIT} == True ]] ; then 
    ${SJOB}
fi

