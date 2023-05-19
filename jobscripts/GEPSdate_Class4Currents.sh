#!/bin/bash
#bash /fs/homeu2/eccc/mrd/ords/rpnenv/dpe000/Class4_Currents/jobscripts/GEPSdate_Class4Currents.sh --date=CCYYMMDD

USAGE="USAGE:  GEPSdate_Class4Currents.sh -d=CCYYMMDD"

SUBMIT=False
FILTER=True
ENAN=False
FCEX=GEPS_STO2X
for i in "$@"
do
case $i in
    -A|--Ensemble_Analysis|--EA)
    ENAN=True
    shift # past argument=value
    ;;
    -F=*|--forecast=*|--fcst=*)
    FCEX="${i#*=}"
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
if [[ -z ${ENSM} ]] ; then
    echo "NEED ENSM"
    echo ${USAGE}
    exit 99 
fi

WDIR=/fs/homeu2/eccc/mrd/ords/rpnenv/dpe000/Class4_Currents
cd ${WDIR}

if [[ ${ENAN} == False ]]; then 
    GEPS=GEPS
else
    GEPS=${FCEX}
fi
BJOB=${WDIR}/JOBS/${GEPS}_Class4Currents.${DATE}.${ENSM}.${ENAN:0:1}.${FILTER:0:1}.sh
NJOB=${WDIR}/JOBS/${GEPS}_Class4Currents.${NEXT}.${ENSM}.${ENAN:0:1}.${FILTER:0:1}.sh
PJOB=${WDIR}/JOBS/${GEPS}_Class4Currents.${DATE}.${ENSM}.${ENAN:0:1}.${FILTER:0:1}.py
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
import Class4CurrentEA
datestr="${DATE}"
date=datadatefile.convert_strint_date(datestr)
if not ( '${ENSM}' == 'A' ):
    ens_list=[${ENSM}]
    print("PROCESSING DATE", datestr, date, ens_list)
    try:
        t0 = time.time()
        if not ${ENAN}:
            Class4Current.process_geps_obs(date=date, ens_list=ens_list, filter=${FILTER})
        else:
            Class4CurrentEA.process_enan_obs(date=date, ens_list=ens_list, filter=${FILTER}, expt='${FCEX}')
        te = time.time() - t0
        print("PROCESS SUCCESS:  TIME ELAPSED ", te)
    except:
        print(traceback.format_exc())
        sys.stdout.flush()
else:
    print("ASSEMBLING DATE", datestr, date)
    try:
        # NEED TO ADD NON FILTER OPTIONS.
        t0 = time.time()
        if not ${ENAN}:
            Class4Current.assemble_ensemble_date(date, obspre='CLASS4_currents_GEPS_FILT/class4', obssuf='GEPS_orca025_currents', iters=['f1','f2'], nens=21, clobber=True)
        else:
            Class4CurrentEA.assemble_ensembleEA_date(date, obspre='CLASS4_currents_ENAN_FILT/class4', obssuf='ENAN_orca025_currents', iters=['f1','f2'], nens=21, clobber=True, expt='${FCEX}')
        print("ASSEMBLE SUCCESS:  TIME ELAPSED ", te)
    except:
        print(traceback.format_exc())
        sys.stdout.flush()
EOP

if [[ ${SUBMIT} == True ]] ; then 
    ${SJOB}
fi

