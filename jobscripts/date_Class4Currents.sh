#!/bin/bash
# ord_soumet /fs/homeu2/eccc/mrd/ords/rpnenv/dpe000/Class4_Currents/jobscripts/date_Class4Currents.sh -cpus 1 -mpi -cm 64000M -t 10800 -shell=/bin/bash
#bash /fs/homeu2/eccc/mrd/ords/rpnenv/dpe000/Class4_Currents/jobscripts/date_Class4Currents.sh --date=CCYYMMDD

USAGE="USAGE:  date_Class4Currents.sh -d=CCYYMMDD"

SUBMIT=False
SINGLE=False
PLOT=False
CHARLY=True
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
    -s|--submit)
    SUBMIT=True
    shift # past argument=value
    ;;
    -1|--single)
    SINGLE=True
    shift # past argument=value
    ;;
    -p|--plot)
    PLOT=True
    shift # past argument=value
    ;;
    -q|--no_plot)
    PLOT=False
    shift # past argument=value
    ;;
    -o|--orig)
    CHARLY=False
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

WDIR=/fs/homeu2/eccc/mrd/ords/rpnenv/dpe000/Class4_Currents
cd ${WDIR}

BJOB=${WDIR}/JOBS/date_Class4Currents.${DATE}.${CHARLY:0:1}${FILTER:0:1}.sh
NJOB=${WDIR}/JOBS/date_Class4Currents.${NEXT}.${CHARLY:0:1}${FILTER:0:1}.sh
PJOB=${WDIR}/JOBS/date_Class4Currents.${DATE}.${CHARLY:0:1}${FILTER:0:1}.py
SJOB="ord_soumet ${BJOB} -cpus 1 -mpi -cm 64000M -t 21600 -shell=/bin/bash"
CJOB="ord_soumet ${NJOB} -cpus 1 -mpi -cm 64000M -t 21600 -shell=/bin/bash"

cat > ${BJOB} << EOJ
#!/bin/bash
####${SJOB}
echo "STARTING JOB for DATE ${DATE}"
cd ${WDIR}
. r.load.dot rpn/utils/20220509
source jobscripts/prepython.sh
#source activate metcarto
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
datefile='DATES/todo.date'
date=datadatefile.convert_strint_date(datestr)
print("PROCESSING DATE", datestr, date)
try:
    Class4Current.process_obsfile(date=date, TEST_SINGLE=${SINGLE}, Plot=${PLOT}, CHARLY=${CHARLY}, filter=${FILTER})
except:
    print(traceback.format_exc())
    sys.stdout.flush()
    if ( not ${SINGLE} ):
        datadatefile.add_to_boolean_file(int(datestr), False, file=datefile)
else:
    if ( not ${SINGLE} ):
        datadatefile.add_to_boolean_file(int(datestr), True, file=datefile)
EOP

if [[ ${SUBMIT} == True ]] ; then 
    ${SJOB}
fi

