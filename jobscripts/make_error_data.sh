#!/bin/bash
# ord_soumet /fs/homeu1/eccc/mrd/ords/rpnenv/dpe000/Class4_Currents/jobscripts/make_error_data.sh -cpus 1 -mpi -cm 64000M -t 10800 -shell=/bin/bash
#bash /fs/homeu1/eccc/mrd/ords/rpnenv/dpe000/Class4_Currents/jobscripts/make_error_data.sh --date=CCYYMMDD

USAGE="USAGE:  date_Class4Currents.sh -d=CCYYMMDD"

FILTER=True
THREE=False
DATE="????????"
IIT=0

for i in "$@"
do
case $i in
    -d=*|--date=*)
    DATE="${i#*=}"
    shift # past argument=value
    ;;
    -e=*|--expt=*|--model=*)
    EXPT="${i#*=}"
    shift # past argument=value
    ;;
    -i=*|--intr=*)
    IIT="${i#*=}"
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
    -3|--three)
    THREE=True
    shift # past argument=value
    ;;
    *)
    echo ${USAGE}
          # unknown option
    ;;
esac
done

WDIR=/fs/homeu1/eccc/mrd/ords/rpnenv/dpe000/Class4_Currents
cd ${WDIR}

TIME=$(date +%s)
BJOB=${WDIR}/JOBS/makeerrors.${EXPT}.${FILTER:0:1}${THREE:0:1}${IIT}.${TIME}.sh
PJOB=${WDIR}/JOBS/makeerrors.${EXPT}.${FILTER:0:1}${THREE:0:1}${IIT}.${TIME}.py
SJOB="ord_soumet ${BJOB} -cpus 1 -mpi -cm 16000M -t 21600 -shell=/bin/bash"
cat > ${BJOB} << EOJ
#!/bin/bash
####${SJOB}
echo "STARTING JOB for DATE ${DATE}"
cd ${WDIR}
source jobscripts/preconda.sh
source activate metcarto
python ${PJOB}
echo "FINISHED JOB for ${EXPT}"
EOJ

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
Class4Current.write_mean_errors_model(datestr, '${EXPT}', $IIT, filter=$FILTER,three=THREE)
EOP

${SJOB}
