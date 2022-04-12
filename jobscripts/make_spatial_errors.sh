#!/bin/bash
# ord_soumet /fs/homeu1/eccc/mrd/ords/rpnenv/dpe000/Class4_Currents/jobscripts/make_spatial_errors.sh -cpus 1 -mpi -cm 64000M -t 10800 -shell=/bin/bash
#bash /fs/homeu1/eccc/mrd/ords/rpnenv/dpe000/Class4_Currents/jobscripts/make_spatial_erros.sh

USAGE="USAGE:  date_Class4Currents.sh --begin=CCYYMMDD --end=CCYYMMDD"

FILTER=True
SPEEDA=False
IIT=0
START_DATE=20190401
FINAL_DATE=20220301

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
    -m=*|--expt=*|--model=*)
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
    -a|--angle|-s|--speed)
    SPEEDA=True
    shift # past argument=value
    ;;
    -v|--vector)
    SPEEDA=False
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

BJOB=${WDIR}/JOBS/makemaps.${EXPT}.${START_DATE}_${FINAL_DATE}.${SPEEDA:0:1}${FILTER:0:1}${IIT}.sh
PJOB=${WDIR}/JOBS/makemaps.${EXPT}.${START_DATE}_${FINAL_DATE}.${SPEEDA:0:1}${FILTER:0:1}${IIT}.py
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
import numpy as np
import datadatefile
import Class4Current
datestart_str="${START_DATE}"
datefinal_str="${FINAL_DATE}"
datestart=Class4Current.check_date(datestart_str, outtype=datetime.datetime)
datefinal=Class4Current.check_date(datefinal_str, outtype=datetime.datetime)
daterange=[datestart, datefinal]

expt='${EXPT}'
filter=${FILTER}
speed=${SPEEDA}
if ( expt == 'GIOPS' ):
    indir='CLASS4_currents_CCMEP_'
    suffix='GIOPS_orca025_currents.'
    outdir='EPLOTS/${EXPT}.'
    if ( filter ): 
      indir=indir+'FILT'
      suffix=suffix+'f'+str(${IIT})
      addout='f'+str(${IIT})
    if ( not filter): 
      indir=indir+'UFIL'
      suffix=suffix+str(${IIT})
      addout='u'+str(${IIT})
    outdir=outdir+addout+'_'
      
if ( expt == 'PSY' ):
    indir='CLASS4_currents_CHARLY'
    suffix='PSY4V3R1_orca12_currents'
    outdir='EPLOTS/${EXPT}.'
    if ( not filter ):
      addout='u'
    if ( filter ):
      suffix=suffix+'-filtr'    
      addout='f'
    outdir=outdir+addout+'_'

Class4Current.make_spatial_plots(indir, suffix, daterange, outdir=outdir, 
                       calculate_error=True, magnitude=speed, plot_daily=True, plot_time_mean=True, 
		       levels=np.arange(-1.0, 1.1, 0.1))

EOP

${SJOB}
