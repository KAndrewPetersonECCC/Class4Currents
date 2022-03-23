#!/bin/bash
#bash /fs/homeu1/eccc/mrd/ords/rpnenv/dpe000/Class4_Currents/jobscripts/date_Class4Currents.sh --start=CCYYMMDD --final=CCYYMMDD

WDIR=/fs/homeu1/eccc/mrd/ords/rpnenv/dpe000/Class4_Currents
USAGE="USAGE:  date_Class4Currents.sh -s=CCYYMMDD -f=CCYYMMDD"

IGNORE=False
PLOT="--p"
CHARLY=""
FILTER="--filter"
for i in "$@"
do
case $i in
    -s=*|--start=*)
    STARTDATE="${i#*=}"
    shift # past argument=value
    ;;
    -f=*|--final=*)
    FINALDATE="${i#*=}"
    shift # past argument=value
    ;;
    -p|--plot)
    PLOT="-p"
    shift # past argument=value
    ;;
    -q|--no_plot)
    PLOT="-q"
    shift # past argument=value
    ;;
    -o|--orig)
    CHARLY="-o"
    shift # past argument=value
    ;;
    --filter)
    FILTER="--filter"
    shift # past argument=value
    ;;
    --not_filter)
    FILTER="--not_filter"
    shift # past argument=value
    ;;
    -i|--ignore)
    IGNORE=True
    shift # past argument=value
    ;;
    *)
    echo ${USAGE}
          # unknown option
    ;;
esac
done

if [[ -z ${STARTDATE} ]] ; then
    echo "NEED STARTDATE"
    echo ${USAGE}
    exit 99 
fi
if [[ ${#STARTDATE} -ne 8 ]] ; then
  echo "NEED STARTDATE CCYYMMDD"
  exit 99
fi
if [[ -z ${FINALDATE} ]] ; then
    echo "NEED FINALDATE"
    echo ${USAGE}
    exit 99 
fi
if [[ ${#FINALDATE} -ne 8 ]] ; then
  echo "NEED FINALDATE CCYYMMDD"
  exit 99
fi

DATE=${STARTDATE}
THISTIME=$(date -u --date ${DATE} +%s)
FINALTIME=$(date -u --date ${FINALDATE} +%s)
ZERO=True

while [[ ${THISTIME} -le ${FINALTIME} ]] ; do 
    NEXT=$(date -u --date "${DATE} + 1 day" +%Y%m%d)
    NEXTTIME=$(date -u --date ${NEXT} +%s)
    if [[ ${ZERO} == True ]];  then
        echo "Submitting Zeroth Job"
	echo "jobscripts/date_Class4Currents.sh  --date=${DATE} ${PLOT} ${CHARLY} ${FILTER} --next=${NEXT} -s"
	bash ${WDIR}/jobscripts/date_Class4Currents.sh  --date=${DATE} ${PLOT} ${CHARLY} ${FILTER} --next=${NEXT} -s
	ZERO=False
    else
        echo "Preparing Other Jobs"
        if [[ ${NEXTTIME} -le ${FINALTIME} ]]; then 
	    echo "jobscripts/date_Class4Currents.sh  --date=${DATE} ${PLOT} ${CHARLY} ${FILTER} --next=${NEXT}"
            bash ${WDIR}/jobscripts/date_Class4Currents.sh  --date=${DATE} ${PLOT} ${CHARLY} ${FILTER} --next=${NEXT}
	else
	    echo "jobscripts/date_Class4Currents.sh  --date=${DATE} ${PLOT} ${CHARLY} ${FILTER}"
	    bash ${WDIR}/jobscripts/date_Class4Currents.sh  --date=${DATE} ${PLOT} ${CHARLY} ${FILTER}
	fi
    fi
     
    DATE=$(date -u --date "${DATE} + 1 day" +%Y%m%d)
    THISTIME=$(date -u --date ${DATE} +%s)
    echo ${DATE} ${THISTIME}
    
done
