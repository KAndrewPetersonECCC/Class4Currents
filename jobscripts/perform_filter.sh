#!/bin/bash
# ord_soumet /home/dpe000/Class4_Currents/jobscripts/perform_filter.sh -cpus 1 -mpi -cm 64000M -t 3600 -shell=/bin/bash

# LOAD editfst
. r.load.dot rpn/utils/19.5

USAGE="USAGE:  perform_filter.sh -f=filter_file -i=input_standard_file -o=output_standard_file --clobber (overwrite input)"
CLOBBER=False

for i in "$@"
do
case $i in
    -f=*|--filter=*)
    FILTER="${i#*=}"
    shift # past argument=value
    ;;
    -i=*|--input=*)
    INPUT="${i#*=}"
    shift # past argument=value
    ;;
    -o=*|--output=*)
    OUTPUT="${i#*=}"
    shift # past argument=value
    ;;
    -O|--clobber|--overwrite)
    CLOBBER=True
    shift # past argument=value
    ;;
    *)
    echo ${USAGE}
          # unknown option
    ;;
esac
done

if [[ -z ${FILTER} ]] ; then
    echo "NEED FILTER FILE"
    echo ${USAGE}
    exit 101
fi

if [[ -z ${INPUT} ]] ; then
    echo "NEED INPUT FILE"
    echo ${USAGE}
    exit 102
fi

if [[ -z ${OUTPUT} ]] ; then
    OUTPUT=${INPUT}.tmp
    echo "DEFAULT OUTPUT: ${OUTPUT}"
fi

if [[ -f ${OUTPUT} ]] ; then
    echo "Removing pre-existing output:  EDITFST appends by default" 
    rm ${OUTPUT}
fi

bash ${FILTER} ${INPUT} ${OUTPUT}

if [[ ${CLOBBER} == True ]] ; then
   mv ${OUTPUT} ${INPUT}
fi
