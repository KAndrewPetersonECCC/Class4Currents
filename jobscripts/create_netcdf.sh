#!/bin/bash
# /fs/homeu2/eccc/mrd/ords/rpnenv/dpe000/Class4_Currents/jobscripts/create_netcdf.sh -i=input_file -o=output_file

for i in "$@"
do
case $i in
    -i=*|--input=*)
    INPUT="${i#*=}"
    shift # past argument=value
    ;;
    -o=*|--output=*)
    OUTPUT="${i#*=}"
    shift # past argument=value
    ;;
    *)
    echo ${USAGE}
          # unknown option
    ;;
esac
done

if [[ -z ${INPUT} ]]; then 
    echo "NEED INPUT FILE -i=input_file"
    exit 99
fi
if [[ -z ${OUTPUT} ]]; then 
    echo "NEED OUTPUT FILE -o=output_file"
    exit 99
fi

if [[ ! -f ${INPUT} ]]; then
    echo "NEED INPUT FILE to exist -- ${INPUT}"
    exit 99
fi
if [[ -f ${OUTPUT} ]] ; then
   echo "WARNING.  FILE WILL BE CLOBBERED -- ${OUTPUT}"
   exit 99
fi

. ssmuse-sh -x eccc/crd/ccmr/EC-CAS/fstd2nc/0.20220829.0
. r.load.dot eccc/mrd/rpn/MIG/ENV/rpnpy/2.1-u2.5  main/opt/intelcomp/inteloneapi-2022.1.2/intelcomp+mpi+mkl

fstd2nc $INPUT $OUTPUT
