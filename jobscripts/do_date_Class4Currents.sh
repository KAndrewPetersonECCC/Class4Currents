#!/bin/bash -x
#bash /fs/homeu2/eccc/mrd/ords/rpnenv/dpe000/Class4_Currents/jobscripts/date_Class4Currents.sh --start=CCYYMMDD --final=CCYYMMDD

WDIR=/fs/homeu2/eccc/mrd/ords/rpnenv/dpe000/Class4_Currents
USAGE="USAGE:  date_Class4Currents.sh -s=CCYYMMDD -f=CCYYMMDD"

NDAYS=0
ORIG=True
GEPS=False
GIOP=False
ENAN=False
ANAL=False
ASSG=False
ASEA=False
FCEX=GEPS_STO2X
EXPT=""
DDIR="DREW"
IGNORE=False
PLOT="-q"
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
    -n=*|--ndays=*)
    NDAYS="${i#*=}"
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
    -u|--not_filter)
    FILTER="--not_filter"
    shift # past argument=value
    ;;
    -g|--geps)
    GEPS=True
    ORIG=False
    shift # past argument=value
    ;;
    -z|--giops|--GIOPS)
    GIOP=True
    ORIG=False
    shift # past argument=value
    ;;
    -A|--ENAN|--Ensemble_Analysis|--EA)
    ENAN=True
    ORIG=False
    shift # past argument=value
    ;;
    --anal|--only_analysis|--analysis)
    ANAL=True
    ORIG=False
    shift # past argument=value
    ;;
    -F=*|--forecast=*|--fcst=*)
    FCEX="${i#*=}"
    shift # past argument=value
    ;;
    -E=*|--experiment=*|--expt=*)
    EXPT="${i#*=}"
    shift # past argument=value
    ;;
    -D=*|--data_dir=*|--ddir=*)
    DDIR="${i#*=}"
    shift # past argument=value
    ;;
    -a|--assemble|--assemble_geps)
    ASSG=True
    ORIG=False
    shift # past argument=value
    ;;
    -b|--assembleEA|--assemble_ensemble|--EB)
    ASEA=True
    ORIG=False
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
    echo ${NDAYS} ${FINALDATE}
if [[ ${NDAYS} -gt 0 ]]; then 
    FINALDATE=$(date -u --date "${STARTDATE} + ${NDAYS} days" +%Y%m%d)
    echo ${NDAYS} ${FINALDATE}
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
        if [[ ${ORIG} == True ]]; then
	        echo "jobscripts/date_Class4Currents.sh  --date=${DATE} ${PLOT} ${CHARLY} ${FILTER} --next=${NEXT} -s"
	        bash ${WDIR}/jobscripts/date_Class4Currents.sh  --date=${DATE} ${PLOT} ${CHARLY} ${FILTER} --next=${NEXT} -s
	    fi
        if [[ ${GIOP} == True ]] ; then 
	        echo "jobscripts/giops_Class4Currents.sh  --date=${DATE} ${CHARLY} ${FILTER} --next=${NEXT} -s"
            bash ${WDIR}/jobscripts/giops_Class4Currents.sh  --date=${DATE} ${CHARLY} ${FILTER} --next=${NEXT} -s
        fi
	    if [[ ${GEPS} == True ]] ; then 
	        for ie in 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 ; do 
	            echo "jobscripts/GEPS_Class4Currents.sh  --date=${DATE} -e=${ie} ${FILTER} --next=${NEXT} -s"
	            bash ${WDIR}/jobscripts/GEPSdate_Class4Currents.sh  --date=${DATE} -e=${ie} ${FILTER} --next=${NEXT} -s
            done
	    fi
	    if [[ ${ENAN} == True ]] ; then 
	        for ie in 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 ; do 
	            echo "jobscripts/GEPS_Class4Currents.sh  --EA --fcst=${FCEX} --date=${DATE} -e=${ie} ${FILTER} --next=${NEXT} -s"
	            bash ${WDIR}/jobscripts/GEPSdate_Class4Currents.sh  --EA --fcst=${FCEX} --date=${DATE} -e=${ie} ${FILTER} --next=${NEXT} -s
            done
	    fi
	    if [[ ${ANAL} == True ]] ; then 
	        echo "jobscripts/ANALdate_Class4Currents.sh  --expt=${EXPT} --date=${DATE} ${FILTER} --next=${NEXT} --ddir=${DDIR} -s"
	        bash ${WDIR}/jobscripts/ANALdate_Class4Currents.sh  --expt=${EXPT} --date=${DATE} ${FILTER} --next=${NEXT} --ddir=${DDIR} -s
	    fi
	    if [[ ${ASSG} == True ]] ; then 
	        echo "jobscripts/GEPS_Class4Currents.sh  --date=${DATE} -e=A ${FILTER} --next=${NEXT} -s"
	        bash ${WDIR}/jobscripts/GEPSdate_Class4Currents.sh  --date=${DATE} -e=A ${FILTER} --next=${NEXT} -s
	    fi	       
	    if [[ ${ASEA} == True ]] ; then 
	        echo "jobscripts/GEPS_Class4Currents.sh  --EA --fcst=${FCEX} --date=${DATE} -e=A ${FILTER} --next=${NEXT} -s"
	        bash ${WDIR}/jobscripts/GEPSdate_Class4Currents.sh  --EA --fcst=${FCEX} --date=${DATE} -e=A ${FILTER} --next=${NEXT} -s
	    fi	       
	    ZERO=False
    else
        echo "Preparing Other Jobs"
        if [[ ${NEXTTIME} -le ${FINALTIME} ]]; then 
            echo "${GEPS} == False  ${ASSG} == False ${ENAN} == False  ${ASEA} == False ${ANAL} == False"
	        #if [[ ${GEPS} == False && ${ASSG} == False  && ${ENAN} == False && ${ASEA} == False && ${ANAL} == False ]] ; then 
            if [[ ${ORIG} == True ]]; then
	            echo "jobscripts/date_Class4Currents.sh  --date=${DATE} ${PLOT} ${CHARLY} ${FILTER} --next=${NEXT}"
                bash ${WDIR}/jobscripts/date_Class4Currents.sh  --date=${DATE} ${PLOT} ${CHARLY} ${FILTER} --next=${NEXT}
            fi
            if [[ ${GIOP} == True ]] ; then 
	            echo "jobscripts/giops_Class4Currents.sh  --date=${DATE} ${CHARLY} ${FILTER} --next=${NEXT}"
                bash ${WDIR}/jobscripts/giops_Class4Currents.sh  --date=${DATE} ${CHARLY} ${FILTER} --next=${NEXT}
            fi
 	        if [[ ${GEPS} == True ]] ; then 
	            for ie in 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 ; do 
	                echo "jobscripts/GEPS_Class4Currents.sh  --date=${DATE} -e=${ie} ${FILTER} --next=${NEXT}"
	                bash ${WDIR}/jobscripts/GEPSdate_Class4Currents.sh  --date=${DATE} -e=${ie} ${FILTER} --next=${NEXT}
                done
            fi
 	        if [[ ${ENAN} == True ]] ; then 
	            for ie in 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 ; do 
	                echo "jobscripts/GEPS_Class4Currents.sh  --EA --fcst=${FCEX} --date=${DATE} -e=${ie} ${FILTER} --next=${NEXT}"
	                bash ${WDIR}/jobscripts/GEPSdate_Class4Currents.sh  --EA --fcst=${FCEX} --date=${DATE} -e=${ie} ${FILTER} --next=${NEXT}
                done
            fi
	        if [[ ${ANAL} == True ]] ; then 
	            echo "jobscripts/ANALdate_Class4Currents.sh  --expt=${EXPT} --date=${DATE} ${FILTER} --next=${NEXT} --ddir=${DDIR} "
	            bash ${WDIR}/jobscripts/ANALdate_Class4Currents.sh  --expt=${EXPT} --date=${DATE} ${FILTER} --next=${NEXT} --ddir=${DDIR}
	        fi
	        if [[ ${ASSG} == True ]] ; then 
	            echo "jobscripts/GEPS_Class4Currents.sh  --date=${DATE} -e=A ${FILTER} --next=${NEXT}"
	            bash ${WDIR}/jobscripts/GEPSdate_Class4Currents.sh  --date=${DATE} -e=A ${FILTER} --next=${NEXT} 
	        fi	       
	        if [[ ${ASEA} == True ]] ; then 
	            echo "jobscripts/GEPS_Class4Currents.sh  --EA --fcst=${FCEX} --date=${DATE} -e=A ${FILTER} --next=${NEXT}"
	            bash ${WDIR}/jobscripts/GEPSdate_Class4Currents.sh  --EA --fcst=${FCEX} --date=${DATE} -e=A ${FILTER} --next=${NEXT} 
	        fi	       
	    else
	        if [[ ${GEPS} == False && ${ASSG} == False && ${ENAN} == False && ${ASEA} == False && ${ANAL} == False ]] ; then 
	            echo "jobscripts/date_Class4Currents.sh  --date=${DATE} ${PLOT} ${CHARLY} ${FILTER}"
	            bash ${WDIR}/jobscripts/date_Class4Currents.sh  --date=${DATE} ${PLOT} ${CHARLY} ${FILTER}
            fi
 	        if [[ ${GEPS} == True ]] ; then 
	            for ie in 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 ; do 
	                echo "jobscripts/GEPS_Class4Currents.sh  --date=${DATE} -e=${ie} ${FILTER}"
	                bash ${WDIR}/jobscripts/GEPSdate_Class4Currents.sh  --date=${DATE} -e=${ie} ${FILTER}
                done
	        fi
 	        if [[ ${ENAN} == True ]] ; then 
	            for ie in 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 ; do 
	                echo "jobscripts/GEPSdate_Class4Currents.sh  --EA --fcst=${FCEX} --date=${DATE} -e=${ie} ${FILTER}"
	                bash ${WDIR}/jobscripts/GEPSdate_Class4Currents.sh   --EA --fcst=${FCEX} --date=${DATE} -e=${ie} ${FILTER}
                done
	        fi
	        if [[ ${ANAL} == True ]] ; then 
	            echo "jobscripts/ANALdate_Class4Currents.sh  --expt=${EXPT} --date=${DATE} ${FILTER} --ddir=${DDIR} "
	            bash ${WDIR}/jobscripts/ANALdate_Class4Currents.sh  --expt=${EXPT} --date=${DATE} ${FILTER} --ddir=${DDIR}
	        fi
 	        if [[ ${ASSG} == True ]] ; then 
	            echo "jobscripts/GEPSdate_Class4Currents.sh  --date=${DATE} -e=A ${FILTER}"
	            bash ${WDIR}/jobscripts/GEPSdate_Class4Currents.sh  --date=${DATE} -e=A ${FILTER} 
            fi
 	        if [[ ${ASEA} == True ]] ; then 
	           echo "jobscripts/GEPSdate_Class4Currents.sh  --EA --fcst=${FCEX} --date=${DATE} -e=A ${FILTER}"
	           bash ${WDIR}/jobscripts/GEPSdate_Class4Currents.sh  --EA --fcst=${FCEX} --date=${DATE} -e=A ${FILTER} 
            fi
	    fi
    fi
     
    DATE=$(date -u --date "${DATE} + 1 day" +%Y%m%d)
    THISTIME=$(date -u --date ${DATE} +%s)
    echo ${DATE} ${THISTIME}
    
done
