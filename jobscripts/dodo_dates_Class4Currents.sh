#!/bin/bash
#bash /fs/homeu1/eccc/mrd/ords/rpnenv/dpe000/Class4_Currents/jobscripts/date_Class4Currents.sh --start=CCYYMMDD --final=CCYYMMDD

WDIR=/fs/homeu1/eccc/mrd/ords/rpnenv/dpe000/Class4_Currents
cd ${WDIR}

FILTER="--filter"
PLOT="--plot"

FILTER="--not_filter"
PLOT="--no_plot"
ssh eccc-ppp4 bash -l ${WDIR}/jobscripts/do_date_Class4Currents.sh -s=20190401 -f=20190430 ${PLOT} ${FILTER} 
ssh eccc-ppp4 bash -l ${WDIR}/jobscripts/do_date_Class4Currents.sh -s=20190601 -f=20190630 ${PLOT} ${FILTER} 
ssh eccc-ppp4 bash -l ${WDIR}/jobscripts/do_date_Class4Currents.sh -s=20190801 -f=20190831 ${PLOT} ${FILTER} 
ssh eccc-ppp4 bash -l ${WDIR}/jobscripts/do_date_Class4Currents.sh -s=20191001 -f=20191031 ${PLOT} ${FILTER} 
ssh eccc-ppp4 bash -l ${WDIR}/jobscripts/do_date_Class4Currents.sh -s=20191201 -f=20191231 ${PLOT} ${FILTER} 
ssh eccc-ppp4 bash -l ${WDIR}/jobscripts/do_date_Class4Currents.sh -s=20200101 -f=20200131 ${PLOT} ${FILTER} 
ssh eccc-ppp4 bash -l ${WDIR}/jobscripts/do_date_Class4Currents.sh -s=20200201 -f=20200229 ${PLOT} ${FILTER} 
ssh eccc-ppp4 bash -l ${WDIR}/jobscripts/do_date_Class4Currents.sh -s=20200401 -f=20200430 ${PLOT} ${FILTER} 
ssh eccc-ppp4 bash -l ${WDIR}/jobscripts/do_date_Class4Currents.sh -s=20200601 -f=20200630 ${PLOT} ${FILTER} 
ssh eccc-ppp4 bash -l ${WDIR}/jobscripts/do_date_Class4Currents.sh -s=20200801 -f=20200831 ${PLOT} ${FILTER} 
ssh eccc-ppp4 bash -l ${WDIR}/jobscripts/do_date_Class4Currents.sh -s=20201001 -f=20201031 ${PLOT} ${FILTER} 
ssh eccc-ppp4 bash -l ${WDIR}/jobscripts/do_date_Class4Currents.sh -s=20201201 -f=20201231 ${PLOT} ${FILTER} 
ssh eccc-ppp4 bash -l ${WDIR}/jobscripts/do_date_Class4Currents.sh -s=20210201 -f=20210228 ${PLOT} ${FILTER} 
ssh eccc-ppp4 bash -l ${WDIR}/jobscripts/do_date_Class4Currents.sh -s=20210401 -f=20210430 ${PLOT} ${FILTER} 
ssh eccc-ppp4 bash -l ${WDIR}/jobscripts/do_date_Class4Currents.sh -s=20210801 -f=20210831 ${PLOT} ${FILTER} 
ssh eccc-ppp4 bash -l ${WDIR}/jobscripts/do_date_Class4Currents.sh -s=20211001 -f=20211031 ${PLOT} ${FILTER} 
ssh eccc-ppp4 bash -l ${WDIR}/jobscripts/do_date_Class4Currents.sh -s=20211201 -f=20211231 ${PLOT} ${FILTER} 
ssh eccc-ppp4 bash -l ${WDIR}/jobscripts/do_date_Class4Currents.sh -s=20220201 -f=20220228 ${PLOT} ${FILTER} 

ssh eccc-ppp4 bash -l ${WDIR}/jobscripts/do_date_Class4Currents.sh -s=20190501 -f=20190531 ${PLOT} ${FILTER} 
ssh eccc-ppp3 bash -l ${WDIR}/jobscripts/do_date_Class4Currents.sh -s=20190701 -f=20190731 ${PLOT} ${FILTER}
ssh eccc-ppp3 bash -l ${WDIR}/jobscripts/do_date_Class4Currents.sh -s=20190901 -f=20190930 ${PLOT} ${FILTER}
ssh eccc-ppp3 bash -l ${WDIR}/jobscripts/do_date_Class4Currents.sh -s=20191101 -f=20191130 ${PLOT} ${FILTER}
ssh eccc-ppp3 bash -l ${WDIR}/jobscripts/do_date_Class4Currents.sh -s=20200101 -f=20200131 ${PLOT} ${FILTER}
ssh eccc-ppp3 bash -l ${WDIR}/jobscripts/do_date_Class4Currents.sh -s=20200301 -f=20200331 ${PLOT} ${FILTER}
ssh eccc-ppp3 bash -l ${WDIR}/jobscripts/do_date_Class4Currents.sh -s=20200501 -f=20200531 ${PLOT} ${FILTER}
ssh eccc-ppp3 bash -l ${WDIR}/jobscripts/do_date_Class4Currents.sh -s=20200701 -f=20200731 ${PLOT} ${FILTER}
ssh eccc-ppp3 bash -l ${WDIR}/jobscripts/do_date_Class4Currents.sh -s=20200901 -f=20200930 ${PLOT} ${FILTER}
ssh eccc-ppp3 bash -l ${WDIR}/jobscripts/do_date_Class4Currents.sh -s=20201101 -f=20201130 ${PLOT} ${FILTER}
ssh eccc-ppp3 bash -l ${WDIR}/jobscripts/do_date_Class4Currents.sh -s=20210101 -f=20210131 ${PLOT} ${FILTER}
ssh eccc-ppp3 bash -l ${WDIR}/jobscripts/do_date_Class4Currents.sh -s=20210301 -f=20210331 ${PLOT} ${FILTER}
ssh eccc-ppp3 bash -l ${WDIR}/jobscripts/do_date_Class4Currents.sh -s=20210501 -f=20210531 ${PLOT} ${FILTER}
ssh eccc-ppp3 bash -l ${WDIR}/jobscripts/do_date_Class4Currents.sh -s=20210701 -f=20210731 ${PLOT} ${FILTER}
ssh eccc-ppp3 bash -l ${WDIR}/jobscripts/do_date_Class4Currents.sh -s=20210901 -f=20210930 ${PLOT} ${FILTER}
ssh eccc-ppp3 bash -l ${WDIR}/jobscripts/do_date_Class4Currents.sh -s=20211101 -f=20211130 ${PLOT} ${FILTER}
ssh eccc-ppp3 bash -l ${WDIR}/jobscripts/do_date_Class4Currents.sh -s=20220101 -f=20220131 ${PLOT} ${FILTER}

bash /fs/homeu2/eccc/mrd/ords/rpnenv/dpe000/Class4_Currents/jobscripts/do_date_Class4Currents.sh --start=2022 --final=${year}${MN}31 --filter --anal --expt=${expt} --ddir=SynObs
