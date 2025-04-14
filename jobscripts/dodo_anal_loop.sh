
bash jobscripts/do_analysis_loop.sh -R=GX340-U2-CONTROLE-LR -E=dev-3.4.0-gx-7dIAU -s=20210105 -f=20220204 -L="7d IAU" --submit
bash jobscripts/do_analysis_loop.sh -R=GX340-U2-CONTROLE-LR -E=dev-3.4.0-gx-7dIAU -s=20190402 -f=20200331 -L="7d IAU" --submit
bash jobscripts/do_analysis_loop.sh -R=gx_dev-3.4.0_DA_y2 -E=gx_dev-3.4.0_DA_y2_MDT2 -s=20190402 -f=20200331 -L="MDT" --submit
bash jobscripts/do_analysis_loop.sh -R=gx_dev-3.4.0_DA_y2 -E=gx_dev-3.4.0_DA_y2_core -s=20190402 -f=20200331 -L="core" --submit
bash jobscripts/do_analysis_loop.sh -R=GX340-U2-CONTROLE-LR -E=gx_dev-3.4.0_DA_y2_core -s=20190402 -f=20200331 -L="core core" --submit

bash jobscripts/do_analysis_loop.sh -E=CLASS4_currents_ENAN_FILT -R=CLASS4_currents_ENAN_FILT -p=ENAN_orca025_currents.f2 -q=ENAN_orca025_currents.f2.000 -s=20210609 -f=20220524 -L="ENSEMBLE"  -C="CONTROL" --ensemble="[2,2]" --submit
bash jobscripts/do_analysis_loop.sh -E=CLASS4_currents_ENAN_UFIL -R=CLASS4_currents_ENAN_UFIL -p=ENAN_orca025_currents.2 -q=ENAN_orca025_currents.2.000  -s=20210609 -f=20220524 -L="UNFILTERED" -C="CONTROL"  --ensemble="[2,2]" --submit
bash jobscripts/do_analysis_loop.sh -E=CLASS4_currents_ENAN_UFIL -R=CLASS4_currents_ENAN_FILT -p=ENAN_orca025_currents.2 -q=ENAN_orca025_currents.f2  -s=20210609 -f=20220524 -L="UNFILTERED"  -C="FILTERED"  --ensemble="[2,2]" --submit

bash jobscripts/do_analysis_loop.sh -E=CLASS4_currents_ENAN_30_FILT -R=CLASS4_currents_ENAN_FILT -p=ENAN_30_orca025_currents.f2.000 -q=ENAN_orca025_currents.f2.000 -s=20210609 -f=20220524 -L="30X-SHAPIRO"  -C="CONTROL" --ensemble="[2,2]" --submit
bash jobscripts/do_analysis_loop.sh -E=CLASS4_currents_ENAN_FILT -R=CLASS4_currents_ENAN_FILT -p=ENAN_orca025_currents.f2 -q=ENAN_orca025_currents.f2.000 -s=20210609 -f=20220524 -L="ENSEMBLE"  -C="CONTROL" --ensemble="[2,2]" --submit
bash jobscripts/do_analysis_loop.sh -R=CLASS4_currents_ENAN_30_FILT -E=CLASS4_currents_ENAN_FILT -q=ENAN_30_orca025_currents.f2.000 -p=ENAN_orca025_currents.f2 -s=20210609 -f=20220524 -C="30X-SHAPIRO"  -L="ENSEMBLE" --ensemble="[2,2]" --submit

bash jobscripts/do_analysis_loop.sh -E=CLASS4_currents_ENAN_30_FILT -R=CLASS4_currents_ENAN_FILT -p=ENAN_30_orca025_currents.f2.000 -q=ENAN_orca025_currents.f2.000 -s=20211002 -f=20220524 -L="30X-SHAPIRO"  --C="CONTROL" --ensemble="[2,2]" --submit
bash jobscripts/do_analysis_loop.sh -E=CLASS4_currents_ENAN_FILT -R=CLASS4_currents_ENAN_FILT -p=ENAN_orca025_currents.f2 -q=ENAN_orca025_currents.f2.000 -s=20211002 -f=20220524 -L="ENSEMBLE"  --C="CONTROL" --ensemble="[2,2]" --submit
bash jobscripts/do_analysis_loop.sh -R=CLASS4_currents_ENAN_30_FILT -E=CLASS4_currents_ENAN_FILT -q=ENAN_30_orca025_currents.f2.000 -p=ENAN_orca025_currents.f2 -s=20211002 -f=20220524 -C="30X-SHAPIRO"  -L="ENSEMBLE" --ensemble="[2,2]" --submit

bash jobscripts/do_analysis_loop.sh -R=Oper -E=HalfArgo -s=0 -f=0 -L="HalfArgo" -C="Oper" --submit
bash jobscripts/do_analysis_loop.sh -R=Oper -E=Free -s=0 -f=0 -L="Free" -C="Oper" --submit

bash jobscripts/do_analysis_loop.sh -R=CNTL -E=NoArgo -s=0 -f=0 -L="NoArgo" -C="CNTL" --submit

bash jobscripts/do_analysis_loop.sh -R=IC4GX340EH22-CTR -E=GX35FCH22V2 -s=20211001 -f=20220331 -L="IC4-H22" -C="IC3" --submit
bash jobscripts/do_analysis_loop.sh -R=IC4GX340EH22-CTR -E=GX35FCE22V2 -s=20220601 -f=20220831 -L="IC4-E22" -C="IC3" --submit

#for expt in NoArgoV2 Oper HalfArgoV2 NoAltV2 NoSSTV2 SSTonlyV2 NoInsituV2 NoMoorV2 Free ; do 
for expt in NoArgoV2 Oper HalfArgoV2 NoAltV2 SSTonlyV3 NoSSTV2 NoInsituV3 NoMoorV2 Free ; do 
    bash jobscripts/do_analysis_loop.sh -R=CNTLV2 -E=${expt} -s=20200101 -f=20221231 --submit
done
for expt in NoArgoV2 Oper HalfArgoV2 NoAltV2 SSTonlyV3 NoSSTV2 NoInsituV3 NoMoorV2 CNTLV2 ; do 
    bash jobscripts/do_analysis_loop.sh -R=Free -E=${expt} -s=20200101 -f=20221231 --submit
done
