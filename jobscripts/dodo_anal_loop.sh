
bash jobscripts/do_analysis_loop.sh -R=GX340-U2-CONTROLE-LR -E=dev-3.4.0-gx-7dIAU -s=20210105 -f=20220204 -L="7d IAU" --submit
bash jobscripts/do_analysis_loop.sh -R=GX340-U2-CONTROLE-LR -E=dev-3.4.0-gx-7dIAU -s=20190402 -f=20200331 -L="7d IAU" --submit
bash jobscripts/do_analysis_loop.sh -R=gx_dev-3.4.0_DA_y2 -E=gx_dev-3.4.0_DA_y2_MDT2 -s=20190402 -f=20200331 -L="MDT" --submit
bash jobscripts/do_analysis_loop.sh -R=gx_dev-3.4.0_DA_y2 -E=gx_dev-3.4.0_DA_y2_core -s=20190402 -f=20200331 -L="core" --submit
bash jobscripts/do_analysis_loop.sh -R=GX340-U2-CONTROLE-LR -E=gx_dev-3.4.0_DA_y2_core -s=20190402 -f=20200331 -L="core core" --submit
