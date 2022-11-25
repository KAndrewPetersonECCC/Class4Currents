#!/bin/bash
# ord_soumet /fs/homeu2/eccc/mrd/ords/rpnenv/dpe000/Class4_Currents/jobscripts/test_corinne.sh -cpus 20 -mpi -cm 8000M -t 18000 -shell=/bin/bash

cd /fs/homeu2/eccc/mrd/ords/rpnenv/dpe000/Class4_Currents
source jobscripts/prepython.sh
python python/test_corinne.py

