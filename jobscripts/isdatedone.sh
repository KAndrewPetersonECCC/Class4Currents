#!/bin/bash

date=$1
datefile=/fs/homeu1/eccc/mrd/ords/rpnenv/dpe000/Class4_Currents/DATES/todo.date

source /fs/homeu1/eccc/mrd/ords/rpnenv/dpe000/Class4_Currents/jobscripts/preconda.sh > /dev/null 2>&1
source activate metcarto

python << EOP
import os
import sys
sys.path.insert(0, '/home/dpe000/Class4_Currents/python')
import datadatefile

datefile='${datefile}'
intdate=${date}
TorF=datadatefile.get_from_boolean_file(intdate, file=datefile)
print(TorF)
EOP
