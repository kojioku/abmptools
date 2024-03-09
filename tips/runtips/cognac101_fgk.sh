#!/bin/bash

# info:
# 1rack = 384nodes
# 1node = 48 physical cores
# whole: 330racks -> 384nodes * 330rack -> 126720 nodes -> 6082560 cores

# ---user input --- #
node=1
proc_per_node=1
jobtime="00:10:00"  # "hour:minutes:seconds"
BINARY_NAME=/vol0003/hp190133/data/programs/OCTA84-fgk/ENGINES/bin/A64FX/cognac101
rscgrp='small'
group='hp190133'
# ----user input end ---- #

fhead=${1%.*}
totalproc=`echo "$node * $proc_per_node" | bc`
OMP_NUM_THREADS=`echo "48 / $proc_per_node" | bc`

OMP_STACKSIZE='256M'

FILE_NAME=${fhead}
NUM_CORE=_${node}n-${OMP_NUM_THREADS}t

IN_NAME=$1
OUT_NAME=${FILE_NAME}${NUM_CORE}.bdf
ERR_NAME=${FILE_NAME}${NUM_CORE}.err

# --- print section ---

echo """#!/bin/bash

#------- pjsub option -------#
#PJM -L \"rscgrp=${rscgrp}\"
#PJM -L \"node=$node\"
#PJM --mpi \"proc=$totalproc,max-proc-per-node=$proc_per_node\"
#PJM -L \"elapse=$jobtime\"
#PJM -g \"${group}\"
#PJM -j

export OMP_NUM_THREADS=${OMP_NUM_THREADS}
export OMP_STACKSIZE=${OMP_STACKSIZE}

export OCTA84_HOME=/vol0003/hp190133/data/programs/OCTA84
export PATH=\$OCTA84_HOME/GOURMET:\$PATH
. /vol0003/hp190133/data/programs/OCTA84/GOURMET/gourmetterm -

#------- Program execution -------#
$BINARY_NAME -I ${IN_NAME} -O ${OUT_NAME} -n $OMP_NUM_THREADS
""" > ${fhead}${NUM_CORE}.sh

# --- run pjsub ---
echo "--- run ${fhead}${NUM_CORE}.sh job (pjsub)---"
pjsub ${fhead}${NUM_CORE}.sh

