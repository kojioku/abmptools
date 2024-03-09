#!/bin/bash

# info:
# 1rack = 384nodes
# 1node = 48 physical cores
# whole: 330racks -> 384nodes * 330rack -> 126720 nodes -> 6082560 cores

# ---user input --- #
node=1
proc_per_node=1
jobtime="10:00:00"  # "hour:minutes:seconds"
BINARY_NAME='/path/to/cognac101'
rscgrp='small'
group='hpxxxxxx'
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
#PJM -x PJM_LLIO_GFSCACHE=/vol0002:/vol0003:/vol0004:/vol0005:/vol0006
#PJM -j

export OMP_NUM_THREADS=${OMP_NUM_THREADS}
export OMP_STACKSIZE=${OMP_STACKSIZE}

export OCTA84_HOME=${OCTA84_HOME}
export PATH=\$OCTA84_HOME/GOURMET:\$PATH
. ${OCTA84_HOME}/GOURMET/gourmetterm -

# export UDF_DEF_PATH="/path/to/OCTA84/ENGINES/udf"

#------- Program execution -------#
$BINARY_NAME -I ${IN_NAME} -O ${OUT_NAME} -n $OMP_NUM_THREADS
""" > ${fhead}${NUM_CORE}.sh

# --- run pjsub ---
echo "--- run ${fhead}${NUM_CORE}.sh job (pjsub)---"
pjsub ${fhead}${NUM_CORE}.sh

