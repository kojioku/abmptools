#!/bin/bash

# info:
# 1rack = 384nodes
# 1node = 48 physical cores
# whole: 330racks -> 384nodes * 330rack -> 126720 nodes -> 6082560 cores

# ---user input --- #
node=192
proc_per_node=1
jobtime="04:00:00"  # "hour:minutes:seconds"
ABINIT_DIR=/data/hp210026/Expt-2-Ver2-Rev5-wo-fint
BINARY_NAME=abinitmp_smp
rscgrp='small'
group='hp210026'
# ----user input end ---- #

fhead=${1%.*}
totalproc=`echo "$node * $proc_per_node" | bc`
OMP_NUM_THREADS=`echo "48 / $proc_per_node" | bc`

OMP_STACKSIZE='256M'

FILE_NAME=${fhead}
NUM_CORE=_${node}n-${proc_per_node}p-${OMP_NUM_THREADS}t

AJF_NAME=${FILE_NAME}.ajf
OUT_NAME=${FILE_NAME}${NUM_CORE}.log
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

module switch lang/tcsds-1.2.31

export OMP_NUM_THREADS=${OMP_NUM_THREADS}
export OMP_STACKSIZE=${OMP_STACKSIZE}

#------- Program execution -------#
# ${ABINIT_DIR}/mkinp_openver1rev20.py < ${AJF_NAME} > ${FILE_NAME}.inp
mpiexec -stdin ${FILE_NAME}.ajf -stdout ${OUT_NAME} -stderr ${ERR_NAME} ${ABINIT_DIR}/${BINARY_NAME}
""" > ${fhead}${NUM_CORE}.sh

# --- run pjsub ---
echo "--- run ${fhead}${NUM_CORE}.sh job (pjsub)---"
pjsub ${fhead}${NUM_CORE}.sh

