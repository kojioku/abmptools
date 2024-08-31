#!/bin/bash

# info:
# 1rack = 384nodes
# 1node = 48 physical cores
# whole: 330racks -> 384nodes * 330rack -> 126720 nodes -> 6082560 cores

# ---user input --- #
node=12
proc_per_node=2
jobtime="04:00:00"  # "hour:minutes:seconds"
ABINIT_DIR=/vol0003/hp190133/data/programs/ABINIT-MP/open/spack2-4
BINARY_NAME=abinitmp_smp
rscgrp='small'
group='hp190133'
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
#PJM -L "freq=2200,eco_state=2"
#PJM -x PJM_LLIO_GFSCACHE=/vol0004:/vol0003
#PJM -j

. /vol0004/apps/oss/spack/share/spack/setup-env.sh
spack load abinitmp@2-4

export OMP_NUM_THREADS=${OMP_NUM_THREADS}
export OMP_STACKSIZE=${OMP_STACKSIZE}

#------- Program execution -------#
${ABINIT_DIR}/mkinp_openver1rev20.py < ${AJF_NAME} > ${FILE_NAME}.inp
mpiexec -stdin ${FILE_NAME}.inp -stdout-proc ${OUT_NAME} -stderr-proc ${ERR_NAME} ${BINARY_NAME}
mv -f ${OUT_NAME}.1.0 ${OUT_NAME}
cat ${ERR_NAME}.1.* > ${ERR_NAME} 2>/dev/null
rm -f ${ERR_NAME}.1.*
""" > ${fhead}${NUM_CORE}.sh

# --- run pjsub ---
echo "--- run ${fhead}${NUM_CORE}.sh job (pjsub)---"
pjsub ${fhead}${NUM_CORE}.sh
