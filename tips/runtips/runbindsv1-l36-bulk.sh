#!/bin/bash

# info:
# 1rack = 384nodes
# 1node = 48 physical cores
# whole: 330racks -> 384nodes * 330rack -> 126720 nodes -> 6082560 cores

# ---user input --- #
inname=$1
replacestr='xxx'
node=16
proc_per_node=2
ndigit=1 # zero padding
jobtime="05:00:00"  # "hour:minutes:seconds"
ABINIT_DIR='/vol0003/hp190133/data/programs/ABINIT-MP/open/rev23bindv1-20230712/lang-tcsds-1.2.39'
BINARY_NAME='abinitmp_smp'
rscgrp='small'
group='hp190133'
lang='lang/tcsds-1.2.39'
# ----user input end ---- #

totalproc=`echo "$node * $proc_per_node" | bc`
OMP_NUM_THREADS=`echo "48 / $proc_per_node" | bc`

OMP_STACKSIZE='8G'

# innameのreplacestr以降を削除
FHead=`echo $inname | sed -e "s/${replacestr}.*//"`
FTail=`echo $inname | sed -e "s/.*${replacestr}//" | sed -e "s/\..*//"`


# --- print section ---

echo """#!/bin/bash

#------- pjsub option -------#
#PJM -L \"rscgrp=${rscgrp}\"
#PJM -L \"node=$node\"
#PJM --mpi \"proc=$totalproc,max-proc-per-node=$proc_per_node\"
#PJM -L \"elapse=$jobtime\"
#PJM -g \"${group}\"
#PJM -L "freq=2200,eco_state=2"
#PJM -j

module switch ${lang}

ABINIT_DIR=${ABINIT_DIR}
BINARY_NAME=${BINARY_NAME}

export OMP_NUM_THREADS=${OMP_NUM_THREADS}
export OMP_STACKSIZE=${OMP_STACKSIZE}

num=\`printf \"%0${ndigit}d\" \${PJM_BULKNUM}\`
FILE_NAME=${FHead}\${num}${FTail}
NUM_CORE=_${node}n-${proc_per_node}p-${OMP_NUM_THREADS}t

AJF_NAME=\${FILE_NAME}.ajf
OUT_NAME=\${FILE_NAME}\${NUM_CORE}.log
ERR_NAME=\${FILE_NAME}\${NUM_CORE}.err


#------- Program execution -------#
\${ABINIT_DIR}/mkinp_bindsv1.py < \${AJF_NAME} > \${FILE_NAME}.inp
mpiexec -stdin \${FILE_NAME}.inp -stdout-proc \${OUT_NAME} -stderr-proc \${ERR_NAME} \${ABINIT_DIR}/\${BINARY_NAME}

mv -f \${OUT_NAME}.1.0 \${OUT_NAME}
cat \${ERR_NAME}.1.* > \${ERR_NAME} 2>/dev/null
rm -f \${ERR_NAME}.1.*
""" > ${FHead}${NUM_CORE}bulk.sh

# --- run pjsub ---
echo "--- run ${FHead}${NUM_CORE}bulk.sh job (pjsub)---"
echo "Please enter pjsub bulk job command."
echo """e.g.) pjsub --bulk --sparam "1-10" ${FHead}${NUM_CORE}bulk.sh"""
# pjsub --bulk --sparam "1-10" ${fhead}${NUM_CORE}.sh
