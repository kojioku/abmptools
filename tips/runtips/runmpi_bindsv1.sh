#!/bin/bash
if [ "$#" -eq 0 ]; then
    echo "Usage: $0 [arguments](.ajf)"
    exit 1
fi

# ---user input --- #
node=4 #4
jobtime="20:00:00"  # "hour:minutes:seconds"
ABINIT_DIR=/sqfs/work/K23A03/programs/ABINIT-MP/bindv1-20230712/BaseCPU2023
BINARY_NAME=abinitmp
MKINP=mkinp_bindsv1.py
class='SQUID'
group='K23A03'
lang='BaseCPU/2023'

fhead=${1%.*}
totalproc=`echo "$node * 76" | bc`

# ----user input end ---- #

FILE_NAME=${fhead}
NUM_CORE=_${node}n-flatmpi

AJF_NAME=${FILE_NAME}.ajf
OUT_NAME=${FILE_NAME}${NUM_CORE}.log
ERR_NAME=${FILE_NAME}${NUM_CORE}.err

echo """ #!/bin/bash
#PBS -q ${class}
#PBS -l cpunum_job=76
#PBS --group=${group}
#PBS -l elapstim_req=${jobtime}
#PBS -b ${node}
#PBS -T intmpi
 
#PBS -v NUM_PROCS=${totalproc}
 
#PBS -v ABINIT_DIR=${ABINIT_DIR}
#PBS -v BINARY_NAME=${BINARY_NAME}
#PBS -v MKINP=${MKINP}
#PBS -v FILE_NAME=${FILE_NAME}
#PBS -v OUT_NAME=${FILE_NAME}${NUM_CORE}
ulimit –s unlimited
 
module load ${lang}

cd \${PBS_O_WORKDIR}
python \${ABINIT_DIR}/\${MKINP} < \${FILE_NAME}.ajf > \${FILE_NAME}.inp
mpirun \${NQSV_MPIOPTS} -np \${NUM_PROCS} \${ABINIT_DIR}/\${BINARY_NAME} < \${FILE_NAME}.inp > \${OUT_NAME}.log
""" > ${FILE_NAME}${NUM_CORE}.sh

# --- run qsub ---
echo "--- run ${FILE_NAME}${NUM_CORE}.sh job (qsub)---"
qsub ${FILE_NAME}${NUM_CORE}.sh

