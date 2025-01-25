#!/bin/bash

#------- pjsub option -------#
#PJM -g hp190133
#PJM -L "rscgrp=small"
#PJM -L "node=16"
#PJM -L "elapse=72:00:00"
#PJM --mpi "proc=192,max-proc-per-node=12"
#PJM -x PJM_LLIO_GFSCACHE=/vol0004
#PJM -S

export OMP_NUM_THREADS=4

. /vol0004/apps/oss/spack/share/spack/setup-env.sh
spack load /p6ktu7p
export LD_LIBRARY_PATH=/lib64:${LD_LIBRARY_PATH}

llio_transfer `which gmx_mpi`

NAME="dpdFEL_PVPVATPGS_WAT70_forReverse_init_out_1-3del_mod_rec201"
# OPTION="-stdout-proc ./output.%j/stdout -stderr-proc ./output.%j/stderr"
mpi='mpiexec'

f1=step1-min.mdp
f2=step2-heat.mdp
f3=step3-npt1.mdp
f4=step4-npt2.mdp
f5=step5-npt3.mdp
# f4=step4-equil.mdp
# f4=step5-prod.mdp

# opt 全原子最適化
gmx grompp -f $f1 -c ${NAME}.gro -r ${NAME}.gro -p ${NAME}.top -o ${NAME}_min.tpr -maxwarn 1
${mpi} -np ${PJM_MPI_PROC} gmx_mpi mdrun -deffnm ${NAME}_min

# # heat system
gmx grompp -f $f2 -c ${NAME}_min.gro -r ${NAME}_min.gro -p ${NAME}.top -o ${NAME}_out2.tpr -maxwarn 1
${mpi} -np ${PJM_MPI_PROC} gmx_mpi mdrun  -deffnm ${NAME}_out2
# 
# # opt density
gmx grompp -f $f3 -c ${NAME}_out2.gro -t ${NAME}_out2.cpt -r ${NAME}_out2.gro -p ${NAME}.top -o ${NAME}_out3.tpr -maxwarn 1
${mpi} -np ${PJM_MPI_PROC} gmx_mpi mdrun  -deffnm ${NAME}_out3
# 
# # equil or prod
gmx grompp -f $f4 -c ${NAME}_out3.gro -t ${NAME}_out3.cpt -r ${NAME}_out3.gro -p ${NAME}.top -o ${NAME}_out4.tpr -maxwarn 1
${mpi} -np ${PJM_MPI_PROC} gmx_mpi mdrun  -deffnm ${NAME}_out4
#
# # equil or prod
gmx grompp -f $f5 -c ${NAME}_out4.gro -t ${NAME}_out4.cpt -r ${NAME}_out4.gro -p ${NAME}.top -o ${NAME}_out5.tpr -maxwarn 1
${mpi} -np ${PJM_MPI_PROC} gmx_mpi mdrun  -deffnm ${NAME}_out5

# prod
# ${mpi} ${OPTION} gmx_mpi grompp -f $f5 -c ${NAME}_out4.gro -r ${NAME}_out4.gro -p ${NAME}.top -o ${NAME}_out5.tpr -maxwarn 1
# ${mpi} -np ${np} ${OPTION} gmx_mpi mdrun -ntomp ${NTOMP} -deffnm ${NAME}_out5
#
#  -cpi ${NAME}_pr.cpt
