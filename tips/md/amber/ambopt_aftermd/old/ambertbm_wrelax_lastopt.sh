#!/bin/sh
# Date: 2020/4/30
# Author: Koji Okuwaki
# ver2.4: fix .in file number
# ver2.3: fix prodction run ps bug
# ver2.2: modify heat interval
# ver2.1: Calculations up to equilibration is applied to the MOE.
# ver2.0: fix for setting cutoff, density, temp, and mask
set -e

# -- user setting --
fbase=$1
np=28

# -------------------------

mpi="mpirun -np $np"
ambermin="pmemd.MPI"
prmtop=$fbase.z.prmtop

t0=.a.0
rstrt=rstrt
n=1

#$ -l f_node=1                 ## * [資源タイプ名]=[使用ノード数個数]
#$ -l h_rt=0:30:00             ## * 経過時間 -> Max 24時間までです。
#$ -p -3                       ## * ジョブの優先度 *デフォルトは-5です。

# カレントディレクトリに移動
#$ -cwd

# Moduleコマンドの初期化
. /etc/profile.d/modules.sh

# CUDA 環境の読み込み
module load amber/16up10_cuda

# Intel Compiler環境の読み込み
module load intel/18.0.1.163
# module load intel/19.0.0.117

# Intel MPI 環境の読み込み * コンパイル時の設定環境を指定
module load intel-mpi/18.1.163
# module load intel-mpi/19.0.117


# AMBER configuration.
# Set the amber executable in $amber
# - subshells will _MOE_AMBER_EXE to override $amber
# - Default to $AMBERHOME which we assume has been set up properly.
# - CUDA_VISIBLE_DEVICES is required for GPU usage.
# - insert mpirun -np # if needed for mpi use


function run_cp_md(){
\cp -f $3 inpcrd
$mpi $2 -O -ref inpcrd -r restrt_tmp
# pmemd -O -ref inpcrd -r restrt_tmp
\cp -f restrt_tmp restrt
ambpdb -c restrt > pdb
mv -f mdin ${1}.mdin
mv -f mdout ${1}.mdout
mv -f restrt ${1}.restrt
# mv -f mdcrd ${1}.mdcrd
mv -f pdb ${1}.pdb
mv -f restrt_tmp restrt
}

function runminlast() {
echo start $1
cat << EOF > mdin
Minimization1
&cntrl
  imin=1, maxcyc=100000, ncyc=3000, drms=0.1,
  ntr=1, restraint_wt=10000, restraintmask='!@H=',

ig=-1,
vlimit=-1,
iwrap=1,

/

EOF

run_cp_md ${1%.*} $ambermin $1
}

# ntwf option was removed for all md mdin because of error (GPU)
# heat (default 100ps)

# -- main --
# -- setup --
\cp $prmtop prmtop
# gencfg
# -- minimization --

# cp Equilibration.restrt $fbase.Equilibration.restrt

trajs=`ls tgt*.rst`
for traj in $trajs
do
    runminlast ${traj}
done

# cp Minimization3.restrt $fbase.a.0.min3.rstrt

exit 0;

