#!/bin/bash
# Author: Koji Okuwaki

# OpenPBS(qsub)ディレクティブ(#PBS始まり行)
## gpu設定
## 長時間実行する場合は、select=ngpus=1:gpu_id=gpu0のようにgpuを指定してください
#PBS -l select=ncpus=10:ngpus=1 

# TSUBAME(qsub)ディレクティブ(#$始まり行)
#$ -l f_node=1                 ## * [資源タイプ名]=[使用ノード数個数]
#$ -l h_rt=24:00:00             ## * 経過時間 -> Max 24時間までです。
#$ -p -5                      ## * ジョブの優先度 *デフォルトは-5です。
# カレントディレクトリに移動
#$ -cwd

set -e

# optimize program using AMBER program after MD
# method: position constrain (iberry) and restraint (restraint)
# input: all .rst files in directory
# output: optimize structure (1.Hopt -> 2.SCopt(w/restraint) -> 3.Hopt)

## カレントディレクトリに移動(OpenPBSの場合)
if [ ""$PBS_O_WORKDIR != "" ]; then cd $PBS_O_WORKDIR; fi

# load amber
module load amber

# -- user setting --
# np=4
#mpi="mpirun -np $np"
mpi=""
ambermin="pmemd.cuda"
# ambermin="pmemd"
# prmtop=$1
prmtop=`ls *.prmtop`
backbone="@N,CA,C,O3',C3',C4',C5',O5',P"
trajs=`ls tgt*.rst`
# ---------

# prepare ps
annealps=5000

# run-step

# md param
cutval=12.0
pval=1.01
dtval=0.0005
tempval=310
restmask="@N,CA,C,O3',C3',C4',C5',O5',P"
# -------------------------

# mpi="mpirun -np $np"
mpi=""
# amber20はcuda版でもMinimize可能
ambermin="pmemd.cuda"
amber="pmemd.cuda"

rstrt=rstrt
n=1

# Moduleコマンドの初期化
. /etc/profile.d/modules.sh

# amberの読み込み
module load amber/20

function run_cp_md(){
\cp -f $3 inpcrd
$mpi $2 -O -ref inpcrd -r restrt_tmp
\cp -f restrt_tmp restrt
ambpdb -c restrt > pdb
mv -f mdin ${1}.mdin
mv -f mdout ${1}.mdout
mv -f restrt ${1}.restrt
# mv -f mdcrd ${1}.mdcrd
mv -f pdb ${1}.pdb
mv -f restrt_tmp restrt
}

function runanneal(){
infile=$1
outfile=$2
echo start anneal
cat << EOF > mdin
Heating
&cntrl
  nstlim=`echo "$annealps/$dtval" | bc`,
  ntr=1, restraint_wt=3.0, restraintmask="${restmask}",

cut=${cutval},
ntt=3, gamma_ln=1.0, temp0=${tempval},
taup=2.0, pres0=${pval},
dt=${dtval},
ntpr=1000, ntwr=1000, ntwx=1000, ntwv=1000, ntwe=1000,
ig=-1,
vlimit=-1,
iwrap=1,

nmropt=1,
/
&wt type='TEMP0', istep1=0, istep2=`echo "$annealps/$dtval" | bc`, value1=${tempval}, value2=0., /
&wt type='END', /

EOF

run_cp_md $outfile $ambermin $infile
}


function runminlast1() {
infile=$1
outfile=$2
echo "step1 Hopt"
cat << EOF > mdin
Minimization1
&cntrl
  imin=1, maxcyc=100000, ncyc=3000, drms=0.1,
  ntr=1, restraint_wt=10000, restraintmask='!@H=',
  ig=-1,
  vlimit=-1,
/
EOF

run_cp_md $outfile $ambermin $infile
}


function runminlast2() {
infile=$1
outfile=$2
echo "step2(backbone:constrain SideChain:restraint)"
cat << EOF > mdin
Minimization2
&cntrl
  imin=1, maxcyc=100000, ncyc=3000, drms=0.5,
  ntr=1, restraint_wt=10000.0, restraintmask="${backbone}",
  ig=-1,
  vlimit=-1,
/
EOF

run_cp_md $outfile $ambermin $infile
}


function runminlast3() {
infile=$1
outfile=$2
echo "step3 Hopt"
cat << EOF > mdin
Minimization3
&cntrl
  imin=1, maxcyc=100000, ncyc=3000, drms=0.1,
  ntr=1, restraint_wt=10000, restraintmask='!@H=',
  ig=-1,
  vlimit=-1,
/
EOF

run_cp_md $outfile $ambermin $infile
}

# -- main section --
\cp $prmtop prmtop
for traj in $trajs
# trajs (*.rst)
do
    echo start $traj
    runanneal ${traj} ${traj%.*}.anneal

    runminlast1 ${traj%.*}.anneal.restrt ${traj%.*}.minlast1
    # in: xxx.rst, out: xxx_minlast1.restrt
    runminlast2 ${traj%.*}.minlast1.restrt ${traj%.*}.minlast2
    # in: (xxx.rst -> xxx.minlast1.restrt), out: xxx.minlast2.restrt
    runminlast3 ${traj%.*}.minlast2.restrt ${traj%.*}.minlast3
    # in: (xxx.rst -> xxx.minlast2.restrt), out: xxx.minlast3.restrt

done

exit 0;

