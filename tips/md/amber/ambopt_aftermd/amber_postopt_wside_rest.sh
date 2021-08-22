#!/bin/bash
# Author: Koji Okuwaki
set -e

# optimize program using AMBER program after MD
# method: position constrain (iberry) and restraint (restraint)
# input: all .rst files in directory
# output: optimize structure (1.Hopt -> 2.SCopt(w/restraint) -> 3.Hopt)

# -- user setting --
np=8
mpi="mpirun -np $np"
ambermin="pmemd.MPI"
prmtop=$1
backbone="@N,CA,C,O3',C3',C4',C5',O5',P"
trajs=`ls tgt*.rst`
# ---------

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
    runminlast1 ${traj} ${traj%.*}.minlast1
    # in: xxx.rst, out: xxx_minlast1.restrt
    runminlast2 ${traj%.*}.minlast1.restrt ${traj%.*}.minlast2
    # in: (xxx.rst -> xxx.minlast1.restrt), out: xxx.minlast2.restrt
    runminlast3 ${traj%.*}.minlast2.restrt ${traj%.*}.minlast3
    # in: (xxx.rst -> xxx.minlast2.restrt), out: xxx.minlast3.restrt

done

exit 0;

