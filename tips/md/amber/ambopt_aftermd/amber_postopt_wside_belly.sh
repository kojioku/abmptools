#!/bin/bash
# Author: Koji Okuwaki
set -e

# optimize program using AMBER program after MD
# method: position constrain (iberry) and restraint (restraint)

# input: all .rst files in directory
# output: optimize structure (1.Hopt -> 2.SCopt(w/restraint) -> 3.Hopt)

# -- user setting --
np=4
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
echo "step1 Hopt"
cat << EOF > mdin
Minimization1
&cntrl
  imin=1, maxcyc=100000, ncyc=3000,
  ibelly=1, !constrain atom positions
  bellymask='!@H=',
  ig=-1,
  vlimit=-1,
  iwrap=1,
/
EOF

run_cp_md ${1%.*}_minlast1 $ambermin $1
}


function runminlast2() {
echo "step2(backbone:constrain SideChain:restraint)"
cat << EOF > mdin
Minimization2
&cntrl
  imin=1, maxcyc=100000, ncyc=3000
  ibelly=1, !constrain atom positions
  bellymask="${backbone}",
  ntr=1, restraint_wt=100, restraintmask="!${backbone}",
  ig=-1,
  vlimit=-1,
  iwrap=1,
/
EOF

run_cp_md ${1%.*}_minlast2 $ambermin $1
}


function runminlast3() {
echo "step3 Hopt"
cat << EOF > mdin
Minimization3
&cntrl
  imin=1, maxcyc=100000, ncyc=3000,
  ibelly=1, !constrain atom positions
  bellymask='!@H=',
  ig=-1,
  vlimit=-1,
  iwrap=1,
/
EOF

run_cp_md ${1%.*}_minlast3 $ambermin $1
}



# -- main section --
\cp $prmtop prmtop
for traj in $trajs
do
    echo start $traj
    runminlast1 ${traj}
    runminlast2 ${traj}
    runminlast3 ${traj}
done

exit 0;

