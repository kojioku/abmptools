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
np=2

# -------------------------

mpi="mpirun -np $np"
ambermin="pmemd.MPI"
prmtop=$1
restmask="@N,CA,C,O3',C3',C4',C5',O5',P"


rstrt=rstrt
n=1

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
  imin=1, maxcyc=100000, ncyc=3000
  ntr=1, restraint_wt=10000, restraintmask="${restmask}",

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

