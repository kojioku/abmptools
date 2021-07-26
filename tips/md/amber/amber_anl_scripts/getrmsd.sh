#!/bin/bash

# 2020 5/31 v2.4
# update 2020/05/12 Koji Okuwaki

# first frame: initial (0ps) structure
# e.g.) heat 100ps, equil 100ps, prod 10000ps, sampletime 1.0ps
# -> 0 means structure at 0ps(initial)
# -> 200 means structure at 200ps (end of equil)
# -> 10200 means structure at 10200ps (last)
# v2.4: add method to exclude minimization.mdcrd

# for tsubame
module load amber 2> /dev/null

#--user setting--
# caputure time info(ps)

# backbone atone mask
backbone="@C,CA,N"
# backbone="@O3',C3',C4',C5',O5',P"
#--user setting end--

ref=`ls *.a.0.coor`
name=`ls *.a.0.coor`
head=${name%%.*}

mkdir $dir 2>/dev/null

# init=`ls *.a.0.coor`
trajs=`ls *.mdcrd`
prmtop=$head.z.prmtop


#-- print section --
# echo etime: $etime
# echo finterval: $finterval
echo reference coordfile: $ref

sleep 3

echo "parm $prmtop" > cpptraj.in
for traj in $trajs
do
    if [ "$traj" == "Minimization1.mdcrd" ] || [ "$traj" == "Minimization2.mdcrd" ] || [ "$traj" == "Minimization3.mdcrd" ]; then
        echo "$traj skip"
    else
        echo "trajin $traj"
        echo "trajin $traj" >> cpptraj.in
    fi
done
echo "parminfo" >> cpptraj.in
echo "reference $ref [ref_data]" >> cpptraj.in
echo "rmsd rtest $backbone ref [ref_data] $backbone out $head.agr" >> cpptraj.in
echo "run" >> cpptraj.in
cpptraj < cpptraj.in


echo "$head.agr was generated."

