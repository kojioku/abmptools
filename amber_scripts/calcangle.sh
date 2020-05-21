#!/bin/bash
# 2020/5/20 v2.0
module load amber 2> /dev/null
#--user setting--
atom1=1
atom2=2
atom3=3
#--user setting end--

prmtop=`ls *.prmtop`
trajs=`ls *0.mdcrd`

echo "parm $prmtop" > cpptraj.in
for traj in $trajs
do
    echo "trajin $traj"
    echo "trajin $traj" >> cpptraj.in
done
echo "angle angle @$atom1 @$atom2 @$atom3 out angle-a$atom1-a$atom2-a$atom3.txt" >> cpptraj.in
echo 'run' >> cpptraj.in
cpptraj < cpptraj.in

echo "angle-a$atom1-a$atom2-a$atom3.txt was generated"

