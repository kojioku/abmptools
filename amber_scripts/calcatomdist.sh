#!/bin/bash
module load amber 2> /dev/null

#--user setting--
atom1=1
atom2=30
#--user setting end--

prmtop=`ls *.prmtop`
trajs=`ls *0.mdcrd`

echo "parm $prmtop" > cpptraj.in
for traj in $trajs
do
    echo "trajin $traj"
    echo "trajin $traj" >> cpptraj.in
done
echo "distance dist @$atom1 @$atom2 out distance-a$atom1-a$atom2.txt" >> cpptraj.in
echo 'run' >> cpptraj.in
cpptraj < cpptraj.in

echo "distance-a$atom1-a$atom2.txt was generated"

