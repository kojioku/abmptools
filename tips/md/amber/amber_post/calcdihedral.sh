#!/bin/bash
# 2020/5/23 v2.4
# v2.4: add method to exclude minimization.mdcrd

module load amber 2> /dev/null
#--user setting--
atom1=1
atom2=2
atom3=3
atom4=4
#--user setting end--

prmtop=`ls *.prmtop`
trajs=`ls *.mdcrd`

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
echo "dihedral dihedral @$atom1 @$atom2 @$atom3 @$atom4 out dihedral-a$atom1-a$atom2-a$atom3-a$atom4.txt" >> cpptraj.in
echo 'run' >> cpptraj.in
cpptraj < cpptraj.in

echo "dihedral-a$atom1-a$atom2-a$atom3-a$atom4.txt was generated"

