#!/bin/bash
# 2020/5/23 v2.4
# v2.4: add method to exclude minimization.mdcrd

module load amber 2> /dev/null
#--user setting--
res1=1
res2=30
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
echo "distance dist :$res1 :$res2 out distance-res$res1-res$res2.txt" >> cpptraj.in
echo 'run' >> cpptraj.in
cpptraj < cpptraj.in

echo "distance-res$res1-res$res2.txt was generated"

