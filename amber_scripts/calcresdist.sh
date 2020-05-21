#!/bin/bash
# 2020/5/20 v2.0
module load amber 2> /dev/null
#--user setting--
res1=1
res2=30
#--user setting end--

prmtop=`ls *.prmtop`
trajs=`ls *0.mdcrd`

echo "parm $prmtop" > cpptraj.in
for traj in $trajs
do
    echo "trajin $traj"
    echo "trajin $traj" >> cpptraj.in
done
echo "distance dist :$res1 :$res2 out distance-res$res1-res$res2.txt" >> cpptraj.in
echo 'run' >> cpptraj.in
cpptraj < cpptraj.in

echo "distance-res$res1-res$res2.txt was generated"

