#!/bin/bash

# 2020/6/12 v3.3
# first frame: initial (0ps) structure
# e.g.) heat 100ps, equil 100ps, prod 10000ps, sampletime 1.0ps
# -> 0 means structure at 0ps(initial)
# -> 200 means structure at 200ps (end of equil)
# -> 10200 means structure at 10200ps (last)
module load amber 2> /dev/null

# Warning: cpptraj trajout(startframe) have a bug. please check a MODEL number in tgt.pdb
# v3.3: add error handling for cpptraj bug
# v3.2: add method to exclude minimization.mdcrd
# v3.0: add tempolary solution for cpp traj trajout(start) bug

#--user setting--
# caputure time info(ps)

prmtop=$1
traj=$2
centerinfo=":1-381"
dir='gmxpdbs-foropt'
headbuf=${traj%.*}
head=${headbuf##*/}

startframe=1001
endframe=2001
interval=20

stimeps=10000  # please specify the start time of production run
intervalps=200

mkdir $dir 2> /dev/null
echo "parm $prmtop" > cpptraj.in
echo "parminfo $prmtop" >> cpptraj.in
echo "trajin $traj" >> cpptraj.in
echo "autoimage anchor $centerinfo origin" >> cpptraj.in

for i in `seq $startframe $interval $endframe`
do
    newtraj=${head}_${stimeps}ps.pdb
    echo "trajout $newtraj pdb onlyframes $i restart" >> cpptraj.in
    stimeps=`echo "$stimeps + $intervalps" | bc`
done
echo 'run' >> cpptraj.in
cpptraj < cpptraj.in

mv ${head}*.pdb $dir
cp $prmtop $dir
echo "$dir/*.pdb was generated."

today=$(date "+%Y%m%d")
ofile="getcutpdb.${today}.log"

echo `date` >> $ofile
echo $prmtop >> $ofile
echo $traj >> $ofile
echo $startframe, $endframe, $interval >> $ofile
