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
stimeps=100000  # please specify the start time of production run
etimeps=200000
intervalps=1000

centerinfo=":1-100"
maskinfo=":1-3426"
stripdist=4.0
maskflag=true
#--user setting end--

dir='mdopted_centerd_mask'
mkdir $dir 2>/dev/null

prmtop=$1
head=$2

# input filter
trajs_filt=''
for num in `seq $stimeps $intervalps $etimeps`
do
    echo $num
    traj=`ls *_${num}ps.restrt`
    trajs_filt="$trajs_filt $traj"
done
echo "input trajs: $trajs_filt"

#-- print section --
echo stimeps: $stimeps, etimeps: $etimeps, interval: $intervalps
echo centerinfo, $centerinfo
echo stripdist, $stripdist
sleep 3

# input section
echo "parm $prmtop" > cpptraj.in
newtraj="${2}.trj"

# write section
for line in $trajs_filt
do
    echo $line
    echo "trajin $line" >> cpptraj.in
done

echo "parminfo $prmtop" >> cpptraj.in
# echo "autoimage anchor $centerinfo origin" >> cpptraj.in
echo "center $centerinfo" >> cpptraj.in
echo "image" >> cpptraj.in
echo "trajout $newtraj" >> cpptraj.in
echo 'run' >> cpptraj.in
cpptraj < cpptraj.in

# mask section
if "$maskflag"; then
    echo "parm $prmtop" > cpptraj_mask.in
    echo "trajin $newtraj" >> cpptraj_mask.in
    echo "reference $newtraj" >> cpptraj_mask.in
    echo "mask \"($maskinfo<:$stripdist)|:NA|:Na+|:CL|:Cl-\" maskpdb $dir/$head.pdb" >> cpptraj_mask.in
    echo 'run' >> cpptraj_mask.in
    cpptraj < cpptraj_mask.in
fi

num=1
for i in `seq $stimeps $intervalps $etimeps`
do
    echo $num to $i
    # future work: add zero pading
    mv $dir/$head.pdb.$num $dir/${head}-${i}ps.pdb
    num=$((num+1))
done

echo "$dir/mdout.pdb ($sframe,$eframe) was generated."
