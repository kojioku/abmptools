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
stime=1100  # please specify the start time of production run
etime=51100
interval=10000

centerinfo=":1-306"
maskinfo=":1-312"
stripdist=4.0
maskflag=true
#--user setting end--

dir='mdopted_centerd_mask'
mkdir $dir 2>/dev/null

prmtop=`ls *.prmtop`
head=${prmtop%%.*}

# trajs=`ls *.mdcrd`
# echo $trajs

# input filter
trajs_filt=''
for num in `seq $stime $interval $etime`
do
    echo $num
    traj=tgt${num}.restrt
    trajs_filt="$trajs_filt $traj"
done
echo "input trajs: $trajs_filt"


#-- print section --
echo stime: $stime, etime: $etime, getinterval: $interval
echo centerinfo, $centerinfo
echo stripdist, $stripdist

sleep 3

# input section
echo "parm $prmtop" > cpptraj.in
newtraj='opted.trj'

# write section
for line in $trajs_filt
do
    echo $line
    echo "trajin $line" >> cpptraj.in
done

echo "parminfo $prmtop" >> cpptraj.in
echo "autoimage anchor $centerinfo origin" >> cpptraj.in
# echo "trajout $newtraj start $sframe_inprod stop $eframe_inprod offset $ivframe" >> cpptraj.in
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
for i in `seq $stime $interval $etime`
do
    echo $num to $i
    # future work: add zero pading
    mv $dir/$head.pdb.$num $dir/${head}-${i}ps.pdb
    num=$((num+1))
done

echo "$dir/mdout.pdb ($sframe,$eframe) was generated."
