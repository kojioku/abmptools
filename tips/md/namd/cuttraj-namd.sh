#!/bin/bash

# first frame: initial (0ps) structure
# e.g.) heat 100ps, equil 100ps, prod 10000ps, sampletime 1.0ps
# -> 0 means structure at 0ps(initial)
# -> 200 means structure at 200ps (end of equil)
# -> 10200 means structure at 10200ps (last)

# Warning: cpptraj trajout(startframe) have a bug. please check a MODEL number in tgt.pdb

# --- user setting ---
## caputure time info(ps)
stime=0
etime=100200
interval=10
sampletime=1.0  # ps
centerinfo=":1"
# --- user setting end ---

# for tsubame
module load amber 2> /dev/null

# tempolary solution for cpp traj trajout(start) bug
stime_inprod=0
etime_inprod=`echo "$etime - $stime" | bc`
sframe_inprod=0
eframe_inprod=`echo "$etime_inprod / $sampletime" | bc`
# -----
sframe=`echo "$stime / $sampletime" | bc`
eframe=`echo "$etime / $sampletime" | bc`
ivframe=`echo "$interval / $sampletime" | bc`

name=`ls *.a.0.coor`
head=${name%%.*}

dirhead='cuttraj_centerd'
nowtime=`date +"%Y%m%d%I%M"`
dir=${dirhead}.${nowtime}
mkdir ${dir} 2>/dev/null

prmtop=`ls *.psf`
xsc=`ls *.a.0.xsc`
init=`ls *.a.0.coor`
trajs=`ls *.coor.dcd`

# prodinit=`ls ${head}*.${stime}.coor`
# if [ "$prodinit" == "" ]; then
#     prodinit="xxxxx"
#     echo "Error!!: check start time setting"
#     exit 0
# fi

# --- input filter ---
trajs_filt=''
for traj in $trajs
do
    buf=${traj%.*}
    buf2=${buf%.*}
    number=${buf2##*.}
    # echo $number, $buf2
    if [ $number -le $stime ]; then
        echo $traj skip
    else
        trajs_filt="$trajs_filt $traj"
    fi

    if [ $number -ge $etime ]; then
        echo $traj: last
        break
    fi
done
echo "input trajs: $trajs_filt"

# --- parm print section ---
echo stime: $stime, etime: $etime, interval: $interval, sampletime: $sampletime
echo centerinfo, $centerinfo

echo sframe, $sframe, eframe, $eframe, ivframe, $ivframe
sleep 3

# --- input section ---
if [ $stime_inprod -eq 0 ]; then
    initbox=`tail -n 1 $xsc | awk '{printf ("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f%16i\n", $2, $6, $10, 90, 90, 90, 1)}'`
    sed -e "9i $initbox" $init > init.coor
    ls init.coor > inbuf.txt
    ls $trajs_filt >> inbuf.txt
    sframe_inprod=`echo "$sframe_inprod + 1" | bc`
    eframe_inprod=`echo "$eframe_inprod + 1" | bc`
else
    ls $trajs_filt > inbuf.txt
fi

newtraj="${head}.dcd"
# [section1] trajectory cut and center solute
echo "parm $prmtop" > cpptraj.in
while read line
do
    # echo $line
    echo "trajin $line" >> cpptraj.in
done < inbuf.txt

echo "parminfo $parm" >> cpptraj.in
# echo "autoimage anchor $centerinfo origin" >> cpptraj.in
echo "center $centerinfo" >> cpptraj.in
echo "image" >> cpptraj.in
echo "trajout $dir/$newtraj dcd start $sframe stop $eframe offset $ivframe" >> cpptraj.in
#echo "trajout $dir/$newtraj pdb multi start $sframe stop $eframe offset $ivframe" >> cpptraj.in

echo 'run' >> cpptraj.in
cpptraj < cpptraj.in

cp $prmtop $dir
echo "${dir}/${head}.dcd ($stime,$etime) was generated."

# [section2] mask droplet
# if "$maskflag"; then
#     echo "parm $prmtop" > cpptraj_mask.in
#     echo "trajin $newtraj" >> cpptraj_mask.in
#     echo "reference $newtraj" >> cpptraj_mask.in
#     echo "mask \"($maskinfo<:$stripdist)|:NA|:Na+|:CL|:Cl-\" maskpdb $dir/$head.pdb" >> cpptraj_mask.in
#     echo 'run' >> cpptraj_mask.in
#     cpptraj < cpptraj_mask.in
# fi
# 
# num=1
# for i in `seq $stime $interval $etime`
# do
#     echo $num to $i
#     # future work: add zero pading
#     mv $dir/$head.pdb.$num $dir/${head}-${i}ps.pdb
#     num=$((num+1))
# done
# 
# echo "${dir}/mdout.pdb ($stime,$etime) was generated."

