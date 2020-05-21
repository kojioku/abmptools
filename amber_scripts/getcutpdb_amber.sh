#!/bin/bash

# first frame: initial (0ps) structure
# e.g.) heat 100ps, equil 100ps, prod 10000ps, sampletime 1.0ps
# -> 0 means structure at 0ps(initial)
# -> 200 means structure at 200ps (end of equil)
# -> 10200 means structure at 10200ps (last)
module load amber 2> /dev/null

# Warning: cpptraj trajout(startframe) have a bug. please check a MODEL number in tgt.pdb
# v3.0: add tempolary solution for cpp traj trajout(start) bug

#--user setting--
# caputure time info(ps)
stime=1100  # please specify the start time of production run
etime=101100
interval=1000
sampletime=1.0

centerinfo=":1-306"
maskinfo=":1-312"
stripdist=4.0
maskflag=true
#--user setting end--


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

dir='mdout_centerd_mask'
mkdir $dir 2>/dev/null

prmtop=`ls *.prmtop`
init=`ls *.a.0.coor`
prodinit=`ls ${head}*.${stime}.rstrt`

trajs=`ls *0.mdcrd`

# input filter
trajs_filt=''
for traj in $trajs
do
    buf=${traj%.*}
    number=${buf##*.}
    if [ $number -le $sframe ]; then
        echo $traj skip
    else
        trajs_filt="$trajs_filt $traj"
    fi

    if [ $number -ge $eframe ]; then
        echo $traj: last 
        break
    fi

done
echo "input trajs: $trajs_filt"

newtraj='tgt.pdb'

#-- print section --
echo stime: $stime, etime: $etime, getinterval: $getinterval, sampletime: $sampletime
echo centerinfo, $centerinfo
echo stripdist, $stripdist

echo sframe, $sframe, eframe, $eframe, ivframe, $ivframe

sleep 3

if [ $stime_inprod -eq 0 ]; then
    # initbox=`tail -n 1 $xsc | awk '{printf ("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f%16i\n", $2, $6, $10, 90, 90, 90, 1)}'`
    # sed -e "9i $initbox" $init > init.coor
    echo 'add init file'
    ls $prodinit > inbuf.txt
    ls $trajs_filt >> inbuf.txt
    sframe_inprod=`echo "$sframe_inprod + 1" | bc`
    eframe_inprod=`echo "$eframe_inprod + 1" | bc`
else
    ls $trajs_filt > inbuf.txt
fi

echo "parm $prmtop" > cpptraj.in

while read line
do
    echo $line
    echo "trajin $line" >> cpptraj.in
done < inbuf.txt

echo "parminfo $prmtop" >> cpptraj.in
echo "autoimage anchor $centerinfo origin" >> cpptraj.in
echo "trajout $newtraj start $sframe_inprod stop $eframe_inprod offset $ivframe" >> cpptraj.in
echo 'run' >> cpptraj.in
cpptraj < cpptraj.in

if "$maskflag"; then
    echo "parm $prmtop" > cpptraj.in
    echo "trajin $newtraj" >> cpptraj.in
    echo "reference $newtraj" >> cpptraj.in
    echo "mask \"($maskinfo<:$stripdist)|:NA\" maskpdb $dir/$head.pdb" >> cpptraj.in
    echo 'run' >> cpptraj.in
    cpptraj < cpptraj.in
fi

num=1
for i in `seq $sframe $ivframe $eframe`
do
    echo $num to $i
    # future work: add zero pading
    mv $dir/$head.pdb.$num $dir/${head}-${i}ps.pdb
    num=$((num+1))
done

echo "$dir/mdout.pdb ($sframe,$eframe) was generated."

