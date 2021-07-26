#!/bin/bash

# v2.2
# update 2021/03/16 Koji Okuwaki

# first frame: initial (0ps) structure
# e.g.) heat 100ps, equil 100ps, prod 10000ps, sampletime 1.0ps
# -> 0 means structure at 0ps(initial)
# -> 200 means structure at 200ps (end of equil)
# -> 10200 means structure at 10200ps (last)

#--user setting--
# caputure time info(ps)
stime=200
etime=10200
getinterval=100
sampletime=1.0
parm="210306dmebcd.a.0.z.prmtop"
# file interval

maskflag=true
centerinfo=":1"
maskinfo=":1-2"
stripdist=8.0

#--user setting end--

sframe=`echo "$stime / $sampletime" | bc`
eframe=`echo "$etime / $sampletime" | bc`
ivframe=`echo "$getinterval / $sampletime" | bc`

name=`ls *.a.0.coor`
head=${name%%.*}

dirbuf='mdout'
dir='mdout_centerd_mask'
mkdir $dirbuf $dir 2>/dev/null

# psf=`ls *.a.0.psf`
xsc=`ls *.a.0.xsc`
init=`ls *.a.0.coor`
traj=`ls *coor.dcd`

newtraj='tgt.trj'

#-- print section --
echo stime: $stime, etime: $etime, getinterval: $getinterval, sampletime: $sampletime
echo centerinfo, $centerinfo
echo stripdist, $stripdist

echo sframe, $sframe, eframe, $eframe, ivframe, $ivframe

sleep 3

echo "parm $parm" > cpptraj.in

if [ $stime = 0 ]; then
    echo 'read first step'
    initbox=`tail -n 1 $xsc | awk '{printf ("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f%16i\n", $2, $6, $10, 90, 90, 90, 1)}'`
    echo "$initbox"
    sed -e "9i $initbox" $init > init.coor
    echo init.coor > inbuf.txt
    ls $traj >> inbuf.txt
    sframe=`echo "$sframe + 1" | bc`
    eframe=`echo "$eframe + 1" | bc`
else
    echo 'read start'
    ls $traj > inbuf.txt
fi

while read line
do
    echo $line
    echo "loadcrd $line name coord" >> cpptraj.in
done < inbuf.txt

echo 'start trajout'
#trajout section
echo "parminfo $parm" >> cpptraj.in
# echo "autoimage anchor $centerinfo origin" >> cpptraj.in
echo "crdout coord mdout/${head}.pdb multi crdframes $sframe,$eframe,$ivframe" >> cpptraj.in
echo 'run' >> cpptraj.in
cpptraj < cpptraj.in

#mask section
if "$maskflag"; then
    echo 'start cut'
    echo "parm $parm" > cpptraj_mask.in
    num=0
    for i in `seq $sframe $ivframe $eframe`
    do
        echo "trajin mdout/${head}.pdb.$i" >> cpptraj_mask.in
        num=$((num+1))
    done
    echo "autoimage anchor $centerinfo origin" >> cpptraj_mask.in
    echo "reference mdout/${head}.pdb.$sframe" >> cpptraj_mask.in
    echo "mask \"$maskinfo<:$stripdist\" maskpdb $dir/$head.pdb" >> cpptraj_mask.in
    echo 'run' >> cpptraj_mask.in
    cpptraj < cpptraj_mask.in
fi

echo "$dir/mdout.pdb ($sframe,$eframe) was generated."

num=1
for i in `seq $sframe $ivframe $eframe`
do
    mv $dir/$head.pdb.$num $dir/${head}.pdb.${i}
    num=$((num+1))
done
