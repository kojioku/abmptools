#!/bin/bash

# first frame: initial (0ps) structure
# e.g.) heat 100ps, equil 100ps, prod 10000ps, sampletime 1.0ps
# -> 0 means structure at 0ps(initial)
# -> 200 means structure at 200ps (end of equil)
# -> 10200 means structure at 10200ps (last)

#--user setting--
# caputure time info(ps)
stime=0
etime=1000
interval=100
sampletime=1.0

centerinfo=":1-10"
#--user setting end--

sframe=`echo "$stime / $sampletime + 1" | bc`
eframe=`echo "$etime / $sampletime + 1" | bc`
ivframe=`echo "$interval / $sampletime" | bc`

name=`ls *.a.0.coor`
head=${name%%.*}

dir='mdout_centerd'
mkdir $dir 2>/dev/null

prmtop=`ls *.psf`
xsc=`ls *.a.0.xsc`
init=`ls *.a.0.coor`
traj=`ls *.coor.dcd`

if [ $stime == 0 ]; then
    initbox=`tail -n 1 $xsc | awk '{printf ("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f%16i\n", $2, $6, $10, 90, 90, 90, 1)}'`
    sed -e "9i $initbox" $init > init.coor
    ls init.coor > inbuf.txt
    ls $traj >> inbuf.txt
else
    ls $traj > inbuf.txt
fi

echo "parm $prmtop" > cpptraj.in
while read line
do
    # echo $line
    echo "trajin $line" >> cpptraj.in
done < inbuf.txt

echo "parminfo $parm" >> cpptraj.in
echo "autoimage anchor $centerinfo origin" >> cpptraj.in
echo "trajout $dir/${head}.pdb multi start $sframe stop $eframe offset $ivframe" >> cpptraj.in

echo 'run' >> cpptraj.in
cpptraj < cpptraj.in

echo "$dir/mdout.pdb ($stime,$etime) was generated."

