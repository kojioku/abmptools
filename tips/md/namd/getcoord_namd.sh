#!/bin/bash

#--user setting--
stime=200
etime=1200
interval=100
sampletime=1.0
#--user setting end--

sframe=`echo "$stime / $sampletime" | bc`
eframe=`echo "$etime / $sampletime" | bc`
ivframe=`echo "$interval / $sampletime" | bc`
name=`ls *.a.0.coor`
head=${name%%.*}


mkdir mdout >/dev/null

prmtop=`ls *.prmtop`
ls *.coor.dcd > inbuf.txt

echo "parm $prmtop" > cpptraj.in
while read line
do
    echo $line 
    echo "loadcrd $line name coord" >> cpptraj.in
done < inbuf.txt

echo "crdout coord mdout/${head}.pdb multi crdframes $sframe,$eframe,$ivframe" >> cpptraj.in
echo 'run' >> cpptraj.in
cpptraj < cpptraj.in

echo "mdout/mdout.pdb ($stime,$etime) was generated."

