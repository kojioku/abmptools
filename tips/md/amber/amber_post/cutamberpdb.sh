#!/bin/bash

# 2020/6/12 v1.0
# first frame: initial (0ps) structure
# e.g.) heat 100ps, equil 100ps, prod 10000ps, sampletime 1.0ps
# -> 0 means structure at 0ps(initial)
# -> 200 means structure at 200ps (end of equil)
# -> 10200 means structure at 10200ps (last)
module load amber 2> /dev/null

# Warning: cpptraj trajout(startframe) have a bug. please check a MODEL number in tgt.pdb

#--user setting--
# caputure time info(ps)
tgtpdb=$1  # please specify the start time of production run

centerinfo=":57"    # ":1-276"
maskinfo=":1-213"
stripdist=4.0      #3.0
maskflag=true
#--user setting end--

tgtframe=1

name=`ls *.a.0.coor`
head=${name%%.*}

dir='pdbcut_centerd_mask'
mkdir $dir 2>/dev/null

prmtop=`ls *.prmtop`
init=`ls *.a.0.coor`

trajs=`ls *.mdcrd`
# echo $trajs

# input filter
trajs_filt=''
newtraj='tgt.pdb'

#-- print section --
echo "## input info"
echo - tgtname, $tgtpdb
echo - centerinfo, $centerinfo
echo - stripdist, $stripdist
# echo - tgtframe, $tgtframe

sleep 3

ls $trajs_filt > inbuf.txt
echo "parm $prmtop" > cpptraj.in

echo "trajin $tgtpdb" >> cpptraj.in

echo "parminfo $prmtop" >> cpptraj.in
echo "autoimage anchor $centerinfo origin" >> cpptraj.in
echo "trajout $newtraj onlyframes $tgtframe" >> cpptraj.in
echo 'run' >> cpptraj.in
cpptraj < cpptraj.in

if "$maskflag"; then
    echo "parm $prmtop" > cpptraj_mask.in
    echo "trajin $newtraj" >> cpptraj_mask.in
    echo "reference $newtraj" >> cpptraj_mask.in
    echo "mask \"($maskinfo<:$stripdist)|:NA|:Na+|:CL|:Cl-\" maskpdb $dir/$head.pdb" >> cpptraj_mask.in
    echo 'run' >> cpptraj_mask.in
    cpptraj < cpptraj_mask.in
fi

# echo 1 to $tgtframe
# future work: add zero pading

oname=${head}_cent${centerinfo:1}_strip${stripdist}.pdb

mv $dir/${head}.pdb.1 $dir/$oname

echo "$dir/$oname was generated."

