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
# centerinfo=":596"
headbuf=${traj%.*}
head=${headbuf##*/}

startframe=1001
endframe=2001 #2001
interval=20

# stimeps=10000  # please specify the start time of production run
# intervalps=200

newtraj=${head}_int${interval}.xtc

mkdir $dir 2> /dev/null
echo "parm $prmtop" > cpptraj.in
# echo "parminfo $prmtop" >> cpptraj.in
echo "trajin $traj" >> cpptraj.in
# echo "autoimage anchor $centerinfo origin" >> cpptraj.in
# echo "autoimage" >> cpptraj.in
# echo "center $centerinfo mass" >> cpptraj.in
echo "trajout ${head}-first.pdb onlyframes 1" >> cpptraj.in
echo "trajout $newtraj start $startframe stop $endframe offset $interval" >> cpptraj.in
echo 'run' >> cpptraj.in
cpptraj < cpptraj.in

today=$(date "+%Y%m%d")
ofile="trajcut.${today}.log"

echo `date` >> $ofile
echo $prmtop >> $ofile
echo $traj >> $ofile
echo $startframe, $endframe, $interval >> $ofile
