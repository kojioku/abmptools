#!/bin/bash
# 2020 5/31 v2.4
# update 2023/02/18 Koji Okuwaki

# first frame: initial (0ps) structure
# e.g.) heat 100ps, equil 100ps, prod 10000ps, sampletime 1.0ps
# -> 0 means structure at 0ps(initial)
# -> 200 means structure at 200ps (end of equil)
# -> 10200 means structure at 10200ps (last)
# v2.4: add method to exclude minimization.mdcrd

# for tsubame
module load amber 2> /dev/null

# --user setting--
## caputure time info(ps)
centerinfo=":1"
## backbone atom mask
tgt=":2"
# backbone="@C,CA,N"
# backbone="@O3',C3',C4',C5',O5',P"

# e.g.)
# :6 -> residue number 6
# :6-12 -> residue number from 6 to 12
# @C,CA,N -> atom name is C,CA,or N (= protein backbone)
# @O3',C3',C4',C5',O5',P" (nucleotide backbone)

# --user setting end--

ref=`ls *.a.0.coor`
name=`ls *.a.0.coor`
head=${name%%.*}

mkdir $dir 2>/dev/null
nowtime=`date +"%Y%m%d%I%M"`

#-- print section --
# echo etime: $etime
# echo finterval: $finterval
echo reference coordfile: $ref
sleep 3

prmtop=`ls *.psf`
xsc=`ls *.a.0.xsc`
init=`ls *.a.0.coor`
traj=`ls *.coor.dcd`

# if [ $stime == 0 ]; then
    initbox=`tail -n 1 $xsc | awk '{printf ("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f%16i\n", $2, $6, $10, 90, 90, 90, 1)}'`
    sed -e "9i $initbox" $init > init.coor
    ls init.coor > inbuf.txt
    ls $traj >> inbuf.txt
# else
    # ls $traj > inbuf.txt
# fi

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
echo "rmsd rtest $tgt out ${head}.${nowtime}.agr nofit" >> cpptraj.in
# echo "reference $ref [ref_data]" >> cpptraj.in
# echo "rmsd rtest $tgt ref [ref_data] $tgt out ${head}.${nowtime}.agr nofit" >> cpptraj.in
echo "run" >> cpptraj.in
cpptraj < cpptraj.in

echo "$head.$nowtime.agr was generated."

