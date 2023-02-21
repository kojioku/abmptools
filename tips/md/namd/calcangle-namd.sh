#!/bin/bash
# 2020/5/23 v2.4
# update 2023/02/18 Koji Okuwaki

# for tsubame
module load amber 2> /dev/null

## -- backbone atom mask
# ":" means Residue, "@" means Atom 
# See CPPTRAJ manual page "https://amber-md.github.io/cpptraj/CPPTRAJ.xhtml#magicparlabel-223" for details.
# e.g.
# :6 -> residue number 6
# :6-12 -> residue number from 6 to 12
# @1-20 -> atom number from 1 to 20
# @C,CA,N -> atomname is C,CA,or N (= protein backbone)
# @O3',C3',C4',C5',O5',P" -> (nucleotide backbone)

#--user setting--
tgt1='@1'
tgt2='@2'
tgt3='@3'
#--user setting end--

ref=`ls *.a.0.coor`
name=`ls *.a.0.coor`
head=${name%%.*}

mkdir $dir 2>/dev/null
nowtime=`date +"%Y%m%d%I%M"`

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
echo "reference $ref [ref_data]" >> cpptraj.in
# echo "autoimage anchor $centerinfo origin" >> cpptraj.in
echo "center $centerinfo" >> cpptraj.in
echo "image" >> cpptraj.in

# echo "rmsd rtest $backbone ref [ref_data] $backbone out ${head}.${nowtime}.agr nofit" >> cpptraj.in

echo "angle angle $tgt1 $tgt2 $tgt3 out $head.angle.$nowtime.log" >> cpptraj.in
echo 'run' >> cpptraj.in
cpptraj < cpptraj.in

echo "$head.angle.$nowtime.log was generated"


