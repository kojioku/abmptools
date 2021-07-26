#!/bin/bash

# 2020 5/31 v2.4
# update 2020/05/12 Koji Okuwaki

# first frame: initial (0ps) structure
# e.g.) heat 100ps, equil 100ps, prod 10000ps, sampletime 1.0ps
# -> 0 means structure at 0ps(initial)
# -> 200 means structure at 200ps (end of equil)
# -> 10200 means structure at 10200ps (last)
# v2.4: add method to exclude minimization.mdcrd

#--user setting--
# caputure time info(ps)

# backbone atone mask
backbone="@CA" #"@C,CA,N"
# backbone="@O3',C3',C4',C5',O5',P"
#--user setting end--

prmtop=$1
tgttraj=$2
ref=$3
out=$4

mkdir $dir 2>/dev/null

# init=`ls *.a.0.coor`

#-- print section --
echo prmtop: $prmtop
echo tgttraj: $tgttraj
echo reference coordfile: $ref
echo out: $out


sleep 3

echo "parm $prmtop" > cpptraj.in
echo "trajin $tgttraj" >> cpptraj.in
echo "parminfo" >> cpptraj.in
echo "reference $ref [ref_data]" >> cpptraj.in
echo "rmsd rtest $backbone ref [ref_data] $backbone out $out.agr" >> cpptraj.in
echo "run" >> cpptraj.in
cpptraj.OMP < cpptraj.in

echo "$out.agr was generated."

