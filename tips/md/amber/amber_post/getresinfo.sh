#!/bin/bash

prmtop=$1
dstr=`date "+%Y%m%d"`

echo "parm $prmtop" > rescpp.in
echo "resinfo" >> rescpp.in
cpptraj.OMP < rescpp.in > resinfo.log

echo "resinfo.${dstr}.log"
