#!/bin/bash

prmtop=$1
newtraj=$2

maskinfo=":1-921"
stripdist=4.0

dir='gmxpdbs-foropt'
head='test'

mkdir $dir > /dev/null
echo "parm $prmtop" > cpptraj_mask.in
echo "trajin $newtraj" >> cpptraj_mask.in
echo "reference $newtraj" >> cpptraj_mask.in
echo "mask \"($maskinfo<:$stripdist)|:NA|:Na+\" maskpdb $dir/$head.pdb" >> cpptraj_mask.in
echo 'run' >> cpptraj_mask.in
cpptraj.OMP < cpptraj_mask.in



