#!/bin/bash

module load amber 2> /dev/null

if [ $# = 0 ]; then
    echo 'error! please set 1 argument (.pdb)'
    exit
fi

inhead=${1%.*}
nohydpdb=${inhead}-nohyd.pdb 
redpdb=${inhead}-p4a.pdb

nohyddir=$PWD/nohyd
p4adir=$PWD/p4a
leapdir=$PWD/leap
mdrundir=$PWD/mdrun

leapscript="$PWD/leap_pc.sh"

mkdir $nohyddir
mkdir $p4adir
mkdir $leapdir
mkdir $mdrundir

cp $1 $nohyddir

# delete h atom
cd $nohyddir
pdb4amber -i $1  -o $nohydpdb --nohyd
cp $nohydpdb $p4adir

# add h atom (amber reduce)
cd $p4adir
pdb4amber -i $nohydpdb -o $redpdb --reduce
cp $redpdb $leapdir

# leap
cd $leapdir
bash $leapscript $redpdb
cp *coor *prmtop $mdrundir


