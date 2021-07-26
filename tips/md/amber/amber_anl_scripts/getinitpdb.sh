#!/bin/bash
module load amber 2> /dev/null

prmtop=`ls *.prmtop`
coor=`ls *.0.coor`
outname='init.pdb'

ambpdb -p $prmtop -c $coor > $outname

echo "init.pdb was generated."
