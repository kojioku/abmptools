#!/bin/bash

# usage
# cpptraj -p ../fujisaki/sbecd7mzp_priz_tryok.z.prmtop -y sbecd7+mzp_priz_50nsdynamics_namd.c.100.coor.dcd -x output.mdcrd

# parm <topology> 
# trajin <input dcd> 
# trajout output.mdcrd 
function mmm(){
cpptraj -p $prmtop -y $1 -x ambtraj/${1%.coor*}.mdcrd
}

mkdir ambtraj 2>/dev/null

prmtop=`ls *.prmtop`
ls *.coor.dcd > n2abuf.txt

while read line
do
    mmm $line 
done < n2abuf.txt
