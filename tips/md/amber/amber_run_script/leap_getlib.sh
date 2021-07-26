
resname='P4C'

echo """
source leaprc.protein.ff14SB
source leaprc.gaff
$resname = loadmol2 $1
check $resname
loadamberparams $2
saveoff $resname $resname.lib
saveamberparm $resname $resname.prmtop $resname.rst7
""" > leaprc
tleap -f leaprc

