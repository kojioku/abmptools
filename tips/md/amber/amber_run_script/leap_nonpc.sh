#!/bin/bash

echo """
source leaprc.protein.ff14SB
source leaprc.DNA.OL15
source leaprc.RNA.OL3
source leaprc.water.tip3p
source leaprc.gaff2
addPdbResMap { { \"LI+\" \"LI\" } { \"NA+\" \"NA\" } { \"MG2\" \"MG\" } { \"CL-\" \"CL\" } { \"K+\" \"K\" } { \"RB+\" \"RB\" } { \"CS+\" \"CS\" } }
addPdbAtomMap { { \"LI+\" \"LI\" } { \"NA+\" \"NA\" } { \"MG2\" \"MG\" } { \"CL-\" \"CL\" } { \"K+\" \"K\" } { \"RB+\" \"RB\" } { \"CS+\" \"CS\" } }
com = createUnit com
mol = loadpdb $1
com = combine { com mol }
savePdb com ${1%.*}.init.pdb
saveAmberParm com ${1%.*}.z.prmtop ${1%.*}.a.0.coor
quit
""" > leaprc
tleap -f leaprc

# source leaprc.DNA.OL15
# source leaprc.RNA.OL3
# source leaprc.water.tip3p
# source leaprc.gaff2

