#!/bin/bash

echo """
verbosity 2
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
addIons2 com Na+ 0
addIons2 com Cl- 0
savePdb com  ${1%.*}.init.pdb
solvateBox com TIP3PBOX 10
saveAmberParm com ${1%.*}.z.prmtop ${1%.*}.a.0.coor
quit
""" > leaprc
tleap -f leaprc


# --- riken  ---
# verbosity 2
# source leaprc.protein.ff14SB
# source leaprc.DNA.OL15
# source leaprc.RNA.OL3
# source leaprc.water.tip3p
# source leaprc.gaff2
# addPdbResMap { { "LI+" "LI" } { "NA+" "NA" } { "MG2" "MG" } { "CL-" "CL" } { "K+" "K" } { "RB+" "RB" } { "CS+" "CS" } }
# addPdbAtomMap { { "LI+" "LI" } { "NA+" "NA" } { "MG2" "MG" } { "CL-" "CL" } { "K+" "K" } { "RB+" "RB" } { "CS+" "CS" } }
# com = createUnit com
# rec_mol1 = loadPdb ./ffknown/rec_mol1.pdb
# com = combine { com rec_mol1 }
# loadAmberParams ./ffunknown/frcmod/lig_mol2.frcmod
# lig_mol2 = loadMol2 ./ffunknown/frcmod/lig_mol2.mol2
# com = combine { com lig_mol2 }
# mol3 = loadPdb ./ffknown/mol3.pdb
# com = combine { com mol3 }
# addIons2 com Na+ 0
# addIons2 com Cl- 0
# savePdb com amber.pdb
# solvateBox com TIP3PBOX 10
# saveAmberParm com prmtop tleap.restrt
# quit
