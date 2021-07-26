#!/bin/bash

# tgtfrcmod=bnp4c2-aftermoe-single.frcmod
tgtfrcmod=myff/dotap.frcmod
tgtfrcmod2=myff/dppc.frcmod
tgtfrcmod3=myff/dypc.frcmod
tgtfrcmod4=myff/tap.frcmod
solvsize=2
# tgtlib=bnp4c2_antechamber/bnp.lib
# tgtres=BNP
# nmol=20

echo """
# gengeral
verbosity 2
source leaprc.protein.ff14SB
source leaprc.DNA.OL15
source leaprc.RNA.OL3
source leaprc.water.tip3p
source leaprc.gaff2
source leaprc.gaff
addPdbResMap { { \"LI+\" \"LI\" } { \"NA+\" \"NA\" } { \"MG2\" \"MG\" } { \"CL-\" \"CL\" } { \"K+\" \"K\" } { \"RB+\" \"RB\" } { \"CS+\" \"CS\" } }
addPdbAtomMap { { \"LI+\" \"LI\" } { \"NA+\" \"NA\" } { \"MG2\" \"MG\" } { \"CL-\" \"CL\" } { \"K+\" \"K\" } { \"RB+\" \"RB\" } { \"CS+\" \"CS\" } }
# read unparamitrized mol param
loadAmberParams ${tgtfrcmod}
loadAmberParams ${tgtfrcmod2}
loadAmberParams ${tgtfrcmod3}
loadAmberParams ${tgtfrcmod4}
# read mol coord
com = createUnit com
loadoff $2
loadoff $3
mol = loadpdb $1
com = combine { com mol }
addIons2 com Na+ 0
addIons2 com Cl- 0
savePdb com  ${1%.*}.init-dry.pdb
solvateBox com TIP3PBOX ${solvsize}
savePdb com  ${1%.*}.init.pdb
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
