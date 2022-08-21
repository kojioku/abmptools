#!/bin/bash
# prepare
# xxx.top, xxx.itp: gromacs topology file
# xxx.prmtop: amber topology file (use parmed.sh)
# xxx.xtc: gromacs trajectory file
# cuttrajgmx-foropt.sh: file for generating cut-pdb

#step0 create cut-pdb
bash cuttrjgmx-foropt.sh 4j0p-1H8_e10ps.prmtop 4j0p-1H8_md-samp_e10ps.xtc
cd gmxpdbs-foropt

#step1 strip prmtop
echo '
parm 4j0p-1H8_e10ps.prmtop
parminfo *
reference 4j0p-1H8_md-samp_e10ps_10000ps.pdb
parmstrip (:1-381>:6.0)&!:NA+&!:NA&!:Na+&!:CL-&!:CL&!:Cl-
parmwrite out 4j0p-1H8_e10ps-10000ps-strip.prmtop
' > cpptraj_strip.in
cpptraj < cpptraj_strip.in

#step2 mask str
echo '
parm 4j0p-1H8_e10ps.prmtop
trajin 4j0p-1H8_md-samp_e10ps_10000ps.pdb
reference 4j0p-1H8_md-samp_e10ps_10000ps.pdb
mask (:1-381<:6.0)|:NA+|:NA|:Na+|:CL|:Cl- maskpdb 4j0p-10000ps-strip.pdb
run
' > cpptraj_mask.in
cpptraj < cpptraj_mask.in

#step3 pdb2trj
echo '
parm 4j0p-1H8_e10ps-10000ps-strip.prmtop
trajin 4j0p-10000ps-strip.pdb.1
trajout 4j0p-10000ps-strip.trj onlyframes 1 restart
run
quit
' > cpptraj_pdb2trj.in
cpptraj < cpptraj_pdb2trj.in

#step4 acpype
acpype -p 4j0p-1H8_e10ps-10000ps-strip.prmtop -x 4j0p-10000ps-strip.trj

cd 4j0p-1H8_e10ps-10000ps-strip.amb2gmx

#step5 edit vaccum energy minimize setting
echo "
integrator    = steep
emtol         = 1000
emstep        = 0.01
nsteps        = 50000
nstxout       = 100

nstlist       = 1
cutoff-scheme = Verlet
ns_type       = grid
coulombtype   = Cut-off
rcoulomb      = 333.3
vdwtype = Cut-off
rvdw          = 333.3
pbc           = xyz
rlist         = 333.3
" > em.mdp

#step6 add restraint to top
cat posre_4j0p-1H8_e10ps-10000ps-strip.itp >>  4j0p-1H8_e10ps-10000ps-strip_GMX.top

#step7 set big box
gmx editconf -f 4j0p-1H8_e10ps-10000ps-strip_GMX.gro -o 4j0p-1H8_e10ps-10000ps-strip_GMXbox.gro -d 700 -bt cubic

#step8 create tpr
gmx grompp -f em.mdp -c 4j0p-1H8_e10ps-10000ps-strip_GMXbox.gro -r 4j0p-1H8_e10ps-10000ps-strip_GMXbox.gro -p 4j0p-1H8_e10ps-10000ps-strip_GMX.top -o 4j0p-1H8_e10ps-10000ps-strip_GMX.tpr

#step9 run
gmx mdrun -ntomp 2 -v -deffnm 4j0p-1H8_e10ps-10000ps-strip_GMX

#step10 cpptraj
echo '
parm ../4j0p-1H8_e10ps-10000ps-strip.prmtop
trajin 4j0p-1H8_e10ps-10000ps-strip_GMX.trr lastframe
autoimage anchor :1-381 origin
trajout 4j0p-1H8_e10ps-10000ps-strip_GMX.pdb
run
quit
' > cpptraj_genoptpdb.in
cpptraj < cpptraj_genoptpdb.in

echo '
parm ../4j0p-1H8_e10ps-10000ps-strip.prmtop
trajin 4j0p-1H8_e10ps-10000ps-strip_GMX.pdb
reference 4j0p-1H8_e10ps-10000ps-strip_GMX.pdb
mask (:1-381<:4.0)|:NA+|:NA|:Na+|:CL|:Cl- maskpdb 4j0p-10000ps-strip-optmasked-4.0ang.pdb
run
quit
' > cpptraj_genoptpdb2.in
cpptraj < cpptraj_genoptpdb2.in

mv 4j0p-10000ps-strip-optmasked-4.0ang.pdb.1 4j0p-10000ps-strip-optmasked-4.0ang.pdb

