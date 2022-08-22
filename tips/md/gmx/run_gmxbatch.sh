#!/bin/bash
# prepare
# xxx.top, xxx.itp: gromacs topology file
# xxx.prmtop: amber topology file (use parmed.sh)
# xxx.xtc: gromacs trajectory file
# cuttrajgmx-foropt.sh: file for generating cut-pdb

##################
prmtop=$1
traj=$2
centerinfo=":1-381"

startframe=1001
endframe=1101 #2001
interval=20

stimeps=10000  # please specify the start time of production run
intervalps=200
##################

dir='gmxpdbs-foropt'
headbuf=${traj%.*}
head=${headbuf##*/}

#step0 create cut-pdb
#bash cuttrjgmx-foropt.sh 4j0p-1H8_e10ps.prmtop 4j0p-1H8_md-samp_e10ps.xtc


function step0_cut-pdb(){
mkdir $dir 2> /dev/null
echo "parm $prmtop" > cpptraj.in
echo "parminfo $prmtop" >> cpptraj.in
echo "trajin $traj" >> cpptraj.in
echo "autoimage anchor $centerinfo origin" >> cpptraj.in

timeps=${stimeps}
for i in `seq $startframe $interval $endframe`
do
    newtraj=${head}_${timeps}ps.pdb
    echo "trajout $newtraj onlyframes $i restart" >> cpptraj.in
    timeps=`echo "$timeps + $intervalps" | bc`
done
echo 'run' >> cpptraj.in
cpptraj < cpptraj.in

mv ${head}*.pdb $dir
cp $prmtop $dir
echo "$dir/*.pdb was generated."

today=$(date "+%Y%m%d")
ofile="getcutpdb.${today}.log"

echo `date` >> $ofile
echo $prmtop >> $ofile
echo $traj >> $ofile
echo $startframe, $endframe, $interval >> $ofile
}

#step1 strip prmtop
function step1_strip-prmtop(){
cd ${dir}
tgtstr=${head}_${tgtps}ps
echo "
parm ${prmtop}
parminfo *
reference ${tgtstr}.pdb
parmstrip (:1-381>:6.0)&!:NA+&!:NA&!:Na+&!:CL-&!:CL&!:Cl-
parmwrite out ${tgtstr}-strip.prmtop
" > cpptraj_strip.in
cpptraj < cpptraj_strip.in
}

#step2 mask str
function step2_mask-str(){
tgtstr=${head}_${tgtps}ps
echo "
parm ${prmtop}
trajin ${tgtstr}.pdb
reference  ${tgtstr}.pdb
mask (:1-381<:6.0)|:NA+|:NA|:Na+|:CL|:Cl- maskpdb ${tgtstr}-strip.pdb
run
" > cpptraj_mask.in
cpptraj < cpptraj_mask.in
}

#step3 pdb2trj
function step3_pdb2trj(){
tgtstr=${head}_${tgtps}ps
echo "
parm ${tgtstr}-strip.prmtop
trajin ${tgtstr}-strip.pdb.1
trajout ${tgtstr}-strip.trj onlyframes 1 restart
run
quit
" > cpptraj_pdb2trj.in
cpptraj < cpptraj_pdb2trj.in
}

#step4 acpype
function step4_acpype(){
tgtstr=${head}_${tgtps}ps
acpype -p ${tgtstr}-strip.prmtop -x ${tgtstr}-strip.trj
}

#step5 edit vaccum energy minimize setting
function step5_edit-vaccum-emset(){
tgtstr=${head}_${tgtps}ps
cd ${tgtstr}-strip.amb2gmx
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
}

#step6 add restraint to top
function step6-9_emrun(){
tgtstr=${head}_${tgtps}ps
cat posre_${tgtstr}-strip.itp >>  ${tgtstr}-strip_GMX.top

#step7 set big box
gmx editconf -f ${tgtstr}-strip_GMX.gro -o ${tgtstr}-strip_GMXbox.gro -d 700 -bt cubic

#step8 create tpr
gmx grompp -f em.mdp -c ${tgtstr}-strip_GMXbox.gro -r ${tgtstr}-strip_GMXbox.gro -p ${tgtstr}-strip_GMX.top -o ${tgtstr}-strip_GMX.tpr

#step9 run
gmx mdrun -ntomp 2 -v -deffnm ${tgtstr}-strip_GMX
}

#step10 getlastpdb
function step10_getlastpdb(){
tgtstr=${head}_${tgtps}ps
echo "
parm ../${tgtstr}-strip.prmtop
trajin ${tgtstr}-strip_GMX.trr lastframe
autoimage anchor $centerinfo origin
trajout ${tgtstr}-strip_GMX.pdb
run
quit
" > cpptraj_genoptpdb.in
cpptraj < cpptraj_genoptpdb.in

echo "
parm ../${tgtstr}-strip.prmtop
trajin ${tgtstr}-strip_GMX.pdb
reference ${tgtstr}-strip_GMX.pdb
mask (${centerinfo}<:4.0)|:NA+|:NA|:Na+|:CL|:Cl- maskpdb ${tgtstr}-optmasked-4.0ang.pdb
run
quit
" > cpptraj_genoptpdb2.in
cpptraj < cpptraj_genoptpdb2.in

mv ${tgtstr}-optmasked-4.0ang.pdb.1 ${tgtstr}-optmasked-4.0ang.pdb
}


step0_cut-pdb
timeps=${stimeps}
for i in `seq $startframe $interval $endframe`
do
    tgtps=${timeps}
    step1_strip-prmtop
    step2_mask-str
    step3_pdb2trj
    step4_acpype
    step5_edit-vaccum-emset
    step6-9_emrun
    step10_getlastpdb
    timeps=`echo "$timeps + $intervalps" | bc`
    cd ../../
done


