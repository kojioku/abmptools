#!/bin/bash
# prepare
# xxx.top, xxx.itp: gromacs topology file
# xxx.prmtop: amber topology file (use parmed.sh)
# xxx.xtc: gromacs trajectory file
# cuttrajgmx-foropt.sh: file for generating cut-pdb

##################
centerinfo=":1-383,596"
maskinfo=":1-383,596"
stripdist=4.0

startframe=1001
endframe=2001 #2001
interval=20

stimeps=10000  # please specify the start time of production run
intervalps=200

optomp=8
optmpi=2

is_rmNaClion=true
##################

if [ $# = 0 ]; then
    echo error! please input argment -p and -x
    exit 1
fi
while getopts p:x: OPTION
do
    case $OPTION in
        p) prmtop=$OPTARG;;
        x) traj=$OPTARG;;
        *) exit 1
    esac
done

if [ -z "$prmtop" ];then
    echo "option p was NOT given, exit."
    exit 1
fi

if [ -z "$traj" ];then
    echo "option x was NOT given, exit."
    exit 1
fi

echo $prmtop
echo $coord

dir='gmxpdbs-foropt'
headbuf=${traj%.*}
head=${headbuf##*/}

#step0 create cut-pdb
#bash cuttrjgmx-foropt.sh 4j0p-1H8_e10ps.prmtop 4j0p-1H8_md-samp_e10ps.xtc

function step0_cuttrj(){
curtime=$(date "+%Y%m%d-%H%M")
if [ -d $dir ]; then
    echo "Backoff! old $dir is renamed to ${dir}.bak.${curtime}"
    mv $dir $dir.bak.${curtime}
    sleep 1
fi
mkdir $dir 2> /dev/null
echo "parm $prmtop" > cpptraj.in
echo "parminfo $prmtop" >> cpptraj.in
echo "trajin $traj" >> cpptraj.in
echo "autoimage anchor $centerinfo origin" >> cpptraj.in

timeps=${stimeps}
for i in `seq $startframe $interval $endframe`
do
    newtraj=${head}_${timeps}ps.rst
    echo "trajout $newtraj onlyframes $i restart" >> cpptraj.in
    timeps=`echo "$timeps + $intervalps" | bc`
done
echo 'run' >> cpptraj.in
cpptraj < cpptraj.in

mv ${head}*.rst $dir
cp $prmtop $dir
echo "$dir/*.rst was generated."

ofile="getcutrst.${curtime}.log"

echo `date` >> $ofile
echo $prmtop >> $ofile
echo $traj >> $ofile
echo $startframe, $endframe, $interval >> $ofile
}



#step4 acpype
function step4_acpype(){
tgtstr=${head}_${tgtps}ps
acpype -p ${prmtop} -x ${tgtstr}.rst -b ${tgtstr}
}

#step5 edit vaccum energy minimize setting
function step5_edit-vaccum-emset(){

echo "
integrator    = steep
nsteps        = 10000
emtol         = 1000
emstep        = 0.01
nstcomm       = 100
comm-mode     = Linear
comm-grps     = System

nstlog        = 1
nstxout       = 0
nstvout       = 0
nstfout       = 0
nstlist       = 100
ns_type       = grid
pbc           = xyz
rlist         = 1.0
coulombtype   = PME
rcoulomb      = 1.2
rvdw          = 1.2
vdwtype       = cut-off
constraints   = none
cutoff-scheme = Verlet
" > em.mdp
}

#step6 add restraint to top
function step6-9_emrun(){
tgtstr=${head}_${tgtps}ps
cat posre_${tgtstr}.itp >>  ${tgtstr}_GMX.top

#step7 set big box
# gmx editconf -f ${tgtstr}_GMX.gro -o ${tgtstr}_GMXbox.gro -d 700 -bt cubic

#step8 create tpr
gmx grompp -f em.mdp -c ${tgtstr}_GMX.gro -r ${tgtstr}_GMX.gro -p ${tgtstr}_GMX.top -o ${tgtstr}_GMX.tpr

#step9 run
gmx mdrun -ntomp ${optomp} -ntmpi ${optmpi} -v -deffnm ${tgtstr}_GMX


}

#step10 getlastpdb
function step10_getlastpdb(){
echo "0
0" > conv.in
gmx trjconv -s ${tgtstr}_GMX.tpr -f ${tgtstr}_GMX.trr -o ${tgtstr}_GMXjump.trr -pbc nojump -center < conv.in

tgtstr=${head}_${tgtps}ps
echo "
parm ../${prmtop}
trajin ${tgtstr}_GMXjump.trr lastframe
unwrap ${maskinfo}
autoimage anchor ${centerinfo} origin
trajout ${tgtstr}_GMXjump.pdb
run
quit
" > cpptraj_genoptpdb.in
cpptraj < cpptraj_genoptpdb.in

echo "
parm ../${prmtop}
trajin ${tgtstr}_GMXjump.pdb
reference ${tgtstr}_GMXjump.pdb
mask (${maskinfo}<:${stripdist})|:NA+|:NA|:Na+|:CL|:Cl- maskpdb ${tgtstr}-optmasked-${stripdist}ang.pdb
run
quit
" > cpptraj_genoptpdb2.in
cpptraj < cpptraj_genoptpdb2.in

mv ${tgtstr}-optmasked-${stripdist}ang.pdb.1 ${tgtstr}-optmasked-${stripdist}ang.pdb
}

function stepop1_rmnaclion(){
cd ${head}-optedpdb
for f in `ls *.pdb`
do
    cat $f |grep -v "NA"| grep -v "CL" > ${f%.*}-noion.pdb
done
}

# step0_cuttrj
timeps=${stimeps}
for i in `seq $startframe $interval $endframe`
do
    tgtps=${timeps}
    tgtstr=${head}_${tgtps}ps
    cd ${dir}
    # step4_acpype
    cd ${tgtstr}.amb2gmx
    # step5_edit-vaccum-emset
    # step6-9_emrun
    step10_getlastpdb
    timeps=`echo "$timeps + $intervalps" | bc`
    cd ../../
done

#last copy optedpdb
curtime=$(date "+%Y%m%d-%H%M")
if [ -d $dir ]; then
    echo "Backoff! old $dir is renamed to ${dir}.bak.${curtime}"
    mv ${head}-optedpdb ${head}-optedpdb.bak.${curtime}
    sleep 1
fi
mkdir ${head}-optedpdb
cp $dir/*amb2gmx/*ang.pdb ${head}-optedpdb

#step-option1 rmnaclion
if [ "${is_rmNaClion}" ]; then
    echo "### rm NaClion ###"
    stepop1_rmnaclion
    echo "Out: ${head}-optedpdb/xxx-noion.pdb"
    mkdir w_ion
    mv *ang.pdb w_ion/
fi

