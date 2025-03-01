#!/bin/bash

##################
tgtindexname="solute"
minscript=min_aftermd.mdp

retainedions="|:NA+|:NA|:Na+|:CL|:Cl-"
maskinfo="1-324"
stripdist="6.0"

##################

pdb_snum=0
pdb_enum=4
pdb_interval=1

optomp=12
# optmpi=2

##################

if [ $# = 0 ]; then
    echo error! please input argment -p and -x
    exit 1
fi
while getopts n:p:f: OPTION
do
    case $OPTION in
        n) solindex=$OPTARG;;
        p) grotop=$OPTARG;;
        f) traj=$OPTARG;;
        *) exit 1
    esac
done

if [ -z "$grotop" ];then
    echo "option p was NOT given, exit."
    exit 1
fi

if [ -z "$traj" ];then
    echo "option x was NOT given, exit."
    exit 1
fi


if [ -z "$solindex" ];then
    echo "option i was NOT given, exit."
    exit 1
fi

echo $grotop
echo $traj
echo $solindex

headbuf=${traj%.*}
head=${headbuf##*/}

# solindexの作成
# gmx make_ndx -f ${traj%.*}_0_.gro -o index.solute.ndx

# 構造分ループ
function minimize() {
# Enter
for i in `seq $pdb_snum $pdb_interval $pdb_enum`
do
    mkdir ${head}_${i}_fmo
    gmx grompp -f $minscript -c ${head}_${i}.gro -r ${head}_${i}.gro -p ${grotop} -o ${head}_${i}_fmo.tpr -maxwarn 1
    # gmx mdrun -ntomp ${optomp} -ntmpi ${optmpi} -v -deffnm ${head}_${i}_fmo
    gmx mdrun -nt ${optomp} -v -deffnm ${head}_${i}_fmo
    mv ${head}_${i}.* ${head}_${i}_fmo/
    mv ${head}_${i}_fmo.* ${head}_${i}_fmo/
done
}

function genref() {
    # mkdir ${head}_ref
    gmx grompp -f $minscript -c ${head}_0.gro -r ${head}_0.gro -p ${grotop} -o ${head}_ref.tpr -maxwarn 1
    # mv ${head}_ref.* ${head}_ref/
    # mv ${head}_0_fmo.* ${head}_ref/
}

function arrangetraj() {
for i in `seq $pdb_snum $pdb_interval $pdb_enum`
do
    cd ${head}_${i}_fmo
    # 周期境界で割れないように分子を移動
    gmx trjconv -f ${head}_${i}_fmo.gro -s ../${head}_ref.tpr -o ${head}_${i}_fmo_center1.gro -pbc whole << EOF
System
EOF
    # 指定したclusterが離れないようにする
    gmx trjconv -f ${head}_${i}_fmo_center1.gro -s ../${head}_ref.tpr -o ${head}_${i}_fmo_center2.gro -n ../$solindex -pbc cluster << EOF
${tgtindexname}
System
EOF
    # 重心の中心移動
    gmx trjconv -f ${head}_${i}_fmo_center2.gro -s ../${head}_ref.tpr -o ${head}_${i}_fmo_center3.gro -n ../$solindex -pbc mol -ur compact -center << EOF
${tgtindexname}
System
EOF
    # タンパク質（溶質）の向きを揃える
    gmx trjconv -f ${head}_${i}_fmo_center3.gro -s ../${head}_ref.tpr -o ${head}_${i}_fmo_center4.pdb -n ../$solindex -center -fit rot+trans << EOF
${tgtindexname}
${tgtindexname}
System
EOF
    # cpptrajでprmtopにして液滴切り出し
    echo "
    parm ../${grotop%.*}.prmtop
    trajin ${head}_${i}_fmo_center4.pdb
    reference ${head}_${i}_fmo_center4.pdb
    mask (:${maskinfo}<:${stripdist})${retainedions} maskpdb ${head}_${i}_fmo_mask.pdb
    run
    quit
    " > cpptraj_genpdb.in
    cpptraj < cpptraj_genpdb.in
    mv ${head}_${i}_fmo_mask.pdb.1 ${head}_${i}_fmo_mask.pdb
    cd ../
done
}

genref
minimize
arrangetraj

#last copy optedpdb
mkdir ${head}-optedpdb
cp *_fmo/*_fmo_mask.pdb ${head}-optedpdb/
