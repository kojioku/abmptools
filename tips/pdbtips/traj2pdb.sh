#!/bin/bash

##################
sepindexname="!Ion" #!Ion, System
tgtindexname="solute"

skipinterval=1000
# rminfo_1="15"
# rminfo_2="24"
minscript=min_aftermd.mdp

retainedions="|:NA+|:NA|:Na+|:CL|:Cl-"
maskinfo="1-324"
stripdist="6.0"

##################

pdb_snum=0
pdb_enum=10
pdb_interval=1

optomp=4
optmpi=2

##################

if [ $# = 0 ]; then
    echo error! please input argment -p and -x
    exit 1
fi
while getopts i:g:p:x: OPTION
do
    case $OPTION in
        i) solindex=$OPTARG;;
        g) inigro=$OPTARG;;
        p) grotop=$OPTARG;;
        x) traj=$OPTARG;;
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

if [ -z "$inigro" ];then
    echo "option g was NOT given, exit."
    exit 1
fi

if [ -z "$solindex" ];then
    echo "option i was NOT given, exit."
    exit 1
fi

echo $grotop
echo $traj
echo $inigro

dir='gmxpdbs-foropt'
headbuf=${traj%.*}
head=${headbuf##*/}

# 指定した分子種の削除
# echo "!${rminfo_1}
# q" > conv_rm1.in
# gmx make_ndx -f ${traj%.*}.gro < conv_rm1.in

# セパレート
function separate() {
#step0 create .gro for each step
curtime=$(date "+%Y%m%d-%H%M")
if [ -d $dir ]; then
    echo "Backoff! old $dir is renamed to ${dir}.bak.${curtime}"
    mv $dir $dir.bak.${curtime}
    sleep 1
fi

# 作業用ディレクトリの作成
mkdir $dir 2> /dev/null

gmx trjconv -f ${traj} -o ${traj%.*}_.gro -s ${inigro} -skip ${skipinterval} -sep -n index.ndx << EOF
${sepindexname}
EOF

# 各種ファイルの移動とコピー
mv ${traj%.*}_*.gro $dir
cp *.prmtop $dir
cp *.itp $dir
cp $minscript $dir
cp *.top $dir
cp $solindex $dir

# cp ${grotop%.*}_min1.tpr $dir
echo "$dir/*.gro was generated."
}

# 構造分ループ
function minimize() {
# Enter
for i in `seq $pdb_snum $pdb_interval $pdb_enum`
do
    mkdir ${traj%.*}_${i}_fmo
    gmx grompp -f $minscript -c ${traj%.*}_${i}.gro -r ${traj%.*}_${i}.gro -p ${grotop} -o ${traj%.*}_${i}_fmo.tpr -maxwarn 1
    gmx mdrun -ntomp ${optomp} -ntmpi ${optmpi} -v -deffnm ${traj%.*}_${i}_fmo
    mv ${traj%.*}_${i}.* ${traj%.*}_${i}_fmo/
    mv ${traj%.*}_${i}_fmo.* ${traj%.*}_${i}_fmo/
done
}

function arrangetraj() {
for i in `seq $pdb_snum $pdb_interval $pdb_enum`
do
    cd ${traj%.*}_${i}_fmo
    # 周期境界で割れないように分子を移動
    gmx trjconv -f ${traj%.*}_${i}_fmo.gro -s ${traj%.*}_${i}_fmo.tpr -o ${traj%.*}_${i}_fmo_center1.gro -pbc whole << EOF
System
EOF
    # 指定したclusterが離れないようにする
    gmx trjconv -f ${traj%.*}_${i}_fmo_center1.gro -s ${traj%.*}_${i}_fmo.tpr -o ${traj%.*}_${i}_fmo_center2.gro -n ../$solindex -pbc cluster << EOF
${tgtindexname}
System
EOF
    # 重心の中心移動
    gmx trjconv -f ${traj%.*}_${i}_fmo_center2.gro -s ${traj%.*}_${i}_fmo.tpr -o ${traj%.*}_${i}_fmo_center3.gro -n ../$solindex -pbc mol -ur compact -center << EOF
${tgtindexname}
System
EOF
    # タンパク質（溶質）の向きを揃える
    gmx trjconv -f ${traj%.*}_${i}_fmo_center3.gro -s ${traj%.*}_${i}_fmo.tpr -o ${traj%.*}_${i}_fmo_center4.pdb -n ../$solindex -center -fit rot+trans << EOF
${tgtindexname}
${tgtindexname}
System
EOF
    # cpptrajでprmtopにして液滴切り出し
    echo "
    parm ../${grotop%.*}.prmtop
    trajin ${traj%.*}_${i}_fmo_center4.pdb
    reference ${traj%.*}_${i}_fmo_center4.pdb
    mask (:${maskinfo}<:${stripdist})${retainedions} maskpdb ${traj%.*}_${i}_fmo_mask.pdb
    run
    quit
    " > cpptraj_genpdb.in
    cpptraj < cpptraj_genpdb.in
    mv ${traj%.*}_${i}_fmo_mask.pdb.1 ${traj%.*}_${i}_fmo_mask.pdb
    cd ../
done
}

# separate
cd ${dir}
# minimize
# arrangetraj

#last copy optedpdb
mkdir ${head}-optedpdb
cp *_fmo/*_fmo_mask.pdb ${head}-optedpdb/
