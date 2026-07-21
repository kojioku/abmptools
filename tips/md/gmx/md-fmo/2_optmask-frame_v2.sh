#!/bin/bash
#
# 2_optmask-frame.sh の置き換え版 (PBC 再構成を cpptraj autoimage 化)。
#
# 旧 2_optmask-frame.sh の PBC 手順:
#   -pbc whole -> -pbc cluster -> -pbc mol -ur compact -center -> -fit rot+trans
# 本版:
#   -pbc whole -> cpptraj autoimage (anchor 指定) -> mask
#
# なぜ変えるか (詳細は README.md):
#   複合体 (タンパク A + タンパク B + リガンド等) が「1 つの GROMACS
#   moleculetype」になっている系では、gmx trjconv -pbc cluster が
#   反復アルゴリズムの初期配置依存で、正しく組み上がった複合体を
#   壊すことがある (200 ns 系の 20 フレーム中 11 で破綻を確認)。
#   壊れたフレームは複合体が数十 A 離れ、水和殻に穴が空く。
#   cpptraj autoimage は prmtop の分子情報から各フラグメントを別分子と
#   認識し「anchor に対する最近接イメージ」を決定論的に選ぶため破綻しない。
#
# 前提: 0_parmed.sh で ${grotop%.*}.prmtop が生成済みであること。
#
##################
tgtindexname="solute"
minscript=min_aftermd.mdp

# --- autoimage の分子指定 (prmtop の残基番号; 系ごとに要編集) ---
#   anchor : 箱の中心に固定する分子 (通常タンパク片方)
#   fixed  : anchor の最近接イメージへ寄せる分子 (相手タンパク + リガンド)
#   mobile : 箱に詰め直す溶媒
# 例) frag A=res 1-840, frag B=res 841-1184, リガンド=res 1185 の場合:
anchormask=":1-840"
fixedmask=":841-1185"
mobilemask=":WAT"

retainedions="|:NA+|:NA|:Na+|:CL|:Cl-"
maskinfo="1-324"          # 液滴に残す溶質の残基レンジ (系ごとに要編集)
stripdist="6.0"           # 溶質からこの距離 [A] 以内の水を残す

# 全フレームを共通の向きに揃えたい場合は 1 (FMO エネルギーには無影響)
dofit=0
fitmask=":1-840@CA"

##################

pdb_snum=0
pdb_enum=4
pdb_interval=1

optomp=12

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
prmtop=${grotop%.*}.prmtop

function genref() {
    gmx grompp -f $minscript -c ${head}_0.gro -r ${head}_0.gro -p ${grotop} -o ${head}_ref.tpr -maxwarn 1
}

function minimize() {
for i in `seq $pdb_snum $pdb_interval $pdb_enum`
do
    mkdir ${head}_${i}_fmo
    gmx grompp -f $minscript -c ${head}_${i}.gro -r ${head}_${i}.gro -p ${grotop} -o ${head}_${i}_fmo.tpr -maxwarn 1
    gmx mdrun -nt ${optomp} -v -deffnm ${head}_${i}_fmo
    mv ${head}_${i}.* ${head}_${i}_fmo/
    mv ${head}_${i}_fmo.* ${head}_${i}_fmo/
done
}

function arrangetraj() {
for i in `seq $pdb_snum $pdb_interval $pdb_enum`
do
    cd ${head}_${i}_fmo

    # step1: 結合をたどって各フラグメント内部が周期境界で割れないようにする
    #        (決定論的で安全。フラグメント間の相対位置はここでは直らない)
    gmx trjconv -f ${head}_${i}_fmo.gro -s ../${head}_ref.tpr \
        -o ${head}_${i}_fmo_whole.pdb -pbc whole << EOF
System
EOF

    # step2: cpptraj autoimage でフラグメント間の相対位置を決定論的に組み直し、
    #        続けて液滴 (溶質から stripdist 以内の水) を切り出す
    if [ "$dofit" = "1" ]; then
        fitline="rms reference ${fitmask}"
    else
        fitline=""
    fi

    cat > cpptraj_arrange.in <<EOF
parm ../${prmtop}
trajin ${head}_${i}_fmo_whole.pdb
autoimage anchor ${anchormask} fixed ${fixedmask} mobile ${mobilemask}
${fitline}
outtraj ${head}_${i}_fmo_arranged.pdb pdb
mask (:${maskinfo}<:${stripdist})${retainedions} maskpdb ${head}_${i}_fmo_mask.pdb
run
quit
EOF
    cpptraj -i cpptraj_arrange.in

    mv ${head}_${i}_fmo_mask.pdb.1 ${head}_${i}_fmo_mask.pdb
    cd ../
done
}

genref
minimize
arrangetraj

mkdir ${head}-optedpdb
cp *_fmo/*_fmo_mask.pdb ${head}-optedpdb/

# 最後に必ず検証。NG フレームは FMO に使わないこと。
# (--frag は anchor/fixed の残基レンジに合わせて編集)
echo "==================== contact check ===================="
python3 check_contact.py --frag A:1-840 --frag B:841-1184 --frag lig:1185 \
    ${head}-optedpdb/*_fmo_mask.pdb
