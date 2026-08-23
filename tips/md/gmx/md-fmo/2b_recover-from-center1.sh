#!/bin/bash
#
# 既存の *_fmo_center1.pdb から液滴を作り直す (救済用)。
#
# 旧版の 2_optmask-frame.sh (-pbc cluster を挟む 4 段だった頃) を流し終えた後で
# 「複合体が離れている / 水和殻に穴がある」と判明した場合、最小化 (gmx mdrun) を
# やり直さずに液滴だけ作り直せる。
# center1 (= -pbc whole 済み) から cpptraj autoimage で組み直す。
# center2 / center3 / center4 は使わない。
#
# 現行の 2_optmask-frame.sh は最初から autoimage を使うので、これから流す場合に
# このスクリプトは要らない。既に旧版で作った *_fmo_center1.pdb がある場合の救済用。
#
# 必要: cpptraj (AmberTools), python3 + numpy + scipy
#
set -e

##################### 系ごとに編集 #####################
head="<traj-head>"                  # 例: 1_trajsep.sh が作った .gro の共通接頭辞
prmtop="<system>.prmtop"            # 0_parmed.sh の出力
workdir="gmxpdbs-foropt"
outdir="${workdir}/${head}-optedpdb-fixed"

# autoimage の分子指定 (prmtop の残基番号)
anchormask=":1-840"                 # 箱の中心に固定する分子
fixedmask=":841-1185"               # anchor の最近接イメージへ寄せる分子
mobilemask=":WAT"                   # 箱に詰め直す溶媒 (無溶媒系では下の行ごと外す)

maskinfo="1-1185"                   # 液滴に残す溶質の残基レンジ
stripdist="6.0"                     # 溶質からこの距離 [A] 以内の水を残す
retainedions="|:NA+|:NA|:Na+|:CL|:Cl-"

snum=0
enum=19

# 検証用のフラグメント指定 (check_contact.py に渡す)
fragargs="--frag A:1-840 --frag B:841-1184 --frag lig:1185"
########################################################

command -v cpptraj >/dev/null || { echo "cpptraj が見つかりません (amber.sh を source)"; exit 1; }

base=$(cd "$(dirname "$0")" && pwd)
mkdir -p "${base}/${outdir}"

for i in $(seq $snum $enum); do
    d="${base}/${workdir}/${head}_${i}_fmo"
    src="${d}/${head}_${i}_fmo_center1.pdb"

    if [ ! -f "$src" ]; then
        echo "SKIP frame ${i}: ${src} がありません"
        continue
    fi

    cat > "${d}/cpptraj_recover.in" <<EOF
parm ${base}/${prmtop}
trajin ${src}
autoimage anchor ${anchormask} fixed ${fixedmask} mobile ${mobilemask}
mask (:${maskinfo}<:${stripdist})${retainedions} maskpdb ${d}/${head}_${i}_fmo_mask.pdb
run
quit
EOF
    cpptraj -i "${d}/cpptraj_recover.in" > "${d}/cpptraj_recover.log" 2>&1

    # cpptraj の maskpdb 出力は連番 (.1) が付く
    mv "${d}/${head}_${i}_fmo_mask.pdb.1" "${d}/${head}_${i}_fmo_mask.pdb"
    cp "${d}/${head}_${i}_fmo_mask.pdb" "${base}/${outdir}/"
    echo "frame ${i} done"
done

echo
echo "==================== contact check ===================="
python3 "${base}/check_contact.py" $fragargs "${base}/${outdir}"/*_fmo_mask.pdb
