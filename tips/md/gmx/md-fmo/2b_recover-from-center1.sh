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

# 検証する溶質の部分 (残基番号。FMO のフラグメント番号ではない)
fragargs="--residues A:1-840 --residues B:841-1184 --residues lig:1185"
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

# 検証用の python は 2_optmask-frame.sh と同じ探し方をする。
# システムの python3 に scipy が無い環境がある (富岳のログインノードは
# 3.6.8 で scipy 無し)。また amber.sh が PYTHONPATH に入れる Amber の
# py3.9 site-packages は、別の python から見ると割り込みになるので外す。
_clean_pp() {
    local out="" e
    IFS=: read -ra _pp <<<"${PYTHONPATH:-}"
    for e in "${_pp[@]:-}"; do
        [ -z "$e" ] && continue
        case "$e" in
            ${AMBERHOME:-/nonexistent}/*) continue;;
            */amber*/lib/python*/site-packages) continue;;
        esac
        out="${out:+$out:}$e"
    done
    printf '%s' "$out"
}

_py=""
for c in "${PYTHON:-}" "${VIRTUAL_ENV:+$VIRTUAL_ENV/bin/python}" \
         "$HOME/fmoenv/bin/python" python3 python; do
    [ -z "$c" ] && continue
    command -v "$c" >/dev/null 2>&1 || [ -x "$c" ] || continue
    if PYTHONPATH="$(_clean_pp)" "$c" -c 'import numpy, scipy' \
            >/dev/null 2>&1; then
        _py=$c; break
    fi
done

if [ -z "$_py" ]; then
    echo "*** 検証していない ***" >&2
    echo "  numpy と scipy の入った python が無いので check_contact.py を" >&2
    echo "  実行できなかった。液滴は出来ているが未検証。後で流すこと:" >&2
    echo "    python ${base}/check_contact.py $fragargs ${base}/${outdir}/*_fmo_mask.pdb" >&2
    exit 2
fi

PYTHONPATH="$(_clean_pp)" "$_py" "${base}/check_contact.py" $fragargs \
    "${base}/${outdir}"/*_fmo_mask.pdb
