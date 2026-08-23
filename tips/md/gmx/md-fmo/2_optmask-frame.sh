#!/bin/bash
#
# 各フレームを最小化し、PBC を組み直して FMO 用の液滴を切り出す。
#
#   gmx grompp/mdrun (最小化)
#     -> gmx trjconv -pbc whole      (各フラグメント内部の割れを直す)
#     -> cpptraj autoimage           (フラグメント間の相対位置を組み直す)
#     -> cpptraj mask                (溶質から stripdist 以内の水を残す)
#     -> check_contact.py            (接触と水和殻を検証)
#
# PBC 再構成に cpptraj autoimage を使う理由 (詳細は README.md):
#   複合体 (タンパク A + タンパク B + リガンド等) が「1 つの GROMACS
#   moleculetype」になっている系では、gmx trjconv -pbc cluster が
#   反復アルゴリズムの初期配置依存で、正しく組み上がった複合体を
#   壊すことがある (200 ns 系の 20 フレーム中 11 で破綻を確認)。
#   壊れたフレームは複合体が数十 A 離れ、水和殻に穴が空く。
#   cpptraj autoimage は prmtop の分子情報から各フラグメントを別分子と
#   認識し「anchor に対する最近接イメージ」を決定論的に選ぶため破綻しない。
#
#   2026-08 まではこの手順を 2_optmask-frame_v2.sh として別に持ち、本体は
#   -pbc whole -> cluster -> mol/compact -> fit の 4 段だった。破綻する方を
#   残しておく理由が無いので、本体をこちらへ統合した。
#
# 前提: 0_parmed.sh で ${grotop%.*}.prmtop が生成済みであること。
#
# 実行例:
#   bash 2_optmask-frame.sh -n index.solute.ndx -p system.top -f traj.xtc
#   bash 2_optmask-frame.sh -n index.solute.ndx -p system.top -f traj.xtc -t 8
#
##################

set -euo pipefail

tgtindexname="solute"
minscript=min_aftermd.mdp

# --- autoimage の分子指定 (prmtop の残基番号; 系ごとに要編集) ---
#   anchor : 箱の中心に固定する分子 (通常タンパク片方)
#   fixed  : anchor の最近接イメージへ寄せる分子 (相手タンパク + リガンド)
#   mobile : 箱に詰め直す溶媒
#
#   ★ 範囲は prmtop の残基数と突き合わせて自動検査する (validate_mask)。
#     cpptraj は範囲を超えたレンジを無警告で丸めるため、たとえば 324 残基の
#     溶質しかない系に ":1-840" と書くと anchor に溶媒まで入り、autoimage の
#     基準が静かに壊れる。実測: 334 残基の系で ":1-840" は全 1055 原子を選択し、
#     警告は一切出なかった。
#
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

# minimize のスレッド数。-t か環境変数 OPT_THREADS で上書きできる。
#   ログインノードで回すことがあるので既定は控えめにする。共有ノードの
#   コアを埋めると他の利用者の作業を止めるため、大きくするのは計算ノードで
#   ジョブとして流すときだけにすること。
optomp=${OPT_THREADS:-4}

# 検証 (check_contact.py) を走らせる python。numpy と scipy が要る。
#   システムの python3 が古い / scipy が無い環境 (富岳のログインノードは
#   python3 が 3.6.8 で scipy 無し) では venv の python を渡す:
#     PYTHON=/path/to/venv/bin/python bash 2_optmask-frame.sh ...
PYTHON=${PYTHON:-python3}

##################

usage() {
    echo "usage: $0 -n <index.ndx> -p <topol.top> -f <traj> [-t <threads>]" >&2
    exit 1
}

if [ $# = 0 ]; then
    echo "error! please input argment -p and -f" >&2
    usage
fi

solindex=""; grotop=""; traj=""
while getopts n:p:f:t: OPTION
do
    case $OPTION in
        n) solindex=$OPTARG;;
        p) grotop=$OPTARG;;
        f) traj=$OPTARG;;
        t) optomp=$OPTARG;;
        *) usage
    esac
done

for v in grotop traj solindex; do
    if [ -z "${!v}" ]; then
        echo "option for ${v} was NOT given, exit." >&2
        usage
    fi
done

echo "$grotop"
echo "$traj"
echo "$solindex"

headbuf=${traj%.*}
head=${headbuf##*/}
prmtop=${grotop%.*}.prmtop

# ---------------------------------------------------------------- helpers

# 出力が出来ているかを段ごとに確認する。set -e だけでは、直前のコマンドが
# 0 で返りつつ空ファイルを吐いた場合 (例: gmx trjconv に .pdb を渡した場合)
# を捕まえられない。
require() {
    if [ ! -s "$1" ]; then
        echo "ERROR: $1 が出来ていない (${2:-直前の段が失敗}) ので中断する。" >&2
        exit 1
    fi
}

#: 溶媒とみなす残基名 (check_contact.py の WATER_RESNAMES と同じ)。
_SOLVENT_RESNAMES="WAT HOH SOL TIP3 T3P"

# Amber prmtop の %FLAG POINTERS 12 番目 = NRES。10I8 固定長で読む。
nres_of_prmtop() {
    awk '/^%FLAG POINTERS/{f=1; next}
         f && /^%FORMAT/{next}
         f { for (i = 1; i + 7 <= length($0); i += 8) {
                 n++
                 if (n == 12) { v = substr($0, i, 8); gsub(/ /, "", v); print v; exit }
             } }' "$1"
}

# マスクが参照する最大残基番号が prmtop の残基数を超えていないか調べる。
# 超えていても cpptraj はエラーにならず、静かに違うものを選ぶ。
validate_mask() {
    local name=$1 mask=$2 nres=$3 mx
    mx=$(printf '%s' "$mask" | grep -oE '[0-9]+' | sort -n | tail -1 || true)
    [ -z "$mx" ] && return 0
    if [ "$mx" -gt "$nres" ]; then
        echo "ERROR: ${name}=\"${mask}\" は残基 ${mx} を指しているが," >&2
        echo "       ${prmtop} は ${nres} 残基しかない。" >&2
        echo "       cpptraj は範囲外を無警告で丸めるので、このままだと" >&2
        echo "       溶媒まで巻き込んだまま静かに間違った結果が出る。" >&2
        echo "       スクリプト冒頭のマスクを系に合わせて直すこと。" >&2
        exit 1
    fi
}

# anchor / fixed に溶媒が入っていないか調べる。
#
# 範囲チェック (validate_mask) だけでは足りない。溶媒つきの系では NRES が
# 水込みで数万になるので、溶質しか無い前提で書いた ":1-840" のようなマスクが
# 範囲内に収まってしまい素通りする。実測: EGFR (NRES=22215) で ":1-840" は
# 検査を通り、水 516 残基を anchor に巻き込む。
# autoimage の anchor に溶媒が入ると基準が壊れるので、ここで止める。
solvent_in_mask() {
    local name=$1 mask=$2 n
    n=$(printf '%s' "$mask" | grep -oE '[0-9]+(-[0-9]+)?' | awk -v RS=' ' '{print}' |
        awk -F- 'NF==1{print $1"|"$1} NF==2{print $1"|"$2}' |
        awk -F'|' -v labels="$(residue_labels)" -v sol="$_SOLVENT_RESNAMES" '
            BEGIN { split(labels, L, " "); split(sol, S, " ")
                    for (i in S) issol[S[i]] = 1 }
            { for (r = $1; r <= $2; r++) if (issol[L[r]]) c++ }
            END { print c + 0 }')
    if [ "${n:-0}" -gt 0 ]; then
        echo "ERROR: ${name}=\"${mask}\" に溶媒が ${n} 残基入っている。" >&2
        echo "       autoimage の anchor / fixed は溶質だけを指すこと。" >&2
        echo "       (範囲チェックは通る。溶媒つきの系では NRES が水込みで" >&2
        echo "        数万になるため、溶質前提のマスクが範囲内に収まってしまう)" >&2
        echo "       スクリプト冒頭のマスクを系に合わせて直すこと。" >&2
        exit 1
    fi
}

# prmtop の %FLAG RESIDUE_LABEL を空白区切りで返す (1 始まりで添字が残基番号)。
residue_labels() {
    awk '/^%FLAG RESIDUE_LABEL/{f=1; next}
         f && /^%FORMAT/{next}
         f && /^%FLAG/{exit}
         f { for (i = 1; i + 3 <= length($0); i += 4) {
                 v = substr($0, i, 4); gsub(/ /, "", v)
                 if (v != "") printf "%s ", v
             } }' "$prmtop"
}

# ---------------------------------------------------------------- checks

require "$minscript" "minimize 用 mdp が無い"
require "$grotop"    "top が無い"
require "$prmtop"    "0_parmed.sh を先に流して prmtop を作ること"

nres=$(nres_of_prmtop "$prmtop")
if [ -z "$nres" ]; then
    echo "ERROR: ${prmtop} から残基数を読めなかった。" >&2
    exit 1
fi
echo "prmtop residues = ${nres}"

validate_mask anchormask "$anchormask" "$nres"
validate_mask fixedmask  "$fixedmask"  "$nres"
validate_mask maskinfo   "$maskinfo"   "$nres"
[ "$dofit" = "1" ] && validate_mask fitmask "$fitmask" "$nres"

solvent_in_mask anchormask "$anchormask"
solvent_in_mask fixedmask  "$fixedmask"

echo "threads for minimize = ${optomp}"

# ---------------------------------------------------------------- steps

function genref() {
    if [ -s "${head}_ref.tpr" ]; then
        echo "[skip] ${head}_ref.tpr は作成済み"
        return
    fi
    require "${head}_0.gro" "1_trajsep.sh の出力が無い"
    gmx grompp -f "$minscript" -c "${head}_0.gro" -r "${head}_0.gro" \
        -p "${grotop}" -o "${head}_ref.tpr" -maxwarn 1
    require "${head}_ref.tpr"
}

function minimize() {
for i in $(seq $pdb_snum $pdb_interval $pdb_enum)
do
    d=${head}_${i}_fmo
    mkdir -p "$d"

    # 済んでいるフレームは飛ばす (再実行できるようにする)
    if [ -s "${d}/${d}.gro" ]; then
        echo "[skip] ${d} は minimize 済み"
        continue
    fi

    # 1 回目は直下、2 回目以降は ${d}/ に入力がある
    src=${head}_${i}.gro
    [ -f "$src" ] || src=${d}/${head}_${i}.gro
    require "$src" "フレーム ${i} の gro が無い"

    gmx grompp -f "$minscript" -c "$src" -r "$src" -p "${grotop}" \
        -o "${d}.tpr" -maxwarn 1
    require "${d}.tpr"

    # -ntmpi/-ntomp を明示し、OMP_NUM_THREADS もそろえる。
    # -nt だけを渡すと、環境に OMP_NUM_THREADS が居るとき
    #   "The total number of threads requested (12) is not divisible by
    #    the number of OpenMP threads requested (8)"
    # で落ちる。HPC のジョブスクリプトはたいてい OMP_NUM_THREADS を撒くので、
    # 環境任せにしない。
    OMP_NUM_THREADS=${optomp} gmx mdrun -ntmpi 1 -ntomp "${optomp}" \
        -v -deffnm "${d}"
    require "${d}.gro" "gmx mdrun が構造を書けていない"

    mv -f ${head}_${i}.* "$d"/ 2>/dev/null || true
    mv -f ${d}.* "$d"/
done
}

function arrangetraj() {
for i in $(seq $pdb_snum $pdb_interval $pdb_enum)
do
    d=${head}_${i}_fmo
    (
    cd "$d"

    # step1: 結合をたどって各フラグメント内部が周期境界で割れないようにする
    #        (決定論的で安全。フラグメント間の相対位置はここでは直らない)
    gmx trjconv -f "${d}.gro" -s "../${head}_ref.tpr" \
        -o "${d}_whole.pdb" -pbc whole << EOF
System
EOF
    require "${d}_whole.pdb" "gmx trjconv -pbc whole が空を吐いた"

    # step2: cpptraj autoimage でフラグメント間の相対位置を決定論的に組み直し、
    #        続けて液滴 (溶質から stripdist 以内の水) を切り出す
    if [ "$dofit" = "1" ]; then
        fitline="rms reference ${fitmask}"
    else
        fitline=""
    fi

    cat > cpptraj_arrange.in <<EOF
parm ../${prmtop}
trajin ${d}_whole.pdb
autoimage anchor ${anchormask} fixed ${fixedmask} mobile ${mobilemask}
${fitline}
outtraj ${d}_arranged.pdb pdb
mask (:${maskinfo}<:${stripdist})${retainedions} maskpdb ${d}_mask.pdb
run
quit
EOF
    cpptraj -i cpptraj_arrange.in

    # maskpdb はフレーム番号を付けて書く
    if [ -s "${d}_mask.pdb.1" ]; then
        mv -f "${d}_mask.pdb.1" "${d}_mask.pdb"
    fi
    require "${d}_mask.pdb" "cpptraj の autoimage/mask が失敗した"
    )
done
}

genref
minimize
arrangetraj

mkdir -p "${head}-optedpdb"
cp *_fmo/*_fmo_mask.pdb "${head}-optedpdb"/

# 最後に必ず検証。NG フレームは FMO に使わないこと。
# (--frag は anchor/fixed の残基レンジに合わせて編集)
echo "==================== contact check ===================="
"${PYTHON}" check_contact.py --frag A:1-840 --frag B:841-1184 --frag lig:1185 \
    "${head}"-optedpdb/*_fmo_mask.pdb
