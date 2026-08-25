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
#   残しておく理由が無いので本体をこちらへ統合し、_v2 はリポジトリ外へ退避した。
#
# 前提: 0_parmed.sh で ${grotop%.*}.prmtop が生成済みであること。
#
# 途中で落ちても、同じコマンドで再実行すれば続きから進む。判断はフレーム
# ディレクトリの中の完了マーカー (.done.minimize / .done.arrange) だけで行い、
# **マーカーは中身を検めたあとに最後に書く**。ディレクトリがあることも、
# ファイルがあることも、済んだ証拠には使わない —— mkdir は作業の前に走るし、
# 殺されて切り詰められたファイルも「存在する」からである。
# 何がどこまで済んでいるかは --check で確認できる。
#
# 上書きで消すことはしない。作り直しが要る場面では park() で
# <名前>.<日時> にどけてから進む。
#
# 実行例:
#   bash 2_optmask-frame.sh -n index.solute.ndx -p system.top -f traj.xtc \\
#       --anchor :1-323 --fixed :324 --solute 1-324 --frames 0:4:1 -t 4
#
#   系ごとの値はすべてコマンドラインで渡せる。検証 (check_contact.py) に渡す
#   残基レンジは anchor / fixed から自動で作るので、同じレンジを二度書かない。
#   bash 2_optmask-frame.sh ... --check              # 進み具合だけ見る
#   bash 2_optmask-frame.sh ... --redo-from arrange  # arrange からやり直す
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
# 既定値。スクリプトを書き換えなくても --anchor / --fixed 等で上書きできるので、
# 系ごとの指定はコマンドラインで渡すことを勧める (編集箇所が散らばらない)。
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

# 検証 (check_contact.py) に使う python は自動で探す (下の find_python)。
# 明示したいときだけ PYTHON=/path/to/python を渡す。
PYTHON=${PYTHON:-}

##################

usage() {
    cat >&2 <<'USAGE'
usage: 2_optmask-frame.sh -n <index.ndx> -p <topol.top> -f <traj> [options]

  -t <threads>        minimize のスレッド数 (既定 4 / $OPT_THREADS)

系ごとの指定 (スクリプトを書き換えず、ここで渡せる):
  --mdp    <file>     minimize の条件ファイル  既定 min_aftermd.mdp
  --anchor <mask>     箱の中心に固定する分子   例 :1-323
  --fixed  <mask>     anchor に寄せる分子      例 :324
  --mobile <mask>     詰め直す溶媒             既定 :WAT
  --solute <range>    液滴に残す溶質のレンジ   例 1-324
  --strip  <A>        溶質から何 A の水を残すか 既定 6.0
  --frames <s:e[:i]>  処理するフレーム番号     例 0:4:1
  --fit    <mask>     全フレームの向きを揃える 例 :1-323@CA
  --residues <name:lo-hi>  検証する溶質の部分 (残基番号)。
                      既定は anchor / fixed から自動導出するので通常は不要
  --check             何もせず、フレームごとの進み具合だけを表示する
  --redo              マーカーを無視して全部やり直す
  --redo-from <stage> minimize | arrange のどちらかからやり直す

途中で落ちた場合はそのまま同じコマンドで再実行すれば、済んだ段は飛ばして
続きから進む。何が済んでいるかは --check で確認できる。
USAGE
    exit 1
}

if [ $# = 0 ]; then
    echo "error! please input argment -p and -f" >&2
    usage
fi

solindex=""; grotop=""; traj=""
fragspec=""         # 空なら anchor/fixed から導出する (下記 derive_fragspec)
mode="run"          # run | check
redo_from=""        # "" | minimize | arrange
while [ $# -gt 0 ]
do
    case $1 in
        -n) solindex=${2:-}; shift 2;;
        -p) grotop=${2:-};   shift 2;;
        -f) traj=${2:-};     shift 2;;
        -t) optomp=${2:-};   shift 2;;
        --mdp)    minscript=${2:-}; shift 2;;
        --anchor) anchormask=${2:-}; shift 2;;
        --fixed)  fixedmask=${2:-};  shift 2;;
        --mobile) mobilemask=${2:-}; shift 2;;
        --solute) maskinfo=${2:-};   shift 2;;
        --strip)  stripdist=${2:-};  shift 2;;
        --fit)    fitmask=${2:-}; dofit=1; shift 2;;
        --frames)
            # s:e:i (i 省略時は 1)
            IFS=: read -r pdb_snum pdb_enum pdb_interval <<<"${2:-}"
            pdb_interval=${pdb_interval:-1}
            case "${pdb_snum}${pdb_enum}${pdb_interval}" in
                *[!0-9]*|'') echo "--frames は s:e[:i] の形で。" >&2; usage;;
            esac
            shift 2;;
        --residues) fragspec="${fragspec} --residues ${2:-}"; shift 2;;
        --check) mode="check"; shift;;
        --redo)  redo_from="minimize"; shift;;
        --redo-from)
            redo_from=${2:-}
            case "$redo_from" in
                minimize|arrange) ;;
                *) echo "--redo-from は minimize か arrange のどちらか。" >&2; usage;;
            esac
            shift 2;;
        -h|--help) usage;;
        *) echo "unknown option: $1" >&2; usage;;
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

# ---------------------------------------------------------------- python

# Amber の site-packages を PYTHONPATH から外した環境変数を組み立てる。
#
# amber.sh は PYTHONPATH に $AMBERHOME/lib/python3.9/site-packages を入れる。
# これは parmed の CLI が動くために必要 (parmed 本体はそこにしか無い) なので
# 消せないが、別の python から見ると 3.9 用のパッケージが割り込んでくる。
# numpy 2 系の venv で Amber の parmed を掴むと
#   ModuleNotFoundError: No module named 'numpy.compat'
# で落ちる。ここでは Amber 由来の項目だけを落とし、利用者が自分で入れた
# PYTHONPATH は残す。
clean_pythonpath() {
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

#: numpy と scipy が入っている python を探す。
#  順に: PYTHON 指定 -> 有効な venv -> ~/fmoenv -> python3 -> python
#  Amber の PYTHONPATH を外した状態で判定する (外さないと 3.9 用が割り込む)。
find_python() {
    local c
    for c in "$PYTHON" "${VIRTUAL_ENV:+$VIRTUAL_ENV/bin/python}" \
             "$HOME/fmoenv/bin/python" python3 python; do
        [ -z "$c" ] && continue
        command -v "$c" >/dev/null 2>&1 || [ -x "$c" ] || continue
        if PYTHONPATH="$(clean_pythonpath)" "$c" -c 'import numpy, scipy' \
                >/dev/null 2>&1; then
            printf '%s' "$c"
            return 0
        fi
    done
    return 1
}

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

# ------------------------------------------------- 再開のための完了判定
#
# 途中で落ちた実行を再開できるようにする。判定の根拠は 3 つ:
#
#   1. 段ごとの完了マーカー (${d}/.done.<stage>) を **最後に** 書く。
#      ディレクトリの存在も、ファイルの存在も、済んだ証拠にはならない。
#      mkdir は作業の前に走るし、途中で殺されたファイルも「存在する」。
#   2. マーカーを書く前に **中身を検証する**。`-s` (空でない) だけでは、
#      書きかけで切り詰められた .gro や .pdb が完成扱いで素通りする。
#   3. 確定は **一時名 -> rename**。同一ファイルシステム内の rename は
#      原子的なので、殺されても半端な最終ファイルが残らない。
#
# 既存のものを消す操作はしない。上書きが要る場面では park() でどける。

#: gro が最後まで書けているか。2 行目の原子数と実際の行数 (N+3) を突き合わせ、
#: 最終行が箱ベクトル (3 or 9 個の数値) であることまで見る。
gro_is_complete() {
    local f=$1 natoms
    [ -s "$f" ] || return 1
    natoms=$(sed -n '2p' "$f" | tr -d '[:space:]')
    case "$natoms" in ''|*[!0-9]*) return 1;; esac
    awk -v want=$((natoms + 3)) '
        END {
            if (NR != want) exit 1
            if (NF != 3 && NF != 9) exit 1
            for (i = 1; i <= NF; i++)
                if ($i !~ /^-?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$/) exit 1
        }' "$f"
}

#: pdb が最後まで書けているか。原子行があり、かつ END で終わっているか、
#: 終端レコードが無い書き手のために「最終行が完全な原子行」でも通す。
pdb_is_complete() {
    local f=$1
    [ -s "$f" ] || return 1
    grep -qE '^(ATOM|HETATM)' "$f" || return 1
    awk '
        NF { last = $0 }
        END {
            if (last ~ /^END/) exit 0
            # 座標欄 (31-54 桁) まで書けている原子行なら、途中で切れていない
            if (last ~ /^(ATOM|HETATM)/ && length(last) >= 54) exit 0
            exit 1
        }' "$f"
}

#: 段が完了しているか (マーカーの有無だけで判断する)。
stage_done() { [ -f "$1/.done.$2" ]; }

#: 段の完了を記録する。検証を通ったあとにだけ呼ぶこと。
mark_done() {
    printf '%s\n' "$(date '+%Y-%m-%d %H:%M:%S')" > "$1/.done.$2"
}

#: --redo / --redo-from で、この段をやり直すべきか。
redo_stage() {
    case "$redo_from" in
        minimize) return 0;;
        arrange)  [ "$1" = "arrange" ] && return 0 || return 1;;
        *)        return 1;;
    esac
}

#: 既存ファイルを消さずにどける。上書きが要る場面では必ずこれを通す。
park() {
    local f=$1 stamp
    [ -e "$f" ] || return 0
    stamp=$(date '+%Y%m%d-%H%M%S')
    mv "$f" "${f}.${stamp}"
    echo "[park] ${f} -> ${f}.${stamp} (上書きせずどけた)"
}

#: ファイル群を dest へ移す。同名があれば park してから移す (mv -f を使わない)。
move_into() {
    local dest=$1 f; shift
    for f in "$@"; do
        [ -e "$f" ] || continue
        park "${dest}/$(basename "$f")"
        mv "$f" "${dest}/"
    done
}

#: 1 フレームの状態を "済 / 途中 / 未" で返す (--check 用)。
frame_status() {
    local d=$1 stage=$2 artifact=$3 kind=$4
    if stage_done "$d" "$stage"; then
        echo "済"
    elif [ -e "$d/$artifact" ]; then
        if "${kind}_is_complete" "$d/$artifact"; then
            echo "途中(成果物あり・マーカー無し)"
        else
            echo "途中(成果物が壊れている)"
        fi
    elif [ -d "$d" ]; then
        echo "未(ディレクトリのみ)"
    else
        echo "未"
    fi
}

function report() {
    # 桁を揃えない。printf の %-Ns はバイト数で詰めるので、日本語が混ざると
    # 表がずれる。文字数で詰めるには locale に依存した細工が要り、HPC の
    # LANG=C なログインノードで壊れる。区切り記号で読ませる方が確実。
    for i in $(seq $pdb_snum $pdb_interval $pdb_enum)
    do
        d=${head}_${i}_fmo
        echo "frame ${i} | minimize: $(frame_status "$d" minimize "${d}.gro" gro)" \
             "| arrange: $(frame_status "$d" arrange "${d}_mask.pdb" pdb)"
    done
    cat <<'NOTE'

  済     … 完了マーカーがある
  途中   … 前回ここで落ちた。そのまま再実行すれば続きから進む
           (成果物が壊れている場合は park してから作り直す)
  未     … 未着手。ディレクトリだけあるのは mkdir の直後に落ちた状態で、
           これも未着手として扱う
NOTE
}




# check_contact.py の --residues を anchor / fixed から作る。
#
# 検証したいフラグメントは autoimage で組み直したフラグメントそのものなので、
# 同じレンジを二度書く必要はない。二度書くと片方だけ直して気づかない、という
# 事故が起きる (編集箇所が冒頭と末尾で 480 行離れていた)。
#
# ":1-323" -> "prot:1-323"、":324" -> "lig:324" のように名前を付ける。
# カンマ区切りや原子指定を含むマスクは 1 つのレンジに落とせないので、その場合は
# --frag を明示してもらう。
derive_fragspec() {
    local out="" i=0 m name lo hi
    for m in "$anchormask" "$fixedmask"; do
        i=$((i + 1))
        name=$([ $i = 1 ] && echo anchor || echo fixed)
        # ":1-323" / ":324" だけを受け付ける
        case "$m" in
            :[0-9]*-[0-9]*) lo=${m#:}; hi=${lo#*-}; lo=${lo%-*};;
            :[0-9]*)        lo=${m#:}; hi=$lo;;
            *) echo "" ; return 0;;
        esac
        case "$lo$hi" in *[!0-9]*) echo ""; return 0;; esac
        out="${out} --residues ${name}:${lo}-${hi}"
    done
    echo "$out"
}

# ---------------------------------------------------------------- checks

# --residues はマスクから作る。prmtop を読まないので、まだ何も無い段階でも
# 決まる (--check をそこで使えるようにするため、ファイル検査より前に置く)。
if [ -z "$fragspec" ]; then
    fragspec=$(derive_fragspec)
fi

cat <<SETTINGS
--- 設定 ---
  anchor : ${anchormask}
  fixed  : ${fixedmask}
  mobile : ${mobilemask}
  solute : ${maskinfo}  (液滴に残す水: 溶質から ${stripdist} A)
  frames : ${pdb_snum}..${pdb_enum} step ${pdb_interval}
  mdp    : ${minscript}
  verify :${fragspec:- (anchor/fixed から決められないので --residues が要る)}
  threads: ${optomp}
SETTINGS

# --check は何も実行しない。**prmtop も mdp も無い段階で状態を見たい**ので、
# 下のファイル検査より前で返す。設定は上で出しているので、どの値で見ている
# のかは分かる。
if [ "$mode" = "check" ]; then
    report
    exit 0
fi

if [ -z "$fragspec" ]; then
    echo "ERROR: anchormask / fixedmask が単純なレンジ (\":1-323\" 等) では" >&2
    echo "       ないので、検証用の残基レンジを自動で決められない。" >&2
    echo "       --residues <名前>:<開始>-<終了> を明示すること。" >&2
    exit 1
fi

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
    # 基準構造は「処理する最初のフレーム」。0 決め打ちにすると --frames で
    # 途中から始めたときに存在しないファイルを探して止まる。
    local ref=${head}_${pdb_snum}.gro
    [ -f "$ref" ] || ref=${head}_${pdb_snum}_fmo/${head}_${pdb_snum}.gro
    require "$ref" "1_trajsep.sh の出力 (${head}_${pdb_snum}.gro) が無い"
    gmx grompp -f "$minscript" -c "$ref" -r "$ref" \
        -p "${grotop}" -o "${head}_ref.tpr" -maxwarn 1
    require "${head}_ref.tpr"
}

function minimize() {
for i in $(seq $pdb_snum $pdb_interval $pdb_enum)
do
    d=${head}_${i}_fmo
    mkdir -p "$d"

    # 済んだフレームは飛ばす。判断はマーカーだけで行う ——
    # ディレクトリがあることも、.gro があることも、済んだ証拠にはならない。
    if stage_done "$d" minimize && ! redo_stage minimize; then
        echo "[skip] ${d} は minimize 済み"
        continue
    fi

    # マーカーが無いのに .gro がある = 前回ここで落ちたか、マーカー導入前の
    # 実行。中身が完全ならマーカーだけ書いて済ませ、壊れていればどけてやり直す。
    if [ -e "${d}/${d}.gro" ] && ! redo_stage minimize; then
        if gro_is_complete "${d}/${d}.gro"; then
            echo "[adopt] ${d}/${d}.gro は完全なのでマーカーだけ付けて飛ばす"
            mark_done "$d" minimize
            continue
        fi
        echo "[broken] ${d}/${d}.gro は書きかけなのでやり直す"
        park "${d}/${d}.gro"
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
    if ! gro_is_complete "${d}.gro"; then
        echo "ERROR: ${d}.gro が最後まで書けていない (2 行目の原子数と行数が" >&2
        echo "       合わないか、箱ベクトルが無い)。mdrun が途中で死んでいる。" >&2
        exit 1
    fi

    move_into "$d" ${head}_${i}.* 2>/dev/null || true
    move_into "$d" ${d}.*

    # ここまで通ってから記録する。落ちればマーカーは書かれず、次回やり直す。
    mark_done "$d" minimize
done
}

function arrangetraj() {
for i in $(seq $pdb_snum $pdb_interval $pdb_enum)
do
    d=${head}_${i}_fmo

    if stage_done "$d" arrange && ! redo_stage arrange; then
        echo "[skip] ${d} は arrange 済み"
        continue
    fi

    # マーカーが無いのに液滴がある場合。完全なら採用、壊れていればどけて作り直す。
    if [ -e "${d}/${d}_mask.pdb" ] && ! redo_stage arrange; then
        if pdb_is_complete "${d}/${d}_mask.pdb"; then
            echo "[adopt] ${d}/${d}_mask.pdb は完全なのでマーカーだけ付けて飛ばす"
            mark_done "$d" arrange
            continue
        fi
        echo "[broken] ${d}/${d}_mask.pdb は書きかけなので作り直す"
        park "${d}/${d}_mask.pdb"
    fi

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
    # 向き揃え。cpptraj の `rms reference` は **reference を先に読み込んで
    # おかないと使えない**。以前は rms の行だけ出していたので dofit=1 は必ず
    #   Error: Reference index 0 not found.
    #   Error: Could not initialize action [rms]
    # で落ちていた。
    #
    # 基準は最初に処理するフレームの whole.pdb にする。全フレームが同じ 1 つの
    # 構造に合わせられて初めて「共通の向き」になるので、フレームごとに違う基準
    # (自分自身など) を使っては意味がない。最初のフレーム自身は自分に合わせる
    # ことになり、動かない。
    fitline=""
    if [ "$dofit" = "1" ]; then
        # 基準は最初のフレームの **arranged**。whole を基準にすると合わない。
        # whole は autoimage 前で、溶質が周期境界をまたいだままのことがある
        # (実測: frame 42 の whole は基準から 144 A 離れていた)。そこへ
        # autoimage 済みのフレームを合わせると、cpptraj は「割れた形」に
        # 最小二乗で寄せようとして中途半端な向きで止まる。実測で fit 後の
        # CA RMSD が 14 A、最適重ね合わせなら 1.5 A で済むところだった。
        _refdir=${head}_${pdb_snum}_fmo
        _ref=../${_refdir}/${_refdir}_arranged.pdb
        if [ "$i" = "$pdb_snum" ]; then
            # 最初のフレームには基準がまだ無い。自分が基準になるので、
            # ここでは向きを変えない (変えるものが無い)。
            _ref=""
        elif [ ! -s "$_ref" ]; then
            echo "ERROR: 向き揃えの基準 ${_ref} が無い。" >&2
            echo "       最初のフレーム (${pdb_snum}) の arrange が済んでいない。" >&2
            echo "       --redo-from arrange で最初から流し直すこと。" >&2
            exit 1
        fi
        if [ -n "$_ref" ]; then
            # 名前を付けて渡す。`reference <file> [tag]` の角括弧はタグで、
            # `rms [tag] ...` では参照できない (cpptraj は黙って「最初の
            # フレーム」に落ち、自分自身に合わせて RMSD 0 のまま何も動かない)。
            # `name <名前>` を付けて `rms ref <名前>` で指す。
            fitline="reference ${_ref} name fitref
rms ref fitref ${fitmask}"
        fi
    fi

    # 液滴は一時名に書かせ、中身を検めてから rename で確定する。cpptraj に
    # 最終名を直接書かせると、途中で死んだときに半端な ${d}_mask.pdb が残り、
    # 次の実行がそれを完成品として拾ってしまう。
    # (${d}_whole.pdb / ${d}_arranged.pdb は同じ段で作り直す中間物。前者は
    #  gmx が #file# として自動退避する)
    cat > cpptraj_arrange.in <<EOF
parm ../${prmtop}
trajin ${d}_whole.pdb
autoimage anchor ${anchormask} fixed ${fixedmask} mobile ${mobilemask}
${fitline}
outtraj ${d}_arranged.pdb pdb
mask (:${maskinfo}<:${stripdist})${retainedions} maskpdb ${d}_mask.pdb.tmp
run
quit
EOF
    cpptraj -i cpptraj_arrange.in

    # maskpdb はフレーム番号を付けて書く
    tmp=${d}_mask.pdb.tmp.1
    [ -s "$tmp" ] || tmp=${d}_mask.pdb.tmp
    require "$tmp" "cpptraj の autoimage/mask が失敗した"
    if ! pdb_is_complete "$tmp"; then
        echo "ERROR: ${tmp} が最後まで書けていない。cpptraj が途中で死んでいる。" >&2
        exit 1
    fi

    park "${d}_mask.pdb"        # --redo で既存があるとき用。消さずにどける
    mv "$tmp" "${d}_mask.pdb"   # 同一 fs 内の rename = 原子的
    )

    # 液滴が出来てから記録する。落ちればマーカーは書かれず、次回やり直す。
    mark_done "$d" arrange
done
}

genref
minimize
arrangetraj

mkdir -p "${head}-optedpdb"
for i in $(seq $pdb_snum $pdb_interval $pdb_enum)
do
    d=${head}_${i}_fmo
    src=${d}/${d}_mask.pdb
    [ -s "$src" ] || continue
    park "${head}-optedpdb/${d}_mask.pdb"
    cp "$src" "${head}-optedpdb"/
done

# 最後に必ず検証。NG フレームは FMO に使わないこと。
# 検証対象は anchor / fixed から導出済み (derive_fragspec)。
# 渡すのは残基番号であって FMO のフラグメント番号ではない。
echo "==================== contact check ===================="

# numpy / scipy の入った python を探す。見つからない場合は液滴を捨てずに
# 検証だけ保留し、後から流せるコマンドを出す (富岳のシステム python3 は
# 3.6.8 で scipy が無く、python 環境は次のセクションで作るため)。
if ! PYTHON=$(find_python); then
    cat >&2 <<UNVERIFIED
*** 検証していない ***
  numpy と scipy の入った python が見つからないので check_contact.py を
  実行できなかった。液滴は出来ているが、**まだ検証されていない**。
  python 環境を用意したあと、次を実行して確認すること:

    python check_contact.py${fragspec} ${head}-optedpdb/*_fmo_mask.pdb

  (富岳のシステム python3 は 3.6.8 で scipy が無い。venv を activate するか
   PYTHON=/path/to/python を渡す)
UNVERIFIED
    exit 2
fi

PYTHONPATH="$(clean_pythonpath)" "${PYTHON}" check_contact.py ${fragspec} \
    "${head}"-optedpdb/*_fmo_mask.pdb
