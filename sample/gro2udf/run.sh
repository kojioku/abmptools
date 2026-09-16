#!/usr/bin/env bash
# ---------------------------------------------------------------------------
#  gro2udf のサンプルを一通り流す。
#
#      bash run.sh
#
#  3 本走ります。
#
#    1. udf-and-gro   既存 UDF の座標を .gro で差し替える (後方互換モード)
#    2. from-top      .top + .gro から UDF を新規に組む
#    3. from-top      + 軌跡 + エネルギー  ← gmx が要る
#
#  1 と 2 は gmx を使いません。3 だけが gmx を呼びます —— `--edr` が
#  `gmx energy` を、`--trajectory` が `gmx trjconv -pbc nojump` を走らせる
#  ためです。gmx が無ければ 3 は飛ばします。
#
#  **書き出し先は run_out/ だけ** (git 管理外)。各モードの output/ に置いて
#  ある UDF は**参照用のスナップショット**なので、ここでは触りません ——
#  サンプルを流しただけで作業ツリーが汚れないように。
#
#  環境変数:
#      GMX=/path/to/gmx        PATH に無い gmx を使う
#      PYTHON=python3          別の python を使う
# ---------------------------------------------------------------------------
set -euo pipefail
cd "$(dirname "$0")"

GMX=${GMX:-gmx}
PYTHON=${PYTHON:-python}
OUT=run_out
mkdir -p "$OUT"

echo "=== 1. udf-and-gro mode ==================================="
# 出力名は常に <udf の名前>_groout.udf で、カレントに出ます。選べないので
# run_out/ に cd してから呼びます。
(
    cd "$OUT"
    "$PYTHON" -m abmptools.gro2udf ../udf_and_gro_mode/input/test.udf \
        ../udf_and_gro_mode/input/output.gro
)
echo "    -> $OUT/test_groout.udf"
echo

echo "=== 2. --from-top mode (topology + coordinates) ============"
"$PYTHON" -m abmptools.gro2udf --from-top \
    gro_top_mode/input/test.top gro_top_mode/input/output.gro \
    --template gro_top_mode/input/template_empty.udf \
    --out "$OUT/output_fromtop.udf"
echo "    -> $OUT/output_fromtop.udf"
echo

echo "=== 3. --from-top + trajectory + energy ===================="
if ! command -v "$GMX" > /dev/null 2>&1; then
    cat <<MSG
    [skip] gmx が見つかりません ($GMX)。
           この段は gmx を呼びます -- --edr が \`gmx energy\` を、
           --trajectory が \`gmx trjconv -pbc nojump\` を走らせます。
           PATH に出すか GMX=/path/to/gmx を指定してください。
MSG
    exit 0
fi

# --edr は .xvg を **.edr の隣に**書きます (出力先は選べません)。input/ を
# 汚さないように、.edr を run_out/ へ写してからそちらを渡します。
cp gro_top_mode/input/md.edr "$OUT/md.edr"

# --trajectory に .xtc を渡すには MDAnalysis が要ります (abmptools の依存では
# ありません)。無ければ gmx で multi-frame .gro に直してから渡します ——
# .gro の読み取りは pure-Python です。どちらの経路でも中身は同じで、
# -pbc nojump は gro2udf 側が既定で通します。
TRAJ=gro_top_mode/input/md.xtc
if ! "$PYTHON" -c "import MDAnalysis" > /dev/null 2>&1; then
    echo "    MDAnalysis が無いので .xtc を .gro に直します"
    TRAJ=$OUT/md_traj.gro
    echo 0 | "$GMX" trjconv -f gro_top_mode/input/md.xtc \
        -s gro_top_mode/input/output.gro -o "$TRAJ" > /dev/null 2>&1
fi

"$PYTHON" -m abmptools.gro2udf --from-top \
    gro_top_mode/input/test.top gro_top_mode/input/output.gro \
    --template gro_top_mode/input/template_empty.udf \
    --mdp gro_top_mode/input/md.mdp \
    --trajectory "$TRAJ" \
    --edr "$OUT/md.edr" \
    --gmx "$GMX" \
    --out "$OUT/md_full.udf"
echo "    -> $OUT/md_full.udf (11 frame + energy)"
echo
echo "[OK] 詳しくは README.md"
