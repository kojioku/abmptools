#!/usr/bin/env bash
# trajectory_thin_nojump.sh — DEPRECATED
#
# 機能は abmptools.trajectory (Python module、 Windows 互換) に移行済み。
# このスクリプトは互換性のため一時的に残置。
#
# 新しい呼び出し方:
#   python -m abmptools.trajectory thin_nojump --traj prod/prod.xtc --tpr prod/prod.tpr --skip 10
# 詳細: ./README.md または python -m abmptools.trajectory --help
#
# --- 以下、 旧スクリプト本体 (引数互換 wrapper として動作) ---
#
# aggregation 系 AA-MD の prod.xtc を「可視化・解析の基本セット」に変換する。
#   - `-pbc nojump`: peptide / cluster が PBC で box を跨ぐのを連続化
#   - `-skip N`:     frame 数を 1/N に間引き (デフォルト N=10、 100 ns 1000 frame → 100 frame = 1 ns stride)
#
# 設計方針 (sample/formulation 全体の共通 post-process):
#   - input  = <sample>/prod/prod.xtc + prod.tpr
#   - output = <sample>/prod/prod_nojump_skip<N>.xtc  (元 file は保持)
#   - 各 aggregation sample の README からこの script を相対参照で呼ぶ
#
# 使い方:
#   bash trajectory_thin_nojump.sh                                 # cwd の prod/ で実行 (default)
#   bash trajectory_thin_nojump.sh <traj.xtc> <tpr> <out> <skip> <group>
#
# 環境:
#   gmxcudaenv または GROMACS が PATH にあること (`micromamba activate gmxcudaenv`)
#
# Example (octreotide_l_aggregation_100ns/ で):
#   cd /tmp/oct_l_agg100
#   bash $ABMPTOOLS/sample/formulation/_postprocess/trajectory_thin_nojump.sh
#   # → prod/prod_nojump_skip10.xtc (~300 MB、 VMD で開ける)
#
# 参考:
#   - memory `feedback_bdf_frame_count.md`: BDF 化は ~100 frames を目安
#   - 100 ns / 1 ns stride = 100 frame は VMD アニメ + 解析の bandwidth に最適

set -euo pipefail

TRAJ="${1:-prod/prod.xtc}"
TPR="${2:-prod/prod.tpr}"
SKIP="${4:-10}"
OUT="${3:-${TRAJ%.xtc}_nojump_skip${SKIP}.xtc}"
GROUP="${5:-System}"

# 引数チェック
[ -f "$TRAJ" ] || { echo "ERROR: trajectory not found: $TRAJ" >&2; exit 1; }
[ -f "$TPR" ]  || { echo "ERROR: tpr not found: $TPR" >&2; exit 1; }
command -v gmx >/dev/null || {
    echo "ERROR: gmx not found. Activate gmxcudaenv first:" >&2
    echo "  micromamba activate gmxcudaenv" >&2
    exit 1
}

echo "=== trajectory_thin_nojump ==="
echo "  Input    : $TRAJ ($(du -h "$TRAJ" | cut -f1))"
echo "  TPR      : $TPR"
echo "  Output   : $OUT"
echo "  Skip     : $SKIP (every ${SKIP}-th frame kept)"
echo "  Group    : $GROUP"
echo ""

# `-pbc nojump` は first frame を基準に unwrap、 cluster が box を跨いでも連続的に追跡可能。
# `-skip N` で frame 数を 1/N に間引き (output_dt = input_dt × N)。
echo "$GROUP" | gmx trjconv \
    -f "$TRAJ" -s "$TPR" -o "$OUT" \
    -pbc nojump -skip "$SKIP"

echo ""
echo "=== Done ==="
echo "  Output: $OUT ($(du -h "$OUT" | cut -f1))"
echo ""
echo "VMD で開く:"
echo "  vmd $TPR $OUT"
