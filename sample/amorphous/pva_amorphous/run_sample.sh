#!/bin/bash
# Sample: build PVA (poly vinyl alcohol) amorphous box for hbond generic-mode demo.
#
# 10-mer atactic PVA, n_mol=30 packed at density 1.2 g/cm^3 (close to bulk PVA ~1.25),
# OpenFF Sage 2.1.0, GROMACS 5-stage annealing MD.
# T_high=400 K instead of OpenFF default 600 K (per
# feedback_openff_thigh_for_small_organic memory: small organics gas out at 600 K).
set -e
cd "$(dirname "$0")"

ENV_BIN=/home/okuwaki/.local/share/mamba/envs/abmptoolsenv/bin
AMBER_BIN=/home/okuwaki/.local/share/mamba/envs/AmberTools25/bin
PY="$ENV_BIN/python"
PACKMOL="$ENV_BIN/packmol"
SDF="input/pva_10mer.sdf"

# OpenFF AM1-BCC requires antechamber/sqm from AmberTools; abmptoolsenv does
# not bundle it, so we expose AmberTools25's binaries via PATH.
export PATH="$AMBER_BIN:$PATH"

$PY ../../../build_amorphous.py \
    --mol "$SDF" \
    --name PVA10 \
    --n_mol 30 \
    --density 1.2 \
    --temperature 300 \
    --T_high 400 \
    --seed 42 \
    --packmol_path "$PACKMOL" \
    --output_dir . \
    -v

echo ""
echo "Build complete. To run the MD:"
echo "  cd md && bash run_all.sh && bash wrap_pbc.sh"
echo ""
echo "After MD, run hbond generic-mode analysis on the produced BDF:"
echo "  python -m abmptools.hbond md/test_05_output_rec*.bdf \\"
echo "    --out-prefix output/pva_hbond \\"
echo "    --classify-mode generic \\"
echo "    --donor-groups hydroxyl \\"
echo "    --acceptor-groups hydroxyl_O \\"
echo "    --colorize-mode both"
