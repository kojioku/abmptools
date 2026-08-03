#!/usr/bin/env bash
# Build a PVP (polyvinylpyrrolidone) 5-mer amorphous cell.
# 1) construct the oligomer SDF, 2) OpenFF Sage 2.1.0 + AM1-BCC + packmol,
# 3) run md/run_all.sh for the 5-stage GROMACS MD.
set -e
cd "$(dirname "$0")"
ENV_BIN="${ENV_BIN:-$(dirname "$(command -v python)")}"   # 環境変数で上書き可
AMBER_BIN="${AMBER_BIN:-$(dirname "$(command -v antechamber 2>/dev/null || echo /usr/bin/false)")}"
PY="$ENV_BIN/python"
export PATH="$AMBER_BIN:$ENV_BIN:$PATH"
( cd input && "$PY" build_pvp_oligomer.py )
$PY -m abmptools.amorphous \
    --mol input/pvp_5mer.sdf \
    --name PVP \
    --n_mol 20 \
    --density 1.2 \
    --temperature 300 \
    --T_high 500 \
    --seed 42 \
    --charge_method am1bcc \
    --packmol_path "$ENV_BIN/packmol" \
    --output_dir . \
    -v
echo "System built. Now: cd md && MDRUN_OPTS='-ntmpi 1 -ntomp 8' bash run_all.sh"
