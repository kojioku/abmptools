#!/usr/bin/env bash
# Build an indomethacin (IMC) amorphous cell: OpenFF Sage 2.1.0 + AM1-BCC,
# packmol packing, then run md/run_all.sh for the 5-stage GROMACS MD.
set -e
cd "$(dirname "$0")"
ENV_BIN=/home/okuwaki/.local/share/mamba/envs/abmptoolsenv/bin
AMBER_BIN=/home/okuwaki/.local/share/mamba/envs/AmberTools25/bin   # AM1-BCC (sqm)
PY="$ENV_BIN/python"
export PATH="$AMBER_BIN:$ENV_BIN:$PATH"
$PY -m abmptools.amorphous \
    --pubchem_name indometacin \
    --name IMC \
    --n_mol 48 \
    --density 1.3 \
    --temperature 300 \
    --T_high 500 \
    --seed 42 \
    --charge_method am1bcc \
    --packmol_path "$ENV_BIN/packmol" \
    --output_dir . \
    -v
echo "System built. Now: cd md && MDRUN_OPTS='-ntmpi 1 -ntomp 8' bash run_all.sh"
