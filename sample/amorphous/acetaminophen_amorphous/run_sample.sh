#!/bin/bash
# Sample: build acetaminophen (paracetamol) amorphous box for an hbond
# generic-mode demo on a **COOH-free** drug.
#
# Acetaminophen (C8H9NO2) has a secondary amide (N-H donor + C=O acceptor) and a
# phenol -OH (donor + acceptor) but NO carboxyl group, so the imc mode
# (COOH 4-species) does not apply — use --classify-mode generic.
#
# 64 molecules packed at density 1.2 g/cm^3, OpenFF Sage 2.1.0, AM1-BCC charges,
# GROMACS 5-stage annealing MD. T_high=400 K instead of OpenFF default 600 K
# (per feedback_openff_thigh_for_small_organic: small organics gas out at 600 K).
set -e
cd "$(dirname "$0")"

ENV_BIN=/home/okuwaki/.local/share/mamba/envs/abmptoolsenv/bin
AMBER_BIN=/home/okuwaki/.local/share/mamba/envs/AmberTools25/bin
PY="$ENV_BIN/python"
SDF="input/pubchem_name_acetaminophen.sdf"   # 3D SDF fetched from PubChem

# OpenFF AM1-BCC needs antechamber/sqm from AmberTools (not bundled in
# abmptoolsenv), so expose AmberTools25 binaries via PATH.
export PATH="$AMBER_BIN:$ENV_BIN:$PATH"

# The SDF ships with this sample; to re-fetch from PubChem instead use
#   --pubchem_name acetaminophen
$PY -m abmptools.amorphous \
    --mol "$SDF" \
    --name APAP \
    --n_mol 64 \
    --density 1.2 \
    --temperature 300 \
    --T_high 400 \
    --seed 42 \
    --charge_method am1bcc \
    --packmol_path "$ENV_BIN/packmol" \
    --output_dir . \
    -v

echo ""
echo "Build complete. To run the MD + convert to UDF:"
echo "  export OMP_NUM_THREADS=4"
echo "  cd md && MDRUN_OPTS='-ntmpi 1 -ntomp 4' bash run_all.sh"
echo "  python wrap_pbc.py && python build_bdf.py"
echo ""
echo "Then the hbond generic-mode analysis (amide + phenol donors/acceptors):"
echo "  python -m abmptools.hbond md/05_npt_final.udf \\"
echo "    --out-prefix output/apap_hbond \\"
echo "    --classify-mode generic \\"
echo "    --donor-groups amide_donor,hydroxyl \\"
echo "    --acceptor-groups amide_O,hydroxyl_O \\"
echo "    --dt 5.0 --no-colorize --no-copy-uncolored"
