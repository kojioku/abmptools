#!/bin/bash
# Sample: build ketoprofen amorphous using a PubChem 3D SDF.
#
# Uses an externally-prepared 3D conformer (H included, MMFF94-optimized),
# which is preferable when the SMILES -> OpenFF single-conformer path
# produces an unsuitable geometry.
#
# The SDF is bundled at input/ketoprofen_pubchem_cid3825.sdf. If you want
# to re-download it from PubChem (CID 3825), remove the file and run this
# script again.
set -e

cd "$(dirname "$0")"

SDF="input/ketoprofen_pubchem_cid3825.sdf"
if [ ! -f "$SDF" ]; then
    echo "Downloading PubChem 3D SDF (CID 3825)..."
    mkdir -p input
    curl -sL "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/3825/SDF?record_type=3d" \
         -o "$SDF"
fi

python ../../../build_amorphous.py \
    --mol "$SDF" \
    --name ketoprofen \
    --n_mol 50 \
    --density 0.8 \
    --temperature 300 \
    --seed 42 \
    --output_dir . \
    -v

echo ""
echo "Build complete. To run the MD and post-process for VMD:"
echo "  cd md && bash run_all.sh && bash wrap_pbc.sh"
