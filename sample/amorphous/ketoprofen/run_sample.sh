#!/bin/bash
# Sample: build ketoprofen amorphous from SMILES (single-component, 50 mols).
#
# Produces ./input ./build ./md under this directory.
# For an external PubChem 3D SDF variant, see ../ketoprofen_pubchem/.
set -e

cd "$(dirname "$0")"

python ../../../build_amorphous.py \
    --smiles "OC(=O)C(C)c1cccc(C(=O)c2ccccc2)c1" \
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
