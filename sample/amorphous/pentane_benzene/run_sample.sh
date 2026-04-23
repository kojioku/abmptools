#!/bin/bash
# Sample: build amorphous pentane/benzene mixture from SMILES.
#
# Produces ./input ./build ./md under this directory.
set -e

cd "$(dirname "$0")"

python ../../../build_amorphous.py \
    --smiles "CCCCC" "c1ccccc1" \
    --name pentane benzene \
    --n_mol 200 50 \
    --density 0.8 \
    --temperature 300 \
    --seed 42 \
    --output_dir . \
    -v

echo ""
echo "Output files:"
find . -type f -not -name "run_sample.sh" | sort

echo ""
echo "Next: cd md && bash run_all.sh && bash wrap_pbc.sh"
