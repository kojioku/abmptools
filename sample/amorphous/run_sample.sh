#!/bin/bash
# Sample: build amorphous pentane/benzene mixture
set -e

cd "$(dirname "$0")"

python ../../build_amorphous.py \
    --smiles "CCCCC" "c1ccccc1" \
    --name pentane benzene \
    --n_mol 200 50 \
    --density 0.8 \
    --temperature 300 \
    --seed 42 \
    --output_dir ./pentane_benzene \
    -v

echo ""
echo "Output files:"
find ./pentane_benzene -type f | sort
