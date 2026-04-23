#!/bin/bash
# Sample: build a multi-component amorphous using a JSON config.
#
# mixture.json here defines pentane (200) + benzene (50) at density 0.8;
# edit it to change system composition or swap for any other JSON schema
# accepted by `build_amorphous.py --config`.
#
# Produces ./input ./build ./md under this directory.
set -e

cd "$(dirname "$0")"

python ../../../build_amorphous.py --config mixture.json -v

echo ""
echo "Output files:"
find . -type f -not -name "run_sample.sh" -not -name "mixture.json" | sort

echo ""
echo "Next: cd md && bash run_all.sh && bash wrap_pbc.sh"
