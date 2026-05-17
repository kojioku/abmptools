#!/usr/bin/env bash
# Sample CLI invocation for IMC amorphous H-bond analysis.

set -euo pipefail

# Resolve sample dir
HERE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
mkdir -p "$HERE/output"

BDF="$HOME/llm-project/SI/IMC_result450.0_out_rec900.bdf"
if [[ ! -f "$BDF" ]]; then
    echo "Error: BDF not found at $BDF" >&2
    exit 1
fi

python -m abmptools.hbond "$BDF" \
    --out-prefix "$HERE/output/imc_hbond" \
    --criteria luzar-chandler \
    --mol-name IMC

echo ""
echo "Done. Open '$HERE/output/imc_hbond_colored.bdf' in OCTA gourmet."
