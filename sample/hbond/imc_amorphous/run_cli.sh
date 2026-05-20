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
    --mol-name IMC \
    --colorize-mode both

echo ""
echo "Done."
echo "  - imc_hbond_action.bdf : Mol_Name preserved + Python action overlay (recommended)"
echo "  - imc_hbond_colored.bdf: Mol_Name renamed (v1.25 legacy, post-render only)"
echo "  - imc_hbond.bdf        : plain copy (J-OCTA pre-render only, no color)"
