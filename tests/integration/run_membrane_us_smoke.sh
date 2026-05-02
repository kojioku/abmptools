#!/usr/bin/env bash
# tests/integration/run_membrane_us_smoke.sh
# ------------------------------------------
# End-to-end smoke test for abmptools.membrane (AMBER backend).
#
#   - poly-Ala 5-mer  +  POPC 32/leaflet  +  TIP3P  +  Na+/Cl- 0.15 M
#   - 7 windows  (z = -1.5 .. +1.5 nm,  dz = 0.5 nm)
#   - 1 ns/window  (production size; smoke runs build only)
#
# What this verifies (no MD run):
#   1. MembraneUSBuilder.build() completes for the AMBER backend
#   2. All expected files (top, gro, ndx, mdps, run.sh) are produced
#   3. Each generated MDP passes `gmx grompp -maxwarn 1`
#
# Run manually:
#   cd <abmptools-repo>
#   bash tests/integration/run_membrane_us_smoke.sh
#
# Optional:
#   - export RUN_DIR=/path/to/output  to keep build artefacts (default: tmp)
#   - the abmptoolsenv conda env must be on PATH (or set ABMPTOOLSENV)

set -euo pipefail

ENV_DEFAULT="${HOME}/.local/share/mamba/envs/abmptoolsenv"
ENV="${ABMPTOOLSENV:-$ENV_DEFAULT}"
if [ ! -x "$ENV/bin/python3" ]; then
    echo "ERROR: abmptoolsenv not found at $ENV"
    echo "  override with ABMPTOOLSENV=/path/to/env"
    exit 1
fi

export PATH="$ENV/bin:$PATH"
export AMBERHOME="${AMBERHOME:-$ENV}"

RUN_DIR="${RUN_DIR:-$(mktemp -d -t membrane_smoke_XXXXXX)}"
echo "=== Membrane US smoke test ==="
echo "  abmptoolsenv:  $ENV"
echo "  AMBERHOME:     $AMBERHOME"
echo "  RUN_DIR:       $RUN_DIR"

mkdir -p "$RUN_DIR"

# ---------- 1. Build phase ----------
"$ENV/bin/python3" - <<PYEOF
import logging
logging.basicConfig(level=logging.INFO, format='%(name)s | %(message)s')

from abmptools.membrane import (
    MembraneConfig, MembraneUSBuilder,
    LipidSpec, PeptideSpec, USProtocol,
    EquilibrationProtocol, PullingProtocol,
)

us = USProtocol(
    z_min_nm=-1.5, z_max_nm=+1.5, window_spacing_nm=0.5,
    window_nsteps=500_000, window_dt_ps=0.002,
    window_nstxtcout=1000, window_nstenergy=500,
)
eq = EquilibrationProtocol(
    em_steps=5000, em_tol=1000.0,
    nvt_nsteps=50_000, npt_nsteps=500_000,
    dt_ps=0.002, temperature_K=310.0,
)
pull = PullingProtocol(
    pull_rate_nm_per_ps=0.001, pull_force_constant=1000.0,
    nsteps=2_500_000, nstxout_compressed=500,
)
cfg = MembraneConfig(
    backend='amber',
    lipids=[LipidSpec(resname='POPC', n_per_leaflet=32)],
    peptide=PeptideSpec(name='aa5', sequence='AAAAA'),
    output_dir='$RUN_DIR',
    seed=42,
    box_xy_nm=6.0, water_thickness_nm=2.5, distance_to_lipid_nm=1.5,
    equilibration=eq, pulling=pull, umbrella=us,
)
result = MembraneUSBuilder(cfg).build()
print()
print('=== build() OK ===')
PYEOF

# ---------- 2. File existence check ----------
echo
echo "=== Checking output files ==="
required=(
    "$RUN_DIR/build/system.top"
    "$RUN_DIR/build/system.gro"
    "$RUN_DIR/build/system.ndx"
    "$RUN_DIR/equil/em.mdp"
    "$RUN_DIR/equil/nvt.mdp"
    "$RUN_DIR/equil/npt.mdp"
    "$RUN_DIR/pull/pull.mdp"
    "$RUN_DIR/run.sh"
    "$RUN_DIR/input/config.json"
)
for i in $(seq -f '%03g' 0 6); do
    required+=("$RUN_DIR/windows/win_$i/window.mdp")
done

ok=1
for f in "${required[@]}"; do
    if [ -f "$f" ]; then
        size=$(stat -c%s "$f")
        printf "  [OK]  %8d bytes  %s\n" "$size" "$f"
    else
        printf "  [FAIL] missing  %s\n" "$f"
        ok=0
    fi
done
[ $ok -eq 1 ] || { echo "ERROR: missing files"; exit 1; }

# ---------- 3. Validate each MDP via gmx grompp ----------
echo
echo "=== Validating MDPs via gmx grompp ==="
TOP="$RUN_DIR/build/system.top"
GRO="$RUN_DIR/build/system.gro"
NDX="$RUN_DIR/build/system.ndx"
TMP_TPR="$RUN_DIR/grompp_check.tpr"
MDPOUT="$RUN_DIR/grompp_check.mdp_out"

mdp_list=(
    "equil/em.mdp"
    "equil/nvt.mdp"
    "equil/npt.mdp"
    "pull/pull.mdp"
)
for i in $(seq -f '%03g' 0 6); do
    mdp_list+=("windows/win_$i/window.mdp")
done

for mdp in "${mdp_list[@]}"; do
    if "$ENV/bin/gmx" grompp -f "$RUN_DIR/$mdp" \
            -p "$TOP" -c "$GRO" -n "$NDX" \
            -o "$TMP_TPR" -po "$MDPOUT" \
            -maxwarn 1 >/dev/null 2>"$RUN_DIR/grompp_check.err"; then
        echo "  [OK]  $mdp"
    else
        echo "  [FAIL] $mdp"
        sed 's/^/    /' "$RUN_DIR/grompp_check.err"
        exit 1
    fi
done

echo
echo "=== smoke test PASSED ==="
echo "  output: $RUN_DIR"
echo "  to run actual MD:  cd '$RUN_DIR' && bash run.sh"
