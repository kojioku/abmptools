#!/usr/bin/env bash
# tests/integration/run_membrane_us_charmm_smoke.sh
# -------------------------------------------------
# End-to-end smoke test for abmptools.membrane (CHARMM36 backend).
#
#   - poly-Ala 5-mer  +  POPC 32/leaflet  +  TIP3  +  SOD/CLA 0.15 M
#   - 7 windows  (z = -1.5 .. +1.5 nm,  dz = 0.5 nm)
#   - 1 ns/window  (production size; smoke runs build only)
#
# Skips with exit 0 (CI-friendly) when CHARMM36_FF_DIR is unset, since the
# Klauda lab GROMACS port is not bundled (download instructions:
# docs/membrane.md "CHARMM36 GROMACS port の取得").
#
# What this verifies (no MD run):
#   1. MembraneUSBuilder.build() completes for backend="charmm36"
#   2. All expected files (top, gro, ndx, mdps, run.sh) are produced
#   3. Each generated MDP passes `gmx grompp -maxwarn 1`
#
# Run manually:
#   export CHARMM36_FF_DIR=/path/to/charmm36-jul2022.ff
#   bash tests/integration/run_membrane_us_charmm_smoke.sh
#
# Optional:
#   - export RUN_DIR=/path/to/output  to keep build artefacts
#   - export ABMPTOOLSENV=/path/to/env  if not at default location

set -euo pipefail

if [ -z "${CHARMM36_FF_DIR:-}" ]; then
    echo "SKIP: CHARMM36_FF_DIR is unset; CHARMM36 backend smoke test skipped."
    echo "      To run, download the Klauda lab GROMACS port and:"
    echo "        export CHARMM36_FF_DIR=/path/to/charmm36-jul2022.ff"
    echo "      See docs/membrane.md (CHARMM36 GROMACS port の取得)."
    exit 0
fi

if [ ! -d "$CHARMM36_FF_DIR" ]; then
    echo "ERROR: CHARMM36_FF_DIR does not exist: $CHARMM36_FF_DIR"
    exit 1
fi

ENV_DEFAULT="${HOME}/.local/share/mamba/envs/abmptoolsenv"
ENV="${ABMPTOOLSENV:-$ENV_DEFAULT}"
if [ ! -x "$ENV/bin/python3" ]; then
    echo "ERROR: abmptoolsenv not found at $ENV"
    echo "  override with ABMPTOOLSENV=/path/to/env"
    exit 1
fi

export PATH="$ENV/bin:$PATH"
export AMBERHOME="${AMBERHOME:-$ENV}"

RUN_DIR="${RUN_DIR:-$(mktemp -d -t membrane_charmm_smoke_XXXXXX)}"
echo "=== Membrane US CHARMM36 smoke test ==="
echo "  abmptoolsenv:    $ENV"
echo "  AMBERHOME:       $AMBERHOME"
echo "  CHARMM36_FF_DIR: $CHARMM36_FF_DIR"
echo "  RUN_DIR:         $RUN_DIR"

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
    backend='charmm36',
    charmm_ff_dir='$CHARMM36_FF_DIR',
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
echo "=== CHARMM36 smoke test PASSED ==="
echo "  output: $RUN_DIR"
echo "  to run actual MD:  cd '$RUN_DIR' && bash run.sh"
