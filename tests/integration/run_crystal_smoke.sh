#!/usr/bin/env bash
# Integration smoke for `abmptools.crystal pipeline --run-local`.
#
# Skipped on CI (no abinitmp). Run manually with:
#     bash tests/integration/run_crystal_smoke.sh
#
# Requirements:
#   - abmptoolsenv with the `crystal` extras (ase, pyyaml)
#   - `abinitmp` on PATH or ABINIT_DIR/ABINIT_BIN env vars set
#
# Successful run leaves the FMO log in the working dir for inspection.
set -euo pipefail

ABINIT_BIN="${ABINIT_BIN:-abinitmp}"
ABINIT_DIR="${ABINIT_DIR:-}"

if ! command -v "${ABINIT_BIN}" >/dev/null 2>&1 \
    && [ ! -x "${ABINIT_DIR}/${ABINIT_BIN}" ]; then
    echo "SKIP: ${ABINIT_BIN} not found on PATH (set ABINIT_DIR/ABINIT_BIN to override)"
    exit 0
fi

WORK="$(mktemp -d -t crystal-smoke-XXXXXX)"
echo "[smoke] WORK=${WORK}"

# Private inputs (cif + UNK.ajf) are not bundled with the public abmptools
# repo; they live in abmptools-sample. Override ABMPTOOLS_SAMPLE_DIR to
# point at your checkout if it isn't at ~/repos/abmptools-sample.
SAMPLE_DIR="${ABMPTOOLS_SAMPLE_DIR:-${HOME}/repos/abmptools-sample}"
SAMPLE="${SAMPLE_DIR}/sample/csp7_ciftest/crystal_reference/R00001_layer3_hf_631g"

if [ ! -f "${SAMPLE}/XXXI-MMFF-R00001.cif" ] || [ ! -f "${SAMPLE}/UNK.ajf" ]; then
    echo "SKIP: csp7 inputs not found at ${SAMPLE}"
    echo "      Set ABMPTOOLS_SAMPLE_DIR to your abmptools-sample checkout."
    exit 0
fi

cp "${SAMPLE}/XXXI-MMFF-R00001.cif" "${WORK}/"
cp "${SAMPLE}/UNK.ajf" "${WORK}/"

# Smaller layer for smoke (layer=2 -> 8 cells / 4 mol per cell = 32 mol).
cat > "${WORK}/crystal.yaml" <<EOF
project_name: crystal_smoke
output_dir: ./out
inputs:
  - cif: XXXI-MMFF-R00001.cif
    layer: 2
    atoms_in_mol: [32]
cif_engine:
  engine: legacy
fragment:
  cutmode: around
  solutes: [0]
  criteria: 4.0
  molname: [UNK]
  template_ajf: ./UNK.ajf
fmo:
  method: HF
  basis_set: STO-3G
  memory: 2000
  cpfflag: false
  abinit_ver: rev23
  npro: 1
  is_xyz: true
hpc:
  scheduler: local
  abinit_dir: "${ABINIT_DIR}"
  binary_name: "${ABINIT_BIN}"
  nodes: 1
  proc_per_node: 1
  omp_threads: 1
  # Direct invocation (no mpirun) so the smoke runs without an MPI
  # runtime; abinitmp is the flat (non-OMP) build at this path.
  mpi_launcher: ""
  elapse: "00:30:00"
postproc:
  enable: false
EOF

cd "${WORK}"
python -m abmptools.crystal pipeline --config crystal.yaml --run-local

# Verify a log file was produced.
LOG="$(ls "${WORK}"/out/XXXI-MMFF-R00001/cifout/layer2/pdb/for_abmp/*.log | head -1)"
echo "[smoke] log: ${LOG}"
if grep -q "ABINIT-MP" "${LOG}"; then
    echo "[smoke] PASS"
    exit 0
else
    echo "[smoke] FAIL: log missing ABINIT-MP banner"
    head -20 "${LOG}"
    exit 1
fi
