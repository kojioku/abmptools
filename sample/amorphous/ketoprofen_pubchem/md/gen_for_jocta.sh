#!/bin/bash
# Post-processing: export J-OCTA compatible inputs from MD outputs.
#
#   gmx energy             : dump every energy term (1..N) to <stage>_energy.xvg
#   gmx trjconv -pbc nojump: keep molecules continuous across PBC for J-OCTA
#                            (in contrast to wrap_pbc.sh, which uses -pbc mol
#                            for VMD-compatible compact unit-cell rendering)
#
# Run this after run_all.sh finishes.
set -e

STAGE="05_npt_final"
N_ENERGY_TERMS=50

# 1. All energy terms -> <stage>_energy.xvg
if [ -f "${STAGE}.edr" ]; then
    echo "Exporting energy terms to ${STAGE}_energy.xvg ..."
    seq "${N_ENERGY_TERMS}" | gmx energy -f "${STAGE}.edr" -o "${STAGE}_energy.xvg"
fi

# 2. Trajectory with -pbc nojump -> <stage>_nojump.gro
#    Prefer .trr (positions + velocities) when available, fall back to .xtc.
INPUT=""
if [ -f "${STAGE}.trr" ]; then
    INPUT="${STAGE}.trr"
elif [ -f "${STAGE}.xtc" ]; then
    INPUT="${STAGE}.xtc"
fi
if [ -n "${INPUT}" ] && [ -f "${STAGE}.tpr" ]; then
    echo "Exporting nojump trajectory to ${STAGE}_nojump.gro (from ${INPUT}) ..."
    echo 0 | gmx trjconv -f "${INPUT}" -s "${STAGE}.tpr" -pbc nojump -o "${STAGE}_nojump.gro" -n "../build/system.ndx"
fi

echo ""
echo "J-OCTA export complete:"
echo "  ${STAGE}_energy.xvg   (gmx energy)"
echo "  ${STAGE}_nojump.gro   (gmx trjconv -pbc nojump)"
