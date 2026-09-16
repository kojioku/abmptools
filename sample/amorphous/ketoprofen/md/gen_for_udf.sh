#!/bin/bash
# Export what OCTA viewer (GOURMET) / gro2udf need from a finished GROMACS run:
#
#   <stage>_energy.xvg   gmx energy, every term
#   <stage>_nojump.gro   gmx trjconv -pbc nojump (molecules stay continuous
#                        across the boundary; wrap_pbc uses -pbc mol instead,
#                        for VMD's compact cell)
#
# This is a thin wrapper. The work is in abmptools.trajectory, which is not
# tied to the amorphous protocol -- the stage is detected from the directory,
# so a Tg run's output works the same way. The Windows-compatible equivalent
# is gen_for_udf.py (same module, no bash needed).
#
#   bash gen_for_udf.sh                 # detect the stage here
#   bash gen_for_udf.sh prod            # a named stage, e.g. after Tg
#   python -m abmptools.trajectory gen_for_udf --help   # every option
set -e

if [ -n "$1" ]; then
    set -- --stage "$1"
fi
exec "${PYTHON:-python}" -m abmptools.trajectory gen_for_udf "$@"
