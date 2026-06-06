"""Trajectory post-process helpers (cross-platform Python wrapper around gmx trjconv).

これまで sample/ や builder で生成していた bash script
(`trajectory_thin_nojump.sh`、 amorphous の `wrap_pbc.sh`、 `gen_for_udf.sh`)
を、Windows でも動く Python API に統一する。

Examples
--------
>>> from abmptools.trajectory import thin_and_nojump
>>> thin_and_nojump(
...     trajectory="prod/prod.xtc",
...     tpr="prod/prod.tpr",
...     skip=10,
... )
PosixPath('prod/prod_nojump_skip10.xtc')

CLI:

    python -m abmptools.trajectory thin_nojump --traj prod/prod.xtc \\
        --tpr prod/prod.tpr --skip 10
"""

from .postprocess import (
    GmxError,
    nojump,
    run_trjconv,
    thin,
    thin_and_nojump,
    wrap_pbc,
)

__all__ = [
    "GmxError",
    "nojump",
    "run_trjconv",
    "thin",
    "thin_and_nojump",
    "wrap_pbc",
]
