# -*- coding: utf-8 -*-
"""
abmptools.genesis.mmgbsa.gbsa_runner
------------------------------------
Run ``atdyn`` (via ``mpirun``) on the three rendered ``.inp`` files
(``complex.inp`` / ``ligand.inp`` / ``receptor.inp``) and capture
each run's stdout to ``<name>.log``.

Mirrors POC ``3_rungbsa.py`` (``mpirun -np 1 atdyn {name}.inp | tee
{name}.log``). The ``-np 1`` choice is intentional: GBSA single-point
on a small protein-ligand complex doesn't benefit from spatial
decomposition and ``atdyn`` MPI rank > 1 can produce decomposition
artefacts on tiny systems.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List

from ._subprocess import CommandError, run_command, write_text
from .inp_writer import SYSTEM_NAMES

logger = logging.getLogger(__name__)


@dataclass
class GBSARunResult:
    """Outputs of running atdyn on one target's three systems."""
    target_dir: Path
    inp_files: Dict[str, Path] = field(default_factory=dict)
    log_files: Dict[str, Path] = field(default_factory=dict)
    dcd_files: Dict[str, Path] = field(default_factory=dict)
    rst_files: Dict[str, Path] = field(default_factory=dict)


def run_gbsa_one_target(
    target_dir: Path,
    inps: Dict[str, str],
    atdyn_path: str = "atdyn",
    mpirun_path: str = "mpirun",
    mpi_processes: int = 1,
) -> GBSARunResult:
    """Write 3 ``.inp`` files into *target_dir* and run atdyn on each.

    Parameters
    ----------
    target_dir
        Per-target build directory (already populated with
        ``complex.{prmtop,inpcrd}`` etc. by Stage 2).
    inps
        Mapping from name (``"complex"``/``"ligand"``/``"receptor"``)
        to rendered ``.inp`` text. Use :func:`inp_writer.render_three_inps`.
    atdyn_path / mpirun_path / mpi_processes
        Tool paths and rank count.

    Returns
    -------
    GBSARunResult
    """
    target_dir = Path(target_dir)
    target_dir.mkdir(parents=True, exist_ok=True)

    result = GBSARunResult(target_dir=target_dir)
    missing = [n for n in SYSTEM_NAMES if n not in inps]
    if missing:
        raise ValueError(f"inps is missing keys: {missing}")

    for name in SYSTEM_NAMES:
        inp_path = target_dir / f"{name}.inp"
        log_path = target_dir / f"{name}.log"
        write_text(inp_path, inps[name])
        result.inp_files[name] = inp_path

        cmd: List[str] = [
            mpirun_path, "-np", str(mpi_processes),
            atdyn_path, str(inp_path),
        ]
        logger.info("Running atdyn (%s): %s", name, " ".join(cmd))
        run_result = run_command(cmd, cwd=target_dir, capture=True)
        # Mirror POC's "tee {name}.log": persist the captured stdout
        # so analyse can parse [STEP4] later.
        log_path.write_text(run_result.stdout or "")
        result.log_files[name] = log_path
        result.dcd_files[name] = target_dir / f"{name}.dcd"
        result.rst_files[name] = target_dir / f"{name}.rst"

        if not log_path.is_file() or log_path.stat().st_size == 0:
            raise CommandError(
                cmd=" ".join(cmd),
                returncode=1,
                stderr=(
                    f"atdyn ({name}) finished but {log_path} is empty. "
                    "Check stdout capture / mpirun environment."
                ),
            )
    return result


__all__ = [
    "GBSARunResult",
    "run_gbsa_one_target",
]
