# -*- coding: utf-8 -*-
"""
abmptools.core.acpype
---------------------
Run ``acpype`` on a small-molecule PDB to generate AMBER GAFF/GAFF2
parameters + AM1-BCC charges, then locate the ``frcmod`` and ``mol2``
outputs.

Typical invocation::

    acpype -i ligand.pdb -c bcc -k maxcyc=0

acpype creates a ``<basename>.acpype/`` directory containing:

- ``<basename>_AC.frcmod``           -- additional GAFF parameters
- ``<basename>_bcc_gaff2.mol2``      -- atom types + AM1-BCC charges
                                        (or ``_bcc_gaff.mol2`` for ``gaff``)
- additional files (prmtop, inpcrd, top, gro, mdp, ...)

We pick the ``.frcmod`` and ``.mol2`` files matching the configured
atom type (gaff vs gaff2).

This lives in :mod:`abmptools.core` rather than inside one subsystem
because more than one caller needs it (``formulation`` here, plus
out-of-tree consumers). Keeping it neutral avoids a subsystem-to-
subsystem dependency.

License posture: acpype is GPL-3.0; abmptools shells out via subprocess
only (no source modification or linking), which is mere aggregation per
the GPL FAQ and compatible with abmptools' Apache-2.0.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional

from ..cg.peptide._subprocess import CommandError, run_command

logger = logging.getLogger(__name__)


@dataclass
class LigandParameterization:
    """Knobs for the ``acpype`` invocation (GAFF/GAFF2 + AM1-BCC).

    Default: ``-c bcc -k maxcyc=0``. ``maxcyc=0`` skips antechamber's
    pre-optimisation since the input PDB usually has reasonable geometry
    already.
    """
    charge_method: str = "bcc"          # "bcc" | "gas" | "user"
    net_charge: Optional[int] = None    # None: let acpype detect
    extra_keys: str = "maxcyc=0"
    atom_type: str = "gaff2"            # "gaff" | "gaff2"
    skip_if_cached: bool = True

    def __post_init__(self) -> None:
        if self.charge_method not in ("bcc", "gas", "user"):
            raise ValueError(
                f"LigandParameterization.charge_method must be "
                f"bcc/gas/user, got {self.charge_method!r}"
            )
        if self.atom_type not in ("gaff", "gaff2"):
            raise ValueError(
                f"LigandParameterization.atom_type must be gaff or gaff2, "
                f"got {self.atom_type!r}"
            )


@dataclass
class AcpypeResult:
    """Locations of acpype's GAFF outputs.

    Attributes
    ----------
    acpype_dir
        The ``<basename>.acpype/`` directory.
    frcmod
        Path to ``<basename>_AC.frcmod`` (additional parameters).
    mol2
        Path to ``<basename>_bcc_gaff2.mol2`` (or ``_gaff.mol2`` for
        ``atom_type=gaff``). Used by tleap as ``loadmol2``.
    """
    acpype_dir: Path
    frcmod: Path
    mol2: Path


def run_acpype(
    ligand_pdb: Path,
    workdir: Path,
    config: LigandParameterization,
    acpype_path: str = "acpype",
) -> AcpypeResult:
    """Run acpype on *ligand_pdb* inside *workdir* and locate outputs.

    Parameters
    ----------
    ligand_pdb
        Source ligand PDB (one molecule, HETATM records).
    workdir
        Directory where acpype runs (cwd for the subprocess). The
        ``<basename>.acpype/`` directory lands here.
    config
        :class:`LigandParameterization` with charge method etc.
    acpype_path
        Executable name or path; default ``"acpype"``.

    Returns
    -------
    AcpypeResult
        Locations of frcmod + mol2 files.
    """
    ligand_pdb = Path(ligand_pdb)
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)

    basename = ligand_pdb.stem
    acpype_dir = workdir / f"{basename}.acpype"

    # Cache: if acpype already ran here and produced frcmod + mol2,
    # skip the (slow) AM1-BCC charge calculation.
    if config.skip_if_cached and acpype_dir.is_dir():
        try:
            cached = _locate_outputs(acpype_dir, config.atom_type)
            logger.info(
                "acpype outputs already present in %s; skipping rerun.",
                acpype_dir,
            )
            return cached
        except FileNotFoundError:
            logger.info(
                "acpype dir %s exists but is incomplete; rerunning.",
                acpype_dir,
            )

    cmd: List[str] = [acpype_path, "-i", str(ligand_pdb)]
    cmd += ["-c", config.charge_method]
    if config.net_charge is not None:
        cmd += ["-n", str(config.net_charge)]
    if config.atom_type == "gaff":
        cmd += ["-a", "gaff"]
    if config.extra_keys:
        cmd += ["-k", config.extra_keys]
    logger.info("Running acpype: %s", " ".join(cmd))
    run_command(cmd, cwd=workdir, capture=True)

    if not acpype_dir.is_dir():
        raise CommandError(
            cmd=" ".join(cmd),
            returncode=1,
            stderr=(
                f"acpype completed but {acpype_dir} is missing. "
                "Check ligand PDB connectivity / heteroatom records."
            ),
        )
    return _locate_outputs(acpype_dir, config.atom_type)


def _locate_outputs(acpype_dir: Path, atom_type: str) -> AcpypeResult:
    """Locate ``.frcmod`` + ``.mol2`` from a populated acpype dir."""
    frcmods = sorted(acpype_dir.glob("*.frcmod"))
    if not frcmods:
        raise FileNotFoundError(
            f"No .frcmod found in {acpype_dir}; acpype may have failed silently."
        )
    # Convention: pick the first match (typically <basename>_AC.frcmod).
    frcmod = frcmods[0]

    # mol2 with charges -- pattern depends on atom type.
    mol2_pattern = (
        "*_bcc_gaff2.mol2" if atom_type == "gaff2" else "*_bcc_gaff.mol2"
    )
    mol2s = sorted(acpype_dir.glob(mol2_pattern))
    if not mol2s:
        # Some acpype versions use ``_NEW.mol2`` for non-bcc; fall back.
        mol2s = sorted(acpype_dir.glob("*.mol2"))
    if not mol2s:
        raise FileNotFoundError(
            f"No .mol2 (pattern {mol2_pattern}) found in {acpype_dir}."
        )
    mol2 = mol2s[0]

    return AcpypeResult(acpype_dir=acpype_dir, frcmod=frcmod, mol2=mol2)


__all__ = [
    "LigandParameterization",
    "AcpypeResult",
    "run_acpype",
]
