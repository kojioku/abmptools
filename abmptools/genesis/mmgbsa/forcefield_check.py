# -*- coding: utf-8 -*-
"""
abmptools.genesis.mmgbsa.forcefield_check
-----------------------------------------
External tool resolution for the MM/GBSA pipeline.

Required tools (subprocess only, never bundled):

- ``atdyn``    -- GENESIS, GBSA single-point energy (LGPL-3.0+)
- ``tleap``    -- AmberTools, parameterize 3 systems
- ``acpype``   -- ligand GAFF + AM1-BCC (GPL-3.0)
- ``mpirun``   -- OpenMPI / MPICH

Optional Python modules (auto-skip in tests):

- ``Bio.PDB``    -- Biopython, PDB splitter
- ``matplotlib`` -- ΔG_bind bar plot

The :func:`report` function prints a human-readable status block;
returns ``True`` only if all REQUIRED tools resolve.
"""
from __future__ import annotations

import importlib
import logging
import shutil
from typing import List, Optional

# Re-use ToolStatus from cg.peptide to avoid duplicating the dataclass.
from abmptools.cg.peptide.forcefield_check import ToolStatus

from .models import MMGBSABuildConfig

logger = logging.getLogger(__name__)


GENESIS_REPO_URL = "https://github.com/genesis-release-r-ccs/genesis"
GENESIS_LICENSE = "LGPL-3.0-or-later"
ACPYPE_LICENSE = "GPL-3.0"


# ---------------------------------------------------------------------------
# External tool resolution
# ---------------------------------------------------------------------------

def check_external_tools(
    atdyn: str = "atdyn",
    tleap: str = "tleap",
    acpype: str = "acpype",
    mpirun: str = "mpirun",
) -> List[ToolStatus]:
    """Resolve each external CLI tool via :func:`shutil.which`."""
    spec = [
        ("atdyn", atdyn,
         f"GENESIS atdyn -- GBSA single-point energy ({GENESIS_LICENSE})",
         True),
        ("tleap", tleap,
         "AmberTools tleap -- AMBER parameterization (free academic+commercial)",
         True),
        ("acpype", acpype,
         f"acpype -- ligand GAFF/GAFF2 + AM1-BCC ({ACPYPE_LICENSE})",
         True),
        ("mpirun", mpirun,
         "OpenMPI / MPICH -- mpirun -np N atdyn",
         True),
    ]
    out: List[ToolStatus] = []
    for label, exe, purpose, required in spec:
        resolved = shutil.which(exe)
        out.append(ToolStatus(
            name=label, purpose=purpose, found=resolved is not None,
            path=resolved, required=required,
        ))
    return out


def check_python_modules() -> List[ToolStatus]:
    """Verify optional Python modules used by the pipeline.

    Returns
    -------
    List[ToolStatus]
        with ``found=True`` if importable; ``required`` is True for
        Biopython (Stage 1 requires it) and matplotlib (Stage 4 plot).
    """
    out: List[ToolStatus] = []
    for label, module_name, purpose, required in [
        ("Biopython", "Bio.PDB",
         "PDB splitter (Stage 1)", True),
        ("matplotlib", "matplotlib",
         "ΔG_bind bar plot (Stage 4)", True),
    ]:
        try:
            importlib.import_module(module_name)
            found = True
        except ImportError:
            found = False
        out.append(ToolStatus(
            name=label,
            purpose=purpose,
            found=found,
            path=None,
            required=required,
        ))
    return out


# ---------------------------------------------------------------------------
# Human-readable report
# ---------------------------------------------------------------------------

def report(cfg: MMGBSABuildConfig) -> bool:
    """Print a status block for :class:`MMGBSABuildConfig`.

    Returns ``True`` if buildable: all REQUIRED tools + Python modules
    resolve. Falsy return on any missing required dependency.
    """
    n_targets = len(cfg.targets) if cfg.targets else "(folder mode)"
    print(
        f"Configuration: OK (project={cfg.project_name}, "
        f"AMBER {cfg.force_field.leaprc_protein.split('.')[-1]} + "
        f"{cfg.force_field.leaprc_water.split('.')[-1]} + "
        f"{cfg.ligand.atom_type.upper()}, "
        f"{n_targets} targets, "
        f"implicit_solvent={cfg.energy.implicit_solvent})"
    )

    tools = check_external_tools(
        atdyn=cfg.atdyn_path,
        tleap=cfg.tleap_path,
        acpype=cfg.acpype_path,
        mpirun=cfg.mpirun_path,
    )
    print("\nExternal tools:")
    required_ok = True
    for t in tools:
        if t.found:
            print(f"  [OK]    {t.name:<8} ({t.path})")
        elif t.required:
            required_ok = False
            print(f"  [MISS]  {t.name:<8} (required) -- {t.purpose}")
        else:
            print(f"  [WARN]  {t.name:<8} -- {t.purpose}")

    py_modules = check_python_modules()
    print("\nPython modules ([mmgbsa] extras):")
    py_ok = True
    for m in py_modules:
        if m.found:
            print(f"  [OK]    {m.name:<11} -- {m.purpose}")
        elif m.required:
            py_ok = False
            print(f"  [MISS]  {m.name:<11} (required) -- {m.purpose}")
        else:
            print(f"  [WARN]  {m.name:<11} -- {m.purpose}")

    # Targets preview.
    print("\nTargets preview:")
    if cfg.targets:
        for i, t in enumerate(cfg.targets, start=1):
            ch = f"chain {t.chain}" if t.chain else "any chain"
            print(f"  T{i:02d}: {t.pdb}  ligand_resno={t.ligand_resno}  ({ch})")
    elif cfg.input_dir:
        print(f"  folder mode: {cfg.input_dir}/*.pdb (resolved at run time)")
    else:
        print("  (no targets configured)")

    # Force field summary.
    ff = cfg.force_field
    print("\nForce field:")
    print(f"  protein  = {ff.leaprc_protein}")
    print(f"  dna      = {ff.leaprc_dna}")
    print(f"  rna      = {ff.leaprc_rna}")
    print(f"  water    = {ff.leaprc_water}")
    print(f"  ligand   = {ff.leaprc_gaff2} + {ff.leaprc_gaff} "
          f"(acpype: {cfg.ligand.charge_method.upper()} charges, "
          f"atom_type={cfg.ligand.atom_type})")

    # Protocol summary.
    print("\nProtocol:")
    print(
        f"  energy   = {cfg.energy.implicit_solvent} ({cfg.energy.electrostatic}, "
        f"cutoff={cfg.energy.cutoffdist_A:.1f} Å)"
    )
    print(
        f"  minimize = {cfg.minimize.method} × {cfg.minimize.nsteps} step(s)"
    )
    print(
        f"  mpi      = {cfg.mpi_processes} rank(s) × atdyn (per system × "
        f"{len(cfg.targets) if cfg.targets else 'N'} targets × 3 systems)"
    )

    # Build hints if anything missing.
    if not required_ok:
        print(f"\nGENESIS build instructions: {GENESIS_REPO_URL}")
        print(
            "Install acpype via: pip install acpype "
            "(or conda install -c conda-forge acpype)"
        )
    if not py_ok:
        print("\nInstall Python deps: pip install abmptools[mmgbsa]")

    return required_ok and py_ok


__all__ = [
    "ACPYPE_LICENSE",
    "GENESIS_LICENSE",
    "GENESIS_REPO_URL",
    "check_external_tools",
    "check_python_modules",
    "report",
]
