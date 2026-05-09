# -*- coding: utf-8 -*-
"""
abmptools.crystal.forcefield_check
-----------------------------------
Dependency presence check for the crystal-FMO pipeline.

Phase A scope:
    - Optional Python modules (``ase``, ``yaml``) for the upcoming Phase C
      ``--engine ase`` backend and YAML config loader
    - Optional external tools (``abinitmp_smp``, ``abinitmp_omp``,
      ``mkinp_openver1rev20.py``) for ``--run-local`` and HPC job-script
      rendering

All dependencies are *optional* in Phase A: the legacy CLI
(``python -m abmptools.readcif``, ``python -m abmptools.pdb2fmo``, ...)
runs with stdlib + numpy/pandas only.

Phase C will tighten ``required=True`` for ``ase`` / ``yaml`` once the
new ``pipeline`` subcommand is wired in.
"""
from __future__ import annotations

import importlib
import logging
import shutil
from typing import List, Optional

# Re-use ToolStatus from cg.peptide to avoid duplicating the dataclass.
from abmptools.cg.peptide.forcefield_check import ToolStatus

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# External tools
# ---------------------------------------------------------------------------

def check_external_tools(
    abinitmp_smp: str = "abinitmp_smp",
    abinitmp_omp: str = "abinitmp_omp",
    mkinp: str = "mkinp_openver1rev20.py",
) -> List[ToolStatus]:
    """Resolve optional external tools used by Phase C/D pipelines.

    All tools are optional in Phase A — the legacy CLI does not need them.
    """
    spec = [
        ("abinitmp_smp", abinitmp_smp,
         "ABINIT-MP (SMP) — used by --run-local and HPC jobscripts",
         False),
        ("abinitmp_omp", abinitmp_omp,
         "ABINIT-MP (OMP variant) — alternative HPC binary",
         False),
        ("mkinp", mkinp,
         "ABINIT-MP input preprocessor (rev20+) — HPC jobscripts",
         False),
    ]
    out: List[ToolStatus] = []
    for label, exe, purpose, required in spec:
        resolved = shutil.which(exe)
        out.append(ToolStatus(
            name=label, purpose=purpose, found=resolved is not None,
            path=resolved, required=required,
        ))
    return out


# ---------------------------------------------------------------------------
# Python modules
# ---------------------------------------------------------------------------

def check_python_modules() -> List[ToolStatus]:
    """Verify Python modules used by Phase C (``ase`` engine, YAML config).

    Phase A reports them as optional ([WARN] only). Phase C will mark
    ``ase`` as required when ``--engine ase`` is selected.
    """
    out: List[ToolStatus] = []
    for label, module_name, purpose, required in [
        ("ase", "ase",
         "Atomic Simulation Environment — Phase C `--engine ase` backend "
         "(optional in Phase A; legacy engine works without it)",
         False),
        ("pyyaml", "yaml",
         "YAML config loader (optional; JSON config also works)",
         False),
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

def report(config_path: Optional[str] = None) -> bool:
    """Print a status block for the crystal subpackage.

    Phase A: returns True unconditionally (no required deps yet) — the
    return value is reserved for Phase C, where ``required=True`` deps
    can fail the validate command.
    """
    print(
        "abmptools.crystal — Phase A skeleton "
        "(legacy CLI re-exported; ase/YAML wiring deferred to Phase C)"
    )

    if config_path:
        print(f"\nConfig file: {config_path} "
              "(Phase A: parsed by legacy tools individually)")

    tools = check_external_tools()
    print("\nExternal tools (all optional in Phase A):")
    for t in tools:
        if t.found:
            print(f"  [OK]    {t.name:<14} ({t.path})")
        else:
            print(f"  [warn]  {t.name:<14} not on PATH — {t.purpose}")

    py_modules = check_python_modules()
    print("\nPython modules (will be required from Phase C):")
    for m in py_modules:
        if m.found:
            print(f"  [OK]    {m.name:<8} -- {m.purpose}")
        else:
            print(f"  [warn]  {m.name:<8} not importable -- {m.purpose}")

    print(
        "\nPhase A status: legacy CLI ready "
        "(python -m abmptools.{readcif,pdb2fmo,ajf2config,pdbmodify,"
        "getifiepieda})."
    )
    print("Phase C will add: pipeline / expand / fragment / jobs / postproc / "
          "nearest subcommands.")
    return True


__all__ = [
    "check_external_tools",
    "check_python_modules",
    "report",
]
