# -*- coding: utf-8 -*-
"""
abmptools.cg.membrane.forcefield_check
---------------------------------------
External tool & Martini 3 force field file checks for the membrane builder.

Extends :mod:`abmptools.cg.peptide.forcefield_check` with:

- ``insane`` as a required tool (in addition to martinize2 / gmx; tleap stays
  optional via the cg.peptide sub-call).
- ``martini_v3.0.0_phospholipids_v1.itp`` as a required ITP (in addition to
  the 3 from cg.peptide).

The package does **not bundle** any cgmartini distribution file (license
未明記); users download ``martini_v300.zip`` from cgmartini.nl and unzip the 4
required ITPs into their ``ff_dir``.
"""
from __future__ import annotations

import logging
import shutil
from pathlib import Path
from typing import List, Optional

from abmptools.cg.peptide.forcefield_check import (
    MARTINI_CITATION,
    MARTINI_DOWNLOAD_URL,
    FileStatus,
    ToolStatus,
    _status,
)
from abmptools.cg.peptide.forcefield_check import (
    REQUIRED_MARTINI_FILES as PEPTIDE_REQUIRED_MARTINI_FILES,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Required force field files (4 ITPs for membrane v1)
# ---------------------------------------------------------------------------

#: Martini 3 phospholipid topology (POPC, DOPC, POPE, … all in one file).
#: Distributed by cgmartini.nl as part of ``martini_v300.zip``.
PHOSPHOLIPIDS_ITP = "martini_v3.0.0_phospholipids_v1.itp"

#: Required ITP files for the membrane build = cg.peptide's 3 + phospholipids.
REQUIRED_MARTINI_FILES: List[str] = list(PEPTIDE_REQUIRED_MARTINI_FILES) + [
    PHOSPHOLIPIDS_ITP,
]


# ---------------------------------------------------------------------------
# External tools
# ---------------------------------------------------------------------------

def check_external_tools(
    insane: str = "insane",
    martinize2: str = "martinize2",
    gmx: str = "gmx",
    tleap: str = "tleap",
) -> List[ToolStatus]:
    """Check availability of external tools needed for the membrane build.

    ``insane`` (GPL-2.0) is required because we rely on it to assemble the
    Martini bilayer. ``martinize2`` is required for CG mapping (sub-called
    via cg.peptide). ``gmx`` is required. ``tleap`` is optional (cg.peptide
    falls back to extended-backbone when absent).
    """
    spec = [
        ("insane", insane,
         "Martini bilayer assembly (GPL-2.0)",
         True),
        ("martinize2", martinize2,
         "CG mapping via cg.peptide sub-call (vermouth-martinize, Apache-2.0)",
         True),
        ("gmx", gmx,
         "GROMACS -- grompp / genion / mdrun / wham",
         True),
        ("tleap", tleap,
         "AmberTools -- atomistic peptide PDB "
         "(推奨; 不在時は extended backbone fallback)",
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
# Martini 3 force field files
# ---------------------------------------------------------------------------

def check_martini_files(itp_dir: str) -> List[FileStatus]:
    """Check required Martini 3 ITPs in *itp_dir*.

    All 4 ITPs are REQUIRED (no optional files for the membrane build --
    water box auto-generation from cg.peptide is irrelevant here because
    insane handles solvation).
    """
    out: List[FileStatus] = []
    base = Path(itp_dir) if itp_dir else None
    for fname in REQUIRED_MARTINI_FILES:
        out.append(_status(base, fname, optional=False))
    return out


# ---------------------------------------------------------------------------
# Human-readable report
# ---------------------------------------------------------------------------

def report(
    tools: List[ToolStatus],
    files: List[FileStatus],
    itp_dir: str,
) -> bool:
    """Print human-readable status. Returns True if buildable as-is."""
    print("External tools:")
    required_ok = True
    for t in tools:
        if t.found:
            print(f"  [OK]    {t.name:<11} ({t.path})")
        elif t.required:
            required_ok = False
            print(f"  [MISS]  {t.name:<11} (required) -- {t.purpose}")
        else:
            print(f"  [WARN]  {t.name:<11} -- {t.purpose}")

    print(f"\nMartini 3 force field files (in {itp_dir or '<unset>'}):")
    required_files_ok = True
    for f in files:
        if f.found:
            print(f"  [OK]      {f.name}")
        else:
            required_files_ok = False
            print(f"  [MISSING] {f.name}")

    if not required_files_ok or not itp_dir:
        print(f"\nDownload required ITP files from:\n  {MARTINI_DOWNLOAD_URL}")
        print(
            "  (extract martini_v3.0.0.itp / _solvents_v1 / _ions_v1 / "
            "_phospholipids_v1.itp from martini_v300.zip)"
        )
        print(f"Please cite: {MARTINI_CITATION}")

    if not any(t.found for t in tools if t.name == "insane"):
        print(
            "\nInstall insane via:\n  pip install insane  (GPL-2.0, "
            "https://github.com/Tsjerk/Insane)"
        )

    return required_ok and required_files_ok
