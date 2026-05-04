# -*- coding: utf-8 -*-
"""
abmptools.cg.peptide.forcefield_check
--------------------------------------
External tool & Martini 3 force field file checks.

Two responsibilities:

1. ``check_external_tools()``
       Verify ``martinize2`` / ``gmx`` / ``tleap`` are reachable via PATH
       (or via configured paths in :class:`PeptideBuildConfig`).

2. ``check_martini_files(itp_dir)``
       Verify required Martini 3 ``.itp`` / ``.gro`` files exist in a
       user-supplied directory. **The package does not bundle these
       files** because cgmartini.nl does not state explicit
       redistribution terms; users must download them separately.
"""
from __future__ import annotations

import logging
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# External tools
# ---------------------------------------------------------------------------

@dataclass
class ToolStatus:
    name: str
    purpose: str
    found: bool
    path: Optional[str] = None
    required: bool = True


def check_external_tools(
    martinize2: str = "martinize2",
    gmx: str = "gmx",
    tleap: str = "tleap",
) -> List[ToolStatus]:
    """Check availability of external tools.

    ``tleap`` is treated as **optional** -- its absence triggers the
    extended-backbone fallback in :mod:`peptide_atomistic`.
    """
    spec = [
        ("martinize2", martinize2,
         "CG mapping (vermouth-martinize, Apache-2.0)", True),
        ("gmx", gmx,
         "GROMACS -- solvate / genion / grompp / make_ndx", True),
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

#: ITP files expected in ``martini_itp_dir`` for a peptide + ion build.
#: Distributed by cgmartini.nl as part of ``martini_v300.zip``.
REQUIRED_MARTINI_FILES = [
    "martini_v3.0.0.itp",
    "martini_v3.0.0_solvents_v1.itp",
    "martini_v3.0.0_ions_v1.itp",
]

#: Optional Martini 3 water box for ``gmx solvate -cs``. **cgmartini.nl does
#: not distribute this file directly** (only ITP is in martini_v300.zip);
#: users typically generate one via ``gmx insert-molecules`` or take it from
#: a Martini 3 tutorial archive. Required only when
#: :attr:`PeptideBuildConfig.solvent_enabled` is True.
OPTIONAL_MARTINI_FILES = [
    "martini_v3.0.0_water.gro",
]

MARTINI_DOWNLOAD_URL = (
    "https://cgmartini.nl/docs/downloads/force-field-parameters/martini3/"
)

MARTINI_CITATION = (
    "Souza et al. 2021, Nat. Methods 18:382-388 "
    "(doi:10.1038/s41592-021-01098-3)"
)


@dataclass
class FileStatus:
    name: str
    found: bool
    path: Optional[str] = None
    optional: bool = False


def check_martini_files(itp_dir: str) -> List[FileStatus]:
    """Check required + optional Martini 3 files in *itp_dir*.

    An empty *itp_dir* string yields all-missing statuses (helps the
    ``validate`` CLI subcommand display useful guidance).
    """
    out: List[FileStatus] = []
    base = Path(itp_dir) if itp_dir else None
    for fname in REQUIRED_MARTINI_FILES:
        out.append(_status(base, fname, optional=False))
    for fname in OPTIONAL_MARTINI_FILES:
        out.append(_status(base, fname, optional=True))
    return out


def _status(
    base: Optional[Path], fname: str, *, optional: bool,
) -> FileStatus:
    if base is None:
        return FileStatus(name=fname, found=False, optional=optional)
    p = base / fname
    return FileStatus(
        name=fname,
        found=p.exists(),
        path=str(p) if p.exists() else None,
        optional=optional,
    )


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
            tag = "[OK]    " if not f.optional else "[opt-OK]"
            print(f"  {tag} {f.name}")
        elif f.optional:
            print(
                f"  [opt]    {f.name}  (only needed when "
                "solvent_enabled=True; gmx insert-molecules で自作可)"
            )
        else:
            required_files_ok = False
            print(f"  [MISSING] {f.name}")

    if not required_files_ok or not itp_dir:
        print(f"\nDownload required ITP files from:\n  {MARTINI_DOWNLOAD_URL}")
        print(f"Please cite: {MARTINI_CITATION}")

    return required_ok and required_files_ok
