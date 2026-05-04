# -*- coding: utf-8 -*-
"""
abmptools.cg.peptide.martinize_runner
--------------------------------------
Wrapper around vermouth-martinize (``martinize2`` CLI) for Martini 3
peptide CG mapping.

vermouth-martinize is distributed under Apache-2.0 (verified against
https://github.com/marrink-lab/vermouth-martinize). This wrapper only
invokes ``martinize2`` as an external subprocess -- vermouth source
code is **not bundled or modified** here.
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Set, Tuple

from ._subprocess import ensure_dir, run_command

logger = logging.getLogger(__name__)


def run_martinize2(
    input_pdb: Path,
    output_dir: Path,
    name: str,
    *,
    martinize2_path: str = "martinize2",
    elastic_network: bool = False,
    maxwarn: int = 100,
) -> Tuple[Path, Path]:
    """Run martinize2 on an atomistic peptide PDB to produce M3 CG topology.

    Parameters
    ----------
    input_pdb
        All-atom peptide PDB.
    output_dir
        Directory to write outputs into.
    name
        Molecule name; used for output filenames.
    martinize2_path
        martinize2 CLI command (default ``"martinize2"``).
    elastic_network
        Add ``-elastic`` (Martini 3 elastic network for stability).
    maxwarn
        Pass to ``-maxwarn`` (default 100; some sidechain-mapping warnings
        are expected with extended-backbone fallback input).

    Returns
    -------
    (itp_path, cg_pdb_path)
    """
    output_dir = ensure_dir(Path(output_dir)).resolve()
    input_pdb = Path(input_pdb).resolve()

    itp_path = output_dir / f"{name}.itp"
    cg_pdb_path = output_dir / f"{name}_cg.pdb"
    top_path = output_dir / f"{name}_martinize.top"

    n_res = _count_residues_pdb(input_pdb)

    cmd = [
        martinize2_path,
        "-f", str(input_pdb),
        "-o", str(top_path),
        "-x", str(cg_pdb_path),
        "-ff", "martini3001",
        "-ss", "C" * n_res,           # treat all residues as coil
        "-maxwarn", str(maxwarn),
    ]
    if elastic_network:
        cmd.append("-elastic")

    logger.info(
        "Running martinize2 for %s (martini3001, %d residues)",
        name, n_res,
    )
    run_command(cmd, cwd=output_dir)

    _rename_itp(output_dir, name, itp_path)

    if not cg_pdb_path.exists():
        raise RuntimeError(f"martinize2 did not produce {cg_pdb_path}")

    logger.info(
        "martinize2 complete: itp=%s, cg_pdb=%s", itp_path, cg_pdb_path,
    )
    return itp_path, cg_pdb_path


def _count_residues_pdb(pdb_path: Path) -> int:
    """Count residues by unique resid in a PDB file."""
    resids: Set[int] = set()
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                try:
                    resid = int(line[22:26].strip())
                    resids.add(resid)
                except ValueError:
                    pass
    return len(resids)


def _rename_itp(output_dir: Path, name: str, target: Path) -> None:
    """Locate martinize2's generated .itp and rename to ``<name>.itp``.

    martinize2 typically writes ``molecule_0.itp``. As a fallback we look
    for any other ``*.itp`` that isn't the target itself or the
    ``_martinize.top`` write-out.
    """
    candidates = list(output_dir.glob("molecule_*.itp"))
    if not candidates:
        candidates = [
            c for c in output_dir.glob("*.itp")
            if c.name != target.name
            and not c.name.endswith("_martinize.top")
        ]

    if not candidates:
        raise RuntimeError(f"No .itp produced by martinize2 in {output_dir}")

    if len(candidates) > 1:
        logger.warning(
            "Multiple .itp candidates after martinize2: %s; renaming first",
            [c.name for c in candidates],
        )

    src = candidates[0]
    if target.exists():
        target.unlink()  # idempotent re-runs
    src.rename(target)
    logger.debug("Renamed %s -> %s", src, target)
