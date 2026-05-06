# -*- coding: utf-8 -*-
"""
abmptools.genesis.mmgbsa.pdb_splitter
-------------------------------------
Split a protein-ligand complex PDB into ``receptor`` (everything except
the ligand residue) and ``ligand`` (the ligand residue alone) PDBs.

Mirrors POC ``sample/gbsa/legacy_scripts/1_devidepdb.py`` (originally
``1_devidepdb.py`` in ``sample/gbsa/gbsa-input/``):

- ``Bio.PDB.PDBIO`` with custom ``Select`` subclasses
- Output filenames: ``<basename>_receptor_<resno>.pdb`` and
  ``<basename>_ligand_<resno>.pdb`` (POC convention)
- Optional chain filter (case-sensitive 1-character ID)

Biopython is loaded lazily inside :func:`split_pdb` so the rest of the
package imports cleanly without ``biopython`` installed (``[mmgbsa]``
extras).
"""
from __future__ import annotations

import logging
import shutil
from pathlib import Path
from typing import Optional, Tuple

from .models import TargetSpec

logger = logging.getLogger(__name__)


def split_pdb(
    pdb_path: Path,
    ligand_resno: int,
    output_dir: Path,
    chain: Optional[str] = None,
    output_basename: Optional[str] = None,
) -> Tuple[Path, Path]:
    """Split a PDB into receptor + ligand files.

    Parameters
    ----------
    pdb_path
        Source PDB file.
    ligand_resno
        Residue number to extract as the ligand. All atoms with this
        residue id (across the optional chain filter) become the ligand
        PDB; everything else becomes the receptor PDB.
    output_dir
        Directory to write the two PDBs (created if missing). Also
        receives a copy of the original PDB.
    chain
        Optional 1-character chain ID. ``None`` accepts any chain.
    output_basename
        Stem for the output PDBs. Default is the source PDB stem.

    Returns
    -------
    Tuple[Path, Path]
        ``(receptor_pdb, ligand_pdb)`` absolute paths.

    Raises
    ------
    FileNotFoundError
        If *pdb_path* does not exist.
    ValueError
        If the ligand residue is not found in the structure (after
        applying the optional chain filter).
    ImportError
        If Biopython is not installed.
    """
    try:
        from Bio import PDB
    except ImportError as exc:  # pragma: no cover -- explicit error for users
        raise ImportError(
            "abmptools.genesis.mmgbsa.pdb_splitter requires Biopython. "
            "Install via: pip install abmptools[mmgbsa]"
        ) from exc

    pdb_path = Path(pdb_path).resolve()
    output_dir = Path(output_dir)
    if not pdb_path.is_file():
        raise FileNotFoundError(f"PDB file not found: {pdb_path}")
    output_dir.mkdir(parents=True, exist_ok=True)

    basename = output_basename if output_basename else pdb_path.stem

    # Copy the original PDB for downstream tools (e.g. `tleap` references).
    dest_orig = output_dir / pdb_path.name
    if dest_orig.resolve() != pdb_path.resolve():
        shutil.copy2(pdb_path, dest_orig)

    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("structure", str(pdb_path))

    # Verify the ligand residue exists.
    found = _residue_exists(structure, ligand_resno, chain)
    if not found:
        raise ValueError(
            f"Ligand residue {ligand_resno} "
            f"{'(chain ' + chain + ')' if chain else ''} "
            f"not found in {pdb_path.name}"
        )

    io = PDB.PDBIO()
    io.set_structure(structure)

    receptor_path = output_dir / f"{basename}_receptor_{ligand_resno}.pdb"
    ligand_path = output_dir / f"{basename}_ligand_{ligand_resno}.pdb"

    io.save(str(receptor_path), _ResidueExcluder(ligand_resno, chain))
    io.save(str(ligand_path), _ResidueSelector(ligand_resno, chain))
    logger.info(
        "Split %s -> receptor=%s ligand=%s",
        pdb_path.name, receptor_path.name, ligand_path.name,
    )
    return receptor_path, ligand_path


def split_target(
    target: TargetSpec,
    output_dir: Path,
    input_dir: Optional[Path] = None,
) -> Tuple[Path, Path]:
    """Split a :class:`TargetSpec` into receptor + ligand PDBs.

    Resolves :attr:`TargetSpec.pdb` against *input_dir* if it is not an
    absolute path or already exists relative to cwd. Per-target output
    goes to ``output_dir/<basename>/``.

    Returns ``(receptor_pdb, ligand_pdb)`` paths.
    """
    pdb_path = _resolve_pdb_path(target.pdb, input_dir)
    basename = pdb_path.stem
    target_dir = Path(output_dir) / basename
    return split_pdb(
        pdb_path=pdb_path,
        ligand_resno=target.ligand_resno,
        output_dir=target_dir,
        chain=target.chain,
        output_basename=basename,
    )


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _resolve_pdb_path(pdb: str, input_dir: Optional[Path]) -> Path:
    """Resolve a TargetSpec.pdb string to an absolute Path."""
    p = Path(pdb)
    if p.is_absolute():
        return p
    if p.exists():
        return p.resolve()
    if input_dir is not None:
        candidate = Path(input_dir) / pdb
        if candidate.exists():
            return candidate.resolve()
    raise FileNotFoundError(
        f"Cannot resolve PDB {pdb!r}; "
        f"tried abs / cwd / input_dir={input_dir!r}"
    )


def _residue_exists(structure, resno: int, chain: Optional[str]) -> bool:
    """Check whether *resno* (optionally on *chain*) exists in *structure*."""
    for model in structure:
        for ch in model:
            if chain is not None and ch.id != chain:
                continue
            for residue in ch:
                if residue.get_id()[1] == resno:
                    return True
    return False


# ---------------------------------------------------------------------------
# Bio.PDB.Select subclasses (callable per-atom filters)
# ---------------------------------------------------------------------------
#
# These need to be defined at module top level so the deferred Biopython
# import doesn't shadow them. Both are constructed and used immediately
# above; subclassing PDB.Select would require Biopython at import time,
# so we duck-type the interface instead (PDB.PDBIO duck-types Select via
# attribute access, so a plain object with the same method names suffices).

class _ResidueExcluder:
    """``PDB.Select``-compatible filter: accept everything *except* the
    given residue number (optionally restricted to a chain)."""

    def __init__(self, resno: int, chain: Optional[str]) -> None:
        self.resno = resno
        self.chain = chain

    def accept_model(self, model) -> int:
        return 1

    def accept_chain(self, chain) -> int:
        return 1

    def accept_residue(self, residue) -> int:
        if self.chain is not None and residue.get_parent().id != self.chain:
            return 1
        return 0 if residue.get_id()[1] == self.resno else 1

    def accept_atom(self, atom) -> int:
        return 1


class _ResidueSelector:
    """``PDB.Select``-compatible filter: accept *only* the given residue."""

    def __init__(self, resno: int, chain: Optional[str]) -> None:
        self.resno = resno
        self.chain = chain

    def accept_model(self, model) -> int:
        return 1

    def accept_chain(self, chain) -> int:
        if self.chain is not None and chain.id != self.chain:
            return 0
        return 1

    def accept_residue(self, residue) -> int:
        if self.chain is not None and residue.get_parent().id != self.chain:
            return 0
        return 1 if residue.get_id()[1] == self.resno else 0

    def accept_atom(self, atom) -> int:
        return 1


__all__ = [
    "split_pdb",
    "split_target",
]
