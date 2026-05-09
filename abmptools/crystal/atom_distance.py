# -*- coding: utf-8 -*-
"""
abmptools.crystal.atom_distance
--------------------------------
Nearest-atom utilities for FMO post-processing.

Reimplements the relevant subset of the historical
``tips/pdbtips/readatomdistpdb.py`` as a clean Python API. Used by the
``postproc`` stage to annotate IFIE/PIEDA tables with the
``n_neighbors`` closest atoms of each peripheral fragment relative to
the central solute fragment.

The PDB parser here is intentionally minimal -- only the columns
needed for distance computation (``HETATM``/``ATOM`` element symbol,
residue id, x/y/z) are extracted. Full-fidelity PDB IO lives in
:mod:`abmptools.pdb_io`; this module is meant for read-only
nearest-neighbour queries on the for_abmp PDBs that the crystal
pipeline emits.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import List, Tuple


# ---------------------------------------------------------------------------
# Lightweight PDB record
# ---------------------------------------------------------------------------

@dataclass
class _AtomRecord:
    serial: int
    name: str
    res_name: str
    res_seq: int
    x: float
    y: float
    z: float
    element: str


def _parse_pdb(pdb_path: str) -> List[_AtomRecord]:
    """Parse ``HETATM`` / ``ATOM`` records into a flat list."""
    records: List[_AtomRecord] = []
    text = Path(pdb_path).read_text()
    for line in text.splitlines():
        if not (line.startswith("HETATM") or line.startswith("ATOM  ")):
            continue
        try:
            serial = int(line[6:11].strip())
            name = line[12:16].strip()
            res_name = line[17:20].strip()
            res_seq = int(line[22:26].strip())
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            element = line[76:78].strip() if len(line) >= 78 else name[0]
        except (ValueError, IndexError):
            continue
        records.append(_AtomRecord(serial, name, res_name, res_seq, x, y, z, element))
    return records


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

@dataclass
class NearestAtom:
    """One nearest-atom result row.

    Attributes
    ----------
    serial
        PDB atom serial number (1-indexed).
    res_seq
        Residue sequence number (used as a fragment id by the crystal
        pipeline since each fragment maps to one residue).
    element
        Atomic symbol from PDB columns 77-78 (or the first atom-name
        char if columns are missing).
    name
        PDB atom name (e.g. ``"O1"``, ``"C12"``).
    distance
        Distance to the query centre, in Å.
    """
    serial: int
    res_seq: int
    element: str
    name: str
    distance: float


def find_nearest_atoms(
    pdb_path: str,
    center_res_seq: int,
    *,
    n_neighbors: int = 3,
    exclude_self: bool = True,
) -> List[NearestAtom]:
    """Return the ``n_neighbors`` atoms closest to *center_res_seq*.

    The "centre" is the centroid of all atoms with the given residue
    sequence number (typical use: solute fragment id for the crystal
    pipeline). When ``exclude_self`` is True (default), atoms in the
    centre residue are skipped from the candidates.

    Parameters
    ----------
    pdb_path
        Path to a PDB file emitted by the crystal pipeline.
    center_res_seq
        Residue sequence number of the centre (typically ``1`` for
        the solute fragment).
    n_neighbors
        Number of nearest atoms to return.
    exclude_self
        Skip atoms belonging to ``center_res_seq``.

    Returns
    -------
    List[NearestAtom]
        Up to ``n_neighbors`` rows, sorted by ascending distance.

    Raises
    ------
    ValueError
        If no atoms with ``center_res_seq`` are found, or if the file
        contains no parseable records.
    """
    if n_neighbors <= 0:
        raise ValueError(f"n_neighbors must be > 0, got {n_neighbors}")

    records = _parse_pdb(pdb_path)
    if not records:
        raise ValueError(f"no parseable HETATM/ATOM records in {pdb_path}")

    centre_atoms = [r for r in records if r.res_seq == center_res_seq]
    if not centre_atoms:
        raise ValueError(
            f"no atoms with res_seq={center_res_seq} in {pdb_path}"
        )

    cx = sum(a.x for a in centre_atoms) / len(centre_atoms)
    cy = sum(a.y for a in centre_atoms) / len(centre_atoms)
    cz = sum(a.z for a in centre_atoms) / len(centre_atoms)

    candidates: List[Tuple[float, _AtomRecord]] = []
    for r in records:
        if exclude_self and r.res_seq == center_res_seq:
            continue
        d = math.sqrt((r.x - cx) ** 2 + (r.y - cy) ** 2 + (r.z - cz) ** 2)
        candidates.append((d, r))

    candidates.sort(key=lambda t: t[0])
    top = candidates[:n_neighbors]
    return [
        NearestAtom(
            serial=r.serial,
            res_seq=r.res_seq,
            element=r.element,
            name=r.name,
            distance=d,
        )
        for d, r in top
    ]


__all__ = ["NearestAtom", "find_nearest_atoms"]
