# -*- coding: utf-8 -*-
"""Unit tests for ``abmptools.crystal.atom_distance``.

Synthetic PDB fixtures only — no external sample data needed.
"""
from __future__ import annotations

import os

import pytest

from abmptools.crystal.atom_distance import (
    NearestAtom,
    find_nearest_atoms,
)


def _write_synthetic_pdb(path: str) -> None:
    """3 residues, 2 atoms each. Centre = res 1; expected sort by distance."""
    text = (
        # Residue 1 (centre): two atoms at +/- 0.5 along x.
        "HETATM    1  C1  UNK     1      -0.500   0.000   0.000  1.00  0.00           C\n"
        "HETATM    2  C2  UNK     1       0.500   0.000   0.000  1.00  0.00           C\n"
        # Residue 2: two atoms at +2 along y, +3 along y.
        "HETATM    3  N1  UNK     2       0.000   2.000   0.000  1.00  0.00           N\n"
        "HETATM    4  N2  UNK     2       0.000   3.000   0.000  1.00  0.00           N\n"
        # Residue 3: two atoms at +5 along z, +6 along z.
        "HETATM    5  O1  UNK     3       0.000   0.000   5.000  1.00  0.00           O\n"
        "HETATM    6  O2  UNK     3       0.000   0.000   6.000  1.00  0.00           O\n"
        "END\n"
    )
    with open(path, "w") as f:
        f.write(text)


def test_find_nearest_atoms_basic(tmp_path):
    pdb = tmp_path / "synthetic.pdb"
    _write_synthetic_pdb(str(pdb))

    rows = find_nearest_atoms(str(pdb), center_res_seq=1, n_neighbors=3)
    assert len(rows) == 3
    # Centre is at (0,0,0) (centroid of -0.5 and +0.5 along x).
    # Closest non-self: N1 at distance 2.0 (res 2).
    assert rows[0].res_seq == 2
    assert rows[0].element == "N"
    assert rows[0].distance == pytest.approx(2.0, abs=1e-6)
    # 2nd: N2 at 3.0.
    assert rows[1].res_seq == 2
    assert rows[1].distance == pytest.approx(3.0, abs=1e-6)
    # 3rd: O1 at 5.0.
    assert rows[2].res_seq == 3
    assert rows[2].element == "O"
    assert rows[2].distance == pytest.approx(5.0, abs=1e-6)


def test_find_nearest_atoms_n_default(tmp_path):
    pdb = tmp_path / "s.pdb"
    _write_synthetic_pdb(str(pdb))
    rows = find_nearest_atoms(str(pdb), center_res_seq=1)  # default n=3
    assert len(rows) == 3


def test_find_nearest_atoms_include_self(tmp_path):
    pdb = tmp_path / "s.pdb"
    _write_synthetic_pdb(str(pdb))
    rows = find_nearest_atoms(
        str(pdb), center_res_seq=1, n_neighbors=2, exclude_self=False,
    )
    # With exclude_self=False, the two res 1 atoms are at 0.5 each from centroid.
    assert len(rows) == 2
    assert all(r.res_seq == 1 for r in rows)


def test_find_nearest_atoms_unknown_centre(tmp_path):
    pdb = tmp_path / "s.pdb"
    _write_synthetic_pdb(str(pdb))
    with pytest.raises(ValueError, match="no atoms with res_seq=99"):
        find_nearest_atoms(str(pdb), center_res_seq=99)


def test_find_nearest_atoms_bad_n(tmp_path):
    pdb = tmp_path / "s.pdb"
    _write_synthetic_pdb(str(pdb))
    with pytest.raises(ValueError, match="n_neighbors must be"):
        find_nearest_atoms(str(pdb), center_res_seq=1, n_neighbors=0)


def test_find_nearest_atoms_empty_pdb(tmp_path):
    pdb = tmp_path / "empty.pdb"
    pdb.write_text("END\n")
    with pytest.raises(ValueError, match="no parseable"):
        find_nearest_atoms(str(pdb), center_res_seq=1)
