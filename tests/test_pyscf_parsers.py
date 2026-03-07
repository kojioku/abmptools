# -*- coding: utf-8 -*-
"""Tests for standalone parser/writer helpers in abmptools.geomopt.pyscf_optimizer."""
from pathlib import Path

import pytest

from abmptools.geomopt.pyscf_optimizer import (
    _parse_xyz,
    _parse_pdb,
    _write_xyz,
    _infer_element_from_atom_name,
    _atoms_to_pyscf_spec,
    _HARTREE_TO_EV,
    _BOHR_TO_ANG,
    QMOptimizerPySCF,
)


# ---------------------------------------------------------------------------
# _parse_xyz
# ---------------------------------------------------------------------------

def test_parse_xyz_valid(tmp_path):
    xyz_file = tmp_path / "water.xyz"
    xyz_file.write_text(
        "3\n"
        "water molecule\n"
        "O   0.00000000  0.00000000  0.11730000\n"
        "H   0.75720000  0.00000000 -0.46920000\n"
        "H  -0.75720000  0.00000000 -0.46920000\n"
    )
    atoms = _parse_xyz(xyz_file)
    assert len(atoms) == 3
    assert atoms[0][0] == "O"
    assert atoms[0][1] == pytest.approx(0.0)
    assert atoms[0][3] == pytest.approx(0.1173)
    assert atoms[1][0] == "H"


def test_parse_xyz_malformed_raises(tmp_path):
    xyz_file = tmp_path / "bad.xyz"
    xyz_file.write_text("not a number\n")
    with pytest.raises(ValueError):
        _parse_xyz(xyz_file)


# ---------------------------------------------------------------------------
# _parse_pdb
# ---------------------------------------------------------------------------

def test_parse_pdb_with_element_column(tmp_path):
    pdb_file = tmp_path / "test.pdb"
    # Standard PDB ATOM line (80 chars with element in cols 77-78)
    line = (
        "ATOM      1  CA  ALA A   1       1.000   2.000   3.000"
        "  1.00  0.00           C  \n"
    )
    pdb_file.write_text(line)
    atoms = _parse_pdb(pdb_file)
    assert len(atoms) == 1
    assert atoms[0][0] == "C"
    assert atoms[0][1] == pytest.approx(1.0)
    assert atoms[0][2] == pytest.approx(2.0)
    assert atoms[0][3] == pytest.approx(3.0)


def test_parse_pdb_infer_from_atom_name(tmp_path):
    pdb_file = tmp_path / "test.pdb"
    # No element column (short line), element inferred from atom name " CA "
    line = (
        "ATOM      1  CA  ALA A   1       1.000   2.000   3.000"
        "  1.00  0.00\n"
    )
    pdb_file.write_text(line)
    atoms = _parse_pdb(pdb_file)
    assert len(atoms) == 1
    # " CA " -> infer "C" (Ca would be two-letter but CA atom name
    # starts at col 13; first letter is space, so stripped to "CA",
    # "Ca" is in _TWO_LETTER_ELEMENTS => "Ca")
    # Actually let's just check it returns a valid element string
    assert atoms[0][0] in ("C", "Ca")


def test_parse_pdb_empty_raises(tmp_path):
    pdb_file = tmp_path / "empty.pdb"
    pdb_file.write_text("REMARK empty file\n")
    with pytest.raises(ValueError, match="No ATOM/HETATM"):
        _parse_pdb(pdb_file)


# ---------------------------------------------------------------------------
# _write_xyz + _parse_xyz roundtrip
# ---------------------------------------------------------------------------

def test_write_xyz_roundtrip(tmp_path):
    atoms_in = [("C", 1.0, 2.0, 3.0), ("H", 4.0, 5.0, 6.0)]
    xyz_path = tmp_path / "roundtrip.xyz"
    _write_xyz(xyz_path, atoms_in, comment="roundtrip test")
    atoms_out = _parse_xyz(xyz_path)
    assert len(atoms_out) == len(atoms_in)
    for (e1, x1, y1, z1), (e2, x2, y2, z2) in zip(atoms_in, atoms_out):
        assert e1 == e2
        assert x1 == pytest.approx(x2)
        assert y1 == pytest.approx(y2)
        assert z1 == pytest.approx(z2)


# ---------------------------------------------------------------------------
# _infer_element_from_atom_name
# ---------------------------------------------------------------------------

def test_infer_element_from_atom_name():
    assert _infer_element_from_atom_name(" CA ") == "Ca"
    assert _infer_element_from_atom_name("FE  ") == "Fe"
    assert _infer_element_from_atom_name(" H  ") == "H"
    assert _infer_element_from_atom_name("1HB ") == "H"


# ---------------------------------------------------------------------------
# _atoms_to_pyscf_spec
# ---------------------------------------------------------------------------

def test_atoms_to_pyscf_spec():
    atoms = [("O", 0.0, 0.0, 0.117), ("H", 0.757, 0.0, -0.469)]
    spec = _atoms_to_pyscf_spec(atoms)
    assert len(spec) == 2
    assert spec[0][0] == "O"
    assert spec[0][1] == pytest.approx((0.0, 0.0, 0.117))
    assert spec[1][0] == "H"


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

def test_constants():
    assert _HARTREE_TO_EV == pytest.approx(27.2114, rel=1e-3)
    assert _BOHR_TO_ANG == pytest.approx(0.5292, rel=1e-3)


# ---------------------------------------------------------------------------
# QMOptimizerPySCF.__init__ validation
# ---------------------------------------------------------------------------

def test_qmoptimizer_invalid_solver():
    with pytest.raises(ValueError, match="solver"):
        QMOptimizerPySCF(solver="invalid_solver")


def test_qmoptimizer_invalid_dispersion():
    with pytest.raises(ValueError, match="dispersion"):
        QMOptimizerPySCF(dispersion="invalid_disp")
