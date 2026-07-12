"""Tests for residue-name assignment in amorphous molecule preparation."""
from types import SimpleNamespace

from abmptools.amorphous.molecule_prep import (
    _residue_name_for, apply_residue_name,
)


def test_residue_name_for():
    assert _residue_name_for("IMC") == "IMC"
    assert _residue_name_for("PVP") == "PVP"
    # long names truncate to the 5-char GROMACS .gro residue field
    assert _residue_name_for("acetaminophen") == "ACETA"
    # non-alphanumerics stripped, upper-cased
    assert _residue_name_for("2-hydroxy") == "2HYDR"
    assert _residue_name_for("") == "MOL"
    assert _residue_name_for(None) == "MOL"


def test_apply_residue_name_sets_metadata():
    atoms = [SimpleNamespace(metadata={}) for _ in range(4)]
    mol = SimpleNamespace(name="PVP", atoms=atoms)
    apply_residue_name(mol)
    assert all(a.metadata["residue_name"] == "PVP" for a in atoms)


def test_apply_residue_name_defaults_to_mol():
    atoms = [SimpleNamespace(metadata={})]
    apply_residue_name(SimpleNamespace(name="", atoms=atoms))
    assert atoms[0].metadata["residue_name"] == "MOL"
