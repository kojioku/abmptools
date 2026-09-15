# -*- coding: utf-8 -*-
"""Partial charges supplied by the caller are used verbatim.

The reason this exists is polymers. Running AM1-BCC on a whole chain costs
far more than on its monomer, and a chain that came out of a random-walk
build can be strained enough that the charges are wrong rather than merely
expensive. The caller who knows the repeat structure can charge a short
oligomer, expand it unit by unit, and hand the numbers in.

What must not happen is a count mismatch passing quietly: a system built
from the wrong charges runs perfectly well and gives wrong answers.
"""
from __future__ import annotations

import pytest

from abmptools.amorphous.molecule_prep import read_charges


def test_reads_numbers_in_any_layout(tmp_path):
    f = tmp_path / "q.txt"
    f.write_text("# PVA 10-mer, from the 3-mer\n-0.1 0.2\n0.3\n\n-0.4  # tail\n")
    assert read_charges(str(f)) == [-0.1, 0.2, 0.3, -0.4]


def test_an_empty_file_is_an_error(tmp_path):
    f = tmp_path / "q.txt"
    f.write_text("# nothing but a comment\n")
    with pytest.raises(ValueError, match="No charges"):
        read_charges(str(f))


def test_a_count_mismatch_is_refused(tmp_path):
    """Wrong length must raise, not truncate.

    Truncating or cycling would build a system that runs and is wrong,
    which is the one outcome worth an exception.
    """
    openff = pytest.importorskip("openff.toolkit")
    from abmptools.amorphous.molecule_prep import _assign_charges_from_file

    mol = openff.Molecule.from_smiles("CCO")     # 9 atoms
    f = tmp_path / "q.txt"
    f.write_text("0.1 0.2 0.3\n")
    with pytest.raises(ValueError, match="3 charges but"):
        _assign_charges_from_file(mol, str(f))


def test_charges_land_on_the_molecule_in_order(tmp_path):
    openff = pytest.importorskip("openff.toolkit")
    from abmptools.amorphous.molecule_prep import _assign_charges_from_file

    mol = openff.Molecule.from_smiles("CCO")
    q = [0.05] * mol.n_atoms
    q[0] = -0.4
    q[-1] = 0.4 - 0.05 * (mol.n_atoms - 2)
    f = tmp_path / "q.txt"
    f.write_text(" ".join(f"{x:.6f}" for x in q))
    _assign_charges_from_file(mol, str(f))
    got = [float(x.m) for x in mol.partial_charges]
    assert got[0] == pytest.approx(-0.4)
    assert sum(got) == pytest.approx(0.0, abs=1e-6)
