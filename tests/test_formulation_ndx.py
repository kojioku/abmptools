# -*- coding: utf-8 -*-
"""Unit tests for abmptools.formulation.ndx."""
from __future__ import annotations

from pathlib import Path

import pytest

from abmptools.formulation.ndx import (
    ION_RESNAMES,
    PROTEIN_RESNAMES,
    SOLVENT_RESNAMES,
    classify_atoms,
    parse_gro_residues,
    write_index_from_gro,
    write_ndx,
)


def _write_minimal_gro(path: Path) -> None:
    """Write a tiny .gro with 1 ALA + 1 CPR + 1 TCH + 2 SOL + 1 NA + 1 CL."""
    lines = [
        "test box",
        "    9",
        # resnr/resname/atomname/atomnr cols are 5 chars each (15 chars
        # of resnr+resname+atomname+atomnr, then x/y/z)
        "    1ALA      N    1   0.000   0.000   0.000",
        "    2CPR      C    2   1.000   0.000   0.000",
        "    3TCH      C    3   2.000   0.000   0.000",
        "    4SOL     OW    4   3.000   0.000   0.000",
        "    5SOL     OW    5   3.500   0.000   0.000",
        "    6NA      NA    6   4.000   0.000   0.000",
        "    7CL      CL    7   5.000   0.000   0.000",
        "    8FOO      X    8   6.000   0.000   0.000",
        "    9ALA     CA    9   0.100   0.000   0.000",
        "  10.000  10.000  10.000",
    ]
    path.write_text("\n".join(lines) + "\n")


def test_parse_gro_residues(tmp_path):
    g = tmp_path / "x.gro"
    _write_minimal_gro(g)
    atoms = parse_gro_residues(str(g))
    assert len(atoms) == 9
    assert atoms[0] == (1, "ALA", "N")
    assert atoms[1] == (2, "CPR", "C")
    assert atoms[5] == (6, "NA", "NA")


def test_classify_atoms_basic(tmp_path):
    g = tmp_path / "x.gro"
    _write_minimal_gro(g)
    atoms = parse_gro_residues(str(g))
    groups = classify_atoms(
        atoms, enhancer_resnames=["CPR"], bile_salt_resnames=["TCH"],
    )
    assert groups["Peptide"] == [1, 9]  # ALA × 2
    assert groups["Enhancer"] == [2]
    assert groups["BileSalt"] == [3]
    assert groups["Solvent"] == [4, 5]
    assert groups["Ions"] == [6, 7]
    assert groups["System"] == list(range(1, 10))


def test_classify_atoms_unclassified_stays_in_system(tmp_path):
    g = tmp_path / "x.gro"
    _write_minimal_gro(g)
    atoms = parse_gro_residues(str(g))
    groups = classify_atoms(
        atoms, enhancer_resnames=["CPR"], bile_salt_resnames=["TCH"],
    )
    # FOO is unclassified → present in System but absent from named groups
    assert 8 in groups["System"]
    assert 8 not in groups["Peptide"]
    assert 8 not in groups["Enhancer"]


def test_classify_atoms_peptide_solute_union(tmp_path):
    g = tmp_path / "x.gro"
    _write_minimal_gro(g)
    atoms = parse_gro_residues(str(g))
    groups = classify_atoms(
        atoms, enhancer_resnames=["CPR"], bile_salt_resnames=["TCH"],
    )
    assert sorted(groups["Peptide_Solute"]) == sorted(
        groups["Peptide"] + groups["Enhancer"] + groups["BileSalt"]
    )


def test_classify_atoms_non_peptide_is_complement(tmp_path):
    """Non_Peptide is the complement of Peptide_Solute (not just Peptide),
    so it does not overlap with Enhancer / BileSalt — required for
    the 2-group thermostat to avoid tc-grps atom overlap."""
    g = tmp_path / "x.gro"
    _write_minimal_gro(g)
    atoms = parse_gro_residues(str(g))
    groups = classify_atoms(
        atoms, enhancer_resnames=["CPR"], bile_salt_resnames=["TCH"],
    )
    solute = set(groups["Peptide_Solute"])
    assert all(i not in solute for i in groups["Non_Peptide"])
    assert set(groups["Non_Peptide"]) | solute == set(groups["System"])


def test_write_ndx_emits_all_groups(tmp_path):
    g = tmp_path / "x.gro"
    _write_minimal_gro(g)
    ndx = tmp_path / "x.ndx"
    write_index_from_gro(
        gro_path=str(g), ndx_path=str(ndx),
        enhancer_resnames=["CPR"], bile_salt_resnames=["TCH"],
    )
    text = ndx.read_text()
    for name in [
        "[ System ]", "[ Peptide ]", "[ Enhancer ]", "[ BileSalt ]",
        "[ Solvent ]", "[ Ions ]", "[ Peptide_Solute ]", "[ Non_Peptide ]",
    ]:
        assert name in text


def test_write_ndx_empty_group_ok(tmp_path):
    """No BileSalt residues in box → BileSalt group still emitted, empty."""
    groups = {
        "System": [1, 2],
        "Peptide": [1, 2],
        "Enhancer": [],
        "BileSalt": [],
        "Solvent": [],
        "Ions": [],
        "Peptide_Solute": [1, 2],
        "Non_Peptide": [],
    }
    ndx = tmp_path / "x.ndx"
    write_ndx(groups, str(ndx))
    text = ndx.read_text()
    assert "[ BileSalt ]" in text


def test_protein_resnames_include_caps():
    assert "ACE" in PROTEIN_RESNAMES
    assert "NME" in PROTEIN_RESNAMES


def test_solvent_aliases_present():
    assert SOLVENT_RESNAMES >= {"WAT", "HOH", "SOL"}


def test_ion_aliases_present():
    assert "NA" in ION_RESNAMES
    assert "CL" in ION_RESNAMES
