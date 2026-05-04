# -*- coding: utf-8 -*-
"""Tests for abmptools.cg.membrane.system_packer."""
from __future__ import annotations

import pytest

from abmptools.cg.membrane import system_packer
from abmptools.cg.membrane.system_packer import classify_residue_cg


# ---------------------------------------------------------------------------
# classify_residue_cg
# ---------------------------------------------------------------------------

class TestClassifyResidueCG:
    @pytest.mark.parametrize("name", ["POPC", "DOPC", "POPE", "DPPC"])
    def test_lipid(self, name):
        assert classify_residue_cg(name) == "lipid"

    @pytest.mark.parametrize("name", ["LYS", "GLY", "ALA", "TRP", "PHE"])
    def test_peptide(self, name):
        assert classify_residue_cg(name) == "peptide"

    def test_solvent_w(self):
        assert classify_residue_cg("W") == "solvent"

    def test_solvent_wn(self):
        assert classify_residue_cg("WN") == "solvent"

    def test_na_normalized(self):
        assert classify_residue_cg("NA") == "na"

    def test_na_with_plus(self):
        assert classify_residue_cg("NA+") == "na"

    def test_cl_normalized(self):
        assert classify_residue_cg("CL") == "cl"

    def test_cl_with_minus(self):
        assert classify_residue_cg("CL-") == "cl"

    def test_other(self):
        assert classify_residue_cg("XYZ") == "other"

    def test_case_insensitive(self):
        assert classify_residue_cg("popc") == "lipid"
        assert classify_residue_cg("lys") == "peptide"


# ---------------------------------------------------------------------------
# write_ndx_from_gro_cg (integration via fixture .gro)
# ---------------------------------------------------------------------------

def _make_gro(tmp_path, atoms):
    """Build a minimal .gro with the given (resnr, resname, atomname) tuples."""
    lines = ["title\n", f"{len(atoms):5d}\n"]
    for i, (resnr, resname, atomname) in enumerate(atoms, start=1):
        lines.append(
            f"{resnr:5d}"
            f"{resname:<5s}"
            f"{atomname:>5s}"
            f"{i:5d}"
            f"{1.0:8.3f}{2.0:8.3f}{3.0:8.3f}\n"
        )
    lines.append("  10.000  10.000  10.000\n")
    p = tmp_path / "system.gro"
    p.write_text("".join(lines))
    return p


def test_ndx_groups_present(tmp_path):
    gro = _make_gro(tmp_path, [
        (1, "LYS", "BB"),
        (1, "LYS", "SC1"),
        (2, "GLY", "BB"),
        (3, "POPC", "NC3"),
        (3, "POPC", "PO4"),
        (4, "POPC", "NC3"),
        (5, "W", "W"),
        (6, "W", "W"),
        (7, "NA", "NA"),
        (8, "CL", "CL"),
    ])
    ndx = tmp_path / "index.ndx"
    system_packer.write_ndx_from_gro_cg(
        gro_path=str(gro), ndx_path=str(ndx),
    )
    text = ndx.read_text()
    assert "[ System ]" in text
    assert "[ Bilayer ]" in text
    assert "[ Peptide ]" in text
    assert "[ W ]" in text
    assert "[ NA ]" in text
    assert "[ CL ]" in text
    assert "[ Non_Bilayer ]" in text


def test_ndx_bilayer_contains_only_lipid_atoms(tmp_path):
    gro = _make_gro(tmp_path, [
        (1, "LYS", "BB"),    # idx 1 - peptide
        (2, "POPC", "NC3"),  # idx 2 - lipid
        (2, "POPC", "PO4"),  # idx 3 - lipid
        (3, "W", "W"),       # idx 4 - solvent
    ])
    ndx = tmp_path / "index.ndx"
    system_packer.write_ndx_from_gro_cg(
        gro_path=str(gro), ndx_path=str(ndx),
    )
    text = ndx.read_text()

    # Find Bilayer block
    sections = _split_ndx_sections(text)
    assert "Bilayer" in sections
    bilayer_atoms = sections["Bilayer"]
    assert bilayer_atoms == [2, 3]


def test_ndx_non_bilayer_excludes_lipid(tmp_path):
    gro = _make_gro(tmp_path, [
        (1, "LYS", "BB"),    # idx 1 - peptide
        (2, "POPC", "NC3"),  # idx 2 - lipid
        (3, "W", "W"),       # idx 3 - solvent
        (4, "NA", "NA"),     # idx 4 - ion
    ])
    ndx = tmp_path / "index.ndx"
    system_packer.write_ndx_from_gro_cg(
        gro_path=str(gro), ndx_path=str(ndx),
    )
    sections = _split_ndx_sections(ndx.read_text())
    assert sections["Non_Bilayer"] == [1, 3, 4]
    assert sections["System"] == [1, 2, 3, 4]


def test_ndx_handles_unnormalized_ions(tmp_path):
    """Ion residues with charge sign (NA+ / CL-) still classify correctly."""
    gro = _make_gro(tmp_path, [
        (1, "POPC", "NC3"),
        (2, "NA+", "NA+"),
        (3, "CL-", "CL-"),
    ])
    ndx = tmp_path / "index.ndx"
    system_packer.write_ndx_from_gro_cg(
        gro_path=str(gro), ndx_path=str(ndx),
    )
    sections = _split_ndx_sections(ndx.read_text())
    assert sections["NA"] == [2]
    assert sections["CL"] == [3]


def test_ndx_omits_empty_groups(tmp_path):
    """A bilayer-only system has no Peptide group."""
    gro = _make_gro(tmp_path, [
        (1, "POPC", "NC3"),
        (2, "POPC", "NC3"),
        (3, "W", "W"),
    ])
    ndx = tmp_path / "index.ndx"
    system_packer.write_ndx_from_gro_cg(
        gro_path=str(gro), ndx_path=str(ndx),
    )
    text = ndx.read_text()
    assert "[ Peptide ]" not in text
    assert "[ NA ]" not in text
    assert "[ CL ]" not in text
    assert "[ Bilayer ]" in text


# ---------------------------------------------------------------------------
# Re-exports work
# ---------------------------------------------------------------------------

def test_add_ions_cg_is_reexported():
    assert callable(system_packer.add_ions_cg)


def test_packed_system_dataclass_reexported():
    assert system_packer.PackedSystem is not None


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _split_ndx_sections(text):
    """Parse a .ndx text into {group_name: [atom_indices]}."""
    sections = {}
    cur = None
    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith("[") and line.endswith("]"):
            cur = line[1:-1].strip()
            sections[cur] = []
        elif cur is not None:
            sections[cur].extend(int(t) for t in line.split())
    return sections
