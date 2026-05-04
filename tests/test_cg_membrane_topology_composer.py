# -*- coding: utf-8 -*-
"""Tests for abmptools.cg.membrane.topology_composer."""
from __future__ import annotations

from pathlib import Path

import pytest

from abmptools.cg.membrane import topology_composer


# ---------------------------------------------------------------------------
# get_moleculetype_name_from_itp
# ---------------------------------------------------------------------------

def test_get_moleculetype_name_basic(tmp_path):
    itp = tmp_path / "kgg.itp"
    itp.write_text(
        "[ moleculetype ]\n"
        "; name nrexcl\n"
        "molecule_0 1\n"
    )
    assert topology_composer.get_moleculetype_name_from_itp(itp) == "molecule_0"


def test_get_moleculetype_name_with_blank_lines_and_comments(tmp_path):
    itp = tmp_path / "kgg.itp"
    itp.write_text(
        "; preamble comment\n"
        "\n"
        "[ moleculetype ]\n"
        "\n"
        "; name nrexcl\n"
        "Protein_A     1\n"
        "[ atoms ]\n"
    )
    assert topology_composer.get_moleculetype_name_from_itp(itp) == "Protein_A"


def test_get_moleculetype_name_no_block_raises(tmp_path):
    itp = tmp_path / "broken.itp"
    itp.write_text("; just a comment\n[ atoms ]\n")
    with pytest.raises(ValueError, match="moleculetype"):
        topology_composer.get_moleculetype_name_from_itp(itp)


# ---------------------------------------------------------------------------
# compose_topology
# ---------------------------------------------------------------------------

@pytest.fixture
def insane_top_fixture(tmp_path):
    """A fixture mimicking insane's topol.top output."""
    p = tmp_path / "insane_topol.top"
    p.write_text(
        '#include "martini.itp"\n'
        '\n'
        '[ system ]\n'
        '; name\n'
        'INSANE!\n'
        '\n'
        '[ molecules ]\n'
        '; name  number\n'
        'Protein         1\n'
        'POPC           62\n'
        'POPC           63\n'
        'W             788\n'
        'NA+            12\n'
        'CL-            13\n'
    )
    return p


@pytest.fixture
def peptide_itp_fixture(tmp_path):
    p = tmp_path / "kgg.itp"
    p.write_text(
        "[ moleculetype ]\n"
        "molecule_0 1\n"
        "[ atoms ]\n"
    )
    return p


def test_compose_topology_includes_all_4_itps(
    insane_top_fixture, peptide_itp_fixture, tmp_path,
):
    out = tmp_path / "topol.top"
    topology_composer.compose_topology(
        insane_top_fixture, out,
        peptide_itp_path=peptide_itp_fixture,
        peptide_count=1,
    )
    text = out.read_text()
    assert '#include "martini_v3.0.0.itp"' in text
    assert '#include "martini_v3.0.0_solvents_v1.itp"' in text
    assert '#include "martini_v3.0.0_ions_v1.itp"' in text
    assert '#include "martini_v3.0.0_phospholipids_v1.itp"' in text


def test_compose_topology_inserts_peptide_itp_include(
    insane_top_fixture, peptide_itp_fixture, tmp_path,
):
    out = tmp_path / "topol.top"
    topology_composer.compose_topology(
        insane_top_fixture, out,
        peptide_itp_path=peptide_itp_fixture,
        peptide_itp_relpath="molecules/kgg/kgg.itp",
    )
    text = out.read_text()
    assert '#include "molecules/kgg/kgg.itp"' in text


def test_compose_topology_uses_absolute_when_no_relpath(
    insane_top_fixture, peptide_itp_fixture, tmp_path,
):
    out = tmp_path / "topol.top"
    topology_composer.compose_topology(
        insane_top_fixture, out,
        peptide_itp_path=peptide_itp_fixture,
    )
    text = out.read_text()
    assert str(peptide_itp_fixture.resolve()) in text


def test_compose_topology_replaces_protein_with_moleculetype(
    insane_top_fixture, peptide_itp_fixture, tmp_path,
):
    out = tmp_path / "topol.top"
    topology_composer.compose_topology(
        insane_top_fixture, out,
        peptide_itp_path=peptide_itp_fixture,
        peptide_count=1,
    )
    text = out.read_text()
    assert "molecule_0" in text
    # Bare "Protein" line in the molecules block is gone
    mol_block = text.split("[ molecules ]")[1]
    assert "Protein " not in mol_block
    assert "Protein\n" not in mol_block


def test_compose_topology_normalises_ion_names_in_molecules_block(
    insane_top_fixture, peptide_itp_fixture, tmp_path,
):
    out = tmp_path / "topol.top"
    topology_composer.compose_topology(
        insane_top_fixture, out,
        peptide_itp_path=peptide_itp_fixture,
    )
    mol_block = out.read_text().split("[ molecules ]")[1]
    assert "NA+" not in mol_block
    assert "CL-" not in mol_block
    # NA / CL still present (counts preserved)
    assert "NA " in mol_block
    assert "CL " in mol_block


def test_compose_topology_preserves_popc_lines(
    insane_top_fixture, peptide_itp_fixture, tmp_path,
):
    out = tmp_path / "topol.top"
    topology_composer.compose_topology(
        insane_top_fixture, out,
        peptide_itp_path=peptide_itp_fixture,
    )
    text = out.read_text()
    # Both POPC lines (per-leaflet) are kept
    assert text.count("POPC") >= 2


def test_compose_topology_uses_custom_peptide_count(
    insane_top_fixture, peptide_itp_fixture, tmp_path,
):
    out = tmp_path / "topol.top"
    topology_composer.compose_topology(
        insane_top_fixture, out,
        peptide_itp_path=peptide_itp_fixture,
        peptide_count=3,
    )
    mol_block = out.read_text().split("[ molecules ]")[1]
    # Find the molecule_0 line
    for line in mol_block.splitlines():
        if "molecule_0" in line:
            parts = line.split()
            assert parts[1] == "3"
            break
    else:
        pytest.fail("molecule_0 line not found in [ molecules ] block")


def test_compose_topology_replaces_system_title(
    insane_top_fixture, peptide_itp_fixture, tmp_path,
):
    out = tmp_path / "topol.top"
    topology_composer.compose_topology(
        insane_top_fixture, out,
        peptide_itp_path=peptide_itp_fixture,
        system_title="My KGG-POPC system",
    )
    text = out.read_text()
    assert "My KGG-POPC system" in text
    assert "INSANE!" not in text  # original placeholder removed


def test_compose_topology_no_martini_include_prepends(
    tmp_path, peptide_itp_fixture,
):
    insane_top = tmp_path / "weird_top.top"
    insane_top.write_text(
        "; no martini include here\n"
        "[ system ]\n"
        "STUB\n"
        "[ molecules ]\n"
        "Protein 1\n"
    )
    out = tmp_path / "out.top"
    topology_composer.compose_topology(
        insane_top, out,
        peptide_itp_path=peptide_itp_fixture,
    )
    text = out.read_text()
    assert '#include "martini_v3.0.0.itp"' in text
    assert text.index("martini_v3.0.0.itp") < text.index("[ molecules ]")


# ---------------------------------------------------------------------------
# normalize_ion_atom_names_gro
# ---------------------------------------------------------------------------

def test_normalize_ion_names_in_gro(tmp_path):
    """Verify NA+/CL- in resname/atomname columns are rewritten."""
    gro = tmp_path / "input.gro"
    # title, n_atoms, then 5+5+5+5 column atom lines, then box
    # Atom line layout:
    #  - cols 0-5  : resnr  (5d)
    #  - cols 5-10 : resname (5s, left-justified)
    #  - cols 10-15: atomname (5s, right-justified)
    #  - cols 15-20: atomnr (5d)
    #  - then 8.3f x 3 coords
    # Build lines explicitly to ensure column alignment.
    def atom_line(resnr, resname, atomname, atomnr, x, y, z):
        return (
            f"{resnr:5d}"
            f"{resname:<5s}"
            f"{atomname:>5s}"
            f"{atomnr:5d}"
            f"{x:8.3f}{y:8.3f}{z:8.3f}\n"
        )

    lines = [
        "test\n",
        "  3\n",
        atom_line(1, "POPC", "NC3", 1, 1.0, 2.0, 3.0),
        atom_line(2, "NA+", "NA+", 2, 4.0, 5.0, 6.0),
        atom_line(3, "CL-", "CL-", 3, 7.0, 8.0, 9.0),
        "  10.000  10.000  10.000\n",
    ]
    gro.write_text("".join(lines))

    out = tmp_path / "out.gro"
    topology_composer.normalize_ion_atom_names_gro(gro, out_path=out)
    text = out.read_text()

    # Atom lines after normalization: NA+/CL- replaced with NA / CL
    assert "NA+" not in text
    assert "CL-" not in text

    # Column widths preserved (length should equal original)
    in_lines = gro.read_text().splitlines(keepends=True)
    out_lines = out.read_text().splitlines(keepends=True)
    assert len(out_lines) == len(in_lines)
    for ol, il in zip(out_lines, in_lines):
        assert len(ol) == len(il)


def test_normalize_in_place_default(tmp_path):
    gro = tmp_path / "input.gro"
    gro.write_text(
        "test\n"
        "  1\n"
        f"{1:5d}{'NA+':<5s}{'NA+':>5s}{1:5d}{1.0:8.3f}{2.0:8.3f}{3.0:8.3f}\n"
        "  10.000  10.000  10.000\n"
    )
    topology_composer.normalize_ion_atom_names_gro(gro)
    assert "NA+" not in gro.read_text()
    assert "NA " in gro.read_text() or "NA\n" in gro.read_text()


def test_normalize_preserves_non_ion_lines(tmp_path):
    gro = tmp_path / "input.gro"
    gro.write_text(
        "test\n"
        "  2\n"
        f"{1:5d}{'POPC':<5s}{'NC3':>5s}{1:5d}{1.0:8.3f}{2.0:8.3f}{3.0:8.3f}\n"
        f"{2:5d}{'W':<5s}{'W':>5s}{2:5d}{4.0:8.3f}{5.0:8.3f}{6.0:8.3f}\n"
        "  10.000  10.000  10.000\n"
    )
    out = tmp_path / "out.gro"
    topology_composer.normalize_ion_atom_names_gro(gro, out_path=out)
    text = out.read_text()
    assert "POPC" in text
    assert "NC3" in text
    assert "W" in text


def test_normalize_handles_empty_gro_gracefully(tmp_path):
    gro = tmp_path / "tiny.gro"
    gro.write_text("title\n")  # only 1 line
    out = tmp_path / "out.gro"
    topology_composer.normalize_ion_atom_names_gro(gro, out_path=out)
    assert out.read_text() == "title\n"
