# -*- coding: utf-8 -*-
"""Tests for abmptools.cg.peptide.top_writer module."""
from __future__ import annotations

from pathlib import Path

from abmptools.cg.peptide.models import PeptideBuildConfig, PeptideSpec
from abmptools.cg.peptide.top_writer import (
    M3_FF_INCLUDES,
    _get_molecule_name_from_itp,
    write_topol_top,
)


def _stub_itp(path: Path, mol_name: str) -> Path:
    path.write_text(
        f"""\
[ moleculetype ]
; molname    nrexcl
{mol_name}    1

[ atoms ]
; (stub)
"""
    )
    return path


class TestGetMoleculeNameFromItp:
    def test_extracts_name(self, tmp_path):
        itp = _stub_itp(tmp_path / "x.itp", "MyMol")
        assert _get_molecule_name_from_itp(itp) == "MyMol"

    def test_default_when_missing(self, tmp_path):
        assert _get_molecule_name_from_itp(tmp_path / "missing.itp") == "molecule"

    def test_handles_no_space_bracket(self, tmp_path):
        # `[moleculetype]` (no spaces) is also accepted
        path = tmp_path / "x.itp"
        path.write_text("[moleculetype]\nFoo  1\n")
        assert _get_molecule_name_from_itp(path) == "Foo"


class TestWriteTopolTop:
    def test_includes_m3_force_fields(self, tmp_path):
        cfg = PeptideBuildConfig(
            peptides=[PeptideSpec(name="kgg", sequence="KGG", count=5)],
        )
        itp = _stub_itp(tmp_path / "kgg.itp", "Protein_A")
        out = write_topol_top(cfg, {"kgg": itp}, tmp_path / "topol.top")

        text = out.read_text()
        for inc in M3_FF_INCLUDES:
            assert inc in text
        assert '#include "molecules/kgg/kgg.itp"' in text
        assert "[ system ]" in text
        assert "[ molecules ]" in text
        assert "Protein_A  5" in text

    def test_solvent_comment_when_enabled(self, tmp_path):
        cfg = PeptideBuildConfig(
            peptides=[PeptideSpec(name="kgg", sequence="KGG")],
            solvent_enabled=True,
        )
        itp = _stub_itp(tmp_path / "kgg.itp", "P")
        out = write_topol_top(cfg, {"kgg": itp}, tmp_path / "topol.top")
        assert "gmx solvate" in out.read_text()

    def test_no_solvent_comment_when_disabled(self, tmp_path):
        cfg = PeptideBuildConfig(
            peptides=[PeptideSpec(name="kgg", sequence="KGG")],
            solvent_enabled=False,
        )
        itp = _stub_itp(tmp_path / "kgg.itp", "P")
        out = write_topol_top(cfg, {"kgg": itp}, tmp_path / "topol.top")
        assert "gmx solvate" not in out.read_text()

    def test_two_peptide_species(self, tmp_path):
        cfg = PeptideBuildConfig(
            peptides=[
                PeptideSpec(name="a", sequence="A", count=3),
                PeptideSpec(name="b", sequence="G", count=7),
            ],
        )
        itp_a = _stub_itp(tmp_path / "a.itp", "MolA")
        itp_b = _stub_itp(tmp_path / "b.itp", "MolB")
        out = write_topol_top(
            cfg, {"a": itp_a, "b": itp_b}, tmp_path / "topol.top",
        )
        text = out.read_text()
        assert "MolA  3" in text
        assert "MolB  7" in text
