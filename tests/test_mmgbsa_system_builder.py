# -*- coding: utf-8 -*-
"""Tests for abmptools.genesis.mmgbsa.{ligand_parameterize,system_builder}."""
from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import pytest

from abmptools.genesis.mmgbsa.ligand_parameterize import (
    AcpypeResult,
    _locate_outputs,
    run_acpype,
)
from abmptools.genesis.mmgbsa.models import ForceFieldSet, LigandParameterization
from abmptools.genesis.mmgbsa.system_builder import (
    build_three_systems,
    render_complex_leaprc,
    render_ligand_leaprc,
    render_receptor_leaprc,
)


# ---------------------------------------------------------------------------
# leaprc rendering
# ---------------------------------------------------------------------------

class TestRenderLeaprc:
    def test_complex_has_six_source_lines(self, tmp_path):
        ff = ForceFieldSet()
        text = render_complex_leaprc(
            ff=ff,
            receptor_pdb=tmp_path / "rec.pdb",
            frcmod=tmp_path / "lig.frcmod",
            mol2=tmp_path / "lig.mol2",
            prmtop=tmp_path / "complex.prmtop",
            inpcrd=tmp_path / "complex.inpcrd",
            out_pdb=tmp_path / "complex.pdb",
        )
        # POC's 6 source lines (protein, dna, rna, water, gaff2, gaff).
        assert "source leaprc.protein.ff14SB" in text
        assert "source leaprc.DNA.OL15" in text
        assert "source leaprc.RNA.OL3" in text
        assert "source leaprc.water.tip3p" in text
        assert "source leaprc.gaff2" in text
        assert "source leaprc.gaff" in text

    def test_complex_has_addpdb_maps(self, tmp_path):
        text = render_complex_leaprc(
            ff=ForceFieldSet(),
            receptor_pdb=tmp_path / "rec.pdb",
            frcmod=tmp_path / "f.frcmod",
            mol2=tmp_path / "m.mol2",
            prmtop=tmp_path / "c.prmtop",
            inpcrd=tmp_path / "c.inpcrd",
            out_pdb=tmp_path / "c.pdb",
        )
        # POC addPdbResMap: ions normalised.
        assert 'addPdbResMap' in text
        assert '{ "NA+" "NA" }' in text
        assert '{ "CL-" "CL" }' in text
        assert 'addPdbAtomMap' in text

    def test_complex_combines_protein_and_ligand(self, tmp_path):
        text = render_complex_leaprc(
            ff=ForceFieldSet(),
            receptor_pdb=tmp_path / "rec.pdb",
            frcmod=tmp_path / "f.frcmod",
            mol2=tmp_path / "m.mol2",
            prmtop=tmp_path / "c.prmtop",
            inpcrd=tmp_path / "c.inpcrd",
            out_pdb=tmp_path / "c.pdb",
        )
        assert "loadAmberParams" in text
        assert "loadmol2" in text
        assert "complex = combine { protein ligand }" in text
        assert "saveAmberParm complex" in text

    def test_ligand_only(self, tmp_path):
        text = render_ligand_leaprc(
            ff=ForceFieldSet(),
            frcmod=tmp_path / "f.frcmod",
            mol2=tmp_path / "m.mol2",
            prmtop=tmp_path / "lig.prmtop",
            inpcrd=tmp_path / "lig.inpcrd",
            out_pdb=tmp_path / "lig.pdb",
        )
        # No "loadpdb" -- ligand has only mol2.
        assert "loadpdb" not in text
        assert "loadmol2" in text
        assert "saveAmberParm ligand" in text

    def test_receptor_only(self, tmp_path):
        text = render_receptor_leaprc(
            ff=ForceFieldSet(),
            receptor_pdb=tmp_path / "rec.pdb",
            prmtop=tmp_path / "rec.prmtop",
            inpcrd=tmp_path / "rec.inpcrd",
            out_pdb=tmp_path / "rec.pdb",
        )
        # No mol2 / loadAmberParams -- receptor is protein-only.
        assert "loadmol2" not in text
        assert "loadAmberParams" not in text
        assert "loadpdb" in text
        assert "saveAmberParm protein" in text

    def test_extra_lines_appended(self, tmp_path):
        ff = ForceFieldSet(extra_lines=["loadAmberParams custom.frcmod"])
        text = render_receptor_leaprc(
            ff=ff,
            receptor_pdb=tmp_path / "rec.pdb",
            prmtop=tmp_path / "r.prmtop",
            inpcrd=tmp_path / "r.inpcrd",
            out_pdb=tmp_path / "r.pdb",
        )
        assert "loadAmberParams custom.frcmod" in text


# ---------------------------------------------------------------------------
# build_three_systems (mocked tleap)
# ---------------------------------------------------------------------------

class TestBuildThreeSystems:
    def test_full_flow_mocked(self, tmp_path):
        workdir = tmp_path / "work"
        receptor_pdb = tmp_path / "receptor.pdb"
        receptor_pdb.write_text("# fake\n")
        frcmod = tmp_path / "lig.frcmod"
        frcmod.write_text("# frcmod\n")
        mol2 = tmp_path / "lig.mol2"
        mol2.write_text("# mol2\n")

        def fake_run(cmd, cwd=None, capture=True, **kwargs):
            cwd = Path(cwd)
            # tleap creates the prmtop/inpcrd referenced in its leaprc.
            # We synthesise complex/ligand/receptor based on which
            # leaprc is being run.
            leaprc = cmd[2]  # ['tleap', '-f', '<leaprc>']
            if leaprc.endswith("leaprc_complex"):
                (cwd / "complex.prmtop").write_text("# complex prmtop\n")
                (cwd / "complex.inpcrd").write_text("# complex inpcrd\n")
                (cwd / "complex.pdb").write_text("# complex pdb\n")
            elif leaprc.endswith("leaprc_ligand"):
                (cwd / "ligand.prmtop").write_text("# ligand prmtop\n")
                (cwd / "ligand.inpcrd").write_text("# ligand inpcrd\n")
                (cwd / "ligand.pdb").write_text("# ligand pdb\n")
            elif leaprc.endswith("leaprc_receptor"):
                (cwd / "receptor.prmtop").write_text("# receptor prmtop\n")
                (cwd / "receptor.inpcrd").write_text("# receptor inpcrd\n")
                (cwd / "receptor.pdb").write_text("# receptor pdb\n")
            from subprocess import CompletedProcess
            return CompletedProcess(cmd, 0, "tleap done\n", "")

        with patch(
            "abmptools.genesis.mmgbsa.system_builder.run_command",
            side_effect=fake_run,
        ):
            result = build_three_systems(
                receptor_pdb=receptor_pdb,
                frcmod=frcmod,
                mol2=mol2,
                workdir=workdir,
                ff=ForceFieldSet(),
            )

        # All 3 systems present + their leaprc + log.
        for sys in (result.complex, result.ligand, result.receptor):
            assert sys.leaprc.exists()
            assert sys.prmtop.exists()
            assert sys.inpcrd.exists()
            assert sys.pdb.exists()
            assert sys.log.exists()

    def test_missing_prmtop_raises(self, tmp_path):
        workdir = tmp_path / "work"

        def fake_run_no_outputs(cmd, cwd=None, capture=True, **kwargs):
            from subprocess import CompletedProcess
            return CompletedProcess(cmd, 0, "ok", "")

        from abmptools.genesis.mmgbsa._subprocess import CommandError
        with patch(
            "abmptools.genesis.mmgbsa.system_builder.run_command",
            side_effect=fake_run_no_outputs,
        ):
            with pytest.raises(CommandError, match="missing"):
                build_three_systems(
                    receptor_pdb=tmp_path / "rec.pdb",
                    frcmod=tmp_path / "f",
                    mol2=tmp_path / "m",
                    workdir=workdir,
                    ff=ForceFieldSet(),
                )


# ---------------------------------------------------------------------------
# acpype runner (mocked)
# ---------------------------------------------------------------------------

class TestRunAcpype:
    def test_full_flow_mocked(self, tmp_path):
        ligand_pdb = tmp_path / "ligand.pdb"
        ligand_pdb.write_text("ATOM\n")
        workdir = tmp_path / "work"

        def fake_run(cmd, cwd=None, capture=True, **kwargs):
            cwd = Path(cwd)
            # acpype creates <basename>.acpype/ with the canonical files.
            acpype_dir = cwd / "ligand.acpype"
            acpype_dir.mkdir()
            (acpype_dir / "ligand_AC.frcmod").write_text("# frcmod\n")
            (acpype_dir / "ligand_bcc_gaff2.mol2").write_text("# mol2\n")
            from subprocess import CompletedProcess
            return CompletedProcess(cmd, 0, "ok", "")

        with patch(
            "abmptools.genesis.mmgbsa.ligand_parameterize.run_command",
            side_effect=fake_run,
        ):
            result = run_acpype(
                ligand_pdb=ligand_pdb,
                workdir=workdir,
                config=LigandParameterization(),
            )

        assert isinstance(result, AcpypeResult)
        assert result.acpype_dir.is_dir()
        assert result.frcmod.name == "ligand_AC.frcmod"
        assert result.mol2.name == "ligand_bcc_gaff2.mol2"

    def test_skip_if_cached(self, tmp_path):
        ligand_pdb = tmp_path / "ligand.pdb"
        ligand_pdb.write_text("ATOM\n")
        workdir = tmp_path / "work"
        acpype_dir = workdir / "ligand.acpype"
        acpype_dir.mkdir(parents=True)
        (acpype_dir / "ligand_AC.frcmod").write_text("# pre-existing\n")
        (acpype_dir / "ligand_bcc_gaff2.mol2").write_text("# pre-existing\n")

        # If skip_if_cached=True (default) and outputs already exist,
        # run_command must not be called.
        with patch(
            "abmptools.genesis.mmgbsa.ligand_parameterize.run_command",
        ) as m_run:
            result = run_acpype(
                ligand_pdb=ligand_pdb,
                workdir=workdir,
                config=LigandParameterization(skip_if_cached=True),
            )
        m_run.assert_not_called()
        assert result.frcmod.name == "ligand_AC.frcmod"

    def test_gaff_atom_type_picks_correct_mol2(self, tmp_path):
        # Test _locate_outputs prefers _bcc_gaff.mol2 when atom_type=gaff.
        acpype_dir = tmp_path / "ligand.acpype"
        acpype_dir.mkdir()
        (acpype_dir / "ligand_AC.frcmod").write_text("# frcmod\n")
        (acpype_dir / "ligand_bcc_gaff.mol2").write_text("# gaff mol2\n")
        (acpype_dir / "ligand_bcc_gaff2.mol2").write_text("# gaff2 mol2\n")

        result_gaff = _locate_outputs(acpype_dir, atom_type="gaff")
        assert result_gaff.mol2.name == "ligand_bcc_gaff.mol2"

        result_gaff2 = _locate_outputs(acpype_dir, atom_type="gaff2")
        assert result_gaff2.mol2.name == "ligand_bcc_gaff2.mol2"

    def test_missing_acpype_dir_raises(self, tmp_path):
        ligand_pdb = tmp_path / "ligand.pdb"
        ligand_pdb.write_text("ATOM\n")
        workdir = tmp_path / "work"

        def fake_run_no_outputs(cmd, cwd=None, capture=True, **kwargs):
            from subprocess import CompletedProcess
            return CompletedProcess(cmd, 0, "ok", "")

        from abmptools.genesis.mmgbsa._subprocess import CommandError
        with patch(
            "abmptools.genesis.mmgbsa.ligand_parameterize.run_command",
            side_effect=fake_run_no_outputs,
        ):
            with pytest.raises(CommandError, match="missing"):
                run_acpype(
                    ligand_pdb=ligand_pdb,
                    workdir=workdir,
                    config=LigandParameterization(skip_if_cached=False),
                )
