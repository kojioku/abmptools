# -*- coding: utf-8 -*-
"""Tests for abmptools.genesis.grest.system_builder."""
from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import pytest

from abmptools.genesis.grest.models import (
    GrestBuildConfig,
    ReplicaTemperatureSpec,
    RESTSelectionSpec,
)
from abmptools.genesis.grest.system_builder import (
    build_amber_system,
    parse_box_dimensions,
    parse_n_protein_residues,
    write_tleap_script,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def base_cfg() -> GrestBuildConfig:
    return GrestBuildConfig(
        input_pdb="/tmp/protein.pdb",
        rest_selection=RESTSelectionSpec(
            mode="explicit", residues=["1-10"]
        ),
        replica_temperatures=ReplicaTemperatureSpec(
            mode="manual",
            temperatures=[300.0, 318.11, 337.11, 357.10],
        ),
    )


# ---------------------------------------------------------------------------
# write_tleap_script
# ---------------------------------------------------------------------------

class TestWriteTleapScript:
    def test_basic_directives(self, base_cfg, tmp_path):
        script = tmp_path / "system.tleap"
        write_tleap_script(
            cfg=base_cfg,
            input_pdb=tmp_path / "protein.pdb",
            prmtop=tmp_path / "system.prmtop",
            coor=tmp_path / "system.coor",
            ref_pdb=tmp_path / "system_ref.pdb",
            out_path=script,
        )
        text = script.read_text()
        assert "source leaprc.protein.ff19SB" in text
        assert "source leaprc.water.tip3p" in text
        assert f"system = loadpdb {tmp_path}/protein.pdb" in text
        assert "solvateBox system TIP3PBOX 10.00" in text
        assert "addions system Na+ 0" in text
        assert "addions system Cl- 0" in text
        assert "saveAmberParm system" in text
        assert "savepdb system" in text
        assert "quit" in text

    def test_neutralize_off(self, base_cfg, tmp_path):
        base_cfg.neutralize = False
        script = tmp_path / "system.tleap"
        write_tleap_script(
            cfg=base_cfg,
            input_pdb=tmp_path / "p.pdb",
            prmtop=tmp_path / "s.prmtop",
            coor=tmp_path / "s.coor",
            ref_pdb=tmp_path / "s_ref.pdb",
            out_path=script,
        )
        text = script.read_text()
        assert "addions" not in text

    def test_extra_ff_sources(self, base_cfg, tmp_path):
        base_cfg.ff_extra = ["leaprc.lipid21", "leaprc.gaff2"]
        script = tmp_path / "system.tleap"
        write_tleap_script(
            cfg=base_cfg,
            input_pdb=tmp_path / "p.pdb",
            prmtop=tmp_path / "s.prmtop",
            coor=tmp_path / "s.coor",
            ref_pdb=tmp_path / "s_ref.pdb",
            out_path=script,
        )
        text = script.read_text()
        assert "source leaprc.lipid21" in text
        assert "source leaprc.gaff2" in text

    def test_salt_concentration_emits_comment(self, base_cfg, tmp_path):
        # v1.20.0 emits a TODO comment rather than a precise ion count.
        base_cfg.salt_concentration_M = 0.15
        script = tmp_path / "system.tleap"
        write_tleap_script(
            cfg=base_cfg,
            input_pdb=tmp_path / "p.pdb",
            prmtop=tmp_path / "s.prmtop",
            coor=tmp_path / "s.coor",
            ref_pdb=tmp_path / "s_ref.pdb",
            out_path=script,
        )
        text = script.read_text()
        assert "salt_concentration_M = 0.15" in text


# ---------------------------------------------------------------------------
# parse_box_dimensions
# ---------------------------------------------------------------------------

class TestParseBoxDimensions:
    def test_parses_orthorhombic(self):
        prmtop = (
            "%VERSION  VERSION_STAMP = V0001.000\n"
            "%FLAG TITLE\n"
            "%FORMAT(20a4)\n"
            "test\n"
            "%FLAG BOX_DIMENSIONS\n"
            "%FORMAT(5E16.8)\n"
            "  9.00000000E+01  8.10568060E+01  8.35178610E+01  9.29056340E+01\n"
            "%FLAG RESIDUE_LABEL\n"
        )
        box = parse_box_dimensions(prmtop)
        assert box[0] == pytest.approx(81.056806)
        assert box[1] == pytest.approx(83.517861)
        assert box[2] == pytest.approx(92.905634)

    def test_missing_flag_raises(self):
        prmtop = "%FLAG TITLE\n%FORMAT(20a4)\ntest\n"
        with pytest.raises(ValueError, match="BOX_DIMENSIONS"):
            parse_box_dimensions(prmtop)

    def test_truncated_block_raises(self):
        prmtop = (
            "%FLAG BOX_DIMENSIONS\n"
            "%FORMAT(5E16.8)\n"
            "  9.00000000E+01  8.10000000E+01\n"  # only 2 values
            "%FLAG TITLE\n"
        )
        with pytest.raises(ValueError, match="numeric tokens"):
            parse_box_dimensions(prmtop)


# ---------------------------------------------------------------------------
# parse_n_protein_residues
# ---------------------------------------------------------------------------

class TestParseNProteinResidues:
    def test_basic_count(self):
        # 4 protein + 3 WAT + 1 Cl- => count=4
        prmtop = (
            "%FLAG RESIDUE_LABEL\n"
            "%FORMAT(20a4)\n"
            "LYS GLY GLY LEU WAT WAT WAT Cl- \n"
            "%FLAG NEXT\n"
        )
        n = parse_n_protein_residues(prmtop)
        assert n == 4

    def test_real_amber_format(self):
        # 20 chars per residue label, 5 protein residues + 2 water on second line.
        line1 = "ALA " + "GLY " + "VAL " + "LEU " + "ILE "  # 5 protein
        line2 = "WAT " + "WAT " * 19  # 20 water labels (filling block)
        prmtop = (
            "%FLAG RESIDUE_LABEL\n"
            "%FORMAT(20a4)\n"
            f"{line1}\n"
            f"{line2}\n"
            "%FLAG NEXT_FLAG\n"
        )
        n = parse_n_protein_residues(prmtop)
        assert n == 5

    def test_missing_flag_raises(self):
        with pytest.raises(ValueError, match="RESIDUE_LABEL"):
            parse_n_protein_residues("")

    def test_excludes_NA_and_CL_uppercase(self):
        # Some prmtops use 'NA' / 'CL' instead of 'Na+' / 'Cl-'.
        prmtop = (
            "%FLAG RESIDUE_LABEL\n"
            "%FORMAT(20a4)\n"
            "ARG LYS NA  CL  WAT \n"
            "%FLAG NEXT\n"
        )
        n = parse_n_protein_residues(prmtop)
        assert n == 2


# ---------------------------------------------------------------------------
# build_amber_system (mocked tleap)
# ---------------------------------------------------------------------------

class TestBuildAmberSystem:
    def test_full_flow_mocked(self, base_cfg, tmp_path):
        workdir = tmp_path / "build"

        # Fake tleap: write prmtop with BOX_DIMENSIONS + RESIDUE_LABEL,
        # plus coor + ref_pdb. Returns success.
        def fake_run(cmd, cwd=None, capture=True, **kwargs):
            cwd = Path(cwd)
            (cwd / "system.prmtop").write_text(
                "%VERSION  VERSION_STAMP = V0001.000\n"
                "%FLAG BOX_DIMENSIONS\n"
                "%FORMAT(5E16.8)\n"
                "  9.0E+01  8.0E+01  8.0E+01  8.5E+01\n"
                "%FLAG RESIDUE_LABEL\n"
                "%FORMAT(20a4)\n"
                "LYS GLY GLY WAT WAT WAT \n"
                "%FLAG NEXT\n"
            )
            (cwd / "system.coor").write_text("# fake coor")
            (cwd / "system_ref.pdb").write_text("# fake ref pdb")
            from subprocess import CompletedProcess
            return CompletedProcess(cmd, 0, "tleap done", "")

        with patch(
            "abmptools.genesis.grest.system_builder.run_command",
            side_effect=fake_run,
        ):
            result = build_amber_system(
                base_cfg, workdir,
                input_pdb=tmp_path / "fake_protein.pdb",
            )

        # Output files exist.
        assert result.prmtop.exists()
        assert result.coor.exists()
        assert result.ref_pdb.exists()
        assert result.leap_log.exists()
        # Box & residue count parsed.
        assert result.box_size_A[0] == pytest.approx(80.0)
        assert result.box_size_A[2] == pytest.approx(85.0)
        assert result.n_protein_residues == 3   # LYS GLY GLY
        # tleap script + log were emitted.
        assert (workdir / "system.tleap").exists()

    def test_missing_prmtop_after_tleap_raises(self, base_cfg, tmp_path):
        workdir = tmp_path / "build"

        def fake_run_no_prmtop(cmd, cwd=None, capture=True, **kwargs):
            from subprocess import CompletedProcess
            return CompletedProcess(cmd, 0, "tleap appeared to finish", "")

        from abmptools.genesis.grest._subprocess import CommandError
        with patch(
            "abmptools.genesis.grest.system_builder.run_command",
            side_effect=fake_run_no_prmtop,
        ):
            with pytest.raises(CommandError, match="missing"):
                build_amber_system(
                    base_cfg, workdir,
                    input_pdb=tmp_path / "fake.pdb",
                )
