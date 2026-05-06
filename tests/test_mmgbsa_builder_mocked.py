# -*- coding: utf-8 -*-
"""End-to-end tests for abmptools.genesis.mmgbsa.builder (subprocess mocked)."""
from __future__ import annotations

import shutil
from pathlib import Path
from unittest.mock import patch

import pytest

# Skip module if Biopython is unavailable (Stage 1 dependency).
pytest.importorskip("Bio.PDB", reason="Biopython required")

from abmptools.genesis.mmgbsa.builder import (
    MMGBSAOrchestrator,
    synthesize_targets_from_folder,
)
from abmptools.genesis.mmgbsa.models import (
    MMGBSABuildConfig,
    TargetSpec,
)


# Copy of the synthetic PDB used in the splitter tests.
SYNTHETIC_PDB = """\
HEADER    SYNTHETIC FOR TESTING
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 20.00           N
ATOM      2  CA  ALA A   1       1.500   0.000   0.000  1.00 20.00           C
ATOM      3  C   ALA A   1       2.000   1.400   0.000  1.00 20.00           C
ATOM      4  O   ALA A   1       1.300   2.400   0.000  1.00 20.00           O
HETATM    5  C1  LIG A 201       8.000   8.000   0.000  1.00 20.00           C
HETATM    6  C2  LIG A 201       9.000   8.500   0.000  1.00 20.00           C
HETATM    7  N1  LIG A 201      10.000   8.000   0.000  1.00 20.00           N
HETATM    8  O1  LIG A 201      10.500   9.000   0.000  1.00 20.00           O
TER
END
"""


# Real fixture log path (from POC). Tests using ``run`` need to inject
# fake atdyn output that contains a [STEP4] block; we copy from this
# fixture for deterministic reuse.
FIXTURE_LOG_DIR = (
    Path(__file__).parent / "fixtures" / "gbsa_logs" / "3772L-rename"
)


@pytest.fixture
def synthetic_input_dir(tmp_path: Path) -> Path:
    """Stage one or two PDB inputs under tmp_path/input/."""
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    (input_dir / "synth1.pdb").write_text(SYNTHETIC_PDB)
    (input_dir / "synth2.pdb").write_text(SYNTHETIC_PDB)
    return input_dir


def _build_cfg(input_dir: Path, output_dir: Path, n_targets: int = 1):
    targets = [
        TargetSpec(pdb=f"synth{i+1}.pdb", ligand_resno=201)
        for i in range(n_targets)
    ]
    return MMGBSABuildConfig(
        targets=targets,
        input_dir=str(input_dir),
        output_dir=str(output_dir),
        project_name="mocked_run",
    )


def _make_fake_run(fixture_log_dir: Path):
    """Return a side_effect mock for run_command across all tools."""
    def fake(cmd, cwd=None, capture=True, **kwargs):
        cwd = Path(cwd)
        argv0 = Path(cmd[0]).name
        if argv0 == "acpype":
            # "acpype -i ligand.pdb -c bcc ..." -- create the .acpype dir.
            ligand_pdb = Path(cmd[2])
            basename = ligand_pdb.stem
            acpype_dir = cwd / f"{basename}.acpype"
            acpype_dir.mkdir(exist_ok=True)
            (acpype_dir / f"{basename}_AC.frcmod").write_text("# fake frcmod\n")
            (acpype_dir / f"{basename}_bcc_gaff2.mol2").write_text("# fake mol2\n")
        elif argv0 == "tleap":
            # "tleap -f leaprc_<name>" -- create <name>.{prmtop,inpcrd,pdb}.
            leaprc = Path(cmd[2])
            name = leaprc.name.replace("leaprc_", "")
            (cwd / f"{name}.prmtop").write_text(f"# fake {name} prmtop\n")
            (cwd / f"{name}.inpcrd").write_text(f"# fake {name} inpcrd\n")
            (cwd / f"{name}.pdb").write_text(f"# fake {name} pdb\n")
        elif argv0 == "mpirun":
            # "mpirun -np 1 atdyn <name>.inp" -- emit STEP4 log via stdout.
            inp_path = Path(cmd[-1])
            name = inp_path.stem
            log_src = fixture_log_dir / f"{name}.log"
            stdout = log_src.read_text() if log_src.is_file() else (
                "[STEP4] Compute Single Point Energy for Molecules\n"
                "header\nheader\nheader\n"
                "    0   -100.0    0.0   0.0   0.0   0.0   0.0   0.0  -50.0\n"
            )
            from subprocess import CompletedProcess
            return CompletedProcess(cmd, 0, stdout, "")
        from subprocess import CompletedProcess
        return CompletedProcess(cmd, 0, "ok\n", "")
    return fake


# ---------------------------------------------------------------------------
# MMGBSAOrchestrator end-to-end
# ---------------------------------------------------------------------------

class TestMMGBSAOrchestratorEndToEnd:
    def test_run_produces_expected_layout(
        self, tmp_path, synthetic_input_dir,
    ):
        cfg = _build_cfg(
            input_dir=synthetic_input_dir,
            output_dir=tmp_path / "out",
            n_targets=1,
        )
        fake = _make_fake_run(FIXTURE_LOG_DIR)
        with patch(
            "abmptools.genesis.mmgbsa.system_builder.run_command",
            side_effect=fake,
        ), patch(
            "abmptools.genesis.mmgbsa.ligand_parameterize.run_command",
            side_effect=fake,
        ), patch(
            "abmptools.genesis.mmgbsa.gbsa_runner.run_command",
            side_effect=fake,
        ):
            orch = MMGBSAOrchestrator(cfg)
            result = orch.run()

        out = Path(result["output_dir"])
        # config dump exists.
        assert (out / "input" / "config.json").exists()
        # Per-target dir populated.
        target_dir = out / "synth1"
        assert (target_dir / "synth1_receptor_201.pdb").exists()
        assert (target_dir / "synth1_ligand_201.pdb").exists()
        assert (target_dir / "complex.prmtop").exists()
        assert (target_dir / "ligand.prmtop").exists()
        assert (target_dir / "receptor.prmtop").exists()
        assert (target_dir / "complex.log").exists()
        assert (target_dir / "ligand.log").exists()
        assert (target_dir / "receptor.log").exists()
        assert (target_dir / "complex.inp").exists()
        # POC compatibility.
        assert (out / "dirnames.in").exists()
        # Analysis outputs.
        assert (out / "analysis" / "analysis_results.csv").exists()
        assert (out / "analysis" / "dg_bind_plot.png").exists()

    def test_result_dict_keys(self, tmp_path, synthetic_input_dir):
        cfg = _build_cfg(synthetic_input_dir, tmp_path / "out", n_targets=1)
        fake = _make_fake_run(FIXTURE_LOG_DIR)
        with patch(
            "abmptools.genesis.mmgbsa.system_builder.run_command",
            side_effect=fake,
        ), patch(
            "abmptools.genesis.mmgbsa.ligand_parameterize.run_command",
            side_effect=fake,
        ), patch(
            "abmptools.genesis.mmgbsa.gbsa_runner.run_command",
            side_effect=fake,
        ):
            result = MMGBSAOrchestrator(cfg).run()
        for key in [
            "output_dir", "config_json", "n_targets", "n_succeeded",
            "n_failed", "failures", "splits", "csv_path", "png_path",
            "targets",
        ]:
            assert key in result
        assert result["n_targets"] == 1
        assert result["n_succeeded"] == 1
        assert result["n_failed"] == 0
        assert result["targets"][0]["dg_bind"] == pytest.approx(
            -38.2700, abs=1e-2
        )

    def test_two_targets_succeed(self, tmp_path, synthetic_input_dir):
        # Both targets reuse the same fixture logs (3772L). The
        # ΔG_bind value will be identical for each.
        cfg = _build_cfg(synthetic_input_dir, tmp_path / "out", n_targets=2)
        fake = _make_fake_run(FIXTURE_LOG_DIR)
        with patch(
            "abmptools.genesis.mmgbsa.system_builder.run_command",
            side_effect=fake,
        ), patch(
            "abmptools.genesis.mmgbsa.ligand_parameterize.run_command",
            side_effect=fake,
        ), patch(
            "abmptools.genesis.mmgbsa.gbsa_runner.run_command",
            side_effect=fake,
        ):
            result = MMGBSAOrchestrator(cfg).run()
        assert result["n_succeeded"] == 2
        for t in result["targets"]:
            assert t["dg_bind"] == pytest.approx(-38.2700, abs=1e-2)


# ---------------------------------------------------------------------------
# Folder-mode helper
# ---------------------------------------------------------------------------

class TestSynthesizeTargetsFromFolder:
    def test_enumerates_pdbs(self, synthetic_input_dir):
        targets = synthesize_targets_from_folder(
            synthetic_input_dir, ligand_resno=201, chain="A"
        )
        assert len(targets) == 2
        for t in targets:
            assert t.ligand_resno == 201
            assert t.chain == "A"
            assert t.pdb.endswith(".pdb")

    def test_missing_input_dir_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError, match="not found"):
            synthesize_targets_from_folder(tmp_path / "missing", 201)

    def test_no_pdbs_raises(self, tmp_path):
        empty = tmp_path / "empty"
        empty.mkdir()
        with pytest.raises(FileNotFoundError, match="No"):
            synthesize_targets_from_folder(empty, 201)
