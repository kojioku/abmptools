# -*- coding: utf-8 -*-
"""Tests for abmptools.genesis.mmgbsa.pdb_splitter.

Biopython is required; tests auto-skip on environments without it
(``[mmgbsa]`` extras), matching abmptools' policy of optional extras
rather than hard core dependencies.
"""
from __future__ import annotations

from pathlib import Path

import pytest

# Skip the entire module if Biopython is unavailable.
pytest.importorskip("Bio.PDB", reason="Biopython required for pdb_splitter tests")

from abmptools.genesis.mmgbsa.models import TargetSpec
from abmptools.genesis.mmgbsa.pdb_splitter import (
    _resolve_pdb_path,
    split_pdb,
    split_target,
)


# ---------------------------------------------------------------------------
# Fixtures: build minimal PDB strings inline
# ---------------------------------------------------------------------------

# 5-atom protein residue (ALA) chain A + 4-atom ligand residue chain A,
# resno 201 (matching POC convention).
SYNTHETIC_PDB = """\
HEADER    SYNTHETIC FOR TESTING
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 20.00           N
ATOM      2  CA  ALA A   1       1.500   0.000   0.000  1.00 20.00           C
ATOM      3  C   ALA A   1       2.000   1.400   0.000  1.00 20.00           C
ATOM      4  O   ALA A   1       1.300   2.400   0.000  1.00 20.00           O
ATOM      5  CB  ALA A   1       2.000  -0.700  -1.200  1.00 20.00           C
ATOM      6  N   ALA A   2       3.300   1.600   0.000  1.00 20.00           N
ATOM      7  CA  ALA A   2       4.000   2.900   0.000  1.00 20.00           C
ATOM      8  C   ALA A   2       5.500   2.700   0.000  1.00 20.00           C
ATOM      9  O   ALA A   2       5.900   1.500   0.000  1.00 20.00           O
ATOM     10  CB  ALA A   2       3.700   3.700  -1.200  1.00 20.00           C
HETATM   11  C1  LIG A 201       8.000   8.000   0.000  1.00 20.00           C
HETATM   12  C2  LIG A 201       9.000   8.500   0.000  1.00 20.00           C
HETATM   13  N1  LIG A 201      10.000   8.000   0.000  1.00 20.00           N
HETATM   14  O1  LIG A 201      10.500   9.000   0.000  1.00 20.00           O
TER
END
"""

# Multi-chain variant: ligand 201 in chain A and another LIG residue
# (resno 201) also in chain B -- chain filter must distinguish.
MULTICHAIN_PDB = """\
ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C
HETATM    2  C1  LIG A 201       3.000   0.000   0.000  1.00 20.00           C
ATOM      3  CA  GLY B   1      10.000   0.000   0.000  1.00 20.00           C
HETATM    4  C1  LIG B 201      13.000   0.000   0.000  1.00 20.00           C
TER
END
"""


@pytest.fixture
def synthetic_pdb(tmp_path: Path) -> Path:
    p = tmp_path / "synth.pdb"
    p.write_text(SYNTHETIC_PDB)
    return p


@pytest.fixture
def multichain_pdb(tmp_path: Path) -> Path:
    p = tmp_path / "multi.pdb"
    p.write_text(MULTICHAIN_PDB)
    return p


# ---------------------------------------------------------------------------
# split_pdb
# ---------------------------------------------------------------------------

class TestSplitPdb:
    def test_basic_split(self, synthetic_pdb, tmp_path):
        out = tmp_path / "out"
        receptor, ligand = split_pdb(synthetic_pdb, ligand_resno=201, output_dir=out)

        assert receptor.exists()
        assert ligand.exists()
        assert receptor.name == "synth_receptor_201.pdb"
        assert ligand.name == "synth_ligand_201.pdb"
        # Original PDB copied into out.
        assert (out / "synth.pdb").exists()

    def test_receptor_excludes_ligand_residue(self, synthetic_pdb, tmp_path):
        receptor, ligand = split_pdb(
            synthetic_pdb, ligand_resno=201, output_dir=tmp_path / "out"
        )
        receptor_text = receptor.read_text()
        assert "LIG" not in receptor_text  # ligand residue absent
        assert "ALA" in receptor_text       # protein retained

    def test_ligand_only_keeps_target_residue(self, synthetic_pdb, tmp_path):
        receptor, ligand = split_pdb(
            synthetic_pdb, ligand_resno=201, output_dir=tmp_path / "out"
        )
        ligand_text = ligand.read_text()
        assert "LIG" in ligand_text
        assert "ALA" not in ligand_text

    def test_chain_filter_picks_correct_residue(self, multichain_pdb, tmp_path):
        out = tmp_path / "out"
        receptor, ligand = split_pdb(
            multichain_pdb, ligand_resno=201, output_dir=out, chain="A"
        )
        ligand_text = ligand.read_text()
        # Only chain-A LIG should be in the ligand file.
        # The chain B LIG must not leak into the ligand PDB.
        # (Bio.PDBIO emits a TER line referencing the residue too, so
        # we filter to HETATM/ATOM rows for a robust count.)
        atom_lines = [
            ln for ln in ligand_text.splitlines()
            if ln.startswith(("HETATM", "ATOM"))
        ]
        assert len(atom_lines) == 1
        assert "LIG A 201" in atom_lines[0]
        assert " B " not in ligand_text  # no chain-B leakage anywhere

    def test_chain_filter_receptor_keeps_other_chains(self, multichain_pdb, tmp_path):
        receptor, ligand = split_pdb(
            multichain_pdb,
            ligand_resno=201,
            output_dir=tmp_path / "out",
            chain="A",
        )
        receptor_text = receptor.read_text()
        # Chain A protein + chain B everything (LIG B 201 too) survive.
        assert "ALA" in receptor_text
        assert "GLY" in receptor_text

    def test_missing_pdb_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError, match="not found"):
            split_pdb(tmp_path / "nonexistent.pdb", 201, tmp_path / "out")

    def test_missing_residue_raises(self, synthetic_pdb, tmp_path):
        with pytest.raises(ValueError, match="not found"):
            split_pdb(synthetic_pdb, ligand_resno=9999, output_dir=tmp_path / "out")

    def test_custom_basename(self, synthetic_pdb, tmp_path):
        out = tmp_path / "out"
        receptor, ligand = split_pdb(
            synthetic_pdb,
            ligand_resno=201,
            output_dir=out,
            output_basename="9MM52-rename",
        )
        assert receptor.name == "9MM52-rename_receptor_201.pdb"
        assert ligand.name == "9MM52-rename_ligand_201.pdb"

    def test_output_dir_created(self, synthetic_pdb, tmp_path):
        # Nested non-existent dir.
        out = tmp_path / "deep" / "nest" / "out"
        receptor, _ = split_pdb(synthetic_pdb, 201, out)
        assert receptor.exists()
        assert out.is_dir()


# ---------------------------------------------------------------------------
# split_target (TargetSpec wrapper)
# ---------------------------------------------------------------------------

class TestSplitTarget:
    def test_basic(self, synthetic_pdb, tmp_path):
        target = TargetSpec(pdb=str(synthetic_pdb), ligand_resno=201)
        out_dir = tmp_path / "out"
        receptor, ligand = split_target(target, output_dir=out_dir)

        # Per-target subdir created automatically (POC convention).
        assert receptor.parent == out_dir / "synth"
        assert receptor.name == "synth_receptor_201.pdb"
        assert ligand.name == "synth_ligand_201.pdb"

    def test_resolves_pdb_via_input_dir(self, tmp_path):
        # Stage the PDB inside an "input" dir, reference it by basename only.
        input_dir = tmp_path / "input"
        input_dir.mkdir()
        pdb = input_dir / "abc.pdb"
        pdb.write_text(SYNTHETIC_PDB)

        target = TargetSpec(pdb="abc.pdb", ligand_resno=201)
        out_dir = tmp_path / "out"
        receptor, ligand = split_target(target, output_dir=out_dir, input_dir=input_dir)
        assert receptor.parent == out_dir / "abc"
        assert receptor.exists()


# ---------------------------------------------------------------------------
# _resolve_pdb_path
# ---------------------------------------------------------------------------

class TestResolvePdbPath:
    def test_absolute(self, tmp_path):
        p = tmp_path / "abs.pdb"
        p.write_text("ATOM\n")
        assert _resolve_pdb_path(str(p), input_dir=None) == p

    def test_via_input_dir(self, tmp_path):
        input_dir = tmp_path / "input"
        input_dir.mkdir()
        (input_dir / "named.pdb").write_text("ATOM\n")
        result = _resolve_pdb_path("named.pdb", input_dir=input_dir)
        assert result.name == "named.pdb"
        assert result.parent.name == "input"

    def test_unresolvable_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError, match="resolve"):
            _resolve_pdb_path("missing.pdb", input_dir=tmp_path)
