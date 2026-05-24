# -*- coding: utf-8 -*-
"""Unit tests for abmptools.formulation.packer."""
from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import pytest

from abmptools.formulation import packer as fpk


def _write_dummy_pdb(p: Path, n_atoms: int = 3) -> None:
    lines = ["REMARK dummy"]
    for i in range(n_atoms):
        lines.append(
            f"HETATM{i+1:>5d}  C   LIG A   1    "
            f"{i:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00           C"
        )
    lines.append("END")
    p.write_text("\n".join(lines) + "\n")


def test_pack_mixed_solution_length_mismatch():
    with pytest.raises(ValueError, match="length mismatch"):
        fpk.pack_mixed_solution(
            component_pdbs=["a.pdb", "b.pdb"], counts=[1],
            box_size_nm=10.0, out_pdb="/tmp/out.pdb",
        )


def test_pack_mixed_solution_calls_run_packmol(tmp_path):
    a = tmp_path / "a.pdb"
    b = tmp_path / "b.pdb"
    _write_dummy_pdb(a)
    _write_dummy_pdb(b)
    out_pdb = tmp_path / "mix.pdb"

    with patch.object(fpk, "run_packmol") as mock_run:
        # simulate packmol writing the output
        mock_run.side_effect = lambda **kw: Path(kw["output_pdb"]).write_text("MIX")
        result = fpk.pack_mixed_solution(
            component_pdbs=[str(a), str(b)],
            counts=[2, 4],
            box_size_nm=10.0,
            out_pdb=str(out_pdb),
            tolerance_A=2.0,
            inner_box_margin_nm=0.5,
            seed=42,
        )
        assert result == str(out_pdb.resolve())
        kwargs = mock_run.call_args.kwargs
        # Inner box should be (10 - 0.5) = 9.5 nm
        assert abs(kwargs["box_size_nm"] - 9.5) < 1e-9
        assert kwargs["counts"] == [2, 4]
        assert kwargs["tolerance"] == 2.0
        assert kwargs["seed"] == 42
        assert kwargs["packmol_path"] == "packmol"


def test_pack_mixed_solution_resolves_paths(tmp_path):
    a = tmp_path / "a.pdb"
    _write_dummy_pdb(a)
    out_pdb = tmp_path / "out" / "mix.pdb"

    with patch.object(fpk, "run_packmol") as mock_run:
        mock_run.side_effect = lambda **kw: Path(kw["output_pdb"]).write_text("OK")
        fpk.pack_mixed_solution(
            component_pdbs=["a.pdb"], counts=[1],
            box_size_nm=5.0, out_pdb=str(out_pdb),
        )
        kwargs = mock_run.call_args.kwargs
        assert kwargs["pdb_paths"][0].endswith("a.pdb")
        assert Path(kwargs["pdb_paths"][0]).is_absolute()


def test_pack_mixed_solution_inner_edge_floor(tmp_path):
    """If margin >= box, fall back to a minimal 0.5 nm edge instead of negative."""
    a = tmp_path / "a.pdb"
    _write_dummy_pdb(a)
    with patch.object(fpk, "run_packmol") as mock_run:
        mock_run.side_effect = lambda **kw: Path(kw["output_pdb"]).write_text("OK")
        fpk.pack_mixed_solution(
            component_pdbs=[str(a)], counts=[1],
            box_size_nm=0.3, out_pdb=str(tmp_path / "out.pdb"),
            inner_box_margin_nm=1.0,
        )
        kwargs = mock_run.call_args.kwargs
        assert kwargs["box_size_nm"] >= 0.5
