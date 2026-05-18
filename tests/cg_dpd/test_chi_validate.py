# -*- coding: utf-8 -*-
"""tests/cg_dpd/test_chi_validate.py — G2 (a→chi 逆変換) + validate() + verify CLI."""
from __future__ import annotations

import math
import subprocess
import sys
from pathlib import Path

import pytest


# --- G2: a → chi 逆変換 ----------------------------------------------------

def test_to_chi_values_from_a_mode():
    """a モードからの chi 逆変換 (Groot-Warren 逆向き)."""
    from abmptools.cg.dpd import AijMatrix
    aij = AijMatrix(
        segments=["A", "B"], pairs=[("A","B", 28.27)], mode="a", aii=25.0,
    )
    chi_vals = aij.to_chi_values()
    # chi = (28.27 - 25.0) * 0.306 = 1.00062
    assert math.isclose(chi_vals[0][2], (28.27 - 25.0) * 0.306, rel_tol=1e-6)


def test_to_chi_values_identity_in_chi_mode():
    """chi モードならそのまま返す."""
    from abmptools.cg.dpd import AijMatrix
    aij = AijMatrix(
        segments=["A", "B"], pairs=[("A","B", -4.108)], mode="chi", aii=25.0,
    )
    assert aij.to_chi_values() == [("A","B", -4.108)]


def test_a_chi_round_trip():
    """a → chi → a の round-trip で値が一致 (Groot-Warren で可逆)."""
    from abmptools.cg.dpd import AijMatrix
    orig = AijMatrix(
        segments=["A","B","W"],
        pairs=[("A","A",25.0),("B","B",25.0),("W","W",25.0),
               ("A","B",28.0),("A","W",35.0),("B","W",30.0)],
        mode="a", aii=25.0,
    )
    chi_vals = orig.to_chi_values()
    # chi 形式の AijMatrix を再構築
    reconst = AijMatrix(segments=orig.segments, pairs=chi_vals, mode="chi", aii=orig.aii)
    a_vals_back = reconst.to_a_values()
    for (i1, j1, v1), (i2, j2, v2) in zip(orig.pairs, a_vals_back):
        assert (i1, j1) == (i2, j2)
        assert math.isclose(v1, v2, abs_tol=1e-9)


def test_write_aij_force_chi(tmp_path):
    """write_aij(out_mode='chi') で a → chi 強制変換."""
    from abmptools.cg.dpd import AijMatrix, read_aij, write_aij
    src = AijMatrix(
        segments=["A","B"], pairs=[("A","A",25.0),("A","B",28.27)],
        mode="a", aii=25.0,
    )
    p = tmp_path / "chi_out.dat"
    write_aij(p, src, out_mode="chi")
    text = p.read_text()
    assert text.startswith("# fcews aij.dat")
    assert "chi = [" in text
    assert "aij = [" not in text
    # 再読み込みで mode='chi' になる
    re = read_aij(p)
    assert re.mode == "chi"
    # chi=(28.27-25)*0.306=1.00062
    assert math.isclose(
        next(v for (i,j,v) in re.pairs if (i,j)==("A","B")),
        (28.27 - 25.0) * 0.306, rel_tol=1e-6,
    )


def test_write_aij_force_a_from_chi(tmp_path):
    """write_aij(out_mode='a') で chi → a 強制変換."""
    from abmptools.cg.dpd import AijMatrix, read_aij, write_aij
    src = AijMatrix(
        segments=["A","B"], pairs=[("A","B", 1.0)], mode="chi", aii=25.0,
    )
    p = tmp_path / "a_out.dat"
    write_aij(p, src, out_mode="a")
    text = p.read_text()
    assert "aij = [" in text
    assert "chi = [" not in text
    re = read_aij(p)
    assert re.mode == "a"
    # a = 25 + 1/0.306 ≈ 28.268
    assert math.isclose(
        next(v for (i,j,v) in re.pairs if (i,j)==("A","B")),
        25.0 + 1.0 / 0.306, rel_tol=1e-6,
    )


def test_write_aij_invalid_out_mode(tmp_path):
    from abmptools.cg.dpd import AijMatrix, write_aij
    src = AijMatrix(segments=["A"], pairs=[("A","A",25.0)], mode="a")
    with pytest.raises(ValueError, match="out_mode must be"):
        write_aij(tmp_path / "x.dat", src, out_mode="invalid")


# --- validate() ------------------------------------------------------------

def test_validate_ok(cholesterol_cg, sample_aij_a_mode):
    """particle_names と aij.dat の segments が一致なら警告なし."""
    from abmptools.cg.dpd import CGDpdBuilder
    builder = CGDpdBuilder.from_files(
        monomer=cholesterol_cg["monomer"], aij=sample_aij_a_mode,
        calc_sett=cholesterol_cg["calc_sett"],
    )
    warnings = builder.validate()
    assert warnings == [], f"unexpected warnings: {warnings}"


def test_validate_missing_segments(cholesterol_cg, tmp_path):
    """monomer の particle 名が aij.dat にない場合は warning."""
    from abmptools.cg.dpd import AijMatrix, CGDpdBuilder, write_aij
    # P0, P1 だけの aij (cholesterol P2..P4 が欠落)
    aij = AijMatrix(
        segments=["P0","P1"],
        pairs=[("P0","P0",25.0),("P1","P1",25.0),("P0","P1",30.0)],
        mode="a", aii=25.0,
    )
    p = tmp_path / "subset.dat"
    write_aij(p, aij)
    builder = CGDpdBuilder.from_files(
        monomer=cholesterol_cg["monomer"], aij=p,
        calc_sett=cholesterol_cg["calc_sett"],
    )
    warnings = builder.validate()
    assert len(warnings) >= 1
    assert "NOT in aij.dat" in warnings[0]
    assert "P2" in warnings[0] and "P3" in warnings[0] and "P4" in warnings[0]


def test_validate_extra_aij_segments(cholesterol_cg, tmp_path):
    """aij.dat に monomer 未使用の segment があれば warning."""
    from abmptools.cg.dpd import AijMatrix, CGDpdBuilder, write_aij
    # P0..P4 に加えて X / Y が aij にだけある
    pairs = [(f"P{i}", f"P{i}", 25.0) for i in range(5)]
    pairs += [("X","X",25.0), ("Y","Y",25.0)]
    pairs += [("P0","X",30.0)]
    aij = AijMatrix(
        segments=["P0","P1","P2","P3","P4","X","Y"],
        pairs=pairs, mode="a", aii=25.0,
    )
    p = tmp_path / "extra.dat"
    write_aij(p, aij)
    builder = CGDpdBuilder.from_files(
        monomer=cholesterol_cg["monomer"], aij=p,
        calc_sett=cholesterol_cg["calc_sett"],
    )
    warnings = builder.validate()
    assert any("NOT used by any monomer" in w for w in warnings)
    assert any("X" in w and "Y" in w for w in warnings)


# --- CLI: verify subcommand ------------------------------------------------

def test_cli_verify_ok(tmp_path, cholesterol_cg, sample_aij_a_mode):
    """CLI verify: 整合性 OK 時は return 0 + OK メッセージ."""
    result = subprocess.run(
        [
            sys.executable, "-m", "abmptools.cg.dpd", "verify",
            "--monomer", str(cholesterol_cg["monomer"]),
            "--aij", str(sample_aij_a_mode),
            "--calc-sett", str(cholesterol_cg["calc_sett"]),
        ],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, f"stderr: {result.stderr}"
    assert "OK" in result.stdout


def test_cli_verify_warning(tmp_path, cholesterol_cg):
    """CLI verify: warning ありなら return 1."""
    from abmptools.cg.dpd import AijMatrix, write_aij
    aij = AijMatrix(
        segments=["P0"], pairs=[("P0","P0",25.0)], mode="a", aii=25.0,
    )
    p = tmp_path / "minimal.dat"
    write_aij(p, aij)
    result = subprocess.run(
        [
            sys.executable, "-m", "abmptools.cg.dpd", "verify",
            "--monomer", str(cholesterol_cg["monomer"]),
            "--aij", str(p),
            "--calc-sett", str(cholesterol_cg["calc_sett"]),
        ],
        capture_output=True, text=True,
    )
    assert result.returncode == 1
    assert "WARNING" in result.stdout
