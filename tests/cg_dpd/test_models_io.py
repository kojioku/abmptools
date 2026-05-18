# -*- coding: utf-8 -*-
"""tests/cg_dpd/test_models_io.py — データクラスと I/O parser (aij / monomer / calc_sett)."""
from __future__ import annotations

import math
from pathlib import Path

import pytest


def test_imports():
    from abmptools.cg.dpd import (
        AijMatrix, CalcSett, DpdSystemSpec, MonomerSpec,
        read_aij, write_aij, aij_to_dict,
        read_monomer, read_monomers_dict, assign_particle_names,
        read_calc_sett,
        patch_dpm, propagate_virtual_mom, write_message_txt,
        write_dpd_udf, CGDpdBuilder,
    )


# --- AijMatrix --------------------------------------------------------------

def test_aij_matrix_a_mode_to_a_values():
    from abmptools.cg.dpd import AijMatrix
    aij = AijMatrix(
        segments=["A", "B"], pairs=[("A","B", 30.0)], mode="a", aii=25.0,
    )
    vals = aij.to_a_values()
    assert vals == [("A","B", 30.0)]


def test_aij_matrix_chi_mode_converts():
    """chi=-4.108, aii=25.0 → a = 25 + (-4.108)/0.306 ≈ 11.58."""
    from abmptools.cg.dpd import AijMatrix
    aij = AijMatrix(
        segments=["A", "B"], pairs=[("A","B", -4.108)], mode="chi", aii=25.0,
    )
    vals = aij.to_a_values()
    assert len(vals) == 1
    assert math.isclose(vals[0][2], 25.0 - 4.108 / 0.306, rel_tol=1e-6)


def test_aij_matrix_get_symmetric():
    from abmptools.cg.dpd import AijMatrix
    aij = AijMatrix(
        segments=["A", "B"], pairs=[("A","B", 30.0)], mode="a",
    )
    assert aij.get("A", "B") == 30.0
    assert aij.get("B", "A") == 30.0  # symmetric lookup
    assert aij.get("A", "C") is None


def test_aij_matrix_invalid_mode():
    from abmptools.cg.dpd import AijMatrix
    with pytest.raises(ValueError):
        AijMatrix(segments=["A"], pairs=[], mode="invalid")


# --- aij_io: round-trip -----------------------------------------------------

def test_aij_io_round_trip_a_mode(sample_aij_a_mode):
    from abmptools.cg.dpd import read_aij
    aij = read_aij(sample_aij_a_mode)
    assert aij.mode == "a"
    assert len(aij.pairs) == 15
    assert ("P0","P0",25.0) in aij.pairs
    assert ("P0","P4",40.0) in aij.pairs


def test_aij_io_round_trip_chi_mode(sample_aij_chi_mode):
    from abmptools.cg.dpd import read_aij
    aij = read_aij(sample_aij_chi_mode)
    assert aij.mode == "chi"
    assert len(aij.pairs) == 3
    assert aij.aii == 25.0
    # chi → a 変換
    a_vals = aij.to_a_values()
    assert math.isclose(a_vals[0][2], 25.0 - 4.108 / 0.306, rel_tol=1e-6)


def test_aij_io_missing_aij_chi(tmp_path):
    from abmptools.cg.dpd import read_aij
    p = tmp_path / "empty.dat"
    p.write_text("# no aij or chi defined\n")
    with pytest.raises(ValueError, match="neither 'aij' nor 'chi'"):
        read_aij(p)


def test_aij_io_missing_file(tmp_path):
    from abmptools.cg.dpd import read_aij
    with pytest.raises(FileNotFoundError):
        read_aij(tmp_path / "nonexistent.dat")


def test_aij_to_dict_symmetric(sample_aij_a_mode):
    from abmptools.cg.dpd import read_aij, aij_to_dict
    d = aij_to_dict(read_aij(sample_aij_a_mode))
    assert d[("P0","P4")] == 40.0
    assert d[("P4","P0")] == 40.0  # symmetric


# --- monomer_io -------------------------------------------------------------

def test_monomer_io_cholesterol(cholesterol_cg):
    from abmptools.cg.dpd import read_monomer
    mono = read_monomer(cholesterol_cg["monomer"])
    assert mono.name == "chol"
    assert mono.n_particles == 5
    assert len(mono.bond12) == 4
    assert len(mono.bond13_150) == 2
    assert len(mono.bond14_150) == 1
    assert len(mono.angle13) == 3
    # particle_names は P0..P4 仮ラベル
    assert mono.particle_names == ["P0","P1","P2","P3","P4"]


def test_monomer_io_assign_particle_names(cholesterol_cg):
    from abmptools.cg.dpd import read_monomer, assign_particle_names
    mono = read_monomer(cholesterol_cg["monomer"])
    new_names = ["A","B","C","D","T"]
    renamed = assign_particle_names(mono, new_names)
    assert renamed.particle_names == new_names
    # 元 spec は変更なし (immutable replace)
    assert mono.particle_names == ["P0","P1","P2","P3","P4"]


def test_monomer_io_wrong_length(cholesterol_cg):
    from abmptools.cg.dpd import read_monomer, assign_particle_names
    mono = read_monomer(cholesterol_cg["monomer"])
    with pytest.raises(ValueError, match="name list length"):
        assign_particle_names(mono, ["A","B"])  # length mismatch


# --- calc_sett_io ----------------------------------------------------------

def test_calc_sett_io_cholesterol(cholesterol_cg):
    from abmptools.cg.dpd import read_calc_sett
    sett = read_calc_sett(cholesterol_cg["calc_sett"])
    assert sett.monomer_file == "chol_monomer"
    assert sett.aij_file == "aij.dat"
    assert sett.box_size == (12.0, 12.0, 12.0)
    assert sett.density == 3.0
    assert sett.dt == 0.05
    assert sett.aii_val == 25.0
