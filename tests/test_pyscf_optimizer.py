# -*- coding: utf-8 -*-
"""Tests for QMOptimizerPySCF's dispersion handling and provenance.

These cover the two ways a D3 run can quietly turn into a plain DFT run:
the provider is looked up under a name it does not export, and the written
xyz still claims the dispersion-corrected level of theory afterwards.
"""
import sys
import types
from pathlib import Path

import pytest

from abmptools.geomopt.pyscf_optimizer import (
    QMOptimizerPySCF,
    _read_xyz_comment,
)


# ---------------------------------------------------------------------------
# _read_xyz_comment
# ---------------------------------------------------------------------------

def test_read_xyz_comment_returns_second_line(tmp_path):
    f = tmp_path / "m.xyz"
    f.write_text(
        "1\n"
        "@@S_CCl4@@ ClC(Cl)(Cl)Cl | delta=17.8 V=97.1\n"
        "C 0.0 0.0 0.0\n"
    )
    assert _read_xyz_comment(f) == "@@S_CCl4@@ ClC(Cl)(Cl)Cl | delta=17.8 V=97.1"


def test_read_xyz_comment_non_xyz_returns_empty(tmp_path):
    f = tmp_path / "m.pdb"
    f.write_text("REMARK something\nATOM ...\n")
    assert _read_xyz_comment(f) == ""


def test_read_xyz_comment_short_file_returns_empty(tmp_path):
    f = tmp_path / "m.xyz"
    f.write_text("0\n")
    assert _read_xyz_comment(f) == ""


def test_read_xyz_comment_missing_file_returns_empty(tmp_path):
    assert _read_xyz_comment(tmp_path / "absent.xyz") == ""


# ---------------------------------------------------------------------------
# _apply_dispersion
# ---------------------------------------------------------------------------

def _fake_sdftd3(recorder):
    """Build a stand-in for the ``dftd3.pyscf`` module of simple-dftd3.

    It exports only what the real package exports -- notably ``energy`` and
    *not* ``DFTD3Model``, which is the shape that used to defeat the lookup.
    """
    mod = types.ModuleType("dftd3.pyscf")

    def energy(mf, method=None, version=None, **kw):
        recorder.update(method=method, version=version)
        mf.with_dftd3 = object()
        return mf

    mod.energy = energy
    return mod


def test_dispersion_none_is_not_applied():
    opt = QMOptimizerPySCF(dispersion="none")
    mf = object()
    assert opt._apply_dispersion(mf) is mf
    assert opt.dispersion_applied is False


def test_dispersion_uses_simple_dftd3_energy_entry_point(monkeypatch):
    rec = {}
    pkg = types.ModuleType("dftd3")
    sub = _fake_sdftd3(rec)
    pkg.pyscf = sub
    monkeypatch.setitem(sys.modules, "dftd3", pkg)
    monkeypatch.setitem(sys.modules, "dftd3.pyscf", sub)

    opt = QMOptimizerPySCF(functional="B3LYP", dispersion="d3bj")
    mf = types.SimpleNamespace()
    out = opt._apply_dispersion(mf)

    assert opt.dispersion_applied is True
    assert hasattr(out, "with_dftd3")
    assert rec == {"method": "B3LYP", "version": "d3bj"}


def test_dispersion_d3_maps_to_zero_damping(monkeypatch):
    rec = {}
    sub = _fake_sdftd3(rec)
    pkg = types.ModuleType("dftd3")
    pkg.pyscf = sub
    monkeypatch.setitem(sys.modules, "dftd3", pkg)
    monkeypatch.setitem(sys.modules, "dftd3.pyscf", sub)

    QMOptimizerPySCF(dispersion="d3")._apply_dispersion(types.SimpleNamespace())
    assert rec["version"] == "d3zero"


def test_dispersion_absent_provider_leaves_flag_false(monkeypatch):
    # Neither provider importable: the call must still return the plain mf.
    for name in ("dftd3", "dftd3.pyscf", "pyscf.dftd3"):
        monkeypatch.setitem(sys.modules, name, None)

    opt = QMOptimizerPySCF(dispersion="d3bj")
    mf = types.SimpleNamespace()
    assert opt._apply_dispersion(mf) is mf
    assert opt.dispersion_applied is False


# ---------------------------------------------------------------------------
# _level_of_theory -- what the written file will claim
# ---------------------------------------------------------------------------

def test_level_of_theory_names_d3_only_when_applied():
    opt = QMOptimizerPySCF(functional="B3LYP", basis="def2-SVP",
                           dispersion="d3bj")
    # Before any run, and after a run where the provider was missing.
    assert opt._level_of_theory() == "B3LYP/def2-SVP"
    opt.dispersion_applied = True
    assert opt._level_of_theory() == "B3LYP-D3BJ/def2-SVP"


def test_level_of_theory_plain_for_dispersion_none():
    opt = QMOptimizerPySCF(functional="PBE0", basis="def2-TZVP",
                           dispersion="none")
    opt.dispersion_applied = True  # must be ignored for "none"
    assert opt._level_of_theory() == "PBE0/def2-TZVP"
