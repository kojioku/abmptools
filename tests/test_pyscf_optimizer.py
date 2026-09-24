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


# ---------------------------------------------------------------------------
# _resolve_eps
# ---------------------------------------------------------------------------

def test_resolve_eps_accepts_name_number_and_none():
    from abmptools.geomopt.pyscf_optimizer import _resolve_eps
    assert _resolve_eps(None) is None
    assert _resolve_eps("water") == pytest.approx(78.3553)
    assert _resolve_eps("WATER") == pytest.approx(78.3553)
    assert _resolve_eps(4.71) == pytest.approx(4.71)
    assert _resolve_eps("4.71") == pytest.approx(4.71)


def test_resolve_eps_rejects_unknown_name():
    from abmptools.geomopt.pyscf_optimizer import _resolve_eps
    with pytest.raises(ValueError):
        _resolve_eps("liquid nitrogen")


# ---------------------------------------------------------------------------
# charge / spin carried on the xyz comment line
# ---------------------------------------------------------------------------

def test_charge_spin_from_comment_reads_both():
    from abmptools.geomopt.pyscf_optimizer import _charge_spin_from_comment
    c, s = _charge_spin_from_comment(
        "@@E_diMe-Pho@@ C2H6O4P | charge=-1 | spin=0")
    assert (c, s) == (-1, 0)


def test_charge_spin_from_comment_handles_plus_and_absence():
    from abmptools.geomopt.pyscf_optimizer import _charge_spin_from_comment
    assert _charge_spin_from_comment("@@F_Choline@@ charge=+1")[0] == 1
    assert _charge_spin_from_comment("@@A_propane@@ C3H8") == (None, None)
    assert _charge_spin_from_comment("") == (None, None)


def test_charge_spin_from_comment_ignores_lookalike_words():
    from abmptools.geomopt.pyscf_optimizer import _charge_spin_from_comment
    # "recharge=" must not be read as "charge=" (\b guards the left edge).
    assert _charge_spin_from_comment("note: recharge=3")[0] is None


# ---------------------------------------------------------------------------
# _apply_solvent
# ---------------------------------------------------------------------------

class _FakeSolvent:
    def __init__(self):
        self.eps = None
        self.method = None


class _FakeMF:
    """Stands in for a PySCF mean-field object that supports solvent models."""

    def __init__(self, supports=("PCM", "SMD", "ddCOSMO")):
        self._supports = supports

    def _make(self, tag):
        if tag not in self._supports:
            raise AttributeError(tag)
        out = _FakeMF(self._supports)
        out.with_solvent = _FakeSolvent()
        out.tag = tag
        return out

    def PCM(self):
        return self._make("PCM")

    def SMD(self):
        return self._make("SMD")

    def ddCOSMO(self):
        return self._make("ddCOSMO")


def test_solvent_none_is_a_passthrough():
    opt = QMOptimizerPySCF(solvent="none")
    mf = _FakeMF()
    assert opt._apply_solvent(mf) is mf
    assert opt.solvent_applied is False


def test_solvent_cpcm_sets_variant_and_eps():
    opt = QMOptimizerPySCF(solvent="cpcm", solvent_eps="water")
    out = opt._apply_solvent(_FakeMF())
    assert opt.solvent_applied is True
    assert out.tag == "PCM"
    assert out.with_solvent.method == "C-PCM"
    assert out.with_solvent.eps == pytest.approx(78.3553)


def test_solvent_plain_pcm_keeps_pyscf_default_variant():
    opt = QMOptimizerPySCF(solvent="pcm", solvent_eps=4.71)
    out = opt._apply_solvent(_FakeMF())
    assert out.with_solvent.method is None
    assert out.with_solvent.eps == pytest.approx(4.71)


def test_solvent_smd_and_ddcosmo_route_to_their_own_builders():
    for name, tag in (("smd", "SMD"), ("ddcosmo", "ddCOSMO")):
        opt = QMOptimizerPySCF(solvent=name)
        assert opt._apply_solvent(_FakeMF()).tag == tag


def test_solvent_unavailable_falls_back_to_gas_phase():
    opt = QMOptimizerPySCF(solvent="smd")
    mf = _FakeMF(supports=("PCM",))       # no SMD in this build
    assert opt._apply_solvent(mf) is mf
    assert opt.solvent_applied is False


def test_unknown_solvent_is_rejected_at_construction():
    with pytest.raises(ValueError):
        QMOptimizerPySCF(solvent="handwavium")


def test_level_of_theory_reports_the_solvent_that_attached():
    opt = QMOptimizerPySCF(functional="B3LYP", basis="def2-SVP",
                           dispersion="d3bj", solvent="cpcm",
                           solvent_eps="water")
    opt.dispersion_applied = True
    # Not attached yet -> the label must not claim a solvent.
    assert opt._level_of_theory() == "B3LYP-D3BJ/def2-SVP"
    opt._apply_solvent(_FakeMF())
    assert opt._level_of_theory() == "B3LYP-D3BJ/def2-SVP [CPCM,eps=78.3553]"
