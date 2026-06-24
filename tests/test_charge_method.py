# -*- coding: utf-8 -*-
"""
Tests for the ``charge_method`` pre-assignment path added to enable the
Windows-native (openff-nagl) route.

`assign_partial_charges` is stubbed so the tests don't need openff-nagl
or AmberTools installed — we only verify the dispatch logic.
"""
from __future__ import annotations

from types import SimpleNamespace

import pytest

from abmptools.amorphous.molecule_prep import _maybe_assign_partial_charges


class _FakeMol:
    def __init__(self):
        self.calls = []

    def assign_partial_charges(self, partial_charge_method, toolkit_registry=None):  # noqa: D401
        # nagl 経路は toolkit_registry=NAGLToolkitWrapper() を明示で渡す
        self.calls.append(partial_charge_method)


def _patch_nagl_wrapper(monkeypatch):
    """openff-nagl 非導入環境でも nagl 経路を試せるよう、 NAGLToolkitWrapper を
    instantiation 不要のダミーに差し替える (実 NAGL は CI windows ジョブで検証)。"""
    import openff.toolkit.utils.nagl_wrapper as nw
    monkeypatch.setattr(nw, "NAGLToolkitWrapper", lambda: object())


def test_am1bcc_and_empty_are_noops():
    """AM1-BCC is handled by Interchange later, so we must not touch the molecule here."""
    mol = _FakeMol()
    _maybe_assign_partial_charges(mol, "", "any_model")
    _maybe_assign_partial_charges(mol, "am1bcc", "any_model")
    assert mol.calls == []


def test_nagl_uses_given_model_name(monkeypatch):
    _patch_nagl_wrapper(monkeypatch)
    mol = _FakeMol()
    _maybe_assign_partial_charges(mol, "nagl", "openff-gnn-am1bcc-0.1.0-rc.3.pt")
    assert mol.calls == ["openff-gnn-am1bcc-0.1.0-rc.3.pt"]


def test_gasteiger_uses_openff_method_name():
    mol = _FakeMol()
    _maybe_assign_partial_charges(mol, "gasteiger", "ignored")
    assert mol.calls == ["gasteiger"]


def test_unknown_method_raises_value_error():
    mol = _FakeMol()
    with pytest.raises(ValueError, match="Unknown charge_method"):
        _maybe_assign_partial_charges(mol, "mulliken", "ignored")


def test_nagl_backend_failure_is_wrapped_with_install_hint(monkeypatch):
    _patch_nagl_wrapper(monkeypatch)

    class _BrokenMol:
        def assign_partial_charges(self, partial_charge_method, toolkit_registry=None):
            raise RuntimeError("underlying toolkit missing")

    with pytest.raises(RuntimeError, match="openff-nagl"):
        _maybe_assign_partial_charges(_BrokenMol(), "nagl", "some_model.pt")


def test_build_config_defaults_keep_am1bcc():
    """BuildConfig default must preserve existing behaviour (AM1-BCC via Interchange)."""
    from abmptools.amorphous.models import BuildConfig, ComponentSpec

    cfg = BuildConfig(components=[ComponentSpec(smiles="C")])
    assert cfg.charge_method == ""
    assert cfg.nagl_model == "openff-gnn-am1bcc-0.1.0-rc.3.pt"


def test_build_config_json_roundtrip_preserves_charge_method(tmp_path):
    from abmptools.amorphous.models import BuildConfig, ComponentSpec

    cfg = BuildConfig(
        components=[ComponentSpec(smiles="C", n_mol=10)],
        charge_method="nagl",
        nagl_model="custom.pt",
    )
    path = tmp_path / "cfg.json"
    cfg.to_json(str(path))
    restored = BuildConfig.from_json(str(path))
    assert restored.charge_method == "nagl"
    assert restored.nagl_model == "custom.pt"
