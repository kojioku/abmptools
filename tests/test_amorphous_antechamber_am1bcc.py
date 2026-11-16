# -*- coding: utf-8 -*-
"""AM1-BCC on Windows through antechamber.exe (OpenFF cannot start antechamber.bat and its
``-pf yes`` runs ``rm``): finding the program, the commands, and that the builder then keeps
the charges (use_precomputed_charges)."""
from __future__ import annotations

import os
import shutil
import sys
from pathlib import Path

import pytest

from abmptools.amorphous import antechamber_am1bcc as W
from abmptools.amorphous import molecule_prep as MP


def _windows_build(root: Path) -> Path:
    (root / "bin" / "wrapped_progs").mkdir(parents=True)
    (root / "bin" / "antechamber.bat").write_text("@echo off\r\n", encoding="ascii")
    exe = root / "bin" / "wrapped_progs" / "antechamber.exe"
    exe.write_bytes(b"MZ")
    return exe


def test_the_exe_behind_amberhome_is_started_with_amberhome_and_bin_on_path(tmp_path, monkeypatch):
    exe = _windows_build(tmp_path / "amber")
    monkeypatch.setenv("AMBERHOME", str(tmp_path / "amber"))
    argv, env = W.antechamber_command()
    assert argv == [str(exe)]
    assert env["AMBERHOME"] == str(tmp_path / "amber")
    assert env["PATH"].split(os.pathsep)[0] == str(tmp_path / "amber" / "bin")


def test_the_exe_beside_a_bat_on_path_is_found(tmp_path, monkeypatch):
    exe = _windows_build(tmp_path / "amber")
    monkeypatch.delenv("AMBERHOME", raising=False)
    monkeypatch.setattr(W.shutil, "which",
                        lambda name: str(tmp_path / "amber" / "bin" / "antechamber.bat"))
    argv, _ = W.antechamber_command()
    assert argv == [str(exe)]


def test_nothing_found_is_none(tmp_path, monkeypatch):
    monkeypatch.delenv("AMBERHOME", raising=False)
    monkeypatch.setattr(W.shutil, "which", lambda name: None)
    assert W.antechamber_command() is None


def test_which_methods_leave_charges_on_the_molecule(monkeypatch):
    monkeypatch.setattr(W, "needed", lambda: False)                 # Linux / macOS
    assert MP.charges_preassigned("nagl") and MP.charges_preassigned("gasteiger")
    assert not MP.charges_preassigned("") and not MP.charges_preassigned("am1bcc")
    monkeypatch.setattr(W, "needed", lambda: True)                  # Windows
    assert MP.charges_preassigned("") and MP.charges_preassigned("am1bcc")


def test_on_windows_am1bcc_is_assigned_before_interchange(monkeypatch):
    seen = []
    monkeypatch.setattr(W, "needed", lambda: True)
    monkeypatch.setattr(W, "assign_am1bcc", lambda mol: seen.append(mol))
    MP._maybe_assign_partial_charges("MOL", "", "model")
    MP._maybe_assign_partial_charges("MOL2", "am1bcc", "model")
    assert seen == ["MOL", "MOL2"]
    monkeypatch.setattr(W, "needed", lambda: False)
    MP._maybe_assign_partial_charges("MOL3", "am1bcc", "model")
    assert seen == ["MOL", "MOL2"]


@pytest.mark.skipif(shutil.which("antechamber") is None, reason="AmberTools not on PATH")
def test_the_charges_equal_openffs_own_am1bcc():
    """Same conformer generation and the same two antechamber commands: same charges."""
    np = pytest.importorskip("numpy")
    off = pytest.importorskip("openff.toolkit")
    from openff.toolkit.utils.ambertools_wrapper import AmberToolsToolkitWrapper
    smi = "CCCCCCCCCC(=O)[O-]"
    a = off.Molecule.from_smiles(smi)
    b = off.Molecule.from_smiles(smi)
    a.assign_partial_charges("am1bcc", toolkit_registry=AmberToolsToolkitWrapper())
    W.assign_am1bcc(b)
    assert np.abs(a.partial_charges.m - b.partial_charges.m).max() == 0.0
    assert abs(b.partial_charges.m.sum() + 1.0) < 1e-9
