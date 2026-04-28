# -*- coding: utf-8 -*-
"""
Tests for abmptools.amorphous.parameterizer.create_interchange.

OpenFF Toolkit と Interchange は heavy dependency なので stub で
forcefield_name の単一/複数受付を検証する。
"""
from __future__ import annotations

import sys
import types
from unittest.mock import MagicMock, patch

import pytest


def _patch_openff(captured):
    """Patch openff.toolkit + openff.interchange + openff.units + numpy."""

    # ------ openff.toolkit ------
    openff_pkg = types.ModuleType("openff")
    toolkit_mod = types.ModuleType("openff.toolkit")

    class _ForceField:
        def __init__(self, *names):
            captured["ff_args"] = names

    class _Topology:
        @classmethod
        def from_pdb(cls, pdb_path, unique_molecules):
            captured["pdb"] = pdb_path
            t = cls()
            t.box_vectors = None
            return t

    class _Molecule:
        pass

    toolkit_mod.ForceField = _ForceField
    toolkit_mod.Topology = _Topology
    toolkit_mod.Molecule = _Molecule

    # ------ openff.interchange ------
    interchange_mod = types.ModuleType("openff.interchange")

    class _Interchange:
        @classmethod
        def from_smirnoff(cls, force_field, topology):
            captured["force_field"] = force_field
            captured["topology"] = topology
            return cls()

    interchange_mod.Interchange = _Interchange

    # ------ openff.units ------
    units_mod = types.ModuleType("openff.units")

    class _Unit:
        def __getattr__(self, name):
            # numpy 乗算と互換にするため、unit attribute は数値 1.0 を返す
            # (実 OpenFF unit は Pint Quantity だがここではスケールのみ)
            return 1.0

    units_mod.unit = _Unit()

    openff_pkg.toolkit = toolkit_mod
    openff_pkg.interchange = interchange_mod
    openff_pkg.units = units_mod

    # ------ openmm (just exists) ------
    openmm_mod = types.ModuleType("openmm")

    return {
        "openff": openff_pkg,
        "openff.toolkit": toolkit_mod,
        "openff.interchange": interchange_mod,
        "openff.units": units_mod,
        "openmm": openmm_mod,
    }


@pytest.fixture
def patched_openff():
    captured = {}
    mods = _patch_openff(captured)
    with patch.dict(sys.modules, mods):
        yield captured


def test_create_interchange_default_forcefield_passes_single_name(
    patched_openff,
):
    """default では openff_unconstrained-2.1.0.offxml の単一 FF が渡る。"""
    from abmptools.amorphous.parameterizer import create_interchange

    create_interchange(
        molecules=[],
        counts=[1],
        box_size_nm=2.5,
        mixture_pdb="/tmp/x.pdb",
    )
    assert patched_openff["ff_args"] == (
        "openff_unconstrained-2.1.0.offxml",
    )


def test_create_interchange_string_forcefield(patched_openff):
    """str を渡すと 1 引数として ForceField に渡る。"""
    from abmptools.amorphous.parameterizer import create_interchange

    create_interchange(
        molecules=[],
        counts=[1],
        box_size_nm=2.5,
        mixture_pdb="/tmp/x.pdb",
        forcefield_name="openff-2.0.0.offxml",
    )
    assert patched_openff["ff_args"] == ("openff-2.0.0.offxml",)


def test_create_interchange_list_forcefield_for_water_override(
    patched_openff,
):
    """list[str] を渡すと stacked FF (water override) として渡る。

    典型: openff organic + tip3p で純 water 系の repulsive σ/ε 問題
    (1.0→0.26 g/cm³) を回避できる。
    """
    from abmptools.amorphous.parameterizer import create_interchange

    create_interchange(
        molecules=[],
        counts=[1, 1],
        box_size_nm=2.5,
        mixture_pdb="/tmp/x.pdb",
        forcefield_name=[
            "openff_unconstrained-2.1.0.offxml",
            "tip3p.offxml",
        ],
    )
    assert patched_openff["ff_args"] == (
        "openff_unconstrained-2.1.0.offxml",
        "tip3p.offxml",
    )


def test_create_interchange_tuple_forcefield(patched_openff):
    """tuple も list と同じく stacked FF として扱う。"""
    from abmptools.amorphous.parameterizer import create_interchange

    create_interchange(
        molecules=[],
        counts=[1],
        box_size_nm=2.5,
        mixture_pdb="/tmp/x.pdb",
        forcefield_name=("openff-2.0.0.offxml", "spce.offxml"),
    )
    assert patched_openff["ff_args"] == (
        "openff-2.0.0.offxml",
        "spce.offxml",
    )


def test_create_interchange_empty_list_raises(patched_openff):
    """空 list を渡すと ValueError。"""
    from abmptools.amorphous.parameterizer import create_interchange

    with pytest.raises(ValueError, match="non-empty"):
        create_interchange(
            molecules=[],
            counts=[1],
            box_size_nm=2.5,
            mixture_pdb="/tmp/x.pdb",
            forcefield_name=[],
        )
