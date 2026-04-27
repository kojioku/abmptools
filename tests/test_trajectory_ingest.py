# -*- coding: utf-8 -*-
"""
Tests for abmptools.amorphous.trajectory_ingest.frames_from_xtc.

MDAnalysis is not a hard dependency of abmptools, so these tests stub it
out and only verify the unit conversion + range-slicing contract.
"""
from __future__ import annotations

import sys
import types
from unittest.mock import patch

import numpy as np
import pytest

from abmptools.amorphous.trajectory_ingest import frames_from_xtc
from abmptools.gro2udf.top_model import GROFrameData


class _StubAtomGroup:
    def __init__(self, positions_A):
        self.positions = np.asarray(positions_A)
        # `make_whole` operates on fragments; for the stub we expose
        # the whole atom group as a single fragment.
        self.fragments = [self]


class _StubTimestep:
    def __init__(self, frame, time):
        self.frame = frame
        self.time = time


class _StubTrajectory:
    def __init__(self, frames):
        self._frames = frames

    def __iter__(self):
        return iter(self._frames)

    def __len__(self):
        return len(self._frames)


class _StubUniverse:
    """Replays two frames with known positions / cell dimensions."""

    def __init__(self, *_args, **_kw):
        # two frames; positions in Å for each atom
        self._positions_by_frame = [
            np.array([[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]]),
            np.array([[0.0, 0.0, 0.0], [20.0, 0.0, 0.0]]),
        ]
        self._current = 0
        self.trajectory = _StubTrajectory([
            _StubTimestep(frame=0, time=0.0),
            _StubTimestep(frame=1, time=1.0),
        ])
        self.dimensions = np.array([50.0, 50.0, 50.0, 90.0, 90.0, 90.0])
        self.atoms = _StubAtomGroup(self._positions_by_frame[0])
        self._added_bonds = None  # captured by add_bonds for inspection
        # Advance atoms.positions alongside trajectory iteration.
        self.trajectory = self._iterating(
            self.trajectory, self._positions_by_frame, self)

    @staticmethod
    def _iterating(inner, positions, universe):
        ts_list = list(inner)

        class _Wrapped:
            def __iter__(self_inner):
                for idx, ts in enumerate(ts_list):
                    universe.atoms = _StubAtomGroup(positions[idx])
                    yield ts

            def __len__(self_inner):
                return len(ts_list)

        return _Wrapped()

    def add_bonds(self, bond_list):
        """Record bond list so the test can verify ``top_path`` parsing."""
        self._added_bonds = list(bond_list)


def _patch_mda():
    mod = types.ModuleType("MDAnalysis")
    mod.Universe = _StubUniverse  # type: ignore[attr-defined]
    # MDAnalysis.lib.mdamath.make_whole stub (recorded for test inspection)
    lib = types.ModuleType("MDAnalysis.lib")
    mdamath = types.ModuleType("MDAnalysis.lib.mdamath")
    _called = {"make_whole": 0}

    def _make_whole(_fragment):
        _called["make_whole"] += 1

    mdamath.make_whole = _make_whole  # type: ignore[attr-defined]
    lib.mdamath = mdamath  # type: ignore[attr-defined]
    mod.lib = lib  # type: ignore[attr-defined]
    mod._make_whole_calls = _called  # type: ignore[attr-defined]
    return mod


def test_frames_from_xtc_converts_angstrom_to_nm():
    with patch.dict(sys.modules, {"MDAnalysis": _patch_mda()}):
        frames = frames_from_xtc("fake.tpr", "fake.xtc")
    assert len(frames) == 2
    # Å / 10 → nm
    assert frames[0].coord_list[1] == pytest.approx([1.0, 0.0, 0.0])
    assert frames[1].coord_list[1] == pytest.approx([2.0, 0.0, 0.0])
    assert frames[0].cell == pytest.approx([5.0, 5.0, 5.0])


def test_frames_from_xtc_slices_range():
    with patch.dict(sys.modules, {"MDAnalysis": _patch_mda()}):
        frames = frames_from_xtc("fake.tpr", "fake.xtc", start=1, end=1)
    assert len(frames) == 1
    assert frames[0].step == 1


def test_frames_from_xtc_returns_groframedata_instances():
    with patch.dict(sys.modules, {"MDAnalysis": _patch_mda()}):
        frames = frames_from_xtc("fake.tpr", "fake.xtc")
    for f in frames:
        assert isinstance(f, GROFrameData)


def test_frames_from_xtc_without_mdanalysis_raises():
    with patch.dict(sys.modules, {"MDAnalysis": None}):
        sys.modules.pop("MDAnalysis", None)
        orig_import = __builtins__["__import__"] if isinstance(
            __builtins__, dict) else __builtins__.__import__

        def _forbid(name, *a, **kw):
            if name == "MDAnalysis":
                raise ImportError("simulated missing MDAnalysis")
            return orig_import(name, *a, **kw)

        with patch("builtins.__import__", side_effect=_forbid):
            with pytest.raises(RuntimeError, match="MDAnalysis is required"):
                frames_from_xtc("fake.tpr", "fake.xtc")


# ---------------------------------------------------------------------------
# top_path → bonds → make_whole
# ---------------------------------------------------------------------------

def _write_minimal_top(tmp_path):
    """Tiny .top with 2 mol types (3- and 4-atom) and 2 instances each."""
    top = tmp_path / "system.top"
    top.write_text(
        "[ defaults ]\n"
        "; nbfunc        comb-rule\n"
        "1               2               yes             0.5     0.8333\n"
        "\n"
        "[ atomtypes ]\n"
        ";name  bond_type  mass  charge  ptype  sigma  epsilon\n"
        "c      c   12.011  0.0  A  0.34  0.36\n"
        "h      h    1.008  0.0  A  0.26  0.06\n"
        "\n"
        "[ moleculetype ]\n"
        "; Name nrexcl\n"
        "MolA   3\n"
        "\n"
        "[ atoms ]\n"
        ";   nr  type  resnr  residu  atom  cgnr  charge   mass\n"
        "    1   c     1      MolA    C1    1    0.0      12.011\n"
        "    2   h     1      MolA    H1    1    0.0       1.008\n"
        "    3   h     1      MolA    H2    1    0.0       1.008\n"
        "\n"
        "[ bonds ]\n"
        "; ai aj funct\n"
        "  1  2  1\n"
        "  1  3  1\n"
        "\n"
        "[ moleculetype ]\n"
        "; Name nrexcl\n"
        "MolB   3\n"
        "\n"
        "[ atoms ]\n"
        "    1   c     1      MolB    C1    1    0.0      12.011\n"
        "    2   c     1      MolB    C2    1    0.0      12.011\n"
        "    3   h     1      MolB    H1    1    0.0       1.008\n"
        "    4   h     1      MolB    H2    1    0.0       1.008\n"
        "\n"
        "[ bonds ]\n"
        "  1  2  1\n"
        "  1  3  1\n"
        "  2  4  1\n"
        "\n"
        "[ system ]\n"
        "Test\n"
        "\n"
        "[ molecules ]\n"
        "MolA   2\n"
        "MolB   1\n"
    )
    return str(top)


def test_bonds_from_top_offsets_correctly(tmp_path):
    """Verify _bonds_from_top accumulates atom offsets across instances."""
    from abmptools.amorphous.trajectory_ingest import _bonds_from_top

    top_path = _write_minimal_top(tmp_path)
    bonds = _bonds_from_top(top_path)
    # Expected (0-based):
    #   instance 0 (MolA, atoms 0..2):    (0,1), (0,2)
    #   instance 1 (MolA, atoms 3..5):    (3,4), (3,5)
    #   instance 2 (MolB, atoms 6..9):    (6,7), (6,8), (7,9)
    expected = [(0, 1), (0, 2), (3, 4), (3, 5), (6, 7), (6, 8), (7, 9)]
    assert bonds == expected


def test_frames_from_xtc_with_top_path_calls_make_whole(tmp_path):
    """When top_path is given, bonds get added and make_whole is invoked."""
    top_path = _write_minimal_top(tmp_path)
    mda_stub = _patch_mda()
    with patch.dict(sys.modules, {
        "MDAnalysis": mda_stub,
        "MDAnalysis.lib": mda_stub.lib,
        "MDAnalysis.lib.mdamath": mda_stub.lib.mdamath,
    }):
        frames = frames_from_xtc(
            "fake.gro", "fake.xtc", top_path=top_path,
        )
    assert len(frames) == 2
    # _StubUniverse exposes 1 fragment; make_whole called once per frame
    assert mda_stub._make_whole_calls["make_whole"] == 2


def test_frames_from_xtc_without_top_path_skips_make_whole(tmp_path):
    """When top_path is omitted, no bonds added and no make_whole calls."""
    mda_stub = _patch_mda()
    with patch.dict(sys.modules, {
        "MDAnalysis": mda_stub,
        "MDAnalysis.lib": mda_stub.lib,
        "MDAnalysis.lib.mdamath": mda_stub.lib.mdamath,
    }):
        frames = frames_from_xtc("fake.gro", "fake.xtc")
    assert len(frames) == 2
    assert mda_stub._make_whole_calls["make_whole"] == 0
