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


def _patch_mda():
    mod = types.ModuleType("MDAnalysis")
    mod.Universe = _StubUniverse  # type: ignore[attr-defined]
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
