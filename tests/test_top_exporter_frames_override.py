# -*- coding: utf-8 -*-
"""
Test that ``TopExporter.export_model(..., frames=...)`` uses the frame
override instead of ``model.frames``.

UDFManager is not a hard dependency in CI, so this test stubs it and
verifies only the frame-selection contract (not UDF text output).
"""
from __future__ import annotations

import sys
import types
from unittest.mock import patch

import pytest

from abmptools.gro2udf.top_model import GROFrameData, TopModel


class _StubUDFManager:
    """Minimal stand-in that swallows every put / get call."""

    def __init__(self, path):
        self._path = path
        self._records = 0

    def totalRecord(self):
        return 0

    def eraseRecord(self, _a, _b):
        pass

    def jump(self, _rec):
        pass

    def newRecord(self, *_a, **_kw):
        self._records += 1

    def put(self, *_a, **_kw):
        pass

    def get(self, *_a, **_kw):
        return 0

    def size(self, *_a, **_kw):
        return 0

    def write(self, *_a, **_kw):
        pass


def _patch_udfmanager():
    mod = types.ModuleType("UDFManager")
    mod.UDFManager = _StubUDFManager  # type: ignore[attr-defined]
    return mod


def _minimal_topmodel(builtin_frames):
    return TopModel(
        comb_rule=2,
        fudge_lj=0.5,
        fudge_qq=0.8333,
        atom_type_specs=[],
        bond_type_specs=[],
        angle_type_specs=[],
        torsion_type_specs=[],
        mass_dict={},
        mol_type_names=[],
        mol_specs=[],
        mol_instance_list=[],
        frames=builtin_frames,
    )


def _collect_frames_written(model, frames_arg, tmp_path):
    from abmptools.gro2udf.top_exporter import TopExporter

    exporter = TopExporter()
    written = []

    def _capture(self_, uobj, model_, frame):
        written.append(frame)

    template = tmp_path / "tmpl.udf"
    template.write_text("")
    out = tmp_path / "out.udf"

    with patch.dict(sys.modules, {"UDFManager": _patch_udfmanager()}), \
         patch.object(TopExporter, "_append_structure", _capture), \
         patch.object(TopExporter, "_write_set_of_molecules", lambda *a, **k: None), \
         patch.object(TopExporter, "_set_default_condition", lambda *a, **k: None), \
         patch.object(TopExporter, "_write_molecular_attributes", lambda *a, **k: None), \
         patch.object(TopExporter, "_write_interactions", lambda *a, **k: None), \
         patch("abmptools.gro2udf.top_exporter.shutil.copy", lambda *a, **k: None):
        exporter.export_model(model, str(template), str(out), frames=frames_arg)
    return written


def test_export_model_uses_model_frames_by_default(tmp_path):
    f0 = GROFrameData(step=0, time=0.0, coord_list=[[0, 0, 0]], cell=[1, 1, 1])
    f1 = GROFrameData(step=1, time=1.0, coord_list=[[1, 1, 1]], cell=[1, 1, 1])
    model = _minimal_topmodel([f0, f1])
    written = _collect_frames_written(model, frames_arg=None, tmp_path=tmp_path)
    assert written == [f0, f1]


def test_export_model_overrides_with_explicit_frames(tmp_path):
    f_default = GROFrameData(
        step=0, time=0.0, coord_list=[[0, 0, 0]], cell=[1, 1, 1])
    f_override1 = GROFrameData(
        step=10, time=10.0, coord_list=[[5, 5, 5]], cell=[2, 2, 2])
    f_override2 = GROFrameData(
        step=20, time=20.0, coord_list=[[6, 6, 6]], cell=[2, 2, 2])
    model = _minimal_topmodel([f_default])
    written = _collect_frames_written(
        model, frames_arg=[f_override1, f_override2], tmp_path=tmp_path)
    assert written == [f_override1, f_override2]
    # model.frames left untouched
    assert model.frames == [f_default]


def test_export_model_empty_frame_list_writes_no_structure(tmp_path):
    model = _minimal_topmodel([GROFrameData(
        step=0, time=0.0, coord_list=[[0, 0, 0]], cell=[1, 1, 1])])
    written = _collect_frames_written(model, frames_arg=[], tmp_path=tmp_path)
    assert written == []
