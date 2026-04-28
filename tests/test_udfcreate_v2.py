# -*- coding: utf-8 -*-
"""Tests for abmptools.udfcreate_v2 (Phase 2c-E clean writer scaffold).

UDFManager isn't a hard dependency, so we stub it for these tests and
verify the put-call sequence rather than running real OCTA tooling.
"""
from __future__ import annotations

import os
import sys
import types
from unittest.mock import patch

import pytest


class _StubUDFManager:
    """Records put calls so tests can assert on what got written."""
    def __init__(self, path):
        self.path = path
        self.calls = []
        self.jumped = None
        self.written = False
        # tests can read these back to verify
        self.values = {}

    def jump(self, rec):
        self.jumped = rec

    def put(self, value, location, *args):
        self.calls.append((value, location, args))
        self.values[location] = value

    def get(self, location, *args):
        return self.values.get(location)

    def write(self):
        self.written = True


@pytest.fixture
def patched_udfmanager():
    mod = types.ModuleType("UDFManager")
    mod.UDFManager = _StubUDFManager  # type: ignore[attr-defined]
    with patch.dict(sys.modules, {"UDFManager": mod}):
        yield mod


# ---------------------------------------------------------------------------
# set_simulation_time
# ---------------------------------------------------------------------------

def test_set_simulation_time_writes_three_paths():
    from abmptools.udfcreate_v2 import set_simulation_time
    u = _StubUDFManager("/tmp/x.udf")
    set_simulation_time(u, dt=0.001, totalstep=20000, outstep=200)
    locs = [c[1] for c in u.calls]
    assert "Simulation_Conditions.Dynamics_Conditions.Time.delta_T" in locs
    assert "Simulation_Conditions.Dynamics_Conditions.Time.Total_Steps" in locs
    assert ("Simulation_Conditions.Dynamics_Conditions.Time."
            "Output_Interval_Steps") in locs
    assert u.values[
        "Simulation_Conditions.Dynamics_Conditions.Time.delta_T"
    ] == 0.001
    assert u.values[
        "Simulation_Conditions.Dynamics_Conditions.Time.Total_Steps"
    ] == 20000
    assert u.values[
        "Simulation_Conditions.Dynamics_Conditions.Time.Output_Interval_Steps"
    ] == 200


def test_set_simulation_time_passes_unit_for_dt():
    """delta_T is written with [ps] tag (UDFManager unit-aware put)."""
    from abmptools.udfcreate_v2 import set_simulation_time
    u = _StubUDFManager("/tmp/x.udf")
    set_simulation_time(u, dt=0.002, totalstep=1, outstep=1)
    dt_call = next(
        c for c in u.calls
        if c[1] == "Simulation_Conditions.Dynamics_Conditions.Time.delta_T"
    )
    assert dt_call[2] == ("[ps]",)


def test_set_simulation_time_coerces_int_steps():
    from abmptools.udfcreate_v2 import set_simulation_time
    u = _StubUDFManager("/tmp/x.udf")
    # pass floats — function should coerce to int
    set_simulation_time(u, dt=0.001, totalstep=20000.0, outstep=200.0)
    assert u.values[
        "Simulation_Conditions.Dynamics_Conditions.Time.Total_Steps"
    ] == 20000
    assert isinstance(u.values[
        "Simulation_Conditions.Dynamics_Conditions.Time.Total_Steps"
    ], int)


# ---------------------------------------------------------------------------
# set_initial_cell
# ---------------------------------------------------------------------------

def test_set_initial_cell_writes_six_paths():
    from abmptools.udfcreate_v2 import set_initial_cell
    u = _StubUDFManager("/tmp/x.udf")
    set_initial_cell(u, cellsize_nm=[3.0, 3.0, 3.0])
    base = "Initial_Structure.Initial_Unit_Cell.Cell_Size"
    for axis in ("a", "b", "c"):
        assert u.values[f"{base}.{axis}"] == 3.0
    for ang in ("alpha", "beta", "gamma"):
        assert u.values[f"{base}.{ang}"] == 90.0


def test_set_initial_cell_supports_anisotropic():
    from abmptools.udfcreate_v2 import set_initial_cell
    u = _StubUDFManager("/tmp/x.udf")
    set_initial_cell(u, cellsize_nm=[2.5, 3.0, 3.5])
    base = "Initial_Structure.Initial_Unit_Cell.Cell_Size"
    assert u.values[f"{base}.a"] == 2.5
    assert u.values[f"{base}.b"] == 3.0
    assert u.values[f"{base}.c"] == 3.5


def test_set_initial_cell_supports_custom_angles():
    from abmptools.udfcreate_v2 import set_initial_cell
    u = _StubUDFManager("/tmp/x.udf")
    set_initial_cell(u, cellsize_nm=[3, 3, 3], angles_deg=[60.0, 70.0, 80.0])
    base = "Initial_Structure.Initial_Unit_Cell.Cell_Size"
    assert u.values[f"{base}.alpha"] == 60.0
    assert u.values[f"{base}.beta"] == 70.0
    assert u.values[f"{base}.gamma"] == 80.0


def test_set_initial_cell_rejects_wrong_shape():
    from abmptools.udfcreate_v2 import set_initial_cell
    u = _StubUDFManager("/tmp/x.udf")
    with pytest.raises(ValueError, match="cellsize_nm"):
        set_initial_cell(u, cellsize_nm=[3.0, 3.0])
    with pytest.raises(ValueError, match="angles_deg"):
        set_initial_cell(u, cellsize_nm=[3, 3, 3], angles_deg=[90, 90])


# ---------------------------------------------------------------------------
# write_skeleton_udf — end-to-end with the stub UDFManager
# ---------------------------------------------------------------------------

def test_write_skeleton_udf_copies_template_then_writes(
    patched_udfmanager, tmp_path
):
    """Drives copy + put + write through the public entry point."""
    from abmptools.udfcreate_v2 import write_skeleton_udf
    out = str(tmp_path / "skeleton.udf")
    written = write_skeleton_udf(
        out,
        dt=0.001, totalstep=5000, outstep=100,
        cellsize_nm=[2.0, 2.5, 3.0],
    )
    assert os.path.isfile(written)
    # template was copied (not zero size)
    assert os.path.getsize(written) > 0


def test_write_skeleton_udf_missing_template_raises(
    patched_udfmanager, tmp_path
):
    from abmptools.udfcreate_v2 import write_skeleton_udf
    with pytest.raises(FileNotFoundError, match="template"):
        write_skeleton_udf(
            str(tmp_path / "x.udf"),
            template_path=str(tmp_path / "nope.udf"),
            cellsize_nm=[1, 1, 1],
        )


def test_write_skeleton_udf_without_udfmanager_raises(tmp_path):
    """Lazy import of UDFManager — if missing, raises RuntimeError."""
    from abmptools.udfcreate_v2 import write_skeleton_udf
    with patch.dict(sys.modules, {"UDFManager": None}):
        sys.modules.pop("UDFManager", None)
        orig_import = __builtins__["__import__"] if isinstance(
            __builtins__, dict) else __builtins__.__import__

        def _forbid(name, *a, **kw):
            if name == "UDFManager":
                raise ImportError("simulated missing UDFManager")
            return orig_import(name, *a, **kw)

        with patch("builtins.__import__", side_effect=_forbid):
            with pytest.raises(RuntimeError, match="UDFManager is required"):
                write_skeleton_udf(
                    str(tmp_path / "x.udf"), cellsize_nm=[1, 1, 1],
                )
