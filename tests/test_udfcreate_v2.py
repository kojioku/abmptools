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


# ---------------------------------------------------------------------------
# D-2: Molecular_Attributes
# ---------------------------------------------------------------------------

class TestSetAtomTypes:
    def test_writes_name_and_mass_per_index(self):
        from abmptools.gro2udf.top_model import AtomTypeSpec
        from abmptools.udfcreate_v2 import set_atom_types
        u = _StubUDFManager("/tmp/x.udf")
        atom_types = [
            AtomTypeSpec(name="c3", mass=12.011, sigma=0.34, epsilon=0.36),
            AtomTypeSpec(name="hc", mass=1.008,  sigma=0.26, epsilon=0.07),
        ]
        set_atom_types(u, atom_types)
        # 2 names + 2 masses
        names = [
            c for c in u.calls
            if c[1] == "Molecular_Attributes.Atom_Type[].Name"
        ]
        masses = [
            c for c in u.calls
            if c[1] == "Molecular_Attributes.Atom_Type[].Mass"
        ]
        assert len(names) == 2 and len(masses) == 2
        # index ordering preserved
        assert names[0][0] == "c3" and names[1][0] == "hc"
        assert masses[0][0] == 12.011

    def test_empty_list_writes_nothing(self):
        from abmptools.udfcreate_v2 import set_atom_types
        u = _StubUDFManager("/tmp/x.udf")
        set_atom_types(u, [])
        assert u.calls == []


class TestSetBondPotentials:
    def test_writes_harmonic_with_units(self):
        from abmptools.gro2udf.top_model import BondTypeSpec
        from abmptools.udfcreate_v2 import set_bond_potentials
        u = _StubUDFManager("/tmp/x.udf")
        bonds = [
            BondTypeSpec(name="c3-c3-0", name1="c3", name2="c3",
                         funct=1, r0=0.151, kb=2.5e5),
        ]
        set_bond_potentials(u, bonds)
        # 4 puts per bond: Name, Potential_Type, R0, Harmonic.K
        assert len(u.calls) == 4
        # R0 has [nm] unit (last positional arg)
        r0_call = next(c for c in u.calls
                       if c[1] == "Molecular_Attributes.Bond_Potential[].R0")
        assert "[nm]" in r0_call[2]
        # Harmonic.K has [kJ/mol/nm^2]
        k_call = next(c for c in u.calls
                      if c[1] == "Molecular_Attributes.Bond_Potential[].Harmonic.K")
        assert "[kJ/mol/nm^2]" in k_call[2]


class TestSetAnglePotentials:
    def test_supplements_theta0_to_cognac_convention(self):
        """COGNAC stores 180 - theta0_gromacs."""
        from abmptools.gro2udf.top_model import AngleTypeSpec
        from abmptools.udfcreate_v2 import set_angle_potentials
        u = _StubUDFManager("/tmp/x.udf")
        angles = [
            AngleTypeSpec(name="c3-c3-c3-0",
                          name1="c3", name2="c3", name3="c3",
                          funct=1, theta0=109.5, k=300.0),
        ]
        set_angle_potentials(u, angles)
        theta_call = next(
            c for c in u.calls
            if c[1] == "Molecular_Attributes.Angle_Potential[].theta0"
        )
        assert theta_call[0] == pytest.approx(70.5)  # 180 - 109.5


class TestSetTorsionPotentials:
    def test_amber_funct1(self):
        from abmptools.gro2udf.top_model import TorsionTypeSpec
        from abmptools.udfcreate_v2 import set_torsion_potentials
        u = _StubUDFManager("/tmp/x.udf")
        ts = [TorsionTypeSpec(
            name="c3-c3-c3-c3-0",
            name1="c3", name2="c3", name3="c3", name4="c3",
            funct=1, improper=False,
            params=[180.0, 0.5, 3],  # PHASE, PK, PN
        )]
        set_torsion_potentials(u, ts)
        # locate Potential_Type
        pt = next(c for c in u.calls
                  if c[1] == "Molecular_Attributes.Torsion_Potential[].Potential_Type")
        assert pt[0] == "Amber"
        pk = next(c for c in u.calls
                  if c[1] == "Molecular_Attributes.Torsion_Potential[].Amber.PK")
        assert pk[0] == 0.5
        assert "[kJ/mol]" in pk[2]

    def test_cosine_polynomial_funct3_trims_zero_tail(self):
        from abmptools.gro2udf.top_model import TorsionTypeSpec
        from abmptools.udfcreate_v2 import set_torsion_potentials
        u = _StubUDFManager("/tmp/x.udf")
        ts = [TorsionTypeSpec(
            name="rb",
            name1="x", name2="x", name3="x", name4="x",
            funct=3, improper=False,
            params=[1.0, 2.0, 3.0, 0.0, 0.0, 0.0],
        )]
        set_torsion_potentials(u, ts)
        n_call = next(
            c for c in u.calls
            if c[1] == "Molecular_Attributes.Torsion_Potential[].Cosine_Polynomial.N"
        )
        # trailing zeros trimmed → N=3 not 6
        assert n_call[0] == 3
        p_calls = [
            c for c in u.calls
            if c[1] == "Molecular_Attributes.Torsion_Potential[].Cosine_Polynomial.p[]"
        ]
        assert len(p_calls) == 3


class TestSetMolecularAttributes:
    def test_drives_all_four_writers(self):
        from abmptools.gro2udf.top_model import (
            AtomTypeSpec, BondTypeSpec, AngleTypeSpec, TorsionTypeSpec,
        )
        from abmptools.udfcreate_v2 import set_molecular_attributes
        u = _StubUDFManager("/tmp/x.udf")
        set_molecular_attributes(
            u,
            atom_types=[AtomTypeSpec(name="c3", mass=12.011,
                                     sigma=0.34, epsilon=0.36)],
            bond_types=[BondTypeSpec(name="c3-c3-0",
                                     name1="c3", name2="c3", funct=1,
                                     r0=0.15, kb=2e5)],
            angle_types=[AngleTypeSpec(name="ang-0", name1="c3",
                                        name2="c3", name3="c3",
                                        funct=1, theta0=109.5, k=300.0)],
            torsion_types=[TorsionTypeSpec(name="t-0",
                                           name1="c3", name2="c3",
                                           name3="c3", name4="c3",
                                           funct=1, improper=False,
                                           params=[180.0, 0.5, 3])],
        )
        # at minimum: atom-type Name + bond Name + angle Name + torsion Name
        loc_set = {c[1] for c in u.calls}
        assert "Molecular_Attributes.Atom_Type[].Name" in loc_set
        assert "Molecular_Attributes.Bond_Potential[].Name" in loc_set
        assert "Molecular_Attributes.Angle_Potential[].Name" in loc_set
        assert "Molecular_Attributes.Torsion_Potential[].Name" in loc_set


# ---------------------------------------------------------------------------
# Legacy-list adapters
# ---------------------------------------------------------------------------

class TestAtomTypesFromLegacy:
    def test_extracts_name_mass_sigma_epsilon(self):
        from abmptools.udfcreate_v2 import atom_types_from_legacy
        # gen_udf-style: [name, mass, ?, epsilon, sigma_raw, ...]
        paramfile = [
            ["c3", 12.011, 0.0, 0.36, 0.34],
            ["hc", 1.008,  0.0, 0.07, 0.26],
        ]
        out = atom_types_from_legacy(paramfile)
        assert len(out) == 2
        assert out[0].name == "c3" and out[0].mass == 12.011
        assert out[0].sigma == 0.34 and out[0].epsilon == 0.36


class TestBondTypesFromLegacy:
    def test_doubles_k_half(self):
        """gen_udf stores k_bond/2 in column 2; v2 should restore K."""
        from abmptools.udfcreate_v2 import bond_types_from_legacy
        atom_list = [["GlobC1", 0, "c3"], ["GlobC2", 0, "c3"]]
        paramfile = [["c3", "c3", 1.25e5, 0.151]]  # k_half=125000
        out = bond_types_from_legacy(paramfile, atom_list)
        assert len(out) == 1
        assert out[0].kb == pytest.approx(2.5e5)
        assert out[0].r0 == 0.151
        assert out[0].name1 == "c3"


class TestAngleTypesFromLegacy:
    def test_keeps_theta0_raw(self):
        """legacy angleparam has theta0 in deg; we keep raw and let the
        writer apply 180-supplement convention."""
        from abmptools.udfcreate_v2 import angle_types_from_legacy
        atom_list = [["G1", 0, "c3"]]
        paramfile = [["c3", "c3", "c3", 300.0, 109.5]]
        out = angle_types_from_legacy(paramfile, atom_list)
        assert out[0].theta0 == 109.5  # not yet supplemented
        assert out[0].k == 300.0


class TestTorsionTypesFromLegacy:
    def test_amber_layout(self):
        from abmptools.udfcreate_v2 import torsion_types_from_legacy
        atom_list = [["G1", 0, "c3"]]
        # [name1, name2, name3, name4, funct, improper, n_params, *params]
        paramfile = [
            ["c3", "c3", "c3", "c3", 1, False, 3, 180.0, 0.5, 3],
        ]
        out = torsion_types_from_legacy(paramfile, atom_list)
        assert len(out) == 1
        assert out[0].funct == 1
        assert out[0].params == [180.0, 0.5, 3.0]


# ---------------------------------------------------------------------------
# write_skeleton_udf still works after the additions (no regression)
# ---------------------------------------------------------------------------

def test_write_skeleton_udf_still_works(patched_udfmanager, tmp_path):
    """D-1 entry point continues to function with D-2 imports added."""
    from abmptools.udfcreate_v2 import write_skeleton_udf
    out = write_skeleton_udf(
        str(tmp_path / "skeleton.udf"),
        dt=0.001, totalstep=1000, outstep=10,
        cellsize_nm=[2.0, 2.0, 2.0],
    )
    assert os.path.isfile(out)
