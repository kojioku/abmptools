# -*- coding: utf-8 -*-
"""
udf_writer.py
-------------
Writes one GRO frame's worth of structural data into the current UDF record.

Contract (gro2udf context):
  AtomPosition.mol_id   = 0-based UDF mol index  (index into mol[])
  AtomPosition.atom_id  = 0-based atom index within mol (index into atom[])
  AtomPosition.x/y/z    = position  [nm]
  AtomPosition.vx/vy/vz = velocity  [nm/ps] -> converted to m/s here

UDF paths written per frame:
  Steps
  Time                                          (value in [ps])
  Structure.Position.mol[m].atom[a].{x,y,z}    [nm]
  Structure.Velocity.mol[m].atom[a].{x,y,z}    [m/s]
  Structure.Unit_Cell.Cell_Size.{a,b,c}         [nm]
  Structure.Unit_Cell.Cell_Size.{alpha,beta,gamma}
  Structure.Unit_Cell.Shear_Strain

Each UDF put first tries with an explicit unit string; on failure it falls
back to the no-unit form.  This replicates the try/except pattern from the
original gro2udf.py.
"""
from __future__ import annotations

from typing import List

from ..core.system_model import AtomPosition, CellGeometry

_VELOCITY_UNIT = 1000   # nm/ps → m/s  (1 nm/ps = 1000 m/s)
_SHEAR_STRAIN  = 0.0    # constant, same as original


class UDFWriter:
    """Writes positional / cell data for one GRO frame into the UDF record.

    Caller is responsible for calling ``udf.newRecord()`` before invoking
    ``write_frame()``.
    """

    def write_frame(
        self,
        udf,
        positions: List[AtomPosition],
        cell: CellGeometry,
        step: int,
        time: float,    # [ps]
    ) -> None:
        """Write Steps, Time, Position, Velocity and Cell to the current record."""
        # --- Steps and Time ---
        udf.put(step, "Steps")
        try:
            udf.put(time, "Time", "[ps]")
        except Exception:
            udf.put(time, "Time")

        # --- Atom positions and velocities ---
        for ap in positions:
            m = ap.mol_id       # 0-based mol index
            a = ap.atom_id      # 0-based atom index within mol
            vx = ap.vx * _VELOCITY_UNIT
            vy = ap.vy * _VELOCITY_UNIT
            vz = ap.vz * _VELOCITY_UNIT
            try:
                udf.put(ap.x, "Structure.Position.mol[].atom[].x", [m, a], "[nm]")
                udf.put(ap.y, "Structure.Position.mol[].atom[].y", [m, a], "[nm]")
                udf.put(ap.z, "Structure.Position.mol[].atom[].z", [m, a], "[nm]")
                udf.put(vx,   "Structure.Velocity.mol[].atom[].x", [m, a], "[m/s]")
                udf.put(vy,   "Structure.Velocity.mol[].atom[].y", [m, a], "[m/s]")
                udf.put(vz,   "Structure.Velocity.mol[].atom[].z", [m, a], "[m/s]")
            except Exception:
                udf.put(ap.x, "Structure.Position.mol[].atom[].x", [m, a])
                udf.put(ap.y, "Structure.Position.mol[].atom[].y", [m, a])
                udf.put(ap.z, "Structure.Position.mol[].atom[].z", [m, a])
                udf.put(vx,   "Structure.Velocity.mol[].atom[].x", [m, a])
                udf.put(vy,   "Structure.Velocity.mol[].atom[].y", [m, a])
                udf.put(vz,   "Structure.Velocity.mol[].atom[].z", [m, a])

        # --- Unit cell ---
        try:
            udf.put(cell.a, "Structure.Unit_Cell.Cell_Size.a", "[nm]")
            udf.put(cell.b, "Structure.Unit_Cell.Cell_Size.b", "[nm]")
            udf.put(cell.c, "Structure.Unit_Cell.Cell_Size.c", "[nm]")
        except Exception:
            udf.put(cell.a, "Structure.Unit_Cell.Cell_Size.a")
            udf.put(cell.b, "Structure.Unit_Cell.Cell_Size.b")
            udf.put(cell.c, "Structure.Unit_Cell.Cell_Size.c")

        udf.put(cell.alpha, "Structure.Unit_Cell.Cell_Size.alpha")
        udf.put(cell.beta,  "Structure.Unit_Cell.Cell_Size.beta")
        udf.put(cell.gamma, "Structure.Unit_Cell.Cell_Size.gamma")
        udf.put(_SHEAR_STRAIN, "Structure.Unit_Cell.Shear_Strain")
