# -*- coding: utf-8 -*-
"""
abmptools.udfcreate_v2
----------------------
Phase 2c-E successor of :mod:`abmptools.udfcreate` (gen_udf path) — emits
COGNAC UDF using a static template + structured ``UDFManager.put`` calls
instead of Python-string concatenation.

Why
~~~
The legacy ``gen_udf`` builds a UDF by concatenating ~700 lines of
hand-written ``put*`` text helpers. This is fragile (recall the
self-angle / duplicate-angle bug fixed via post-process clean_top.py)
and difficult to test piecewise. ``TopExporter`` (gro2udf direction) has
shown that a template-UDF + UDFManager.put approach is much cleaner —
this module brings the same approach to the udf_create direction.

Status (2026-04-28, D-1 scaffolding)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Migrated so far:
  - Simulation_Conditions.Dynamics_Conditions.Time.{delta_T, Total_Steps,
    Output_Interval_Steps}
  - Initial_Structure.Initial_Unit_Cell.Cell_Size.{a, b, c, alpha, beta, gamma}

Remaining (planned, in order of size):
  D-2: Molecular_Attributes (atom_param + lj/bond/angle/torsion params)
  D-3: Set_of_Molecules + Structure.Position
  D-4: switch fcewsmb call sites + deprecate ``udfcreate.gen_udf``

Until D-4 lands, this is an additive module — ``gen_udf`` stays
authoritative and v2 is a per-section migration target.

Public API
~~~~~~~~~~
``set_simulation_time(uobj, dt, totalstep, outstep)``
``set_initial_cell(uobj, cellsize_nm)``
``write_skeleton_udf(out_path, *, template_path=None, dt, totalstep,
                     outstep, cellsize_nm)``
"""
from __future__ import annotations

import os
import shutil
from typing import List, Optional, Sequence


def _default_template_path() -> str:
    """Bundled UDF template (re-uses the gro2udf default_template.udf)."""
    here = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(here, "gro2udf", "default_template.udf")


def set_simulation_time(uobj, dt: float, totalstep: int,
                        outstep: int) -> None:
    """Populate ``Simulation_Conditions.Dynamics_Conditions.Time``.

    Equivalent to the ``{dt, totalstep, outstep}`` triple in the legacy
    ``putsimulationcondition`` text template.

    Parameters
    ----------
    uobj : UDFManager instance, jumped to the common record (``-1``).
    dt : float
        Time step in **picoseconds**.
    totalstep : int
        Total number of MD steps.
    outstep : int
        Output interval in steps.
    """
    uobj.put(dt, "Simulation_Conditions.Dynamics_Conditions.Time.delta_T",
             "[ps]")
    uobj.put(int(totalstep),
             "Simulation_Conditions.Dynamics_Conditions.Time.Total_Steps")
    uobj.put(int(outstep),
             "Simulation_Conditions.Dynamics_Conditions.Time.Output_Interval_Steps")


def set_initial_cell(uobj, cellsize_nm: Sequence[float],
                     angles_deg: Sequence[float] = (90.0, 90.0, 90.0)) -> None:
    """Populate ``Initial_Structure.Initial_Unit_Cell.Cell_Size``.

    Equivalent to the ``{0.0,{a,b,c,90,90,90}0.0}`` tuple in the legacy
    ``putinitialstructure`` text template.

    Parameters
    ----------
    uobj : UDFManager instance, jumped to the common record (``-1``).
    cellsize_nm : sequence of 3 floats
        Box edge lengths ``[a, b, c]`` in **nm**. Note: COGNAC's UDF
        stores cell sizes as plain numbers without an explicit unit, so
        they end up in whatever unit the schema implies (Å in the
        existing flow). We follow the legacy behaviour of writing the
        raw numbers; conversion to Å (if needed) is the caller's
        responsibility.
    angles_deg : sequence of 3 floats, optional
        ``[alpha, beta, gamma]`` in degrees. Defaults to orthorhombic
        ``(90, 90, 90)``.
    """
    if len(cellsize_nm) != 3:
        raise ValueError(
            f"cellsize_nm must have 3 elements, got {len(cellsize_nm)}"
        )
    if len(angles_deg) != 3:
        raise ValueError(
            f"angles_deg must have 3 elements, got {len(angles_deg)}"
        )
    base = "Initial_Structure.Initial_Unit_Cell.Cell_Size"
    uobj.put(float(cellsize_nm[0]), f"{base}.a")
    uobj.put(float(cellsize_nm[1]), f"{base}.b")
    uobj.put(float(cellsize_nm[2]), f"{base}.c")
    uobj.put(float(angles_deg[0]), f"{base}.alpha")
    uobj.put(float(angles_deg[1]), f"{base}.beta")
    uobj.put(float(angles_deg[2]), f"{base}.gamma")


def write_skeleton_udf(out_path: str, *,
                       dt: float = 2.04584957182253e-02,
                       totalstep: int = 10000,
                       outstep: int = 100,
                       cellsize_nm: Sequence[float] = (10.0, 10.0, 10.0),
                       template_path: Optional[str] = None) -> str:
    """Copy *template_path* to *out_path* and populate the v2 sections.

    This is a convenience wrapper that drives ``set_simulation_time`` +
    ``set_initial_cell`` end-to-end. It produces a UDF that is *valid*
    but missing all force-field / molecule data — the remaining
    sections (Molecular_Attributes, Set_of_Molecules, Structure) are
    still emitted by the legacy gen_udf path until D-2/D-3 land.

    Parameters
    ----------
    out_path : str
        Destination UDF (overwritten if it exists).
    dt, totalstep, outstep, cellsize_nm
        See :func:`set_simulation_time` and :func:`set_initial_cell`.
    template_path : str, optional
        Override the bundled template; defaults to the gro2udf
        ``default_template.udf``.

    Returns
    -------
    str
        Absolute path of the written UDF.

    Raises
    ------
    RuntimeError
        When UDFManager is not importable.
    """
    try:
        from UDFManager import UDFManager  # type: ignore
    except ImportError as e:
        raise RuntimeError(
            "UDFManager is required by udfcreate_v2.write_skeleton_udf "
            "(install OCTA / UDFManager Python bindings)."
        ) from e

    template = template_path or _default_template_path()
    if not os.path.isfile(template):
        raise FileNotFoundError(
            f"UDF template not found: {template}"
        )

    out_path = os.path.abspath(out_path)
    shutil.copy(template, out_path)

    uobj = UDFManager(out_path)
    uobj.jump(-1)
    set_simulation_time(uobj, dt, totalstep, outstep)
    set_initial_cell(uobj, cellsize_nm)
    uobj.write()

    return out_path
