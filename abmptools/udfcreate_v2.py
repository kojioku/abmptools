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

Status (2026-04-28, D-2: Molecular_Attributes)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Migrated so far:
  D-1:
    - Simulation_Conditions.Dynamics_Conditions.Time.{delta_T,
      Total_Steps, Output_Interval_Steps}
    - Initial_Structure.Initial_Unit_Cell.Cell_Size.{a,b,c,alpha,beta,gamma}
  D-2:
    - Molecular_Attributes.Atom_Type[].{Name, Mass}
    - Molecular_Attributes.Bond_Potential[].{Name, Potential_Type, R0,
      Harmonic.K}
    - Molecular_Attributes.Angle_Potential[].{Name, Potential_Type,
      theta0, Theta.K}
    - Molecular_Attributes.Torsion_Potential[].{Name, Potential_Type,
      Amber.PK/PN/PHASE/IDIVF or Cosine_Polynomial.K/N/p[]}
    - Plus legacy-list adapters (paramfile / atom_list shapes used by
      gen_udf input)

Remaining (planned, in order of size):
  D-3: Set_of_Molecules + Structure.Position
  D-4: switch fcewsmb call sites + deprecate ``udfcreate.gen_udf``

Until D-4 lands, this is an additive module — ``gen_udf`` stays
authoritative and v2 is a per-section migration target.

Public API
~~~~~~~~~~
``set_simulation_time(uobj, dt, totalstep, outstep)``
``set_initial_cell(uobj, cellsize_nm)``
``set_atom_types(uobj, atom_types)``
``set_bond_potentials(uobj, bond_types)``
``set_angle_potentials(uobj, angle_types)``
``set_torsion_potentials(uobj, torsion_types)``
``set_molecular_attributes(uobj, *, atom_types, bond_types, angle_types,
                           torsion_types)``
``atom_types_from_legacy(paramfile)``
``bond_types_from_legacy(paramfile, atom_list)``
``angle_types_from_legacy(paramfile, atom_list)``
``torsion_types_from_legacy(paramfile, atom_list)``
``write_skeleton_udf(out_path, *, ...)``
"""
from __future__ import annotations

import logging
import math
import os
import shutil
from typing import List, Optional, Sequence

from .gro2udf.top_model import (
    AtomTypeSpec,
    BondTypeSpec,
    AngleTypeSpec,
    TorsionTypeSpec,
)

logger = logging.getLogger(__name__)


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


def set_atom_types(uobj, atom_types: Sequence[AtomTypeSpec]) -> None:
    """Populate ``Molecular_Attributes.Atom_Type[]`` (Name + Mass).

    Equivalent to the ``putatomparam`` text helper.

    Parameters
    ----------
    uobj : UDFManager instance, jumped to common record (-1).
    atom_types : sequence of :class:`AtomTypeSpec`
        Force-field atom types in the order they should appear in the
        UDF (typically order of first appearance across molecules).
    """
    for j, at in enumerate(atom_types):
        uobj.put(at.name,
                 "Molecular_Attributes.Atom_Type[].Name", [j])
        uobj.put(float(at.mass),
                 "Molecular_Attributes.Atom_Type[].Mass", [j])


def set_bond_potentials(uobj, bond_types: Sequence[BondTypeSpec]) -> None:
    """Populate ``Molecular_Attributes.Bond_Potential[]`` (Harmonic).

    Equivalent to the ``putbondparam`` text helper. Currently only
    ``funct=1`` (Harmonic) is supported — the same scope gen_udf has.

    Parameters
    ----------
    uobj : UDFManager instance, jumped to common record (-1).
    bond_types : sequence of :class:`BondTypeSpec`
    """
    for j, bt in enumerate(bond_types):
        uobj.put(bt.name,
                 "Molecular_Attributes.Bond_Potential[].Name", [j])
        uobj.put("Harmonic",
                 "Molecular_Attributes.Bond_Potential[].Potential_Type", [j])
        uobj.put(float(bt.r0),
                 "Molecular_Attributes.Bond_Potential[].R0", [j], "[nm]")
        uobj.put(float(bt.kb),
                 "Molecular_Attributes.Bond_Potential[].Harmonic.K",
                 [j], "[kJ/mol/nm^2]")


def set_angle_potentials(uobj, angle_types: Sequence[AngleTypeSpec]) -> None:
    """Populate ``Molecular_Attributes.Angle_Potential[]`` (Theta).

    Equivalent to ``putangleparam``. ``theta0`` is converted from the
    GROMACS-style absolute angle to the COGNAC supplement convention
    (``180 - theta0_gro``) just like :class:`TopExporter` does.

    Parameters
    ----------
    uobj : UDFManager instance, jumped to common record (-1).
    angle_types : sequence of :class:`AngleTypeSpec`
    """
    for j, at in enumerate(angle_types):
        q0 = 180.0 - float(at.theta0)
        uobj.put(at.name,
                 "Molecular_Attributes.Angle_Potential[].Name", [j])
        uobj.put("Theta",
                 "Molecular_Attributes.Angle_Potential[].Potential_Type", [j])
        uobj.put(q0,
                 "Molecular_Attributes.Angle_Potential[].theta0", [j])
        uobj.put(float(at.k),
                 "Molecular_Attributes.Angle_Potential[].Theta.K",
                 [j], "[kJ/mol/rad^2]")


def set_torsion_potentials(uobj, torsion_types: Sequence[TorsionTypeSpec]) -> None:
    """Populate ``Molecular_Attributes.Torsion_Potential[]``.

    Supports Amber (funct 1/4/9, params=[PHASE, PK, PN]) and
    Ryckaert-Bellemans / Cosine_Polynomial (funct 3, params=[C0..C5]).
    Same convention as :class:`TopExporter`.

    Parameters
    ----------
    uobj : UDFManager instance, jumped to common record (-1).
    torsion_types : sequence of :class:`TorsionTypeSpec`
    """
    for j, tt in enumerate(torsion_types):
        funct = tt.funct
        params = tt.params
        uobj.put(tt.name,
                 "Molecular_Attributes.Torsion_Potential[].Name", [j])
        if funct in (1, 9, 4):
            phase, pk, pn = params[0], params[1], params[2]
            uobj.put("Amber",
                     "Molecular_Attributes.Torsion_Potential[].Potential_Type",
                     [j])
            uobj.put(float(pk),
                     "Molecular_Attributes.Torsion_Potential[].Amber.PK",
                     [j], "[kJ/mol]")
            uobj.put(1,
                     "Molecular_Attributes.Torsion_Potential[].Amber.IDIVF",
                     [j])
            uobj.put(float(pn),
                     "Molecular_Attributes.Torsion_Potential[].Amber.PN",
                     [j])
            uobj.put(float(phase),
                     "Molecular_Attributes.Torsion_Potential[].Amber.PHASE",
                     [j])
            uobj.put(1,
                     "Molecular_Attributes.Torsion_Potential[].Amber.trans_is_0",
                     [j])
        elif funct == 3:
            # trim trailing zero coefficients
            nparams = len(params)
            while nparams > 0 and params[nparams - 1] == 0.0:
                nparams -= 1
            uobj.put("Cosine_Polynomial",
                     "Molecular_Attributes.Torsion_Potential[].Potential_Type",
                     [j])
            uobj.put(int(nparams),
                     "Molecular_Attributes.Torsion_Potential[].Cosine_Polynomial.N",
                     [j])
            uobj.put(1.0,
                     "Molecular_Attributes.Torsion_Potential[].Cosine_Polynomial.K",
                     [j], "[kJ/mol]")
            for k in range(nparams):
                uobj.put(float(params[k]),
                         "Molecular_Attributes.Torsion_Potential[].Cosine_Polynomial.p[]",
                         [j, k])
        else:
            logger.warning(
                "torsion funct=%s is not supported (%s)", funct, tt.name)


def set_molecular_attributes(uobj, *,
                              atom_types: Sequence[AtomTypeSpec],
                              bond_types: Sequence[BondTypeSpec] = (),
                              angle_types: Sequence[AngleTypeSpec] = (),
                              torsion_types: Sequence[TorsionTypeSpec] = ()
                              ) -> None:
    """Convenience wrapper that drives all four set_* writers.

    Replaces the legacy ``putmolecularattributes`` orchestrator.
    """
    set_atom_types(uobj, atom_types)
    set_bond_potentials(uobj, bond_types)
    set_angle_potentials(uobj, angle_types)
    set_torsion_potentials(uobj, torsion_types)


# ---------------------------------------------------------------------------
# Legacy-list adapters
# ---------------------------------------------------------------------------
# Convert from the raw Python lists that gen_udf consumes into the
# dataclass-based form used by udfcreate_v2. Each adapter is the
# "single-section" inverse of the corresponding put* helper in
# abmptools.udfcreate.
#
# Layouts (matching abmptools.udfcreate.udfcreate.put*):
#   atom_list[i]        = [name, ?, type_name, ...]      (FF-typed atom row)
#   ljparam[i]          = [type_name, ?, ?, epsilon, sigma_raw, ...]
#                         where sigma_raw is the GROMACS-style sigma in nm
#   bondparam[i]        = [name1, name2, k_half, R0]    (k_half = K/2 [kJ/mol/nm^2])
#   angleparam[i]       = [name1, name2, name3, k_angle, theta0_gro_deg]
#   torsionparam[i]     = [name1, name2, name3, name4, funct,
#                          improper, n_params, *params]
# ---------------------------------------------------------------------------


def atom_types_from_legacy(paramfile: Sequence) -> List[AtomTypeSpec]:
    """Convert ``putatomparam``/``putljparam`` input into AtomTypeSpec.

    The ``paramfile`` rows used by gen_udf are
    ``[name, ?, ?, epsilon, sigma_raw, ...]`` where ``sigma_raw`` is the
    GROMACS-style σ in nm; ``mass`` is at index 1.
    """
    out: List[AtomTypeSpec] = []
    for row in paramfile:
        name = str(row[0])
        mass = float(row[1])
        epsilon = float(row[3]) if len(row) > 3 else 0.0
        sigma = float(row[4]) if len(row) > 4 else 0.0
        out.append(AtomTypeSpec(name=name, mass=mass,
                                sigma=sigma, epsilon=epsilon))
    return out


def bond_types_from_legacy(paramfile: Sequence,
                           atom_list: Sequence) -> List[BondTypeSpec]:
    """Convert ``putbondparam`` input into BondTypeSpec.

    ``putbondparam`` builds ``"name1-name2"`` from ``atom_list[mol[i]][0]``
    (the atom global names rather than type names), so we mirror that
    here. ``k_half`` in the legacy is ``K/2``; we restore K = k_half × 2.
    """
    out: List[BondTypeSpec] = []
    for row in paramfile:
        type1, type2 = str(row[0]), str(row[1])
        k_bond = float(row[2]) * 2.0
        r0 = float(row[3])
        # mirror the gen_udf naming: lookup global atom name from atom_list
        n1, n2 = type1, type2
        for atom in atom_list:
            if atom[2] == type1:
                n1 = atom[0]
                break
        for atom in atom_list:
            if atom[2] == type2:
                n2 = atom[0]
                break
        out.append(BondTypeSpec(name=f"{n1}-{n2}",
                                name1=type1, name2=type2,
                                funct=1, r0=r0, kb=k_bond))
    return out


def angle_types_from_legacy(paramfile: Sequence,
                            atom_list: Sequence) -> List[AngleTypeSpec]:
    """Convert ``putangleparam`` input into AngleTypeSpec.

    ``theta0`` in the legacy paramfile is the GROMACS-style angle (deg)
    *before* the COGNAC 180-supplement conversion. We keep it raw here
    because :func:`set_angle_potentials` re-applies the supplement.
    """
    out: List[AngleTypeSpec] = []
    for row in paramfile:
        t1, t2, t3 = str(row[0]), str(row[1]), str(row[2])
        k = float(row[3])
        theta0 = float(row[4])
        n1, n2, n3 = t1, t2, t3
        for atom in atom_list:
            if atom[2] == t1:
                n1 = atom[0]; break
        for atom in atom_list:
            if atom[2] == t2:
                n2 = atom[0]; break
        for atom in atom_list:
            if atom[2] == t3:
                n3 = atom[0]; break
        out.append(AngleTypeSpec(name=f"{n1}-{n2}-{n3}",
                                 name1=t1, name2=t2, name3=t3,
                                 funct=1, theta0=theta0, k=k))
    return out


def torsion_types_from_legacy(paramfile: Sequence,
                              atom_list: Sequence) -> List[TorsionTypeSpec]:
    """Convert ``puttorsionparam`` input into TorsionTypeSpec.

    Layout assumed: ``[name1, name2, name3, name4, funct, improper,
    n_params, *params]``. The first 4 fields are FF type names; the
    parameter list shape depends on funct.
    """
    out: List[TorsionTypeSpec] = []
    for row in paramfile:
        t1, t2, t3, t4 = (str(row[0]), str(row[1]), str(row[2]), str(row[3]))
        funct = int(row[4])
        improper = bool(row[5])
        n_params = int(row[6])
        params = list(row[7:7 + n_params])
        n1, n2, n3, n4 = t1, t2, t3, t4
        for atom in atom_list:
            if atom[2] == t1:
                n1 = atom[0]; break
        for atom in atom_list:
            if atom[2] == t2:
                n2 = atom[0]; break
        for atom in atom_list:
            if atom[2] == t3:
                n3 = atom[0]; break
        for atom in atom_list:
            if atom[2] == t4:
                n4 = atom[0]; break
        out.append(TorsionTypeSpec(name=f"{n1}-{n2}-{n3}-{n4}",
                                   name1=t1, name2=t2, name3=t3, name4=t4,
                                   funct=funct, improper=improper,
                                   params=[float(p) for p in params]))
    return out


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
