# -*- coding: utf-8 -*-
"""Write an angle potential that has no closed form as a numerical table.

Restricted bending (GROMACS angle funct 10) is

    V(t) = k/2 (cos t - cos t0)^2 / sin^2 t

and the engine has no such function. It does take a **table**, read as a
function of one variable, and that reproduces the potential to whatever
accuracy the mesh is given -- this is not an approximation of the form, only
of the sampling, and the sampling error is measured in the tests.

Two conventions matter and both are read off the engine's own source
(``src/tableangle.cpp``):

* the table variable is ``cos`` of the **exterior** angle,
  ``var = dr0 . dr1 / |dr0||dr1|``, which is ``-cos t`` of the interior angle
  GROMACS uses -- the same supplement that bonded angles already follow;
* the values are energies in the file's own energy unit.

The table file is a small UDF of its own (``Mesh`` + ``FieldValue`` +
``Unit_Parameter``), the layout the engine ships as
``sample/tablepotential/angle_table.udf``.
"""
from __future__ import annotations

import logging
import math
from pathlib import Path
from typing import List, Optional, Sequence

logger = logging.getLogger(__name__)

__all__ = ["RESTRICTED_BENDING", "reb_energy", "table_values", "write_angle_table"]

#: GROMACS angle funct for restricted bending.
RESTRICTED_BENDING = 10

#: Divisions of the ``cos`` axis. 2000 puts the points 0.001 apart, which
#: keeps the interpolation error below a thousandth of a kJ/mol away from
#: the ends (``tests/test_gro2udf_angle_table.py`` measures it).
DEFAULT_DIVISIONS = 2000

#: The potential diverges where ``sin t`` reaches zero, so the two end points
#: are written at this value instead of infinity. It is far above anything a
#: simulation reaches, and the engine never needs a force there because the
#: geometry cannot get past it.
END_VALUE_KJ = 1.0e6


def reb_energy(var: float, k: float, theta0_deg: float) -> float:
    """Restricted bending at ``var = cos(exterior angle)`` [kJ/mol].

    *k* is the force constant from the topology [kJ/mol] and *theta0_deg* the
    equilibrium angle in GROMACS' (interior) convention.
    """
    sin2 = 1.0 - var * var
    if sin2 <= 1.0e-12:
        return END_VALUE_KJ
    cos_t = -var                      # interior angle
    cos_t0 = math.cos(math.radians(theta0_deg))
    value = 0.5 * k * (cos_t - cos_t0) ** 2 / sin2
    return min(value, END_VALUE_KJ)


def table_values(k: float, theta0_deg: float,
                 divisions: int = DEFAULT_DIVISIONS) -> List[float]:
    """The ``divisions + 1`` energies [kJ/mol] over ``var`` from -1 to 1."""
    step = 2.0 / divisions
    return [reb_energy(-1.0 + i * step, k, theta0_deg)
            for i in range(divisions + 1)]


def write_angle_table(path: Path, name: str, k: float, theta0_deg: float,
                      energy_unit_kj: float = 4.184,
                      mass_unit_amu: float = 1.0,
                      length_unit_nm: float = 0.1,
                      divisions: int = DEFAULT_DIVISIONS) -> Path:
    """Write one table UDF and return its path.

    *energy_unit_kj* must be the same unit the main UDF uses, because the
    values are written in it.
    """
    path = Path(path)
    values = table_values(k, theta0_deg, divisions)
    rows = "\n".join("          {:.8g}".format(v / energy_unit_kj) for v in values)
    head = [
        "Table Potential FOR COGNAC",
        '\\include{"cognac_unit.udf"}',
        "\\begin{def}",
        "class MeshAxis:{",
        '\tvalues[]:double\t"Parameter of axis. [0]:min, [1]:max, [2]:divisions"',
        "}",
        "Mesh:{",
        '\tname:KEY\t"Name of mesh"',
        '\ttype:string\t"Type of mesh. Always \'Regular\' in COGNAC"',
        '\taxes[]:MeshAxis\t"Array of the parameter for each axis"',
        "}",
        "FieldValue: {",
        '    values[]:double\t[epsilon]\t"Field data"',
        "}",
        "Unit_Parameter:{",
        '\tName: string\t"Name"',
        '\tComment:string\t"Comment"',
        '\tMass:double\t[amu]\t"Unit mass (atomic mass unit)"',
        '\tEnergy:double\t[kJ/mol]\t"Unit energy(kJ/mol)"',
        '\tLength:double\t[nm]\t"Unit length(nm)"',
        '} "Parameters for unit conversion"',
        "\\end{def}",
        "",
        "\\begin{data}",
        "",
        "Mesh:{",
        '\tname "{}"'.format(name),
        '\ttype "REGULAR"',
        "\t[ devision",
        "\t{{ [ -1.0 1.0 {} ] }}".format(divisions),
        "\t]",
        "}",
        "",
        "FieldValue:{",
        "\t[ value",
        rows,
        "\t]",
        "}",
        "",
        'Unit_Parameter:{{"{}"  "restricted bending, written by gro2udf"  '
        "{:g} {:g} {:g} }}".format(name, mass_unit_amu, energy_unit_kj,
                                   length_unit_nm),
        "\\end{data}",
        "",
    ]
    path.write_text("\n".join(head), encoding="utf-8", newline="\n")
    logger.info("angle table written: %s (%d points, k = %g, theta0 = %g)",
                path, divisions + 1, k, theta0_deg)
    return path
