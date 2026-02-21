# -*- coding: utf-8 -*-
"""
gro_writer.py
-------------
Writes the GROMACS .gro structure file from a SystemModel.
"""
from __future__ import annotations
import math
from ....core.system_model import SystemModel, CellGeometry


class GroWriter:
    """Generates the content of a .gro file."""

    def write(self, model: SystemModel, filepath: str) -> None:
        content = self._build(model)
        with open(filepath, "w") as f:
            f.write(content)

    def _build(self, model: SystemModel) -> str:
        total_atoms = len(model.atom_positions)
        lines = []
        lines.append(model.title)
        lines.append(str(total_atoms))

        for ap in model.atom_positions:
            line = (
                "%5d" % ap.mol_id
                + "{0:<5}".format(ap.mol_name_short)
                + "%5s" % ap.atom_gro_name
                + "%5i" % ap.atom_id
                + "%8.3f" % ap.x
                + "%8.3f" % ap.y
                + "%8.3f" % ap.z
                + "%8.4f" % ap.vx
                + "%8.4f" % ap.vy
                + "%8.4f" % ap.vz
            )
            lines.append(line)

        lines.append(self._cell_line(model.cell))
        return "\n".join(lines) + "\n"

    @staticmethod
    def _cell_line(cell: CellGeometry) -> str:
        if cell.is_rectangular():
            return "{} {} {}".format(cell.a, cell.b, cell.c)
        else:
            # Full triclinic box matrix
            a = cell.a
            b = cell.b
            c = cell.c
            cosa = math.cos(math.radians(cell.alpha))
            cosb = math.cos(math.radians(cell.beta))
            cosg = math.cos(math.radians(cell.gamma))
            sing = math.sin(math.radians(cell.gamma))

            ax = a
            ay = 0.0
            az = 0.0
            bx = b * cosg
            by = b * sing
            bz = 0.0
            lc = c
            cx = lc * cosb
            cy = lc * (cosa - cosb * cosg) / sing
            cz = (lc * lc - cx * cx - cy * cy) ** 0.5

            fmt = "  {:.5f}" * 9
            return fmt.format(ax, by, cz, ay, az, bx, bz, cx, cy)
