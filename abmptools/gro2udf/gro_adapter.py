# -*- coding: utf-8 -*-
"""
gro_adapter.py
--------------
Converts a GROFrame into the SystemModel intermediate data types:
  List[AtomPosition]  and  CellGeometry

Reuses the existing abmptools.udf2gro.system_model dataclasses so that
the gro2udf pipeline shares the same intermediate representation.

Interpretation of AtomPosition fields in this gro2udf context:
  mol_id          -> 0-based UDF mol index       (used as index in mol[])
  atom_id         -> 0-based atom index within mol (used as index in atom[])
  x, y, z         -> position [nm]
  vx, vy, vz      -> velocity [nm/ps]  (UDFWriter converts to m/s when writing)
  mol_name_short  -> empty string (not used by UDFWriter)
  atom_gro_name   -> empty string (not used by UDFWriter)

Molecule boundary detection replicates the original gro2udf.py logic:
the first 5 characters of each atom line encode the residue number; a
change in those 5 characters marks the start of a new UDF molecule.
"""
from __future__ import annotations

import math
from typing import List, Tuple

from ..udf2gro.system_model import AtomPosition, CellGeometry
from .gro_parser import GROFrame


class GROAdapter:
    """Converts one GROFrame to (positions, cell) for UDFWriter."""

    def to_positions_and_cell(
        self,
        frame: GROFrame,
    ) -> Tuple[List[AtomPosition], CellGeometry]:
        positions = self._build_positions(frame)
        cell = self._build_cell(frame)
        return positions, cell

    # ------------------------------------------------------------------
    # Atom positions
    # ------------------------------------------------------------------

    def _build_positions(self, frame: GROFrame) -> List[AtomPosition]:
        """Convert atom lines to AtomPosition list.

        mol_id and atom_id store 0-based UDF indices (not GRO 1-based IDs).
        """
        mol_index = -1      # will be 0 after first increment
        atom_index = -1     # will be 0 after first increment in a new mol
        mol_idx_s = ""      # tracks the residue-col string

        positions: List[AtomPosition] = []

        for atom in frame.atoms:
            if atom.residue_col != mol_idx_s:
                # New molecule detected
                mol_index += 1
                atom_index = 0
                mol_idx_s = atom.residue_col
            else:
                atom_index += 1

            positions.append(AtomPosition(
                mol_id=mol_index,       # 0-based UDF mol index
                mol_name_short="",      # unused in gro2udf direction
                atom_gro_name="",       # unused in gro2udf direction
                atom_id=atom_index,     # 0-based atom index within mol
                x=atom.x,
                y=atom.y,
                z=atom.z,
                vx=atom.vx,             # [nm/ps]; UDFWriter multiplies × 1000
                vy=atom.vy,
                vz=atom.vz,
            ))

        return positions

    # ------------------------------------------------------------------
    # Cell geometry
    # ------------------------------------------------------------------

    @staticmethod
    def _build_cell(frame: GROFrame) -> CellGeometry:
        """Convert GRO box vectors to CellGeometry.

        GROMACS stores triclinic boxes as 9 values (v1x, v2y, v3z, v1y, v1z,
        v2x, v2z, v3x, v3y), i.e. temp[0..8] mapped to:
          a = (temp[0], temp[3], temp[4])
          b = (temp[5], temp[1], temp[6])
          c = (temp[7], temp[8], temp[2])
        This matches the original gro2udf.py normalization.
        """
        box = frame.box_vals

        if len(box) == 9:
            a_vec = [box[0], box[3], box[4]]
            b_vec = [box[5], box[1], box[6]]
            c_vec = [box[7], box[8], box[2]]
            la = _norm(a_vec)
            lb = _norm(b_vec)
            lc = _norm(c_vec)
            alpha = math.degrees(math.acos(_dot(b_vec, c_vec) / (lb * lc)))
            beta  = math.degrees(math.acos(_dot(a_vec, c_vec) / (la * lc)))
            gamma = math.degrees(math.acos(_dot(a_vec, b_vec) / (la * lb)))
            return CellGeometry(a=la, b=lb, c=lc,
                                alpha=alpha, beta=beta, gamma=gamma)

        elif len(box) == 3:
            return CellGeometry(a=box[0], b=box[1], c=box[2])

        else:
            print("warning : cell size ", box)
            return CellGeometry(a=0.0, b=0.0, c=0.0)


# ---------------------------------------------------------------------------
# Module-level helpers (pure functions, no side effects)
# ---------------------------------------------------------------------------

def _norm(v) -> float:
    return math.sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2)


def _dot(v1, v2) -> float:
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]
