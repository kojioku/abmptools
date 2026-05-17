"""
bdf_reader.py
-------------
UDFManager wrapper for COGNAC BDF/UDF trajectory I/O.

Loads:
- topology (atoms, bonds per molecule, atom types)
- cell geometry per frame (orthogonal box assumed)
- atomic positions per frame
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Tuple

import numpy as np


@dataclass
class AtomInfo:
    """Per-atom topology info (constant across frames)."""
    atom_id: int
    atom_name: str
    atom_type: str


@dataclass
class BondInfo:
    """Per-bond topology info."""
    type_name: str
    atom1: int
    atom2: int
    order: float


@dataclass
class MoleculeTopology:
    """One molecule's static topology."""
    mol_name: str
    atoms: List[AtomInfo]
    bonds: List[BondInfo]


@dataclass
class CellBox:
    """Orthogonal cell box (Å)."""
    a: float
    b: float
    c: float

    def as_array(self) -> np.ndarray:
        return np.array([self.a, self.b, self.c], dtype=float)


@dataclass
class TrajectoryFrame:
    """One snapshot: positions for all molecules + cell."""
    record: int
    cell: CellBox
    positions: List[np.ndarray] = field(default_factory=list)


class BDFTrajectory:
    """
    UDFManager-backed trajectory reader.

    Use:
        traj = BDFTrajectory("foo.bdf")
        traj.load_topology()
        for frame in traj.iter_frames():
            ...
    """

    def __init__(self, path: str):
        from UDFManager import UDFManager
        self.path = path
        self._udf = UDFManager(path)
        self.n_records = self._udf.totalRecord()
        self.molecules: List[MoleculeTopology] = []

    def load_topology(self) -> None:
        """Read Set_of_Molecules.molecule[] once."""
        u = self._udf
        n_mol = u.size("Set_of_Molecules.molecule[]")
        if n_mol is None or n_mol <= 0:
            raise RuntimeError(
                f"No molecules found in Set_of_Molecules.molecule[] in {self.path}"
            )

        self.molecules = []
        for i in range(n_mol):
            mol_name = u.get(f"Set_of_Molecules.molecule[{i}].Mol_Name") or "UNK"
            n_atoms = u.size(f"Set_of_Molecules.molecule[{i}].atom[]") or 0
            atoms = []
            for j in range(n_atoms):
                aid = u.get(f"Set_of_Molecules.molecule[{i}].atom[{j}].Atom_ID")
                name = u.get(f"Set_of_Molecules.molecule[{i}].atom[{j}].Atom_Name")
                atype = u.get(f"Set_of_Molecules.molecule[{i}].atom[{j}].Atom_Type_Name")
                atoms.append(AtomInfo(atom_id=aid, atom_name=name, atom_type=atype))

            n_bonds = u.size(f"Set_of_Molecules.molecule[{i}].bond[]") or 0
            bonds = []
            for k in range(n_bonds):
                bond = u.get(f"Set_of_Molecules.molecule[{i}].bond[{k}]")
                bonds.append(BondInfo(
                    type_name=bond[0], atom1=int(bond[1]),
                    atom2=int(bond[2]), order=float(bond[3])
                ))
            self.molecules.append(MoleculeTopology(
                mol_name=mol_name, atoms=atoms, bonds=bonds
            ))

    def get_cell(self, record: int) -> CellBox:
        """Get cell for given record. Falls back to Initial_Structure if absent."""
        u = self._udf
        u.jump(record)
        try:
            cell = u.get("Structure.Unit_Cell.Cell_Size")
            if cell is not None and len(cell) >= 3:
                return CellBox(a=float(cell[0]), b=float(cell[1]), c=float(cell[2]))
        except Exception:
            pass
        init = u.get("Initial_Structure.Initial_Unit_Cell")
        if init is None:
            raise RuntimeError("Cannot read cell geometry")
        cell_arr = init[1]
        return CellBox(a=float(cell_arr[0]), b=float(cell_arr[1]), c=float(cell_arr[2]))

    def get_frame(self, record: int) -> TrajectoryFrame:
        """Read one frame's positions + cell."""
        u = self._udf
        u.jump(record)
        cell = self.get_cell(record)
        positions = []
        for i in range(len(self.molecules)):
            n_atoms = len(self.molecules[i].atoms)
            coords = np.zeros((n_atoms, 3), dtype=float)
            for j in range(n_atoms):
                pos = u.get(f"Structure.Position.mol[{i}].atom[{j}]")
                coords[j] = [float(pos[0]), float(pos[1]), float(pos[2])]
            positions.append(coords)
        return TrajectoryFrame(record=record, cell=cell, positions=positions)

    def iter_frames(self, start: int = 0, end: int = -1):
        """Iterate frames in [start, end). If end<0, iterate all."""
        if end < 0:
            end = self.n_records
        for rec in range(start, end):
            yield self.get_frame(rec)
