# -*- coding: utf-8 -*-
"""
top_model.py
------------
Intermediate representation for GROMACS TOP/ITP + GRO → UDF conversion.

All data needed to write a COGNAC UDF is gathered into :class:`TopModel`
before any UDFManager calls are made.  This cleanly separates the parsing
stage from the writing stage.

Reuses CellGeometry from abmptools.core.system_model.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------

#: Boltzmann constant in COGNAC internal units [amu * Å² / (ps² * K)]
#: k_B = R / N_A  converted from (kg·m²/s²/K) → (amu·Å²/ps²/K)
#: = 1.38064852e-23 J/K / (1.66053906660e-23 J / (amu·Å²/ps²))
#: = 0.83144626
KB_AMU_A2_PS2_K: float = 0.83144626

#: Fallback Ewald R_cutoff [Å] used only when GRO frames are unavailable.
#: Under normal circumstances, R_cutoff is computed by the Deserno & Holm
#: (JCP 1998) formula: R_cutoff = sqrt(11.5) / (sqrt(pi)*(5.5N/V²)^(1/6))
#: (result converted nm → Å by ×10).
EWALD_R_CUTOFF_DEFAULT: float = 10.0

from ..core.system_model import CellGeometry


# ---------------------------------------------------------------------------
# Built-in atomic-mass → element-symbol map (avoids external conf/mass.conf)
# ---------------------------------------------------------------------------
_MASS_TO_ELEMENT: Dict[int, str] = {
    1: "H",   2: "He",  3: "Li",  4: "Be",  5: "B",
    6: "Li",  7: "Li",  9: "Be",  10: "B",  11: "B",
    12: "C",  14: "N",  16: "O",  19: "F",  20: "Ne",
    23: "Na", 24: "Mg", 27: "Al", 28: "Si", 31: "P",
    32: "S",  35: "Cl", 36: "Ar", 39: "K",  40: "Ca",
    45: "Sc", 48: "Ti", 51: "V",  52: "Cr", 55: "Mn",
    56: "Fe", 59: "Co", 58: "Ni", 64: "Cu", 65: "Zn",
    80: "Br", 108: "Ag", 127: "I",
}
# Patch ambiguous entries for common force-field masses
_MASS_TO_ELEMENT[1] = "H"    # 1 amu  → H
_MASS_TO_ELEMENT[12] = "C"   # 12 amu → C
_MASS_TO_ELEMENT[14] = "N"   # 14 amu → N
_MASS_TO_ELEMENT[16] = "O"   # 16 amu → O


def mass_to_element(mass: float) -> str:
    """Return the element symbol for *mass* [amu], or ``"X"`` if unknown."""
    m = round(mass)
    return _MASS_TO_ELEMENT.get(m, "X")


# ---------------------------------------------------------------------------
# Type-level data (global potential definitions)
# ---------------------------------------------------------------------------

@dataclass
class AtomTypeSpec:
    """One entry from ``[ atomtypes ]``."""
    name: str       # e.g. "ca"
    mass: float     # [amu]
    sigma: float    # [nm]
    epsilon: float  # [kJ/mol]


@dataclass
class BondTypeSpec:
    """One unique Harmonic bond potential type."""
    name: str       # e.g. "ca-ca-0"  (getBondName + "-" + type_index)
    name1: str      # atom type 1
    name2: str      # atom type 2
    funct: int      # 1 = Harmonic
    r0: float       # [nm]
    kb: float       # [kJ/mol/nm^2]


@dataclass
class AngleTypeSpec:
    """One unique angle potential type."""
    name: str       # e.g. "ca-ca-ca-0"
    name1: str
    name2: str
    name3: str
    funct: int      # 1 = Theta
    theta0: float   # [deg] (raw value from TOP file)
    k: float        # [kJ/mol/rad^2]


@dataclass
class TorsionTypeSpec:
    """One unique torsion potential type (may be Amber multi-term)."""
    name: str           # base name, e.g. "ca-ca-ca-ca-0" or "ca-ca-ca-ha-4-oopa"
    name1: str
    name2: str
    name3: str
    name4: str
    funct: int          # 1, 3, 4, or 9
    improper: bool
    params: List[float] # raw flat list from TOP


# ---------------------------------------------------------------------------
# Per-molecule topology
# ---------------------------------------------------------------------------

@dataclass
class MolAtomSpec:
    """One atom in a molecule topology."""
    index_1based: int   # 1-based index within the molecule
    atom_name: str      # original GROMACS atom name from [ atoms ], e.g. "ca1"
    element: str        # element symbol from mass mapping, e.g. "C"
    type_name: str      # force-field type, e.g. "ca"
    charge: float       # partial charge [e]
    global_atom_id: int # 0-based sequential global atom ID (assigned by adapter)


@dataclass
class MolBondSpec:
    """Bond entry in a molecule topology."""
    atom1: int          # 1-based within molecule
    atom2: int
    type_index: int     # index into TopModel.bond_type_specs
    potential_name: str # e.g. "ca-ca-0"


@dataclass
class MolAngleSpec:
    """Angle entry in a molecule topology."""
    atom1: int
    atom2: int
    atom3: int
    type_index: int
    potential_name: str


@dataclass
class MolTorsionSpec:
    """Torsion entry in a molecule topology (one row per potential term)."""
    atom1: int
    atom2: int
    atom3: int
    atom4: int
    funct: int
    improper: bool
    n_params: int           # number of params in the type (used for multi-term)
    potential_name: str     # full name, e.g. "ca-ca-ca-ca-0" or "ca-ca-ca-ca-0:1"


@dataclass
class MolSpec:
    """Topology for one unique molecule type."""
    name: str
    atoms: List[MolAtomSpec] = field(default_factory=list)
    bonds: List[MolBondSpec] = field(default_factory=list)
    angles: List[MolAngleSpec] = field(default_factory=list)
    torsions: List[MolTorsionSpec] = field(default_factory=list)


# ---------------------------------------------------------------------------
# GRO frame
# ---------------------------------------------------------------------------

@dataclass
class GROFrameData:
    """One snapshot read from a .gro file."""
    step: int
    time: float             # [ps]
    coord_list: List[List[float]]   # [[x,y,z], ...] in [nm], one per atom
    cell: List[float]               # 3 values [nm] (orthogonal box)


# ---------------------------------------------------------------------------
# Top-level intermediate representation
# ---------------------------------------------------------------------------

@dataclass
class TopModel:
    """
    Complete intermediate representation for GROMACS TOP+GRO → UDF conversion.

    Produced by :class:`~abmptools.gro2udf.top_adapter.TopAdapter`.
    Consumed by :class:`~abmptools.gro2udf.top_exporter.TopExporter`.
    """

    # --- force-field defaults ---
    comb_rule: int          # 2=Lorentz-Berthelot, 3=geometric (OPLS)
    fudge_lj: float
    fudge_qq: float

    # --- global type definitions ---
    atom_type_specs: List[AtomTypeSpec]
    bond_type_specs: List[BondTypeSpec]
    angle_type_specs: List[AngleTypeSpec]
    torsion_type_specs: List[TorsionTypeSpec]
    mass_dict: Dict[str, float]             # type_name -> mass [amu]

    # --- topology ---
    mol_type_names: List[str]               # unique mol type names in order
    mol_specs: List[MolSpec]               # one MolSpec per unique type
    mol_instance_list: List[str]           # flat list of all instances

    # --- structure frames from GRO ---
    frames: List[GROFrameData] = field(default_factory=list)

    # --- simulation parameters (populated from .mdp when available) ---
    #: Total number of atoms in the system (sum over all molecule instances)
    n_atoms_total: int = 0
    #: Reference temperature [K]  (from mdp: ref_t)
    ref_t: float = 300.0
    #: Nose-Hoover relaxation time [ps]  (from mdp: tau_t)
    tau_t: float = 0.1
    #: Ewald real-space cutoff [Å]  (Deserno & Holm formula from GRO box, or default)
    ewald_r_cutoff: float = EWALD_R_CUTOFF_DEFAULT
