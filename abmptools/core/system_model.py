# -*- coding: utf-8 -*-
"""
system_model.py
---------------
Intermediate representation (pure data) for UDF → Gromacs conversion.

Design rule: No Gromacs-specific string syntax here.
All formatting is handled by the Writer layer.
"""
from __future__ import annotations
from dataclasses import dataclass, field
from typing import List, Optional, Dict, Tuple


__all__ = [
    "AtomType",
    "AtomRecord",
    "BondRecord",
    "PairRecord",
    "AngleRecord",
    "DihedralRecord",
    "MoleculeTopology",
    "AtomPosition",
    "CellGeometry",
    "NdxData",
    "SimulationParams",
    "SystemModel",
]


# ---------------------------------------------------------------------------
# Atom-type level
# ---------------------------------------------------------------------------

@dataclass
class AtomType:
    """One row in [ atomtypes ]."""
    name: str           # e.g. "ca", "ha"
    mass: float         # [amu]
    sigma: float        # [nm]  (0.0 if self-interaction not found)
    epsilon: float      # [kJ/mol]  (0.0 if not found)


# ---------------------------------------------------------------------------
# Per-molecule topology
# ---------------------------------------------------------------------------

@dataclass
class AtomRecord:
    """One row in [ atoms ]."""
    index: int          # 1-based serial number within molecule
    type_name: str      # force-field atom type (e.g. "ca")
    gro_name: str       # name in gro/top atom column (first2 chars + hex index)
    charge: float       # [e]


@dataclass
class BondRecord:
    """One row in [ bonds ]."""
    atom1: int          # 1-based
    atom2: int          # 1-based
    funct: str          # "1" (Harmonic) or "7" (FENE_LJ)
    r0: float           # [nm]
    kb: float           # [kJ/mol/nm^2]


@dataclass
class PairRecord:
    """One row in [ pairs ]."""
    atom1: int          # 1-based
    atom2: int          # 1-based


@dataclass
class AngleRecord:
    """One row in [ angles ]."""
    atom1: int
    atom2: int
    atom3: int
    theta0: float       # [deg] (already converted: 180 - theta0_udf)
    k: float            # [kJ/mol/rad^2]


@dataclass
class DihedralRecord:
    """One row in [ dihedrals ].

    funct values:
      "1"  Proper Dih. (single)
      "3"  Ryckaert-Bellemans (C0..C5)
      "4"  Improper Dih. (periodic)
      "9"  Proper Dih. (multiple, Amber)
    params is interpreted by funct:
      funct "1","4","9": [phi_s, k, mult]
      funct "3":         [C0, C1, C2, C3, C4, C5]
    """
    atom1: int
    atom2: int
    atom3: int
    atom4: int
    funct: str
    params: List[float]


@dataclass
class MoleculeTopology:
    """Topology for one molecule species."""
    udf_name: str           # original Mol_Name in UDF
    gro_name: str           # Gromacs name, e.g. "M0000"
    nrexcl: int             # exclusion distance (always 3)
    atoms: List[AtomRecord] = field(default_factory=list)
    bonds: List[BondRecord] = field(default_factory=list)
    pairs: List[PairRecord] = field(default_factory=list)
    angles: List[AngleRecord] = field(default_factory=list)
    dihedrals: List[DihedralRecord] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Per-atom positions (for .gro)
# ---------------------------------------------------------------------------

@dataclass
class AtomPosition:
    """One atom entry in the .gro body."""
    mol_id: int             # 1-based (capped at 99999)
    mol_name_short: str     # ≤5 chars Gromacs mol name
    atom_gro_name: str      # atomname_in_gro (≤5 chars)
    atom_id: int            # 1-based Atom_ID (capped at 99999)
    x: float                # [nm]
    y: float                # [nm]
    z: float                # [nm]
    vx: float               # [nm/ps]
    vy: float               # [nm/ps]
    vz: float               # [nm/ps]


# ---------------------------------------------------------------------------
# Cell geometry
# ---------------------------------------------------------------------------

@dataclass
class CellGeometry:
    """Unit cell.  a,b,c in nm; alpha,beta,gamma in degrees."""
    a: float
    b: float
    c: float
    alpha: float = 90.0
    beta: float  = 90.0
    gamma: float = 90.0

    def is_rectangular(self, thres: float = 1e-5) -> bool:
        return (abs(self.alpha - 90.0) < thres and
                abs(self.beta  - 90.0) < thres and
                abs(self.gamma - 90.0) < thres)


# ---------------------------------------------------------------------------
# NDX data (constraint index file, optional)
# ---------------------------------------------------------------------------

@dataclass
class NdxData:
    """Data needed to write the .ndx file (only when constraints exist)."""
    atom_id_max: int
    mol_groups: Dict[str, List[int]]        # gro_molname -> sorted atom_id list
    constraint_atom_ids: List[int]
    constr_axis: Optional[List[str]]        # e.g. ["YES","YES","NO"]


# ---------------------------------------------------------------------------
# Simulation parameters (for .mdp)
# ---------------------------------------------------------------------------

@dataclass
class SimulationParams:
    """All data needed to produce the .mdp file."""
    title: str
    algorithm: str          # raw UDF value, e.g. "NPT_Andersen_Nose_Hoover"

    # --- time ---
    nsteps: int
    dt: float               # [ps]
    outputinterval: int     # Output_Interval_Steps
    outputinterval2: int    # = outputinterval // 10 if >= 10000 else same

    # --- integrator (derived) ---
    integrator: str         # "md-vv", "md", "sd", etc.
    ld_seed: Optional[int] = None   # for Kremer_Grest (sd)

    # --- velocity ---
    vel_gen: bool = False
    gen_temp: Optional[float] = None   # [K], only if vel_gen

    # --- constraints ---
    rattle_bond: bool = False
    rattle_angle: bool = False

    # --- electrostatics ---
    calcQQ: int = 0
    qq_algorithm: str = ""   # "Ewald", "Cutoff_Coulomb", ...

    # --- cutoffs ---
    lj_cutoff: float = 0.0     # [nm] max LJ cutoff (set by adapter)
    coulomb_cutoff: float = 0.0  # [nm]

    # --- temperature ---
    t_coupl: str = "no"        # "nose-hoover", "berendsen", "no"
    tau_t: float = 0.1         # [ps]
    ref_t: float = 300.0       # [K]

    # --- pressure ---
    p_coupl: str = "no"        # "MTTK", "Parrinello-Rahman", "berendsen", "no"
    pcoupltype: str = "isotropic"
    tau_p: float = 2.0         # [ps]
    ref_p: float = 1.0         # [bar]  (scalar for isotropic)
    ref_p_tensor: Optional[List[float]] = None  # 6-component for anisotropic
    compressibility: float = 4.5e-5
    compressibility_tensor: Optional[List[float]] = None  # 6-comp for anisotropic

    # --- tail correction / vdw ---
    tail_correction: int = 0

    # --- pbc ---
    pbc: str = "xyz"           # "xyz" or "no"
    periodic_mol: bool = False

    # --- deformation (optional) ---
    deform_vel: Optional[List[float]] = None   # 6-component [nm/ps]

    # --- constraint freeze (optional) ---
    freeze_grps: Optional[str] = None
    freeze_dim: Optional[str] = None


# ---------------------------------------------------------------------------
# Top-level system model
# ---------------------------------------------------------------------------

@dataclass
class SystemModel:
    """Complete intermediate representation of a UDF system."""

    # --- metadata ---
    title: str                  # e.g. "test.udf"
    udf_path: str               # full path to input UDF

    # --- force field defaults ---
    comb_rule: int              # 2=Lorentz-Berthelot, 3=geometric (OPLS)
    flags14: int                # 0 or 1
    fudgeLJ: float
    fudgeQQ: float
    calcQQ: int                 # 1=electrostatics ON

    # --- atom types ---
    atom_types: List[AtomType] = field(default_factory=list)

    # --- topology per molecule species ---
    mol_topologies: List[MoleculeTopology] = field(default_factory=list)

    # --- molecule sequence for [ molecules ] section ---
    # List of (gro_name, count) pairs in order of appearance
    mol_sequence: List[Tuple[str, int]] = field(default_factory=list)

    # --- structure (for .gro) ---
    atom_positions: List[AtomPosition] = field(default_factory=list)
    cell: Optional[CellGeometry] = None

    # --- simulation parameters (for .mdp) ---
    sim_params: Optional[SimulationParams] = None

    # --- constraint index (for .ndx, None if no constraints) ---
    ndx_data: Optional[NdxData] = None
