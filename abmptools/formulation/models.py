# -*- coding: utf-8 -*-
"""Dataclasses for AA peptide-formulation mixed-solution builds.

The schema mirrors Hossain et al. 2023 (Nanoscale 15, 19180) but
uses commercially-permissive force fields (ff14SB + GAFF2 + TIP3P
+ Joung-Cheatham).
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple


# ---------------------------------------------------------------------------
# Component specs
# ---------------------------------------------------------------------------


@dataclass
class PeptideSpec:
    """One peptide species in the formulation box.

    Provide exactly one of *sequence* (one-letter AAs, natural only)
    or *pdb_path* (pre-built PDB, e.g. RCSB 2G4M for insulin).

    disulfide_bonds entries use the format ``"CHAIN:RESID:ATOM"``,
    e.g. ``("A:6:SG", "A:11:SG")``. They are emitted as tleap
    ``bond`` directives on each replicate; per-copy chain renaming is
    handled by :mod:`abmptools.formulation.peptide_atomistic`.
    """

    name: str = "peptide"
    sequence: str = ""
    pdb_path: str = ""
    n_copies: int = 1
    disulfide_bonds: List[Tuple[str, str]] = field(default_factory=list)
    cap_n: str = ""  # "ACE" or empty
    cap_c: str = ""  # "NME" or empty
    # "tleap" (default, ff14SB, natural L-AA only) or "gaff"
    # (whole-peptide → acpype/GAFF2 + AM1-BCC, accepts arbitrary
    # PDB including D-amino acids and non-standard residues).
    parameterize_method: str = "tleap"
    resname: str = ""             # used in GAFF mode (3-letter tag)
    net_charge: int = 0           # for acpype -n in GAFF mode

    def __post_init__(self) -> None:
        sources = [self.sequence, self.pdb_path]
        non_empty = [s for s in sources if s]
        if not non_empty:
            raise ValueError(
                f"PeptideSpec[{self.name}] requires exactly one of "
                "'sequence' or 'pdb_path'."
            )
        if len(non_empty) > 1:
            raise ValueError(
                f"PeptideSpec[{self.name}] accepts exactly one of "
                "'sequence' or 'pdb_path'."
            )
        if self.n_copies < 1:
            raise ValueError(
                f"PeptideSpec[{self.name}] n_copies must be >= 1, "
                f"got {self.n_copies}."
            )
        for pair in self.disulfide_bonds:
            if len(pair) != 2:
                raise ValueError(
                    f"PeptideSpec[{self.name}] disulfide_bonds entries "
                    f"must be 2-tuples of 'CHAIN:RESID:ATOM' strings, "
                    f"got {pair!r}."
                )
        if self.parameterize_method not in ("tleap", "gaff"):
            raise ValueError(
                f"PeptideSpec[{self.name}] parameterize_method must be "
                f"'tleap' or 'gaff', got {self.parameterize_method!r}."
            )
        if self.parameterize_method == "gaff":
            if not self.pdb_path:
                raise ValueError(
                    f"PeptideSpec[{self.name}] GAFF mode requires "
                    "pdb_path (whole-peptide PDB to feed to acpype)."
                )
            if not self.resname:
                raise ValueError(
                    f"PeptideSpec[{self.name}] GAFF mode requires "
                    "resname (3-letter tag for acpype mol2/itp)."
                )


@dataclass
class EnhancerSpec:
    """A permeation enhancer (e.g. sodium caprate, SNAC).

    Charged and neutral copies are tracked separately so that the
    fatty-acid pKa-shift convention from Hossain 2023 (50:50
    neutral:charged at intestinal pH) can be modelled. Each form is
    parameterised independently via acpype, since charge differences
    produce different AM1-BCC partial charges.

    *resname* is the 3-letter tag written into the .mol2 / .itp /
    output PDB. It must be unique across all species in a build (the
    builder validates this).
    """

    name: str = "enhancer"
    smiles_neutral: str = ""
    smiles_charged: str = ""
    n_neutral: int = 0
    n_charged: int = 0
    resname: str = ""
    pdb_path_neutral: str = ""  # alternative to smiles
    pdb_path_charged: str = ""

    def __post_init__(self) -> None:
        if not self.resname:
            raise ValueError(f"EnhancerSpec[{self.name}] resname required.")
        if len(self.resname) > 4:
            raise ValueError(
                f"EnhancerSpec[{self.name}] resname must be <= 4 chars."
            )
        if self.n_neutral < 0 or self.n_charged < 0:
            raise ValueError(
                f"EnhancerSpec[{self.name}] n_neutral / n_charged must be >= 0."
            )
        if self.n_neutral + self.n_charged == 0:
            raise ValueError(
                f"EnhancerSpec[{self.name}] n_neutral + n_charged must be > 0."
            )
        if self.n_neutral > 0 and not (self.smiles_neutral or self.pdb_path_neutral):
            raise ValueError(
                f"EnhancerSpec[{self.name}] n_neutral={self.n_neutral} but no "
                "smiles_neutral / pdb_path_neutral given."
            )
        if self.n_charged > 0 and not (self.smiles_charged or self.pdb_path_charged):
            raise ValueError(
                f"EnhancerSpec[{self.name}] n_charged={self.n_charged} but no "
                "smiles_charged / pdb_path_charged given."
            )


@dataclass
class BileSaltSpec:
    """A bile salt species (e.g. taurocholate).

    *net_charge* is passed to acpype's ``-n`` flag.
    """

    name: str = "bile_salt"
    smiles: str = ""
    pdb_path: str = ""
    n_copies: int = 0
    resname: str = ""
    net_charge: int = 0

    def __post_init__(self) -> None:
        if not self.resname:
            raise ValueError(f"BileSaltSpec[{self.name}] resname required.")
        if len(self.resname) > 4:
            raise ValueError(
                f"BileSaltSpec[{self.name}] resname must be <= 4 chars."
            )
        if self.n_copies < 0:
            raise ValueError(
                f"BileSaltSpec[{self.name}] n_copies must be >= 0."
            )
        if self.n_copies > 0 and not (self.smiles or self.pdb_path):
            raise ValueError(
                f"BileSaltSpec[{self.name}] n_copies={self.n_copies} but "
                "neither smiles nor pdb_path given."
            )


# ---------------------------------------------------------------------------
# Composite specs
# ---------------------------------------------------------------------------


@dataclass
class SystemSpec:
    """Composition + geometry of one formulation box.

    *packmol_layout* controls the initial geometry:

    - ``"random"`` (default): all species are packed uniformly into the
      cubic box (Hossain 2023 starting configuration A).
    - ``"cluster_core"``: peptide(s) are placed at the box centre
      (within *cluster_core_radius_nm*), enhancers are packed into
      the surrounding shell (between *cluster_core_radius_nm* and
      *cluster_shell_radius_nm*), and bile salts + remaining volume
      are filled outside. This produces a **pre-formed micelle-like
      aggregate** at t=0, suitable for short-production release-PMF
      smokes (Hossain 2023 starting configuration B).
    """

    peptides: List[PeptideSpec] = field(default_factory=list)
    enhancers: List[EnhancerSpec] = field(default_factory=list)
    bile_salts: List[BileSaltSpec] = field(default_factory=list)
    box_size_nm: float = 10.0
    salt_concentration_M: float = 0.15
    neutralize: bool = True
    packmol_layout: str = "random"        # "random" | "cluster_core"
    cluster_core_radius_nm: float = 1.0   # peptide placement sphere radius
    cluster_shell_radius_nm: float = 3.0  # enhancer outer sphere radius

    def __post_init__(self) -> None:
        if not self.peptides:
            raise ValueError("SystemSpec requires at least one PeptideSpec.")
        if self.box_size_nm <= 0:
            raise ValueError(
                f"SystemSpec box_size_nm must be > 0, got {self.box_size_nm}."
            )
        if self.salt_concentration_M < 0:
            raise ValueError(
                f"SystemSpec salt_concentration_M must be >= 0, "
                f"got {self.salt_concentration_M}."
            )
        resnames: List[str] = []
        for spec in self.enhancers:
            resnames.append(spec.resname)
        for spec in self.bile_salts:
            resnames.append(spec.resname)
        if len(set(resnames)) != len(resnames):
            raise ValueError(
                f"SystemSpec resnames must be unique across enhancers/bile "
                f"salts, got duplicates in {resnames}."
            )
        if self.packmol_layout not in ("random", "cluster_core"):
            raise ValueError(
                f"SystemSpec packmol_layout must be 'random' or "
                f"'cluster_core', got {self.packmol_layout!r}."
            )
        if self.cluster_core_radius_nm <= 0:
            raise ValueError(
                f"SystemSpec cluster_core_radius_nm must be > 0, "
                f"got {self.cluster_core_radius_nm}."
            )
        if self.cluster_shell_radius_nm <= self.cluster_core_radius_nm:
            raise ValueError(
                "SystemSpec cluster_shell_radius_nm must be > "
                "cluster_core_radius_nm."
            )

    @property
    def total_peptide_copies(self) -> int:
        return sum(p.n_copies for p in self.peptides)


@dataclass
class EquilibrationProtocol:
    """em + nvt + npt protocol (Hossain 2023 defaults)."""

    em_steps: int = 10000
    em_tol: float = 1000.0
    nvt_nsteps: int = 50000   # 100 ps at dt=2 fs
    npt_nsteps: int = 50000
    dt_ps: float = 0.002
    temperature_K: float = 310.0       # 37 °C
    pressure_bar: float = 1.0
    tau_t_ps: float = 0.5
    tau_p_ps: float = 5.0
    compressibility: float = 4.5e-5
    nstxout_compressed: int = 5000     # 10 ps stride
    nstenergy: int = 1000
    gen_seed: Optional[int] = None


@dataclass
class ProductionProtocol:
    """Production NPT (Parrinello-Rahman, isotropic)."""

    nsteps: int = 500000              # 1 ns smoke at dt=2 fs
    dt_ps: float = 0.002
    temperature_K: float = 310.0
    pressure_bar: float = 1.0
    tau_t_ps: float = 1.0
    tau_p_ps: float = 5.0
    compressibility: float = 4.5e-5
    nstxout_compressed: int = 5000
    nstenergy: int = 1000


@dataclass
class USProtocol:
    """Umbrella sampling for peptide release from an aggregate."""

    n_windows: int = 12
    window_spacing_nm: float = 0.2
    window_nsteps: int = 6000000      # 12 ns at dt=2 fs
    window_dt_ps: float = 0.002
    force_constant_kj_mol_nm2: float = 1000.0
    pull_rate_nm_per_ps: float = -0.001
    pull_nsteps: int = 5000000        # 10 ns pull
    target_peptide_idx: int = 0
    enhancer_tc_group: str = "peptide"  # "peptide" | "solvent"


# ---------------------------------------------------------------------------
# Top-level build config (JSON round-trippable)
# ---------------------------------------------------------------------------


@dataclass
class FormulationBuildConfig:
    """End-to-end configuration for a formulation build."""

    system: SystemSpec
    equilibration: EquilibrationProtocol = field(
        default_factory=EquilibrationProtocol
    )
    production: ProductionProtocol = field(default_factory=ProductionProtocol)
    release_us: Optional[USProtocol] = None
    seed: Optional[int] = None
    output_dir: str = "."

    # tool paths
    tleap_path: str = "tleap"
    acpype_path: str = "acpype"
    packmol_path: str = "packmol"
    gmx_path: str = "gmx"

    # force-field selection (commercial-permissive defaults)
    amber_protein_ff: str = "leaprc.protein.ff14SB"
    amber_water_ff: str = "leaprc.water.tip3p"
    gaff_version: str = "gaff2"

    # packmol behaviour
    packmol_tolerance_A: float = 2.0
    packmol_inner_box_margin_nm: float = 0.5  # pack to (box - margin), let solvatebox fill

    def to_json(self, path: str) -> None:
        Path(path).write_text(
            json.dumps(_to_jsonable(self), indent=2, ensure_ascii=False)
        )

    @classmethod
    def from_json(cls, path: str) -> "FormulationBuildConfig":
        data = json.loads(Path(path).read_text())
        return cls._from_dict(data)

    @classmethod
    def _from_dict(cls, data: Dict[str, Any]) -> "FormulationBuildConfig":
        sys_dict = data["system"]
        peptides = [
            PeptideSpec(
                **{
                    **p,
                    "disulfide_bonds": [tuple(pair) for pair in p.get("disulfide_bonds", [])],
                }
            )
            for p in sys_dict.get("peptides", [])
        ]
        enhancers = [EnhancerSpec(**e) for e in sys_dict.get("enhancers", [])]
        bile_salts = [BileSaltSpec(**b) for b in sys_dict.get("bile_salts", [])]
        sys_kwargs = {
            k: v for k, v in sys_dict.items()
            if k not in ("peptides", "enhancers", "bile_salts")
        }
        system = SystemSpec(
            peptides=peptides, enhancers=enhancers, bile_salts=bile_salts,
            **sys_kwargs,
        )
        equil = EquilibrationProtocol(**data.get("equilibration", {}))
        prod = ProductionProtocol(**data.get("production", {}))
        us = None
        if data.get("release_us"):
            us = USProtocol(**data["release_us"])
        top = {
            k: v for k, v in data.items()
            if k not in (
                "system", "equilibration", "production", "release_us"
            )
        }
        return cls(
            system=system, equilibration=equil, production=prod,
            release_us=us, **top,
        )


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _to_jsonable(obj: Any) -> Any:
    if hasattr(obj, "__dataclass_fields__"):
        return {k: _to_jsonable(v) for k, v in asdict(obj).items()}
    if isinstance(obj, (list, tuple)):
        return [_to_jsonable(v) for v in obj]
    if isinstance(obj, dict):
        return {k: _to_jsonable(v) for k, v in obj.items()}
    return obj
