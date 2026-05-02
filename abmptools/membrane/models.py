# -*- coding: utf-8 -*-
"""
abmptools.membrane.models
-------------------------
Dataclasses for lipid-bilayer umbrella-sampling builds.

Sister to :mod:`abmptools.amorphous.models`. Same dataclass-based
configuration pattern with JSON round-trip.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import List, Literal, Optional


Backend = Literal["amber", "charmm36"]


@dataclass
class LipidSpec:
    """One lipid species in the bilayer.

    Attributes
    ----------
    resname : str
        Residue name as recognised by the chosen backend.
        AMBER Lipid21 uses head/tail-split residue names internally
        (e.g. ``POPC`` is split into PC + OL + PA), but for user-facing
        config we accept the **whole-lipid** name (``POPC``, ``DOPC``,
        ``DPPC``, ``POPE``, ``POPG``, ``CHL1``, …) and let the backend
        handle the split via packmol-memgen / charmmlipid2amber.
    n_per_leaflet : int
        Number of lipids of this species per leaflet (top / bottom).
        Total in box is ``2 * n_per_leaflet``.
    """
    resname: str
    n_per_leaflet: int


@dataclass
class PeptideSpec:
    """Peptide permeant.

    Provide one of *sequence* (one-letter code) or *pdb_path* (pre-built
    structure). If only *sequence* is given, the backend will build a
    capped, extended-chain initial structure via tleap (AMBER) or
    pdb2gmx (CHARMM36).

    Attributes
    ----------
    name : str
        Short tag used in filenames (e.g. ``"WT_helix"``).
    sequence : str
        One-letter amino-acid code, e.g. ``"AAAAA"`` for poly-Ala 5-mer.
        Empty if *pdb_path* is given.
    pdb_path : str
        Path to a pre-built peptide PDB. Empty if *sequence* is given.
    n_copies : int
        Usually 1 for PMF (single permeant). Multiple copies make z-axis
        crowded; supported but discouraged.
    cap_n : str
        N-terminal cap. ``""`` (free amine), ``"ACE"`` (acetyl).
    cap_c : str
        C-terminal cap. ``""`` (free carboxylate), ``"NME"`` (N-methyl).
    initial_z_offset_nm : float
        Initial peptide centre-of-mass z relative to bilayer centre.
        Default 3.0 nm (well above the upper leaflet) so pulling pulls
        the peptide *into* the membrane.
    """
    name: str = "peptide"
    sequence: str = ""
    pdb_path: str = ""
    n_copies: int = 1
    cap_n: str = "ACE"
    cap_c: str = "NME"
    initial_z_offset_nm: float = 3.0

    def __post_init__(self) -> None:
        if not self.sequence and not self.pdb_path:
            raise ValueError(
                "PeptideSpec requires either 'sequence' or 'pdb_path'."
            )
        if self.sequence and self.pdb_path:
            raise ValueError(
                "PeptideSpec accepts exactly one of 'sequence' or 'pdb_path'."
            )


@dataclass
class IonSpec:
    """Ion / salt configuration."""
    cation: str = "Na+"          # AMBER name; CHARMM uses "SOD"
    anion: str = "Cl-"           # AMBER name; CHARMM uses "CLA"
    salt_concentration_M: float = 0.15
    neutralize: bool = True      # add ions to net-neutralise


@dataclass
class USProtocol:
    """Umbrella-sampling protocol parameters."""
    # Reaction-coordinate range (relative to bilayer centre, nm).
    z_min_nm: float = -3.0
    z_max_nm: float = +3.0
    window_spacing_nm: float = 0.1
    # Per-window settings.
    force_constant_kj_mol_nm2: float = 1000.0
    window_nsteps: int = 10_000_000     # 20 ns at dt=2 fs
    window_dt_ps: float = 0.002
    window_nstxtcout: int = 5000
    window_nstenergy: int = 5000
    # Pull-coord definitions.
    pull_group1: str = "Bilayer"
    pull_group2: str = "Peptide"
    # ``direction`` gives a *signed* projection along ``pull_vec``
    # (peptide_z - bilayer_z when vec=0 0 1) and is compatible with
    # semiisotropic pressure coupling (which scales z dynamically).
    # ``direction-periodic`` is incompatible: GROMACS rejects it with
    # "Can not have dynamic box while using pull geometry
    # 'direction-periodic' (dim z)". ``distance`` returns |z| (always
    # positive) and is unsuitable for two-sided US.
    pull_geometry: str = "direction"
    pull_vec: str = "0 0 1"             # pulling direction (lab frame)
    pull_dim: str = "N N Y"             # only z (for COM calc)

    @property
    def n_windows(self) -> int:
        """Number of windows (inclusive of both endpoints)."""
        span = self.z_max_nm - self.z_min_nm
        return int(round(span / self.window_spacing_nm)) + 1


@dataclass
class EquilibrationProtocol:
    """Pre-pulling equilibration MD protocol (semiisotropic for membrane)."""
    em_steps: int = 50000
    em_tol: float = 1000.0
    nvt_nsteps: int = 100_000          # 0.1 ns (heat to 310 K)
    npt_nsteps: int = 1_000_000        # 2 ns (relax box)
    dt_ps: float = 0.002
    temperature_K: float = 310.0
    pressure_bar: float = 1.0
    tau_t_ps: float = 1.0
    tau_p_ps: float = 5.0
    nstxout_compressed: int = 5000
    nstenergy: int = 1000
    gen_seed: Optional[int] = None


@dataclass
class PullingProtocol:
    """Reaction-coordinate generation (peptide pulled into membrane)."""
    pull_rate_nm_per_ps: float = 0.001    # 1 nm/ns
    pull_force_constant: float = 1000.0
    nsteps: int = 5_000_000               # 10 ns
    nstxout_compressed: int = 500         # finer than equilibration


@dataclass
class MembraneConfig:
    """Top-level configuration for a membrane US build."""
    # --- backend ---
    backend: Backend = "amber"

    # --- system composition ---
    lipids: List[LipidSpec] = field(default_factory=list)
    peptide: Optional[PeptideSpec] = None
    ions: IonSpec = field(default_factory=IonSpec)

    # --- box geometry ---
    box_xy_nm: float = 8.0           # square xy
    water_thickness_nm: float = 2.5  # per side
    distance_to_lipid_nm: float = 1.5  # peptide initial distance from membrane

    # --- protocols ---
    equilibration: EquilibrationProtocol = field(default_factory=EquilibrationProtocol)
    pulling: PullingProtocol = field(default_factory=PullingProtocol)
    umbrella: USProtocol = field(default_factory=USProtocol)

    # --- reproducibility / paths ---
    seed: Optional[int] = None
    output_dir: str = "."

    # --- backend-specific overrides ---
    # AMBER
    amber_protein_ff: str = "leaprc.protein.ff19SB"
    amber_lipid_ff: str = "leaprc.lipid21"
    amber_water_ff: str = "leaprc.water.tip3p"
    # CHARMM36
    charmm_ff_dir: str = ""          # path to charmm36-jul2022.ff (Klauda port)

    # --- external tool paths ---
    packmol_memgen_path: str = "packmol-memgen"
    tleap_path: str = "tleap"
    gmx_path: str = "gmx"

    def __post_init__(self) -> None:
        if self.backend not in ("amber", "charmm36"):
            raise ValueError(
                f"backend must be 'amber' or 'charmm36', got {self.backend!r}"
            )
        if not self.lipids:
            raise ValueError("MembraneConfig.lipids must not be empty.")
        if self.peptide is None:
            raise ValueError(
                "MembraneConfig.peptide is required (no permeant = no PMF)."
            )
        if self.backend == "charmm36" and not self.charmm_ff_dir:
            raise ValueError(
                "backend='charmm36' requires charmm_ff_dir "
                "(path to charmm36-jul2022.ff). Set up the Klauda port "
                "first; CHARMM-GUI auto-generation is not permitted by "
                "this package's commercial-use license rule."
            )

    def to_json(self, path: str) -> None:
        """Save configuration to JSON for reproducibility."""
        Path(path).write_text(
            json.dumps(asdict(self), indent=2, ensure_ascii=False)
        )

    @classmethod
    def from_json(cls, path: str) -> "MembraneConfig":
        """Load configuration from JSON."""
        raw = json.loads(Path(path).read_text())
        lipids = [LipidSpec(**l) for l in raw.pop("lipids", [])]
        pep_raw = raw.pop("peptide", None)
        peptide = PeptideSpec(**pep_raw) if pep_raw else None
        ions = IonSpec(**raw.pop("ions"))
        eq = EquilibrationProtocol(**raw.pop("equilibration"))
        pull = PullingProtocol(**raw.pop("pulling"))
        umb = USProtocol(**raw.pop("umbrella"))
        return cls(
            lipids=lipids,
            peptide=peptide,
            ions=ions,
            equilibration=eq,
            pulling=pull,
            umbrella=umb,
            **raw,
        )
