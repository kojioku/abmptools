# -*- coding: utf-8 -*-
"""
abmptools.cg.membrane.models
-----------------------------
Dataclasses for Martini 3 peptide-membrane US system building.

Schema follows :mod:`abmptools.cg.peptide.models` and
:mod:`abmptools.membrane.models` conventions: ``@dataclass`` with
``to_json/from_json`` round-trip, no Pydantic.
"""
from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import List, Optional

from abmptools.cg.peptide.models import STANDARD_AA


# ---------------------------------------------------------------------------
# Composition
# ---------------------------------------------------------------------------

@dataclass
class LipidMix:
    """A single lipid species in the bilayer.

    For v1, exactly one species is supported (POPC default). Multi-species
    mixtures are deferred to v1.20+: when the time comes this dataclass will
    grow ``ratio`` and the validation in :class:`MembraneCGBuildConfig` will
    accept multiple entries.

    Attributes
    ----------
    resname
        Martini 3 lipid residue name as found in
        ``martini_v3.0.0_phospholipids_v1.itp`` (e.g. ``POPC``).
    n_per_leaflet
        Target lipid count per leaflet. Total in box is ``2 * n_per_leaflet``.
        insane will place upper/lower symmetrically; some lipids may be
        displaced when a peptide is inserted.
    """
    resname: str = "POPC"
    n_per_leaflet: int = 64

    def __post_init__(self) -> None:
        if not self.resname:
            raise ValueError("LipidMix.resname must be non-empty")
        if self.n_per_leaflet < 1:
            raise ValueError(
                f"LipidMix.n_per_leaflet must be >= 1, got {self.n_per_leaflet}"
            )


@dataclass
class PeptideMembraneSpec:
    """Permeant peptide for the PMF.

    Provide one of *sequence* (one-letter code) or both *cg_pdb_path* +
    *cg_itp_path* (pre-built CG topology). When *sequence* is given,
    :class:`abmptools.cg.peptide.PeptideCGBuilder` is sub-called with
    ``solvent_enabled=False, mdp_*=False`` to produce the CG PDB and ITP.

    Attributes
    ----------
    name
        Short tag used in filenames (e.g. ``"kgg"``).
    sequence
        One-letter amino-acid code, e.g. ``"KGG"``. Empty if pre-built CG
        files are given.
    cg_pdb_path
        Path to pre-built CG PDB. Empty if *sequence* is given.
    cg_itp_path
        Path to pre-built CG ITP. Must be supplied alongside *cg_pdb_path*.
    initial_z_offset_nm
        Initial peptide centre-of-mass z relative to bilayer centre, passed
        to ``insane -dm``. Default 3.0 nm (well above upper leaflet).
    """
    name: str = "kgg"
    sequence: str = "KGG"
    cg_pdb_path: str = ""
    cg_itp_path: str = ""
    initial_z_offset_nm: float = 3.0

    def __post_init__(self) -> None:
        if not self.name:
            raise ValueError("PeptideMembraneSpec.name must be non-empty")

        has_seq = bool(self.sequence)
        has_pdb = bool(self.cg_pdb_path) or bool(self.cg_itp_path)
        if not has_seq and not has_pdb:
            raise ValueError(
                "PeptideMembraneSpec requires either 'sequence' or "
                "('cg_pdb_path' + 'cg_itp_path')."
            )
        if has_seq and has_pdb:
            raise ValueError(
                "PeptideMembraneSpec accepts exactly one of 'sequence' or "
                "pre-built CG files; got both."
            )
        if has_pdb and not (self.cg_pdb_path and self.cg_itp_path):
            raise ValueError(
                "When supplying pre-built CG files, both 'cg_pdb_path' and "
                "'cg_itp_path' must be set."
            )

        if has_seq:
            seq = self.sequence.upper()
            invalid = set(seq) - STANDARD_AA
            if invalid:
                raise ValueError(
                    f"Invalid amino acid(s) in sequence {seq!r}: "
                    f"{sorted(invalid)}. Supported: {''.join(sorted(STANDARD_AA))}"
                )
            self.sequence = seq


# ---------------------------------------------------------------------------
# Protocols
# ---------------------------------------------------------------------------

@dataclass
class EquilibrationCGProtocol:
    """CG equilibration protocol: em -> nvt -> npt-semiisotropic.

    Default values target Martini 3 / dt=20 fs / reaction-field (rcoulomb=1.1
    nm), with a 2-group thermostat (Bilayer / Non_Bilayer) and semiisotropic
    P-coupling (tau_p=12 ps standard for Martini).
    """
    em_steps: int = 5000
    em_tol: float = 100.0
    nvt_nsteps: int = 25_000          # 0.5 ns at dt=20 fs
    npt_nsteps: int = 250_000         # 5 ns at dt=20 fs
    dt_ps: float = 0.020
    temperature_K: float = 310.0
    pressure_bar: float = 1.0
    tau_t_ps: float = 1.0
    tau_p_ps: float = 12.0            # Martini standard
    nstxout_compressed: int = 5000
    nstenergy: int = 1000
    gen_seed: Optional[int] = None

    def __post_init__(self) -> None:
        if self.dt_ps <= 0:
            raise ValueError(f"dt_ps must be > 0, got {self.dt_ps}")
        if self.temperature_K <= 0:
            raise ValueError(
                f"temperature_K must be > 0, got {self.temperature_K}"
            )
        for label, ns in (
            ("em_steps", self.em_steps),
            ("nvt_nsteps", self.nvt_nsteps),
            ("npt_nsteps", self.npt_nsteps),
        ):
            if ns < 0:
                raise ValueError(f"{label} must be >= 0, got {ns}")


@dataclass
class PullingCGProtocol:
    """Reaction-coordinate generation (peptide pulled into/through bilayer).

    Defaults: 1 nm/ns rate, k=1000 kJ/mol/nm², 5 ns total (250,000 steps at
    dt=20 fs). At 1 nm/ns this sweeps roughly 5 nm of z, more than enough to
    cover the umbrella range (-1.5 to +1.5 nm).
    """
    pull_rate_nm_per_ps: float = 0.001
    pull_force_constant: float = 1000.0
    nsteps: int = 250_000
    nstxout_compressed: int = 500

    def __post_init__(self) -> None:
        if self.pull_force_constant < 0:
            raise ValueError(
                f"pull_force_constant must be >= 0, got "
                f"{self.pull_force_constant}"
            )
        if self.nsteps < 0:
            raise ValueError(f"nsteps must be >= 0, got {self.nsteps}")


@dataclass
class UmbrellaCGProtocol:
    """Umbrella-sampling protocol parameters.

    Defaults match the AA membrane convention scaled to CG dt=20 fs:
    13 windows × 1 ns × k=1000 kJ/mol/nm² × 0.25 nm spacing.
    Geometry is ``direction`` (lab-frame z projection), compatible with
    semiisotropic P-coupling. ``direction-periodic`` is rejected by GROMACS
    when the box is dynamic.
    """
    z_min_nm: float = -1.5
    z_max_nm: float = +1.5
    window_spacing_nm: float = 0.25
    force_constant_kj_mol_nm2: float = 1000.0
    window_nsteps: int = 50_000        # 1 ns at dt=20 fs
    window_dt_ps: float = 0.020
    window_nstxout_compressed: int = 1000
    window_nstenergy: int = 500
    pull_group1: str = "Bilayer"
    pull_group2: str = "Peptide"
    pull_geometry: str = "direction"
    pull_vec: str = "0 0 1"
    pull_dim: str = "N N Y"

    def __post_init__(self) -> None:
        if self.z_min_nm >= self.z_max_nm:
            raise ValueError(
                f"z_min_nm ({self.z_min_nm}) must be < z_max_nm "
                f"({self.z_max_nm})"
            )
        if self.window_spacing_nm <= 0:
            raise ValueError(
                f"window_spacing_nm must be > 0, got {self.window_spacing_nm}"
            )
        if self.force_constant_kj_mol_nm2 < 0:
            raise ValueError(
                f"force_constant_kj_mol_nm2 must be >= 0, got "
                f"{self.force_constant_kj_mol_nm2}"
            )
        if self.window_dt_ps <= 0:
            raise ValueError(
                f"window_dt_ps must be > 0, got {self.window_dt_ps}"
            )
        if self.window_nsteps < 0:
            raise ValueError(
                f"window_nsteps must be >= 0, got {self.window_nsteps}"
            )

    @property
    def n_windows(self) -> int:
        """Number of windows including both endpoints."""
        span = self.z_max_nm - self.z_min_nm
        return int(round(span / self.window_spacing_nm)) + 1


# ---------------------------------------------------------------------------
# Top-level config
# ---------------------------------------------------------------------------

@dataclass
class MembraneCGBuildConfig:
    """Full configuration for a Martini 3 peptide-membrane US build."""
    lipids: List[LipidMix] = field(default_factory=lambda: [LipidMix()])
    peptide: Optional[PeptideMembraneSpec] = None

    # Box geometry. insane computes xy from APL + lipid count + ``-d``;
    # ``box_z_nm`` is set explicitly so peptide+water+bilayer fits the
    # initial_z_offset on top.
    insane_d_nm: float = 8.0
    box_z_nm: float = 14.0
    insane_pbc: str = "hexagonal"
    insane_extra_args: List[str] = field(default_factory=list)

    # Solvent / ions.
    solvent_type: str = "W"            # Martini W
    nacl_molar: float = 0.15
    neutralize: bool = True

    # Protocols.
    equilibration: EquilibrationCGProtocol = field(
        default_factory=EquilibrationCGProtocol
    )
    pulling: PullingCGProtocol = field(default_factory=PullingCGProtocol)
    umbrella: UmbrellaCGProtocol = field(default_factory=UmbrellaCGProtocol)

    # Force field directory (4 ITPs incl. phospholipids).
    martini_itp_dir: str = ""

    # External tool paths.
    insane_path: str = "insane"
    martinize2_path: str = "martinize2"
    gmx_path: str = "gmx"
    tleap_path: str = "tleap"

    # Reproducibility / paths.
    seed: Optional[int] = None
    output_dir: str = "."

    # Smoke-test / debug toggles.
    skip_pulling: bool = False
    grompp_maxwarn: int = 5

    def __post_init__(self) -> None:
        if not self.lipids:
            raise ValueError("MembraneCGBuildConfig.lipids must not be empty")
        if len(self.lipids) > 1:
            # v1 limitation. When mixture support lands, validate ratios here.
            raise ValueError(
                "v1.19 supports a single lipid species only. "
                "Multi-lipid mixtures are deferred to a future release."
            )
        if self.peptide is None:
            raise ValueError(
                "MembraneCGBuildConfig.peptide is required (no permeant = no PMF)."
            )
        if self.insane_d_nm <= 0:
            raise ValueError(
                f"insane_d_nm must be > 0, got {self.insane_d_nm}"
            )
        if self.box_z_nm <= 0:
            raise ValueError(f"box_z_nm must be > 0, got {self.box_z_nm}")
        if self.nacl_molar < 0:
            raise ValueError(
                f"nacl_molar must be >= 0, got {self.nacl_molar}"
            )
        if self.grompp_maxwarn < 0:
            raise ValueError(
                f"grompp_maxwarn must be >= 0, got {self.grompp_maxwarn}"
            )

    def to_json(self, path: str) -> None:
        """Save configuration to JSON for reproducibility."""
        Path(path).write_text(
            json.dumps(asdict(self), indent=2, ensure_ascii=False)
        )

    @classmethod
    def from_json(cls, path: str) -> "MembraneCGBuildConfig":
        """Load configuration from JSON.

        Handles all 5 nested dataclass fields (lipids list, peptide,
        equilibration, pulling, umbrella).
        """
        raw = json.loads(Path(path).read_text())
        lipids = [LipidMix(**lp) for lp in raw.pop("lipids", [])]
        pep_raw = raw.pop("peptide", None)
        peptide = PeptideMembraneSpec(**pep_raw) if pep_raw else None
        eq = EquilibrationCGProtocol(**raw.pop("equilibration"))
        pull = PullingCGProtocol(**raw.pop("pulling"))
        umb = UmbrellaCGProtocol(**raw.pop("umbrella"))
        return cls(
            lipids=lipids,
            peptide=peptide,
            equilibration=eq,
            pulling=pull,
            umbrella=umb,
            **raw,
        )
