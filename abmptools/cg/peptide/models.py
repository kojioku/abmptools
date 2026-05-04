# -*- coding: utf-8 -*-
"""
abmptools.cg.peptide.models
----------------------------
Dataclasses for Martini 3 peptide CG system building.
"""
from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import List, Optional


STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")


@dataclass
class PeptideSpec:
    """Specification for a single peptide species in the system.

    Sequence is given in 1-letter standard amino acid codes (case-insensitive,
    normalized to uppercase). Non-standard residues raise ``ValueError``.
    """
    name: str = ""
    sequence: str = ""
    count: int = 1

    def __post_init__(self) -> None:
        if not self.name:
            raise ValueError("PeptideSpec.name must be non-empty")
        if not self.sequence:
            raise ValueError("PeptideSpec.sequence must be non-empty")
        seq = self.sequence.upper()
        invalid = set(seq) - STANDARD_AA
        if invalid:
            raise ValueError(
                f"Invalid amino acid(s) in sequence {seq!r}: {sorted(invalid)}. "
                f"Supported: {''.join(sorted(STANDARD_AA))}"
            )
        self.sequence = seq
        if self.count < 1:
            raise ValueError(
                f"PeptideSpec.count must be >= 1, got {self.count}"
            )


@dataclass
class PeptideBuildConfig:
    """Full configuration for a Martini 3 peptide CG build."""
    peptides: List[PeptideSpec] = field(default_factory=list)

    # --- Box ---
    box_size_nm: float = 0.0
    box_lengths_nm: Optional[List[float]] = None  # rectangular [x, y, z]

    # --- Solvent / ions ---
    solvent_enabled: bool = True
    neutralize: bool = True
    nacl_molar: float = 0.15

    # --- MD parameters ---
    temperature: float = 310.0
    seed: Optional[int] = None

    # --- Force field files (user-supplied; not bundled, see forcefield_check) ---
    # Directory containing martini_v3.0.0*.itp + martini_v3.0.0_water.gro.
    # Empty string -> resolved to output_dir/ff at validate/build time.
    martini_itp_dir: str = ""

    # --- External tool paths ---
    martinize2_path: str = "martinize2"
    gmx_path: str = "gmx"
    tleap_path: str = "tleap"

    # --- Output ---
    output_dir: str = "."

    # --- MDP toggles + nsteps ---
    mdp_em: bool = True
    mdp_nvt: bool = True
    mdp_npt: bool = True
    mdp_md: bool = True
    em_steps: int = 50000
    nvt_nsteps: int = 250000
    npt_nsteps: int = 250000
    md_nsteps: int = 5000000
    dt_fs: float = 20.0  # Martini 3 では typical 20 fs

    def __post_init__(self) -> None:
        if self.box_lengths_nm is not None:
            if len(self.box_lengths_nm) != 3:
                raise ValueError(
                    "box_lengths_nm must have 3 values [x, y, z], "
                    f"got {len(self.box_lengths_nm)}"
                )
            if any(L <= 0 for L in self.box_lengths_nm):
                raise ValueError(
                    f"box_lengths_nm values must be positive: {self.box_lengths_nm}"
                )
        if self.box_size_nm < 0:
            raise ValueError(
                f"box_size_nm must be >= 0, got {self.box_size_nm}"
            )
        if self.peptides:
            names = [p.name for p in self.peptides]
            if len(names) != len(set(names)):
                raise ValueError(f"Peptide names must be unique: {names}")
        if self.nacl_molar < 0:
            raise ValueError(
                f"nacl_molar must be >= 0, got {self.nacl_molar}"
            )
        if self.dt_fs <= 0:
            raise ValueError(f"dt_fs must be > 0, got {self.dt_fs}")
        for label, ns in (
            ("em_steps", self.em_steps),
            ("nvt_nsteps", self.nvt_nsteps),
            ("npt_nsteps", self.npt_nsteps),
            ("md_nsteps", self.md_nsteps),
        ):
            if ns < 0:
                raise ValueError(f"{label} must be >= 0, got {ns}")

    def to_json(self, path: str) -> None:
        """Save configuration to JSON for reproducibility."""
        data = asdict(self)
        Path(path).write_text(json.dumps(data, indent=2, ensure_ascii=False))

    @classmethod
    def from_json(cls, path: str) -> "PeptideBuildConfig":
        """Load configuration from JSON."""
        raw = json.loads(Path(path).read_text())
        peps = [PeptideSpec(**p) for p in raw.pop("peptides", [])]
        return cls(peptides=peps, **raw)
