# -*- coding: utf-8 -*-
"""
abmptools.amorphous.models
---------------------------
Dataclasses for amorphous structure building.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import List, Optional


@dataclass
class ComponentSpec:
    """Specification for a single molecular component.

    Provide exactly one of *smiles* or *sdf_path*.
    """
    name: str = ""
    smiles: str = ""
    sdf_path: str = ""
    n_mol: int = 0
    weight_fraction: float = 0.0
    molecular_weight: float = 0.0  # filled at runtime

    def __post_init__(self) -> None:
        if not self.smiles and not self.sdf_path:
            raise ValueError("ComponentSpec requires either 'smiles' or 'sdf_path'.")


@dataclass
class BuildConfig:
    """Full configuration for an amorphous build."""
    components: List[ComponentSpec] = field(default_factory=list)
    total_molecules: int = 0
    box_size_nm: float = 0.0
    density_g_cm3: float = 0.8
    temperature: float = 300.0
    T_high: float = 600.0
    pressure: float = 1.0
    forcefield: str = "openff_unconstrained-2.1.0.offxml"
    packmol_tolerance: float = 2.0
    packmol_path: str = "packmol"
    seed: Optional[int] = None
    output_dir: str = "."

    # --- MDP overrides ---
    em_steps: int = 50000
    em_tol: float = 1000.0
    nvt_high_nsteps: int = 100000
    npt_high_nsteps: int = 200000
    anneal_nsteps: int = 500000
    npt_low_nsteps: int = 500000
    dt: float = 0.001
    tau_t: float = 0.1
    tau_p: float = 2.0
    nstxout_compressed: int = 5000
    nstenergy: int = 1000

    def to_json(self, path: str) -> None:
        """Save configuration to JSON for reproducibility."""
        data = asdict(self)
        Path(path).write_text(json.dumps(data, indent=2, ensure_ascii=False))

    @classmethod
    def from_json(cls, path: str) -> "BuildConfig":
        """Load configuration from JSON."""
        raw = json.loads(Path(path).read_text())
        comps = [ComponentSpec(**c) for c in raw.pop("components", [])]
        return cls(components=comps, **raw)
