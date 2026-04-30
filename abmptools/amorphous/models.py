# -*- coding: utf-8 -*-
"""
abmptools.amorphous.models
---------------------------
Dataclasses for amorphous structure building.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field, asdict
from typing import Any, List, Optional
from pathlib import Path


@dataclass
class ComponentSpec:
    """Specification for a single molecular component.

    Provide exactly one of *smiles*, *sdf_path*, or *pdb_path*.

    *pdb_path* (Phase 9-a): for pre-built oligomer / polymer structures
    that the caller already constructed (e.g. fcewsmb's
    ``getpolymer()`` writes ``polymer/pdb/<name>.pdb`` for ``repeat > 1``
    segments). The PDB is fed to OpenFF Topology directly via
    ``Topology.from_pdb``; the OpenFF Toolkit must be able to assign
    parameters to the resulting molecule (works for simple sp3/sp2
    organics, may fail for complex polymers — fall back to the legacy
    UDF route in that case).
    """
    name: str = ""
    smiles: str = ""
    sdf_path: str = ""
    pdb_path: str = ""
    n_mol: int = 0
    weight_fraction: float = 0.0
    molecular_weight: float = 0.0  # filled at runtime

    def __post_init__(self) -> None:
        sources = [self.smiles, self.sdf_path, self.pdb_path]
        non_empty = [s for s in sources if s]
        if not non_empty:
            raise ValueError(
                "ComponentSpec requires one of 'smiles', 'sdf_path', or 'pdb_path'."
            )
        if len(non_empty) > 1:
            raise ValueError(
                "ComponentSpec accepts exactly one of 'smiles', 'sdf_path', "
                "or 'pdb_path'; got: "
                + ", ".join(
                    name for name, val in
                    [('smiles', self.smiles), ('sdf_path', self.sdf_path),
                     ('pdb_path', self.pdb_path)]
                    if val
                )
            )


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
    # OpenFF force field OFFXML name(s).
    # str: 単一 FF (default、organic 全般を openff_unconstrained-2.1.0 で扱う)
    # List[str]: stacked FFs。後ろが前を override する SMIRKS-overlay 機構で、
    #   typical 用途は water に専用モデルを当てること:
    #   ['openff_unconstrained-2.1.0.offxml', 'tip3p.offxml']
    #   GAFF water の repulsive σ/ε で純 water 系が膨張する問題
    #   (1.0→0.26 g/cm³) の解決策。OpenFF 標準の 'tip3p.offxml' /
    #   'tip3p_fb.offxml' / 'spce.offxml' 等を最後に追加すれば water
    #   分子に SMIRKS マッチして上書きされる。
    forcefield: Any = "openff_unconstrained-2.1.0.offxml"
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
