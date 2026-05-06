# -*- coding: utf-8 -*-
"""
abmptools.genesis.mmgbsa.models
-------------------------------
Dataclasses for GENESIS MM/GBSA single-point ΔG_bind calculations.

Schema follows :mod:`abmptools.genesis.grest.models` and
:mod:`abmptools.cg.membrane.models` conventions: ``@dataclass`` +
``to_json/from_json`` round-trip, no Pydantic.

Six dataclasses:

- :class:`TargetSpec`             -- one protein-ligand binding target
- :class:`ForceFieldSet`          -- tleap ``source`` lines + addPdbResMap
- :class:`LigandParameterization` -- acpype invocation knobs
- :class:`EnergyProtocol`         -- atdyn [ENERGY] + GBSA params
- :class:`MinimizationProtocol`   -- atdyn [MINIMIZE] block (POC: 1-step)
- :class:`MMGBSABuildConfig`      -- top-level config, aggregates above

Cross-field invariant: at least one of ``targets`` (config mode) or
``input_dir`` (folder mode) must be non-empty. Enforced in
``__post_init__``.
"""
from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import List, Optional, Tuple


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

#: POC default ``addPdbResMap`` (counterion canonicalisation).
DEFAULT_ADD_PDB_RES_MAP: List[Tuple[str, str]] = [
    ("LI+", "LI"), ("NA+", "NA"), ("MG2", "MG"), ("CL-", "CL"),
    ("K+", "K"), ("RB+", "RB"), ("CS+", "CS"),
]


# ---------------------------------------------------------------------------
# TargetSpec
# ---------------------------------------------------------------------------

@dataclass
class TargetSpec:
    """One protein-ligand binding target.

    Attributes
    ----------
    pdb
        Path to the PDB file. May be absolute or basename relative to
        :attr:`MMGBSABuildConfig.input_dir`. The basename without ``.pdb``
        is used as the per-target output directory, matching the POC
        convention ``<pdb_basename>/<pdb_basename>_{ligand,receptor}_<resno>.pdb``.
    ligand_resno
        Heteroatom residue id to extract as the ligand. Must be > 0.
    chain
        Optional 1-character chain ID filter; ``None`` accepts any chain.
    name
        Display name. Default is the PDB stem (filled at orchestrator
        level if unset).
    """
    pdb: str = ""
    ligand_resno: int = 0
    chain: Optional[str] = None
    name: Optional[str] = None

    def __post_init__(self) -> None:
        if not self.pdb:
            raise ValueError("TargetSpec.pdb is required.")
        if self.ligand_resno <= 0:
            raise ValueError(
                f"TargetSpec.ligand_resno must be > 0, got {self.ligand_resno}"
            )
        if self.chain is not None and len(self.chain) != 1:
            raise ValueError(
                f"TargetSpec.chain must be a 1-character ID, got {self.chain!r}"
            )


# ---------------------------------------------------------------------------
# ForceFieldSet
# ---------------------------------------------------------------------------

@dataclass
class ForceFieldSet:
    """``tleap`` ``source`` lines + ``addPdbResMap``/``addPdbAtomMap`` overrides.

    Defaults reproduce the POC's 6 ``leaprc`` lines (ff14SB + DNA.OL15 +
    RNA.OL3 + TIP3P + GAFF2 + GAFF). Override per-build for ff14iSB /
    DNA.OL21 / etc.
    """
    leaprc_protein: str = "leaprc.protein.ff14SB"
    leaprc_dna: str = "leaprc.DNA.OL15"
    leaprc_rna: str = "leaprc.RNA.OL3"
    leaprc_water: str = "leaprc.water.tip3p"
    leaprc_gaff2: str = "leaprc.gaff2"
    leaprc_gaff: str = "leaprc.gaff"
    extra_lines: List[str] = field(default_factory=list)
    add_pdb_res_map: List[Tuple[str, str]] = field(
        default_factory=lambda: list(DEFAULT_ADD_PDB_RES_MAP)
    )
    add_pdb_atom_map: List[Tuple[str, str]] = field(
        default_factory=lambda: list(DEFAULT_ADD_PDB_RES_MAP)
    )

    def __post_init__(self) -> None:
        # Normalise tuples (JSON round-trip yields lists).
        self.add_pdb_res_map = [tuple(p) for p in self.add_pdb_res_map]
        self.add_pdb_atom_map = [tuple(p) for p in self.add_pdb_atom_map]
        for label, pairs in (
            ("add_pdb_res_map", self.add_pdb_res_map),
            ("add_pdb_atom_map", self.add_pdb_atom_map),
        ):
            for pair in pairs:
                if len(pair) != 2:
                    raise ValueError(
                        f"ForceFieldSet.{label} entries must be 2-tuples, "
                        f"got {pair!r}"
                    )


# ---------------------------------------------------------------------------
# LigandParameterization
# ---------------------------------------------------------------------------

@dataclass
class LigandParameterization:
    """Knobs for the ``acpype`` invocation (ligand GAFF/GAFF2 + AM1-BCC).

    Default: ``-c bcc -k maxcyc=0`` (POC). ``maxcyc=0`` skips antechamber's
    pre-optimisation since the input PDB usually has reasonable geometry
    already.
    """
    charge_method: str = "bcc"          # "bcc" | "gas" | "user"
    net_charge: Optional[int] = None    # None: let acpype detect
    extra_keys: str = "maxcyc=0"
    atom_type: str = "gaff2"            # "gaff" | "gaff2"
    skip_if_cached: bool = True

    def __post_init__(self) -> None:
        if self.charge_method not in ("bcc", "gas", "user"):
            raise ValueError(
                f"LigandParameterization.charge_method must be "
                f"bcc/gas/user, got {self.charge_method!r}"
            )
        if self.atom_type not in ("gaff", "gaff2"):
            raise ValueError(
                f"LigandParameterization.atom_type must be gaff or gaff2, "
                f"got {self.atom_type!r}"
            )


# ---------------------------------------------------------------------------
# EnergyProtocol
# ---------------------------------------------------------------------------

@dataclass
class EnergyProtocol:
    """Single-point GBSA ``[ENERGY]`` block parameters (POC values).

    GBSA implicit solvent is the load-bearing setting:

    - ``gbsa_salt_cons``: Debye screening (M); 0.15 ≈ physiological
    - ``gbsa_surf_tens``: surface area term coefficient (kcal/(mol·Å²))
    - ``gbsa_vdw_offset``: dielectric boundary offset (Å)

    NOBC requires very large ``cutoffdist_A`` (POC: 99.9 Å).
    """
    electrostatic: str = "CUTOFF"
    switchdist_A: float = 99.9
    cutoffdist_A: float = 99.9
    pairlistdist_A: float = 100.0
    implicit_solvent: str = "GBSA"
    gbsa_salt_cons: float = 0.15
    gbsa_surf_tens: float = 0.0072
    gbsa_vdw_offset: float = 0.09

    def __post_init__(self) -> None:
        if self.implicit_solvent not in ("GBSA", "NONE"):
            raise ValueError(
                f"EnergyProtocol.implicit_solvent must be GBSA or NONE, "
                f"got {self.implicit_solvent!r}"
            )
        if self.electrostatic not in ("CUTOFF", "PME"):
            raise ValueError(
                f"EnergyProtocol.electrostatic must be CUTOFF or PME, "
                f"got {self.electrostatic!r}"
            )
        if self.cutoffdist_A <= 0:
            raise ValueError(
                f"EnergyProtocol.cutoffdist_A must be > 0, got {self.cutoffdist_A}"
            )
        if self.pairlistdist_A < self.cutoffdist_A:
            raise ValueError(
                f"EnergyProtocol.pairlistdist_A ({self.pairlistdist_A}) must "
                f"be >= cutoffdist_A ({self.cutoffdist_A})"
            )


# ---------------------------------------------------------------------------
# MinimizationProtocol
# ---------------------------------------------------------------------------

@dataclass
class MinimizationProtocol:
    """``[MINIMIZE]`` block. POC: 1-step single-point evaluation.

    For ensemble averaging (deferred to v1.22.x) you would replace this
    with multi-frame trajectory analysis. v1.22.0 ships single-point only.
    """
    method: str = "SD"          # "SD" | "LBFGS"
    nsteps: int = 1
    eneout_period: int = 1
    crdout_period: int = 1
    rstout_period: int = 1
    nbupdate_period: int = 1

    def __post_init__(self) -> None:
        if self.method not in ("SD", "LBFGS"):
            raise ValueError(
                f"MinimizationProtocol.method must be SD or LBFGS, "
                f"got {self.method!r}"
            )
        for label, ns in (
            ("nsteps", self.nsteps),
            ("eneout_period", self.eneout_period),
            ("crdout_period", self.crdout_period),
            ("rstout_period", self.rstout_period),
            ("nbupdate_period", self.nbupdate_period),
        ):
            if ns < 1:
                raise ValueError(
                    f"MinimizationProtocol.{label} must be >= 1, got {ns}"
                )


# ---------------------------------------------------------------------------
# MMGBSABuildConfig (top-level)
# ---------------------------------------------------------------------------

@dataclass
class MMGBSABuildConfig:
    """End-to-end MM/GBSA configuration.

    Mode "C" — both inputs are supported:

    - **JSON config (preferred)**: explicit ``targets: [TargetSpec, ...]``.
      Re-runnable, version-controlled.
    - **Folder mode (POC compatibility)**: ``input_dir`` set, ``targets``
      empty. ``MMGBSAOrchestrator`` enumerates ``*.pdb`` in ``input_dir``
      and synthesises ``targets`` lazily before pipeline start. The
      persisted JSON dump records the resolved ``targets`` list.
    """
    targets: List[TargetSpec] = field(default_factory=list)
    input_dir: str = ""
    project_name: str = "mmgbsa_run"

    # Force field + ligand parameterization
    force_field: ForceFieldSet = field(default_factory=ForceFieldSet)
    ligand: LigandParameterization = field(default_factory=LigandParameterization)

    # Protocols
    energy: EnergyProtocol = field(default_factory=EnergyProtocol)
    minimize: MinimizationProtocol = field(default_factory=MinimizationProtocol)

    # MPI / execution (POC: -np 1 — no spatial decomposition needed for
    # single-point GBSA on small systems)
    mpi_processes: int = 1
    mpirun_path: str = "mpirun"

    # External tool paths
    atdyn_path: str = "atdyn"
    tleap_path: str = "tleap"
    acpype_path: str = "acpype"

    # Output / reproducibility
    output_dir: str = "."
    seed: Optional[int] = None
    fail_fast: bool = True

    def __post_init__(self) -> None:
        if not self.targets and not self.input_dir:
            raise ValueError(
                "MMGBSABuildConfig requires either non-empty 'targets' or "
                "non-empty 'input_dir' (folder-mode shortcut)."
            )
        if self.mpi_processes < 1:
            raise ValueError(
                f"MMGBSABuildConfig.mpi_processes must be >= 1, "
                f"got {self.mpi_processes}"
            )
        if not self.project_name:
            raise ValueError("MMGBSABuildConfig.project_name must be non-empty.")

    def to_json(self, path: str) -> None:
        """Save configuration to JSON for reproducibility."""
        Path(path).write_text(
            json.dumps(asdict(self), indent=2, ensure_ascii=False)
        )

    @classmethod
    def from_json(cls, path: str) -> "MMGBSABuildConfig":
        """Load configuration from JSON.

        Round-trips all 5 nested dataclass fields (``targets`` list +
        ``force_field`` + ``ligand`` + ``energy`` + ``minimize``).
        """
        raw = json.loads(Path(path).read_text())
        targets = [TargetSpec(**t) for t in raw.pop("targets", [])]
        force_field = ForceFieldSet(**raw.pop("force_field"))
        ligand = LigandParameterization(**raw.pop("ligand"))
        energy = EnergyProtocol(**raw.pop("energy"))
        minimize = MinimizationProtocol(**raw.pop("minimize"))
        return cls(
            targets=targets,
            force_field=force_field,
            ligand=ligand,
            energy=energy,
            minimize=minimize,
            **raw,
        )


__all__ = [
    "DEFAULT_ADD_PDB_RES_MAP",
    "EnergyProtocol",
    "ForceFieldSet",
    "LigandParameterization",
    "MMGBSABuildConfig",
    "MinimizationProtocol",
    "TargetSpec",
]
