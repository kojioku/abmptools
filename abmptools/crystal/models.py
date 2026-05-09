# -*- coding: utf-8 -*-
"""
abmptools.crystal.models
-------------------------
Dataclasses for the crystal-FMO pipeline.

Schema follows :mod:`abmptools.genesis.mmgbsa.models` and
:mod:`abmptools.cg.peptide.models` conventions: ``@dataclass`` +
``to_json/from_json`` round-trip; YAML support is added on top via
``to_yaml/from_yaml`` because crystal configs benefit from human-friendly
nested formatting.

Seven dataclass leaves + one top-level :class:`CrystalBuildConfig`:

- :class:`CIFInputSpec`        -- one CIF input (cif path, layer, atoms_in_mol)
- :class:`CIFEngineConfig`     -- ``legacy`` vs ``ase`` backend selector
- :class:`FragmentTemplate`    -- ``input_param`` + ``segment_data.dat`` merged
- :class:`FMOMethod`           -- AJF method/basis/memory/abinit_ver/etc.
- :class:`HPCJobSpec`          -- PJM/SLURM/PBS/local jobscript spec
- :class:`PostProcessSpec`     -- IFIE/PIEDA + nearest-atom postprocessing
- :class:`CrystalBuildConfig`  -- top-level, aggregates the above

YAML schema example (also accessible via ``python -m abmptools.crystal example``)::

    project_name: csp7_layer5_around6
    output_dir: ./out
    inputs:
      - cif: XXXI-MMFF-R00001.cif
        layer: 5
        atoms_in_mol: [32]
    cif_engine:
      engine: legacy   # or 'ase'
    fragment:
      cutmode: around
      solutes: [0]
      criteria: 6.0
      molname: [UNK]
      pieda: true
      cmm: true
    fmo:
      method: MP2
      basis_set: 6-31Gdag
      memory: 6000
      is_xyz: true     # direct-coordinate AJF (full precision)
    hpc:
      scheduler: PJM
      queue: small
      group: hp190133
      nodes: 12
      proc_per_node: 2
      elapse: '24:00:00'
    postproc:
      enable: true
      frag_target: '1-10'
      annotate_nearest_atoms: true
"""
from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Union


# ---------------------------------------------------------------------------
# CIFInputSpec
# ---------------------------------------------------------------------------

@dataclass
class CIFInputSpec:
    """One CIF input.

    Attributes
    ----------
    cif
        Path to the CIF file. May be absolute or relative to the config
        file's directory (resolved by :class:`CrystalBuildConfig` when
        a config is loaded).
    layer
        Supercell expansion layer count. ``layer=N`` produces a
        :math:`(2N+1)^3` supercell. Default 5 (= 1331 cells, matches
        csp7 sample).
    atoms_in_mol
        Atoms per asymmetric-unit molecule. Single-element list ``[32]``
        for Z'=1, two-element ``[32, 32]`` for Z'=2 (heterodimer in
        asymmetric unit).
    asymmetric_only
        Skip supercell expansion -- emit only the asymmetric unit.
        Mostly useful for debug/visualisation.
    name
        Optional display name. Defaults to the CIF stem.
    """
    cif: str = ""
    layer: int = 5
    atoms_in_mol: List[int] = field(default_factory=lambda: [32])
    asymmetric_only: bool = False
    name: Optional[str] = None

    def __post_init__(self) -> None:
        if not self.cif:
            raise ValueError("CIFInputSpec.cif is required.")
        if self.layer < 1:
            raise ValueError(
                f"CIFInputSpec.layer must be >= 1, got {self.layer}"
            )
        if not self.atoms_in_mol:
            raise ValueError("CIFInputSpec.atoms_in_mol must be non-empty.")
        for n in self.atoms_in_mol:
            if n <= 0:
                raise ValueError(
                    f"CIFInputSpec.atoms_in_mol entries must be > 0, got {n}"
                )


# ---------------------------------------------------------------------------
# CIFEngineConfig
# ---------------------------------------------------------------------------

@dataclass
class CIFEngineConfig:
    """CIF backend selector.

    Phase C ships two engines:

    - ``legacy`` -- :mod:`abmptools.readcif` self-rolled parser; covers
      the symmetry operations and layer counts hard-coded for csp7
      (P21/n) workflows. Default for backward compatibility.
    - ``ase`` -- :func:`ase.io.read` based; symmetry expansion via
      ``_symmetry_equiv_pos_as_xyz`` parser, supercell via
      ``Atoms.repeat``, molecule detection via
      :func:`ase.neighborlist.neighbor_list` connected components.
    """
    engine: str = "legacy"
    detect_z_prime: bool = True
    use_neighbor_bonds: bool = True
    bond_tolerance: float = 0.4

    def __post_init__(self) -> None:
        if self.engine not in ("legacy", "ase"):
            raise ValueError(
                f"CIFEngineConfig.engine must be 'legacy' or 'ase', "
                f"got {self.engine!r}"
            )
        if self.bond_tolerance < 0:
            raise ValueError(
                f"CIFEngineConfig.bond_tolerance must be >= 0, "
                f"got {self.bond_tolerance}"
            )


# ---------------------------------------------------------------------------
# FragmentTemplate
# ---------------------------------------------------------------------------

@dataclass
class FragmentTemplate:
    """Merged ``input_param`` + ``segment_data.dat`` template.

    Mirrors the keys consumed by :class:`abmptools.setfmo` in the legacy
    pipeline (``cutmode``, ``solutes``, ``criteria``, ...). The
    ``seg_data_path`` field lets callers point at a custom
    ``segment_data.dat`` file; if ``None``, the orchestrator emits a
    minimal default segment block from the CIF input.
    """
    cutmode: str = "around"
    solutes: List[int] = field(default_factory=lambda: [0])
    criteria: Union[float, List[float]] = 6.0
    tgtpos: List[float] = field(default_factory=list)
    molname: List[str] = field(default_factory=lambda: ["UNK"])
    getmode: str = "rfile"
    pieda: bool = True
    cmm: bool = True
    solv_flag: bool = False
    seg_data_path: Optional[str] = None
    template_ajf: Optional[str] = None

    def __post_init__(self) -> None:
        valid_cutmodes = {"sphere", "cube", "around", "neutral", "none"}
        if self.cutmode not in valid_cutmodes:
            raise ValueError(
                f"FragmentTemplate.cutmode must be one of "
                f"{sorted(valid_cutmodes)}, got {self.cutmode!r}"
            )
        if not self.molname:
            raise ValueError("FragmentTemplate.molname must be non-empty.")
        if self.cutmode == "cube":
            if not isinstance(self.criteria, list) or len(self.criteria) != 3:
                raise ValueError(
                    "FragmentTemplate.criteria must be [x, y, z] for "
                    f"cutmode='cube', got {self.criteria!r}"
                )


# ---------------------------------------------------------------------------
# FMOMethod
# ---------------------------------------------------------------------------

@dataclass
class FMOMethod:
    """ABINIT-MP method / basis / output knobs.

    The ``is_xyz`` flag controls direct-coordinate AJF emission
    (``Natom={N}`` + ``&XYZ`` block with full-precision coordinates).
    Default is True because crystal pipelines benefit from avoiding
    the PDB ``%8.3f`` truncation -- this matches the csp7 sample's
    historical output.
    """
    method: str = "MP2"
    basis_set: str = "6-31Gdag"
    cpfflag: bool = True
    abinit_ver: str = "rev23"
    memory: int = 6000
    npro: int = 1
    is_xyz: bool = True
    extra_ajf_lines: List[str] = field(default_factory=list)

    def __post_init__(self) -> None:
        if self.memory <= 0:
            raise ValueError(
                f"FMOMethod.memory must be > 0, got {self.memory}"
            )
        if self.npro <= 0:
            raise ValueError(
                f"FMOMethod.npro must be > 0, got {self.npro}"
            )


# ---------------------------------------------------------------------------
# HPCJobSpec
# ---------------------------------------------------------------------------

@dataclass
class HPCJobSpec:
    """HPC scheduler abstraction (PJM / SLURM / PBS / local).

    For ``scheduler='local'`` the orchestrator emits a plain ``bash``
    script that invokes ``abinitmp`` directly (used by ``--run-local``
    in the smoke tutorial). For ``PJM`` / ``SLURM`` / ``PBS`` the
    rendered script targets the corresponding queue manager. A user
    can supply ``template_override`` (path to a ``string.Template``
    text file) to bypass the bundled templates entirely.
    """
    scheduler: str = "PJM"
    queue: str = "small"
    group: Optional[str] = None
    nodes: int = 12
    proc_per_node: int = 2
    omp_threads: int = 12
    elapse: str = "24:00:00"
    abinit_dir: str = ""
    binary_name: str = "abinitmp_smp"
    mkinp_path: str = ""
    extra_lines: List[str] = field(default_factory=list)
    template_override: Optional[str] = None

    def __post_init__(self) -> None:
        valid_schedulers = {"PJM", "SLURM", "PBS", "local"}
        if self.scheduler not in valid_schedulers:
            raise ValueError(
                f"HPCJobSpec.scheduler must be one of "
                f"{sorted(valid_schedulers)}, got {self.scheduler!r}"
            )
        if self.nodes <= 0:
            raise ValueError(f"HPCJobSpec.nodes must be > 0, got {self.nodes}")
        if self.proc_per_node <= 0:
            raise ValueError(
                f"HPCJobSpec.proc_per_node must be > 0, "
                f"got {self.proc_per_node}"
            )


# ---------------------------------------------------------------------------
# PostProcessSpec
# ---------------------------------------------------------------------------

@dataclass
class PostProcessSpec:
    """IFIE/PIEDA + nearest-atom postprocessing.

    Wraps :func:`abmptools.getifiepieda` (CSV emission) and
    :func:`abmptools.crystal.atom_distance.find_nearest_atoms`
    (per-fragment nearest-atom annotation). Disabled by default because
    postproc requires log files that only exist after ABINIT-MP runs.
    """
    enable: bool = False
    frag_target: Optional[str] = None
    distance: Optional[float] = None
    annotate_nearest_atoms: bool = True
    nearest_atom_count: int = 3
    output_dir: str = "postproc"

    def __post_init__(self) -> None:
        if self.nearest_atom_count <= 0:
            raise ValueError(
                f"PostProcessSpec.nearest_atom_count must be > 0, "
                f"got {self.nearest_atom_count}"
            )


# ---------------------------------------------------------------------------
# CrystalBuildConfig (top-level)
# ---------------------------------------------------------------------------

@dataclass
class CrystalBuildConfig:
    """Top-level crystal-FMO build configuration.

    Attributes are nested dataclasses so JSON / YAML round-trip is
    well-defined. Use :meth:`from_yaml` / :meth:`from_json` to load
    persisted configs and :meth:`to_yaml` / :meth:`to_json` to emit
    them. All paths in :class:`CIFInputSpec.cif` and
    :attr:`output_dir` are interpreted relative to the **config file's
    directory** when loaded via :meth:`from_yaml` / :meth:`from_json`,
    not the current working directory.
    """
    inputs: List[CIFInputSpec] = field(default_factory=list)
    cif_engine: CIFEngineConfig = field(default_factory=CIFEngineConfig)
    fragment: FragmentTemplate = field(default_factory=FragmentTemplate)
    fmo: FMOMethod = field(default_factory=FMOMethod)
    hpc: HPCJobSpec = field(default_factory=HPCJobSpec)
    postproc: PostProcessSpec = field(default_factory=PostProcessSpec)
    output_dir: str = "./crystal_out"
    project_name: str = "crystal_run"
    fail_fast: bool = True
    run_local: bool = False

    def __post_init__(self) -> None:
        if not self.inputs:
            raise ValueError(
                "CrystalBuildConfig.inputs must list at least one CIF input."
            )
        if not self.project_name:
            raise ValueError("CrystalBuildConfig.project_name must be non-empty.")

    # ------------------------------------------------------------------
    # JSON round-trip
    # ------------------------------------------------------------------

    def to_json(self, path: str) -> None:
        """Save configuration to JSON for reproducibility."""
        Path(path).write_text(
            json.dumps(asdict(self), indent=2, ensure_ascii=False)
        )

    @classmethod
    def from_json(cls, path: str) -> "CrystalBuildConfig":
        """Load configuration from JSON.

        Round-trips all 6 nested dataclass fields (``inputs`` list +
        ``cif_engine`` + ``fragment`` + ``fmo`` + ``hpc`` + ``postproc``).
        """
        raw = json.loads(Path(path).read_text())
        return cls._from_dict(raw)

    # ------------------------------------------------------------------
    # YAML round-trip
    # ------------------------------------------------------------------

    def to_yaml(self, path: str) -> None:
        """Save configuration to YAML.

        Requires :mod:`yaml`; raises :class:`ImportError` if pyyaml is
        not installed.
        """
        try:
            import yaml
        except ImportError as exc:
            raise ImportError(
                "pyyaml is required for YAML configs; install via "
                "`mamba install pyyaml` or `pip install abmptools[crystal]`."
            ) from exc
        Path(path).write_text(
            yaml.safe_dump(asdict(self), sort_keys=False, allow_unicode=True)
        )

    @classmethod
    def from_yaml(cls, path: str) -> "CrystalBuildConfig":
        """Load configuration from YAML."""
        try:
            import yaml
        except ImportError as exc:
            raise ImportError(
                "pyyaml is required for YAML configs; install via "
                "`mamba install pyyaml` or `pip install abmptools[crystal]`."
            ) from exc
        raw = yaml.safe_load(Path(path).read_text())
        if raw is None:
            raise ValueError(f"YAML config is empty: {path}")
        return cls._from_dict(raw)

    @classmethod
    def _from_dict(cls, raw: Dict[str, Any]) -> "CrystalBuildConfig":
        raw = dict(raw)
        inputs_raw = raw.pop("inputs", [])
        inputs = [CIFInputSpec(**i) for i in inputs_raw]
        cif_engine = CIFEngineConfig(**raw.pop("cif_engine", {}))
        fragment = FragmentTemplate(**raw.pop("fragment", {}))
        fmo = FMOMethod(**raw.pop("fmo", {}))
        hpc = HPCJobSpec(**raw.pop("hpc", {}))
        postproc = PostProcessSpec(**raw.pop("postproc", {}))
        return cls(
            inputs=inputs,
            cif_engine=cif_engine,
            fragment=fragment,
            fmo=fmo,
            hpc=hpc,
            postproc=postproc,
            **raw,
        )


__all__ = [
    "CIFInputSpec",
    "CIFEngineConfig",
    "FragmentTemplate",
    "FMOMethod",
    "HPCJobSpec",
    "PostProcessSpec",
    "CrystalBuildConfig",
]
