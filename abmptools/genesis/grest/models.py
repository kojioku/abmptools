# -*- coding: utf-8 -*-
"""
abmptools.genesis.grest.models
------------------------------
Dataclasses for GENESIS gREST_SSCR system building.

Schema follows :mod:`abmptools.cg.membrane.models` and
:mod:`abmptools.cg.peptide.models` conventions: ``@dataclass`` with
``to_json/from_json`` round-trip, no Pydantic.

The six dataclasses are:

- :class:`RESTSelectionSpec` -- which residues become REST solute
  (``explicit`` list or ``around`` mask).
- :class:`ReplicaTemperatureSpec` -- temperature ladder
  (``auto`` geometric or ``manual`` list).
- :class:`MinimizationStage` -- atdyn minimisation control parameters.
- :class:`EquilibrationStage` -- spdyn equilibration control parameters.
- :class:`GrestStage` -- spdyn gREST production control parameters.
- :class:`GrestBuildConfig` -- top-level configuration aggregating the
  above plus AMBER force-field selection, MPI mapping, and external
  tool paths.
"""
from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import List, Optional


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

#: GENESIS REST ``param_type`` tokens (Setup_Remd_Solute_Tempering).
#: ALL = whole molecule; B = bonds; A = angles; U = Urey-Bradley;
#: D = dihedrals; RB = Ryckaert-Bellemans; I = impropers; CM = CMAP;
#: CT = contact; C = charge; L = Lennard-Jones.
VALID_PARAM_TYPES = frozenset(
    {"ALL", "B", "A", "U", "D", "RB", "I", "CM", "CT", "C", "L"}
)


# ---------------------------------------------------------------------------
# Selection & temperatures
# ---------------------------------------------------------------------------

@dataclass
class RESTSelectionSpec:
    """Residue selection for the REST solute region.

    Two modes:

    - ``mode="explicit"``: comma-separated residue numbers or
      hyphenated ranges, e.g. ``["1-138"]`` or ``["21", "96", "274-275"]``.
      Stored as ``residues: List[str]`` (string-typed entries; integer
      lists are converted to strings on validation so the JSON form is
      type-stable).
    - ``mode="around"``: cpptraj-style mask on a centre residue plus
      radius in Angstrom, e.g. ``center="rno:96"`` and ``radius_A=5.0``.

    ``param_types`` controls the GENESIS ``Setup_Remd_Solute_Tempering``
    flags (which interactions get tempered). Default is ``["C", "L"]``
    (charge + Lennard-Jones), matching the POC output
    (``CHARGE=T, LJ=T`` in ``Setup_Remd_Solute_Tempering`` block) and
    the SSCR-typical setting in tutorial 12.3 step3.
    """
    mode: str = "explicit"
    residues: List[str] = field(default_factory=list)
    center: str = ""
    radius_A: float = 5.0
    param_types: List[str] = field(default_factory=lambda: ["C", "L"])
    include_heavy_only: bool = False

    def __post_init__(self) -> None:
        if self.mode not in ("explicit", "around"):
            raise ValueError(
                f"RESTSelectionSpec.mode must be 'explicit' or 'around', "
                f"got {self.mode!r}"
            )
        # Normalise residues to List[str].
        self.residues = [str(r) for r in self.residues]
        if self.mode == "explicit":
            if not self.residues:
                raise ValueError(
                    "RESTSelectionSpec(mode='explicit') requires "
                    "non-empty residues."
                )
        if self.mode == "around":
            if not self.center:
                raise ValueError(
                    "RESTSelectionSpec(mode='around') requires "
                    "non-empty center (e.g. 'rno:96')."
                )
            if self.radius_A <= 0:
                raise ValueError(
                    f"RESTSelectionSpec.radius_A must be > 0, got {self.radius_A}"
                )
        # Validate param_types tokens.
        bad = [p for p in self.param_types if p not in VALID_PARAM_TYPES]
        if bad:
            raise ValueError(
                f"Unknown param_type tokens: {bad}. "
                f"Allowed: {sorted(VALID_PARAM_TYPES)}"
            )


@dataclass
class ReplicaTemperatureSpec:
    """Replica temperature ladder.

    Two modes:

    - ``mode="auto"``:
        ``method="geometric"``: ``T_i = T_min * (T_max/T_min)**(i/(N-1))``
        ``method="vanderspoel"``: Patriksson-van der Spoel approximate
        acceptance-ratio formula. **Deferred to v1.20.x**; raises
        ``NotImplementedError`` in v1.20.0.
    - ``mode="manual"``: a fixed list, monotonically increasing,
      e.g. ``[300.0, 318.11, 337.11, 357.10]`` (POC values).
    """
    mode: str = "auto"
    method: str = "geometric"
    T_min: float = 300.0
    T_max: float = 357.10
    n_replicas: int = 4
    temperatures: List[float] = field(default_factory=list)

    def __post_init__(self) -> None:
        if self.mode not in ("auto", "manual"):
            raise ValueError(
                f"ReplicaTemperatureSpec.mode must be 'auto' or 'manual', "
                f"got {self.mode!r}"
            )
        if self.mode == "auto":
            if self.method not in ("geometric", "vanderspoel"):
                raise ValueError(
                    f"ReplicaTemperatureSpec.method must be 'geometric' "
                    f"or 'vanderspoel', got {self.method!r}"
                )
            if self.n_replicas < 2:
                raise ValueError(
                    f"ReplicaTemperatureSpec.n_replicas must be >= 2, "
                    f"got {self.n_replicas}"
                )
            if self.T_min <= 0:
                raise ValueError(
                    f"ReplicaTemperatureSpec.T_min must be > 0, got {self.T_min}"
                )
            if self.T_max <= self.T_min:
                raise ValueError(
                    f"ReplicaTemperatureSpec.T_max ({self.T_max}) must be > "
                    f"T_min ({self.T_min})"
                )
        else:  # manual
            if len(self.temperatures) < 2:
                raise ValueError(
                    "ReplicaTemperatureSpec(mode='manual') requires "
                    ">= 2 temperatures."
                )
            for t in self.temperatures:
                if t <= 0:
                    raise ValueError(
                        f"All temperatures must be > 0, got {t}"
                    )
            if list(self.temperatures) != sorted(self.temperatures):
                raise ValueError(
                    "ReplicaTemperatureSpec.temperatures must be "
                    "monotonically increasing."
                )


# ---------------------------------------------------------------------------
# Stage protocols
# ---------------------------------------------------------------------------

@dataclass
class MinimizationStage:
    """atdyn step -- energy minimisation under PBC + position restraint.

    Defaults follow GENESIS tutorial 12.3 step1.
    """
    method: str = "SD"
    nsteps: int = 10_000
    eneout_period: int = 250
    crdout_period: int = 250
    rstout_period: int = 10_000
    posres_force_kcal_per_mol: float = 1.0
    cutoffdist_A: float = 12.0
    pairlistdist_A: float = 13.5

    def __post_init__(self) -> None:
        if self.method not in ("SD",):
            raise ValueError(
                f"MinimizationStage.method must be 'SD' (only option in "
                f"GENESIS atdyn), got {self.method!r}"
            )
        for label, ns in (
            ("nsteps", self.nsteps),
            ("eneout_period", self.eneout_period),
            ("crdout_period", self.crdout_period),
            ("rstout_period", self.rstout_period),
        ):
            if ns < 0:
                raise ValueError(f"MinimizationStage.{label} must be >= 0, got {ns}")
        if self.posres_force_kcal_per_mol < 0:
            raise ValueError(
                f"MinimizationStage.posres_force_kcal_per_mol must be >= 0, "
                f"got {self.posres_force_kcal_per_mol}"
            )
        if self.pairlistdist_A < self.cutoffdist_A:
            raise ValueError(
                f"MinimizationStage.pairlistdist_A ({self.pairlistdist_A}) "
                f"must be >= cutoffdist_A ({self.cutoffdist_A})"
            )


@dataclass
class EquilibrationStage:
    """spdyn step -- short MD to relax box + heat to target T.

    Defaults follow GENESIS tutorial 12.3 step2 (single NPT phase).
    Hydrogen-mass-repartition (HMR) is enabled by default to permit
    dt = 2 fs without bond-constraint deadlock.
    """
    integrator: str = "VVER"
    timestep_ps: float = 0.002
    nsteps: int = 25_000
    eneout_period: int = 500
    crdout_period: int = 500
    rstout_period: int = 25_000
    ensemble: str = "NPT"
    tpcontrol: str = "BUSSI"
    temperature_K: float = 300.0
    pressure_bar: float = 1.0
    hydrogen_mr: bool = True
    hmr_ratio: float = 3.0
    hmr_ratio_xh1: float = 2.0
    nbupdate_period: int = 10
    cutoffdist_A: float = 8.0
    pairlistdist_A: float = 9.5

    def __post_init__(self) -> None:
        if self.integrator not in ("LEAP", "VVER", "VRES"):
            raise ValueError(
                f"EquilibrationStage.integrator must be one of "
                f"LEAP/VVER/VRES, got {self.integrator!r}"
            )
        if self.timestep_ps <= 0:
            raise ValueError(
                f"EquilibrationStage.timestep_ps must be > 0, got {self.timestep_ps}"
            )
        if self.ensemble not in ("NVE", "NVT", "NPT"):
            raise ValueError(
                f"EquilibrationStage.ensemble must be one of NVE/NVT/NPT, "
                f"got {self.ensemble!r}"
            )
        if self.temperature_K <= 0:
            raise ValueError(
                f"EquilibrationStage.temperature_K must be > 0, got {self.temperature_K}"
            )
        for label, ns in (
            ("nsteps", self.nsteps),
            ("eneout_period", self.eneout_period),
            ("crdout_period", self.crdout_period),
            ("rstout_period", self.rstout_period),
        ):
            if ns < 0:
                raise ValueError(f"EquilibrationStage.{label} must be >= 0, got {ns}")
        if self.pairlistdist_A < self.cutoffdist_A:
            raise ValueError(
                f"EquilibrationStage.pairlistdist_A ({self.pairlistdist_A}) "
                f"must be >= cutoffdist_A ({self.cutoffdist_A})"
            )


@dataclass
class GrestStage:
    """spdyn step -- production gREST_SSCR.

    Defaults follow GENESIS tutorial 12.3 step3:
    VRES r-RESPA / dt = 3.5 fs / NVT / Bussi thermostat /
    nsteps = 3 M (= 10.5 ns) / exchange every 3000 steps (10.5 ps).
    POC used 60 M nsteps; users override per-config.

    Constraints (verified against POC stderr):

    - ``mod(nsteps, 2 * exchange_period * dimension) == 0`` is a
      ``Setup_Remd`` requirement. ``dimension == 1`` for SSCR; the
      check is hardcoded as ``mod(nsteps, 2 * exchange_period) == 0``.
    - ``mod(nsteps, rstout_period) == 0`` is a ``Setup_Remd_Restart``
      requirement.
    """
    integrator: str = "VRES"
    timestep_ps: float = 0.0035
    nsteps: int = 3_000_000
    exchange_period: int = 3_000
    eneout_period: int = 300
    crdout_period: int = 300
    rstout_period: int = 30_000
    nbupdate_period: int = 6
    elec_long_period: int = 2
    thermostat_period: int = 6
    barostat_period: int = 6
    ensemble: str = "NVT"
    tpcontrol: str = "BUSSI"
    temperature_K: float = 300.0
    pressure_bar: float = 1.0
    cutoffdist_A: float = 8.0
    pairlistdist_A: float = 9.5
    analysis_grest: bool = True

    def __post_init__(self) -> None:
        if self.integrator not in ("LEAP", "VVER", "VRES"):
            raise ValueError(
                f"GrestStage.integrator must be one of LEAP/VVER/VRES, "
                f"got {self.integrator!r}"
            )
        if self.timestep_ps <= 0:
            raise ValueError(
                f"GrestStage.timestep_ps must be > 0, got {self.timestep_ps}"
            )
        if self.ensemble not in ("NVE", "NVT", "NPT"):
            raise ValueError(
                f"GrestStage.ensemble must be one of NVE/NVT/NPT, "
                f"got {self.ensemble!r}"
            )
        if self.temperature_K <= 0:
            raise ValueError(
                f"GrestStage.temperature_K must be > 0, got {self.temperature_K}"
            )
        if self.exchange_period <= 0:
            raise ValueError(
                f"GrestStage.exchange_period must be > 0, "
                f"got {self.exchange_period}"
            )
        if self.rstout_period <= 0:
            raise ValueError(
                f"GrestStage.rstout_period must be > 0, "
                f"got {self.rstout_period}"
            )
        if self.nsteps <= 0:
            raise ValueError(f"GrestStage.nsteps must be > 0, got {self.nsteps}")
        # Setup_Remd hard requirement.
        if self.nsteps % (2 * self.exchange_period) != 0:
            raise ValueError(
                f"GrestStage.nsteps ({self.nsteps}) must be a multiple of "
                f"2*exchange_period ({2 * self.exchange_period}). "
                f"GENESIS Setup_Remd rejects otherwise."
            )
        if self.nsteps % self.rstout_period != 0:
            raise ValueError(
                f"GrestStage.nsteps ({self.nsteps}) must be a multiple of "
                f"rstout_period ({self.rstout_period})."
            )
        if self.pairlistdist_A < self.cutoffdist_A:
            raise ValueError(
                f"GrestStage.pairlistdist_A ({self.pairlistdist_A}) "
                f"must be >= cutoffdist_A ({self.cutoffdist_A})"
            )


# ---------------------------------------------------------------------------
# Top-level config
# ---------------------------------------------------------------------------

@dataclass
class GrestBuildConfig:
    """End-to-end gREST_SSCR build configuration.

    Aggregates the input PDB, AMBER force-field selection, REST
    selection, replica ladder, three stage protocols, MPI mapping, and
    external tool paths.

    Cross-field invariant
    ---------------------
    The lowest replica temperature must equal :attr:`GrestStage.temperature_K`
    so that the lowest-T replica is the canonical (un-tempered)
    simulation. Enforced in :meth:`__post_init__`.
    """
    # Input
    input_pdb: str = ""
    project_name: str = "grest_run"

    # Force field (AMBER only in v1.20)
    ff_protein: str = "leaprc.protein.ff19SB"
    ff_water: str = "leaprc.water.tip3p"
    ff_extra: List[str] = field(default_factory=list)
    solvatebox_padding_A: float = 10.0
    neutralize: bool = True
    salt_concentration_M: float = 0.0

    # REST + ladder
    rest_selection: RESTSelectionSpec = field(default_factory=RESTSelectionSpec)
    replica_temperatures: ReplicaTemperatureSpec = field(
        default_factory=ReplicaTemperatureSpec
    )

    # Protocols
    minimize: MinimizationStage = field(default_factory=MinimizationStage)
    equilibrate: EquilibrationStage = field(default_factory=EquilibrationStage)
    grest: GrestStage = field(default_factory=GrestStage)

    # MPI / execution (local mpirun only in v1.20)
    mpi_processes_per_replica: int = 2
    omp_num_threads: int = 2
    mpirun_path: str = "mpirun"

    # External tool paths
    tleap_path: str = "tleap"
    spdyn_path: str = "spdyn"
    atdyn_path: str = "atdyn"
    remd_convert_path: str = "remd_convert"
    cpptraj_path: str = "cpptraj"

    # Output / reproducibility
    output_dir: str = "."
    seed: Optional[int] = None

    def __post_init__(self) -> None:
        if not self.input_pdb:
            raise ValueError("GrestBuildConfig.input_pdb is required.")
        if not self.project_name:
            raise ValueError("GrestBuildConfig.project_name must be non-empty.")
        if self.solvatebox_padding_A <= 0:
            raise ValueError(
                f"GrestBuildConfig.solvatebox_padding_A must be > 0, "
                f"got {self.solvatebox_padding_A}"
            )
        if self.salt_concentration_M < 0:
            raise ValueError(
                f"GrestBuildConfig.salt_concentration_M must be >= 0, "
                f"got {self.salt_concentration_M}"
            )
        if self.mpi_processes_per_replica < 1:
            raise ValueError(
                f"GrestBuildConfig.mpi_processes_per_replica must be >= 1, "
                f"got {self.mpi_processes_per_replica}"
            )
        if self.omp_num_threads < 1:
            raise ValueError(
                f"GrestBuildConfig.omp_num_threads must be >= 1, "
                f"got {self.omp_num_threads}"
            )
        # Cross-field: lowest-T replica must equal grest temperature.
        rt = self.replica_temperatures
        if rt.mode == "manual":
            t_low = rt.temperatures[0]
        else:
            t_low = rt.T_min
        if abs(t_low - self.grest.temperature_K) > 0.01:
            raise ValueError(
                f"GrestStage.temperature_K ({self.grest.temperature_K}) "
                f"must match the lowest replica temperature ({t_low}); "
                f"GENESIS REMD treats the first ladder entry as the "
                f"canonical simulation."
            )

    def to_json(self, path: str) -> None:
        """Save configuration to JSON for reproducibility."""
        Path(path).write_text(
            json.dumps(asdict(self), indent=2, ensure_ascii=False)
        )

    @classmethod
    def from_json(cls, path: str) -> "GrestBuildConfig":
        """Load configuration from JSON.

        Handles all 5 nested dataclass fields (``rest_selection``,
        ``replica_temperatures``, ``minimize``, ``equilibrate``,
        ``grest``).
        """
        raw = json.loads(Path(path).read_text())
        rest = RESTSelectionSpec(**raw.pop("rest_selection"))
        ladder = ReplicaTemperatureSpec(**raw.pop("replica_temperatures"))
        mini = MinimizationStage(**raw.pop("minimize"))
        equil = EquilibrationStage(**raw.pop("equilibrate"))
        grest_stage = GrestStage(**raw.pop("grest"))
        return cls(
            rest_selection=rest,
            replica_temperatures=ladder,
            minimize=mini,
            equilibrate=equil,
            grest=grest_stage,
            **raw,
        )
