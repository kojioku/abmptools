# `abmptools.genesis.grest` — GENESIS gREST_SSCR Builder + Analysis

This document is the subsystem reference for the
`abmptools.genesis.grest` package, introduced in v1.20.0.
It builds and analyses GENESIS gREST_SSCR (generalized
Replica-Exchange with Solute Tempering — Solute Side-Chain
Repartitioning) systems end-to-end from a single JSON config.

For a step-by-step walkthrough see
[`tutorial_grest.md`](./tutorial_grest.md).

---

## 1. Scope and license

**What it does**: takes a protein PDB, generates AMBER `prmtop` +
`coor` via `tleap`, picks REST solute residues (explicit list or
`around` mask), produces a temperature ladder, renders the four
GENESIS `.inp` files (`step1_minimize` / `step2_equilibrate` /
`step3_grest` / `step5_remd_convert`), and emits a `run.sh`
orchestrator. After the user runs MD, the `analyze` subcommand
post-processes the trajectories: replica transition plot, acceptance
ratio, parameter-sorted lowest-T trajectory, and a 1D distance PMF.

**License posture**: GENESIS itself is **LGPL-3.0-or-later**. abmptools
**does not bundle** any GENESIS source or binaries; it shells out via
`subprocess`, which is *mere aggregation* under LGPL §5/§6 and
compatible with abmptools' **Apache-2.0** (v1.23.0+; ≤ v1.22.0 was MIT).
Users build GENESIS themselves from
<https://github.com/genesis-release-r-ccs/genesis>.

AmberTools (`tleap`, `cpptraj`) is free for both academic and
commercial redistribution.

---

## 2. Quick start

```bash
# Install the package in editable mode + grest extras (matplotlib).
cd abmptools/
python -m pip install -e .[grest]

# Generate a default config and inspect it.
python -m abmptools.genesis.grest example > /tmp/cfg.json

# Verify external tools and config self-consistency.
python -m abmptools.genesis.grest validate --config /tmp/cfg.json

# Build the .inp set (writes to ./grest_run/).
python -m abmptools.genesis.grest build --config /tmp/cfg.json -o ./grest_run

# Run MD (user-orchestrated; abmptools merely emits run.sh).
bash ./grest_run/run.sh

# Post-MD analysis.
python -m abmptools.genesis.grest analyze \
    --config /tmp/cfg.json --run-dir ./grest_run \
    --what all --mask1 "rno:96 and an:NZ" --mask2 "rno:274 and an:OD1"
```

---

## 3. Pipeline

### 3.1 Build (5 stages)

```
1. tleap parameterise          tleap         -> system.prmtop + .coor + ref.pdb
2. REST residue resolution     pure / cpptraj -> rest_residues.txt
3. Temperature ladder          pure          -> temperature_ladder.txt
4. GENESIS .inp rendering      pure          -> 4 .inp files in inp/
5. run.sh + HPC scaffolds      pure          -> run.sh + _pjsub / _sbatch
```

All five stages are *idempotent* given the same config + `output_dir`:
re-running overwrites the previous outputs.

### 3.2 Output layout

```
output_dir/
├── input/
│   └── config.json                    # frozen config dump
├── build/
│   ├── system.tleap                   # auto-rendered tleap script
│   ├── system.prmtop                  # AMBER topology
│   ├── system.coor                    # AMBER inpcrd
│   ├── system_ref.pdb                 # reference PDB
│   └── system_tleap.log               # tleap stdout/stderr
├── inp/
│   ├── step1_minimize.inp
│   ├── step2_equilibrate.inp
│   ├── step3_grest.inp                # the gREST_SSCR control file
│   └── step5_remd_convert.inp
├── logs/                              # populated by run.sh
├── run.sh                             # mpirun -np N spdyn orchestrator
├── run_pjsub.sh                       # Fugaku scaffold (comment-only)
└── run_sbatch.sh                      # SLURM scaffold (comment-only)
```

### 3.3 Analysis

`analyze --what` accepts any subset of `{sort, transition,
acceptance, pmf}` (or `all`).

| Task | Output | Implementation |
|---|---|---|
| `sort` | `build/param1.dcd` (lowest-T) | `remd_convert` subprocess |
| `transition` | `analysis/replica_transition.png/.csv` | `.rem` parsing + matplotlib |
| `acceptance` | `analysis/acceptance_ratio.png/.csv` | `REMD>` log parsing + matplotlib |
| `pmf` | `analysis/dist_pmf.png/.xvg` | cpptraj distance + numpy histogram + `-kT log P` |

---

## 4. Public API

### 4.1 Dataclasses

```python
from abmptools.genesis.grest.models import (
    GrestBuildConfig,
    RESTSelectionSpec,
    ReplicaTemperatureSpec,
    MinimizationStage,
    EquilibrationStage,
    GrestStage,
)
```

Cross-field invariant: `GrestStage.temperature_K` must equal the lowest
ladder temperature (the canonical "real" replica). Enforced in
`GrestBuildConfig.__post_init__`.

### 4.2 REST residue selection

Two modes:

- **`mode="explicit"`** — provide a `residues` list of strings:
  ```python
  RESTSelectionSpec(mode="explicit", residues=["1-138"])
  RESTSelectionSpec(mode="explicit", residues=["21,96,274-275"])
  ```
- **`mode="around"`** — provide a `center` mask + `radius_A`:
  ```python
  RESTSelectionSpec(mode="around", center="rno:96", radius_A=5.0)
  ```
  This is resolved at build time via `cpptraj` `<:radius` mask
  (requires `cpptraj` on the PATH).

The `param_types` field controls which interactions get tempered. Default
`["C", "L"]` (CHARGE + LJ) matches the SSCR-typical setting in tutorial
12.3 step3 and the POC log's `Setup_Remd_Solute_Tempering` block.

### 4.3 Temperature ladder

```python
ReplicaTemperatureSpec(mode="manual",
                       temperatures=[300.0, 318.11, 337.11, 357.10])
ReplicaTemperatureSpec(mode="auto", method="geometric",
                       T_min=300.0, T_max=420.0, n_replicas=12)
```

`method="vanderspoel"` is reserved for v1.20.x and currently raises
`NotImplementedError`.

### 4.4 Builder

```python
from abmptools.genesis.grest.builder import GrestBuilder
result = GrestBuilder(cfg).build()
# result is a Dict[str, Any] with keys:
#   output_dir, build_dir, inp_dir, config_json,
#   prmtop, coor, ref_pdb, leap_log, box_size_A,
#   n_protein_residues, rest_residues, rest_selection_string,
#   n_rest_residues, temperature_ladder, n_replicas,
#   inp_files, run_script, run_pjsub, run_sbatch
```

### 4.5 Analysis

```python
from abmptools.genesis.grest import analysis

# (a) parameter sort
analysis.run_remd_convert(cfg, inp_path, workdir)
# (b) replica transition
analysis.plot_replica_transition(rem_files, out_png, out_csv)
# (c) acceptance ratio
analysis.plot_acceptance_ratio(log_path, out_png, burn_in=100)
# (d) 1D distance PMF
distances = analysis.compute_distance_timeseries_cpptraj(
    prmtop, dcd, mask1="rno:96 and an:NZ", mask2="rno:274 and an:OD1",
    workdir=out_dir, cpptraj_path=cfg.cpptraj_path)
analysis.compute_distance_pmf(distances, T_K=300.0, out_png, out_xvg)
```

---

## 5. CLI

```
python -m abmptools.genesis.grest {example, validate, build, analyze}
```

| Subcommand | Required flags | Behaviour |
|---|---|---|
| `example` | (none) | prints a default JSON config (POC-style 4 replicas, 138 explicit residues) |
| `validate` | `--config` | external tool resolution, REST preview, ladder + ratios, MPI mapping |
| `build` | `--config` (`-o` overrides `output_dir`) | runs the 5-stage pipeline |
| `analyze` | `--config`, `--run-dir`, `--what` | runs any subset of `{sort,transition,acceptance,pmf,all}`; `pmf` also requires `--mask1` and `--mask2` |

Every subcommand accepts `-v / --verbose` to flip the `abmptools.genesis.grest`
logger to DEBUG.

---

## 6. Design rationale

### 6.1 Why subprocess (no Python binding) for GENESIS?

GENESIS is LGPL-3.0-or-later. Python bindings would require dynamic
linking (and possibly LGPL-aware redistribution); subprocess is the
simplest "mere aggregation" path that preserves abmptools' Apache-2.0
posture. Coverage of GENESIS' control-file interface (the four `.inp`
templates) is sufficient for the gREST_SSCR workflow.

### 6.2 Why AMBER, not CHARMM?

The POC uses AMBER ff19SB + TIP3P; tutorial 12.3 uses the same. CHARMM
support would require psfgen or CHARMM-GUI input, which adds another
external dependency without a matching POC reference. v1.20.0 ships
AMBER-only; CHARMM may be added later.

### 6.3 Why `param_type1 = C L`?

Per the POC log's `Setup_Remd_Solute_Tempering>` output, only `CHARGE`
and `LJ` were `T` in the validated production run. SSCR (Solute
Side-Chain Repartitioning) targets sidechain torsions implicitly via the
LJ + electrostatic temperning of the selected residues, without needing
explicit dihedral tempering for protein side-chains.

### 6.4 Why `around` mode?

The POC's manual edit history (`LYS21` → `21`, late residue-list
correction) shows that residue picking is a frequent source of error.
`around mode` lets the user say "anything within 5 Å of residue 96"
and have `cpptraj` resolve the actual list at build time. The
`validate` subcommand previews what `explicit` mode would have given
for cross-checking.

### 6.5 Why local `mpirun` only?

Fugaku (PJSUB) and SLURM (sbatch) submission scripts are emitted as
comment-only scaffolds in the build dir, but abmptools does not run
integration tests against any HPC scheduler. Each site has its own
queue/qos/runtime conventions that cannot be portably encoded; the
scaffold gives users a known-working starting point.

### 6.6 Why `nsteps mod (2 * exchange_period) == 0`?

GENESIS's `Setup_Remd` rejects otherwise. Enforced in
`GrestStage.__post_init__` so the user never sees this error from
spdyn at runtime. POC log: `Setup_Remd> mod(nsteps,(2*exchange_period*dimension)) must be zero`.

---

## 7. Caveats and known issues

- **GENESIS build under icx**: `fileio_data_.c` calls `ftello64` /
  `fseeko64` which clang's icx doesn't auto-import. Patch:
  `sed -i 's/ftello64/ftello/g; s/fseeko64/fseeko/g' analysis/src/lib/fileio_data_.c`,
  then `./configure CC=icx && make`.
- **Acceptance ratio inheritance**: when restarting from `step3.rst`,
  the first ~100 exchange events show non-1 acceptance counts inherited
  from the prior run. Default `burn_in=100` excludes them.
- **`convert_ids = `** (empty) in `step5_remd_convert.inp` produces the
  lowest-T trajectory only (POC convention). Set
  `--convert-ids 1 2 3 4` if you need every parameter index.
- **Salt concentration > 0**: v1.20.0 emits a `# salt_concentration_M`
  comment in the tleap script but does not auto-derive an exact ion
  count. For non-zero salt, add `addions system Na+ N` / `addions system
  Cl- N` manually before `saveAmberParm`.
- **`vanderspoel` ladder method**: deferred. Use `manual` or `geometric`
  in v1.20.0.

---

## 8. Default values (selected)

| Field | Default | Rationale |
|---|---|---|
| `ff_protein` | `leaprc.protein.ff19SB` | POC + tutorial 12.3 |
| `ff_water` | `leaprc.water.tip3p` | POC + tutorial 12.3 |
| `solvatebox_padding_A` | 10.0 | POC value |
| `rest_selection.param_types` | `["C", "L"]` | POC `Setup_Remd_Solute_Tempering>` (CHARGE=T LJ=T) |
| `replica_temperatures.{T_min, T_max, n_replicas}` | `300.0, 357.10, 4` | POC ladder |
| `minimize.nsteps` | 10 000 | tutorial 12.3 step1 |
| `minimize.cutoffdist_A` | 12.0 | POC + AMBER convention |
| `equilibrate.timestep_ps` | 0.002 | dt = 2 fs (HMR enables 2 fs without bond-constraint deadlock) |
| `equilibrate.nsteps` | 25 000 | 50 ps NPT relaxation |
| `equilibrate.ensemble` | NPT | density relaxation before gREST NVT |
| `grest.timestep_ps` | 0.0035 | tutorial 12.3 step3 (3.5 fs r-RESPA) |
| `grest.nsteps` | 3 000 000 | tutorial 12.3 step3 (10.5 ns); POC scales to 60 M (210 ns) |
| `grest.exchange_period` | 3 000 | exchange every 10.5 ps at dt=3.5 fs |
| `grest.ensemble` | NVT | tutorial 12.3 step3 |
| `mpi_processes_per_replica` | 2 | POC: 8 total / 4 replicas = 2 |
| `omp_num_threads` | 2 | matches POC's MPI/OpenMP balance |

---

## 9. References

- POC notes: `/home/okuwaki/llm-project/SI/grest-POC.md`
- GENESIS: <https://github.com/genesis-release-r-ccs/genesis>
- Tutorial 12.3 (gREST_SSCR): <https://mdgenesis.org/tutorials/genesis_tutorial_12.3_2022/>
- gREST_SSCR paper: Re, Kim & Sugita, *Int. J. Mol. Sci.* 22:270 (2021),
  <https://www.mdpi.com/1422-0067/22/1/270>
