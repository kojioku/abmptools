# `abmptools.crystal` — Organic-Crystal FMO Pipeline

Subsystem reference for the `abmptools.crystal` package, introduced
in v1.23.0 (Phase A/B skeleton + regression fixture) and completed
in v1.24.0 (Phase C: full orchestrator + 8-subcommand CLI) /
v1.25.0 (Phase D: `--run-local` + ASE PDB writer).

For step-by-step ops see [`tutorial_crystal_fmo.md`](./tutorial_crystal_fmo.md).

---

## 1. Scope

**What it does**: takes one or more small-molecule crystal CIF files
and produces ABINIT-MP FMO calculation inputs (`.ajf`) along with HPC
jobscripts. The pipeline runs five stages:

1. **CIF → supercell** — expand the asymmetric unit into a
   `layer^3`-cell supercell (legacy parser or ASE backend)
2. **PDB → for_abmp/** — fragment cut around a target solute
   molecule (configurable `cutmode`) and emit the FMO `.ajf` (with
   full-precision `&XYZ` block when `is_xyz=True`)
3. **HPC jobscripts** — render PJM/SLURM/PBS/`local` jobscripts and a
   batch submitter (`runbatch.sh`)
4. **Local run (optional)** — invoke `abinitmp` directly per AJF when
   `--run-local` is set (Phase D-1)
5. **Postprocessing (optional)** — IFIE/PIEDA CSV via
   `getifiepieda` + nearest-atom annotation (Phase C-6)

Sibling subpackages: `cg.peptide` / `cg.membrane` / `genesis.grest` /
`genesis.mmgbsa` / `fragmenter` (all in v1.18.0–1.22.0).

---

## 2. Configuration schema (`CrystalBuildConfig`)

YAML / JSON round-trip via `to_yaml` / `from_yaml` / `to_json` /
`from_json`. The schema is one top-level dataclass (`CrystalBuildConfig`)
plus six leaves:

| Field | Type | Notes |
|---|---|---|
| `inputs` | `List[CIFInputSpec]` | One entry per CIF (`cif`, `layer`, `atoms_in_mol`, `name`) |
| `cif_engine` | `CIFEngineConfig` | `engine: 'legacy' \| 'ase'`, `bond_tolerance` |
| `fragment` | `FragmentTemplate` | `cutmode`, `solutes`, `criteria`, `molname`, `template_ajf` (path required) |
| `fmo` | `FMOMethod` | `method`, `basis_set`, `memory`, `is_xyz` (default True for full-precision `&XYZ` block) |
| `hpc` | `HPCJobSpec` | `scheduler: 'PJM' \| 'SLURM' \| 'PBS' \| 'local'`, queue, nodes, etc. |
| `postproc` | `PostProcessSpec` | `enable`, `frag_target`, `nearest_atom_count` |
| `output_dir` / `project_name` / `fail_fast` / `run_local` | top-level | |

Example (`sample/crystal/csp7_smoke/crystal.yaml`):

```yaml
project_name: csp7_r00001_smoke
output_dir: ./out
inputs:
  - cif: XXXI-MMFF-R00001.cif
    layer: 5
    atoms_in_mol: [32]
cif_engine:
  engine: legacy
fragment:
  cutmode: around
  solutes: [0]
  criteria: 6.0
  molname: [UNK]
  template_ajf: ./UNK.ajf
fmo:
  method: MP2
  basis_set: 6-31Gdag
  memory: 6000
  is_xyz: true
hpc:
  scheduler: PJM
  queue: small
  group: hp190133
  nodes: 12
  proc_per_node: 2
  omp_threads: 24
  abinit_dir: /data/hp190133/programs/ABINIT-MP
```

---

## 3. CIF backends: legacy vs ASE

| Aspect | `engine='legacy'` (default) | `engine='ase'` |
|---|---|---|
| Implementation | `abmptools.readcif` self-rolled parser (982 lines) | `ase.io.read` + `ase.spacegroup.crystal` |
| Symmetry operations | Hard-coded table (1/3, 2/3, 1/4, 3/4, 1/6, 5/6 + axes) | `_symmetry_equiv_pos_as_xyz` parser, all standard groups |
| Layer count | Hard-coded if-else for 1, 4, 5 (extends via shifts) | Any `layer >= 1` via `Atoms.repeat((N, N, N))` |
| Output coord origin | Centred (asymmetric unit at lattice midpoint) | Cell origin (asymmetric unit at corner) |
| Atom ordering | Inherits `readcif` order | ASE `Atoms` order (per-element groups) |
| Byte-equivalence with v1.22.0 | **Yes** (Phase B fixture covers R00001/R00002/R00004) | **No** — different supercell origin and atom order; valid AJF emitted but `Natom` and `&XYZ` block size differ |

**When to switch to ASE**: the legacy parser silently fails on
unsupported symmetry operators (e.g. trigonal/hexagonal `1/3+x-y`
with rare denominators). ASE handles all standard space groups via
its general-form xyz parser. ASE also makes `layer` arbitrary —
useful when you want layer 7 or 9 for convergence tests.

**When to stay on legacy**: byte-equivalence with the historical
csp7 1500-structure pipeline matters; layer ∈ {1, 4, 5} suffices;
asymmetric unit origin is the convention you've validated.

---

## 4. `cutmode` semantics

| Mode | What it cuts | When to use |
|---|---|---|
| `around` | Molecules whose any heavy atom is within `criteria` Å of any atom of `solutes[0]` | Default for csp7-style local FMO around one target molecule |
| `sphere` | Molecules whose any atom is within `criteria` Å of `tgtpos` (3D point) | Fixed-point solvation analysis |
| `cube` | Molecules whose any atom lies inside the `[tgtpos[i] ± criteria[i]]` box | Block extraction along principal axes |
| `neutral` | Like `around` but adds counterions to keep total charge = 0 | Charged solute systems |
| `none` | All molecules in the input PDB are kept | Hand-curated input PDBs |

`solutes[0]` is the 1-based residue ID of the target molecule (1 in
csp7's convention since each cocrystal molecule occupies one residue
slot in the supercell PDB).

---

## 5. AJF format: direct (`is_xyz=True`) vs `ReadGeom` (`is_xyz=False`)

The `is_xyz` flag on `FMOMethod` controls the AJF body layout:

- **`is_xyz=True` (default for crystal pipelines)** — the AJF carries
  `Natom={N}` plus a populated `&XYZ` block with full-precision
  coordinates (e.g. `2.891228423347176`). Required when supercell
  coordinates have decimals beyond `%8.3f` (the PDB column-fixed
  format), which is typically the case for ASE-derived or
  high-precision CIF data. `pdb_io.exportardxyzfull` reads coordinates
  from a sibling `<base>.xyz` next to the input PDB, so the XYZ file
  must be present in the working directory when invoking `pdb2fmo`.
- **`is_xyz=False`** — `Natom=0` and `ReadGeom='<base>.pdb'` redirects
  ABINIT-MP to read coordinates from the PDB at runtime. Coordinates
  are limited to PDB precision (`%8.3f`, i.e. 3 decimal Å). Useful
  when working strictly within PDB conventions and full-precision
  coordinates are not required.

The csp7 v1.22.0 baseline used `is_xyz=True` (matched in the Phase B
fixture); the bytes of the `&XYZ` block are byte-equivalent to the
historical output.

---

## 6. HPC job templates

Four schedulers (`PJM` / `SLURM` / `PBS` / `local`), each rendered
from a `string.Template` body in `job_templates.py`:

- **PJM** — Fugaku-style `pjsub` directives, mirrors the legacy
  `runab23q2_fgk.sh`. Includes `module switch lang/tcsds-1.2.31`,
  `mkinp_openver1rev20.py` preprocessor, `mpiexec -stdin` invocation.
- **SLURM** — `#SBATCH` directives, `mpirun -np {total_proc}`.
- **PBS** — PBS Pro / Torque, `select=N:ncpus=...:mpiprocs=...`.
- **`local`** — plain `bash` invocation; used by `--run-local` and
  small smoke runs.

Override the bundled templates via `HPCJobSpec.template_override`
(path to a custom `string.Template` text file) when site-specific
flags are required.

The orchestrator also writes a `runbatch.sh` driver alongside the
per-AJF jobscripts; it loops `pjsub` / `sbatch` / `qsub` over each
script. Running a single submit command `bash runbatch.sh` from the
output directory submits the entire job set.

---

## 7. CLI subcommands

`abmp-crystal {example,validate,expand,fragment,jobs,pipeline,postproc,nearest}`:

| Subcommand | Stages | Use case |
|---|---|---|
| `example` | — | Print the example YAML config to stdout |
| `validate` | — | Check config + optional dependencies (`ase` / `pyyaml` / `abinitmp`) |
| `expand` | 1 | CIF → `cifout/layer<L>/{pdb,xyz}/` only |
| `fragment` | 1+2 | CIF → `for_abmp/*.{ajf,pdb}` |
| `jobs` | 1+2+3 | + HPC jobscripts + `runbatch.sh` |
| `pipeline` | all | Full pipeline; add `--run-local` to invoke `abinitmp` |
| `postproc` | 5 | Stage 5 standalone (logs → IFIE/PIEDA + nearest atoms) |
| `nearest` | — | Standalone nearest-atom lookup on a PDB |

Both YAML (`*.yaml` / `*.yml`) and JSON config files are accepted;
the file extension dispatches to `from_yaml` / `from_json`.

---

## 8. Output directory layout

For `output_dir = ./out` with one CIF input named
`XXXI-MMFF-R00001.cif`:

```
out/
├── crystal_config.resolved.json     ← snapshot of the resolved config
└── XXXI-MMFF-R00001/
    ├── XXXI-MMFF-R00001.cif         ← copied input
    └── cifout/
        └── layer5/
            ├── pdb/
            │   ├── XXXI-MMFF-R00001layer5Zp1.pdb
            │   ├── XXXI-MMFF-R00001layer5Zp1.xyz
            │   ├── input_param         ← drivers (synthesised from config)
            │   ├── segment_data.dat
            │   ├── UNK.ajf             ← copied template
            │   └── for_abmp/
            │       ├── XXXI-MMFF-R00001layer5Zp1-around_ar6.0.ajf  ← FMO input
            │       ├── XXXI-MMFF-R00001layer5Zp1-around_ar6.0.pdb
            │       ├── XXXI-MMFF-R00001layer5Zp1-around_ar6.0.xyz
            │       ├── *_12n-2p-24t.sh ← rendered PJM jobscript
            │       └── runbatch.sh     ← batch submitter
        └── xyz/
            └── XXXI-MMFF-R00001layer5Zp1.xyz
```

When `postproc.enable=true`, the orchestrator also writes
`out/postproc/` containing the IFIE/PIEDA CSVs (`csv/*.csv`) and
optional nearest-atom-annotated CSVs.

---

## 9. Failure modes

- **Unsupported space group (legacy)**: silent silent-fail
  (`getsymcoord` returns the asymmetric unit only). Switch to
  `engine='ase'`.
- **Wrong `atoms_in_mol`**: `detect_molecules` raises `ValueError`
  with the observed sizes; either correct `atoms_in_mol` or adjust
  `cif_engine.bond_tolerance`.
- **Missing `template_ajf`**: `CrystalOrchestrator.run()` raises
  `ValueError` requesting a path to a minimal ABINIT-MP AJF template.
  The csp7 sample (`sample/crystal/csp7_smoke/UNK.ajf`) is the
  canonical template.
- **`abinitmp` binary not found in `--run-local`**: `run_abinit`
  raises `FileNotFoundError` listing the resolution attempts. Set
  `HPCJobSpec.abinit_dir` + `binary_name` or place `abinitmp` on
  PATH.
- **`pdb2fmo` failure** (e.g. atoms_in_mol mismatch propagates):
  raises `RuntimeError` from the in-process call; the underlying
  exception is in the stack trace (no subprocess opacity).

---

## 10. License notes

- Package licence: MIT (planned re-licence to Apache-2.0; see
  `project_abmptools_apache2_relicense_plan.md` in the dev memory).
- ASE: LGPL-2.1+, dynamically imported. No source modification or
  vendoring; mere aggregation.
- abinitmp: external binary; never bundled.

The crystal subpackage adds no GPL-only or copyleft-only dependencies.
