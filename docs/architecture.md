# Architecture

## Design Overview

ABMPTools follows an **inheritance-based layered architecture**:

- A **core I/O layer** provides coordinate reading, molecular math, and fragment management.
- **Data managers** (`CPFManager`, `LOGManager`) handle ABINIT-MP-specific file formats independently.
- **Analysis and setup modules** compose the above layers into CLI-driven workflows.
- A **Fortran extension** accelerates the performance-critical IFIE/PIEDA reading path.

Key design patterns:
- **Class inheritance chain** for progressive capability building.
- **Pandas DataFrames** as the universal internal data representation.
- **argparse-based CLIs** invoked via `python -m abmptools.<module>`.
- **`multiprocessing.Pool`** for embarrassingly parallel workloads.

## Class Inheritance Hierarchy

```
molcalc                          (coordinate math, distance, rotation, PBC)
  └── udf_io                    (OCTA UDF file reading)
        └── udfrm_io            (UDF → PDB/XYZ conversion)

mol_io                           (base molecular I/O — XYZ, coordinate conversion)
  └── abinit_io                 (ABINIT-MP I/O — fragment definitions)
        └── pdb_io              (PDB parsing — residues, atoms, read modes)
              └── anlfmo        (advanced FMO analysis — log parsing, IFIE)

setfmo(pdb_io, udfcreate, udfrm_io)   (FMO setup — multiple inheritance)

CPFManager                       (standalone — CPF parse/write)
LOGManager                       (standalone — ABINIT-MP log parse)
```

**Assumption**: The inheritance hierarchy is inferred from class definitions in each module. `setfmo` uses multiple inheritance combining PDB, UDF creation, and UDF conversion capabilities.

## Core I/O Layer

### `mol_io.py` — Molecular I/O Base

The foundation class. Provides:
- `read_xyz()` — Read XYZ coordinate files.
- `getatoms()` — Extract atom lists from structures.
- `convert_xyzs_pdb()` — XYZ to PDB format conversion.

### `molcalc.py` — Coordinate Mathematics

Standalone utility class for:
- Distance calculation (`getdist()`), center of mass (`getCenter()`).
- Molecular translation (`moveMolTrans()`) and rotation (`moveMolEuler()`).
- Periodic boundary condition handling.
- LAMMPS data/trajectory parsing (`parse_lammps_data()`, `parse_lammps_trajectory()`).
- Uses `multiprocessing.Pool` for distance matrix calculations.

### `abinit_io.py` — ABINIT-MP I/O Base

Extends `mol_io` with:
- Fragment definition management.
- Fragment geometry extraction.
- Support for various ABINIT-MP calculation modes and options.
- Parallel processing for fragment setup.

### `pdb_io.py` — PDB Parser

Extends `abinit_io` with:
- `readpdb()` — Parse PDB files with modes: `TER`, `resnum`, `rfile`.
- Residue and atom property extraction.
- Extended atom numbering support (5+ digit atoms).
- Optional `UDFManager` integration for OCTA MD data.

## Data Managers

### `CPFManager` (`cpfmanager.py`)

Central data hub for CPF files. Key methods:
- `parse(filename)` — Read CPF files (auto-detects version).
- `read_header()`, `read_atominfo()`, `read_fraginfo()`, `read_dimer()` — Section parsers.
- `write(title, filename)` — Generate CPF output.

Stores data in pandas DataFrames:
- `atominfo` — Per-atom properties (coords, charges, residue info).
- `fraginfo` — Fragment-level data.
- `diminfo` — Dimer interaction data (distance, energy components).
- `mominfo` — Monomer/dipole moment information.

Supported CPF versions: 4.201, 7.0, 10, 23 (including `.gz` compressed files).

### `LOGManager` (`logmanager.py`)

Parses ABINIT-MP calculation log files. Key methods:
- `parse(filename)` — Full log parsing.
- `getversion()` — Detect ABINIT-MP version from log header.
- `getcondition()` — Extract calculation conditions.
- `getfraginfo()` — Extract fragment definitions.

Returns: `condition` dict, `fraginfo` dict, `static_data` dict.
Supports nucleic acids and protein-nucleic acid complexes (added v1.14.6).

## Analysis Modules

### `getifiepieda.py` — IFIE/PIEDA Extraction

The primary analysis CLI. Supports multiple modes:

| Mode | Flag | Description |
|------|------|-------------|
| Fragment pair | `--frag i j` | IFIE between specified fragments |
| Distance filter | `--frag i -d r` | Fragments within distance `r` of fragment `i` |
| Molecule-based | `--mol i -d r` | Interactions from molecule `i` within distance `r` |
| Fragment matrix | `--ffmatrix i1-i2 j1-j2` | Full i×j IFIE matrix |
| Time-fragment matrix | `--tfmatrix i-j k-l` | Time-series fragment interaction matrix |
| Multi-sample | `--multi i -t start end interval` | Time-series IFIE with parallel reading |
| In-molecule | `--fraginmol i j MOL k` | Intra-molecular fragment interactions |

Uses the Fortran library (`readifiepiedalib.so`) when available; falls back to Python with `-nof90`.

### `anlfmo.py` — Advanced FMO Analysis

Inherits from `pdb_io`. Provides:
- `readlog()` — Comprehensive log parsing with fragment interaction extraction.
- Time-dependent IFIE calculation.
- Multi-sample analysis with parallel processing.
- Column definitions for CSV output formatting.

### `generate_difie.py` — Dynamic IFIE

Creates time-averaged CPF files from MD trajectory snapshots:
- `getcpfobj()` — Load CPF at specific timestep.
- `getavestddf()` — Compute mean/std across trajectory.
- Outputs a single CPF with `M-` (mean) and `S-` (std) prefixed columns.
- Parallel processing via `-np` flag.

## FMO Setup Pipeline

### `generateajf.py` — AJF Template Generation

Full-featured CLI for ABINIT-MP input file creation:
- Input: PDB file or config dictionary.
- Options: method (HF/MP2), basis set, solvation (PB), PIEDA, CPF output, RESP, DGEMM, dispersion, BSSE.
- Version-aware output for different ABINIT-MP revisions.

### `setfmo.py` — FMO Setup Orchestrator

Multiple-inheritance class combining PDB, UDF, and configuration capabilities:
- `setrfmoparam()` — Configure FMO parameters.
- Cut modes: `sphere`, `cube`, `around`, `neutral`, `none`.
- Solute/solvent specification, ion handling (`remain`/`remove`).
- Fragment assignment for molecular aggregates.

### `pdb2fmo.py` — PDB to FMO Converter

CLI for non-protein systems (molecular aggregates, polymers):
- Reads config file with fragment definitions for one molecule.
- Assigns fragments to all molecules in the system.
- Produces AJF + PDB output set.

## File Conversion Modules

```
LOG  ──→ log2cpf   ──→ CPF
LOG  ──→ log2config ──→ Config dict
AJF  ──→ ajf2config ──→ segment_data.dat
CIF  ──→ readcif    ──→ Cartesian coordinates
CPF  ──→ convertcpf ──→ CPF (different version / filtered)
UDF  ──→ udfrm_io   ──→ PDB / XYZ
AJF  ──→ ajfserial  ──→ Numbered AJF files
PDB  ──→ pdbmodify  ──→ PDB (edited)
```

## Subpackages

In addition to the flat `abmptools/*.py` layer described above, the package ships several subpackages that each encapsulate a self-contained workflow:

- **`abmptools.gro2udf` / `abmptools.udf2gro`** — bidirectional conversion between GROMACS (`.gro / .top / .mdp / .itp`) and OCTA COGNAC UDF. See [`gro2udf.md`](gro2udf.md) / [`udf2gro.md`](udf2gro.md).
- **`abmptools.geomopt`** — pluggable geometry optimization (MACE / OpenFF-OpenMM / PySCF-DFT). Each backend is a lazy-imported class under `geomopt/`. See [`geomopt.md`](geomopt.md) and [`qmopt.md`](qmopt.md).
- **`abmptools.amorphous`** — multi-component amorphous builder (SMILES / SDF / PubChem CID → Packmol packing → OpenFF parameterization with AM1-BCC → GROMACS 5-stage annealing protocol + VMD-friendly PBC wrap script). The `amorphous.pubchem` submodule is a dependency-free (`urllib`-only) wrapper around the PubChem PUG REST API, raising `PubChemNo3DError` when no 3D conformer is available. See [`amorphous.md`](amorphous.md).
- **`abmptools.membrane`** — peptide-bilayer umbrella-sampling builder (packmol-memgen lipid + peptide PDB → AMBER ff19SB+Lipid21 via tleap+parmed → semiisotropic NPT equilibration → z-pulling → per-window MDPs → `gmx wham`). Designed to be commercial-license-clean: CGenFF web server and CHARMM-GUI auto-generation are forbidden by design. See [`membrane.md`](membrane.md).
- **`abmptools.cg.peptide`** — Martini 3 peptide CG system builder (residue sequence + counts + box → `tleap`/extended-backbone fallback for atomistic PDB → `martinize2 -ff martini3001` for CG mapping → `gmx insert-molecules` → `gmx solvate` (Martini W, auto-generated when needed) → `gmx genion` for NaCl 0.15 M → em/nvt/npt/md.mdp + `run.sh`). vermouth-martinize (Apache-2.0) is invoked via `subprocess` only (no source bundled or modified). See [`cg/peptide/README.md`](../abmptools/cg/peptide/README.md). New in v1.18.0.
- **`abmptools.cg.membrane`** — Martini 3 peptide-membrane PMF builder (umbrella sampling). Internally sub-calls `abmptools.cg.peptide` for the CG peptide; embeds it via `insane` (GPL-2.0, subprocess only) into a POPC bilayer; post-processes the topology (4 ITP includes + `Protein → molecule_0` rename + `NA+/CL- → NA/CL` normalization); generates 13-31 windows of NPT-semiisotropic umbrella MDPs (with `pull-group1-pbcatom` injected for the bilayer group); reuses `abmptools.membrane.{pulling,pmf,mdp_us_protocol}` helpers via duck-typed imports for `gmx wham` analysis. CG dt=20 fs gives 30-100× wall-time speedup over AA membrane. See [`cg_membrane.md`](cg_membrane.md) and [`tutorial_cg_membrane_us.md`](tutorial_cg_membrane_us.md). New in v1.19.0.
- **`abmptools.genesis.grest`** — GENESIS gREST_SSCR (generalized Replica-Exchange with Solute Tempering — Solute Side-Chain Repartitioning) builder + analysis. Drives `tleap` for AMBER ff19SB + TIP3P parameterization, resolves REST solute residues via explicit list or `cpptraj` `<:radius` mask, generates a geometric / manual temperature ladder, renders 4 GENESIS `.inp` files (minimize / equilibrate / grest / remd_convert), and emits `mpirun -np N spdyn` driven `run.sh`. Post-MD `analyze` subcommand: replica transition plot, acceptance ratio plot, 1D distance PMF (`-kT log P(r)` from cpptraj distance time-series). GENESIS (LGPL-3.0+) and AmberTools are subprocess-only (no bundling, no linking). See [`grest.md`](grest.md) and [`tutorial_grest.md`](tutorial_grest.md). First module under the `abmptools.genesis/` namespace, new in v1.20.0.
- **`abmptools.fragmenter`** — FMO automatic fragment splitter for small molecules / lipids / polymers. Loads a PDB via RDKit (4-strategy bond perception with proximity / CONECT-only / `sanitize=False` / `obabel` fallback), groups molecules by heavy-atom-only canonical SMILES, walks the heavy-atom graph diameter (2-pass BFS) accumulating MW + side-chain MW until ≥ `target_mw` (default 200 g/mol) and proposes C-C cuts (excluding ring / multi-bond / hetero-adjacent), then applies the cuts via `RWMol` and exports an `log2config`-compatible `segment_data.dat`. UI: `python -m abmptools.fragmenter {suggest,apply,example}` (CLI / SVG+JSON review bundle) or `AutoFragmenter.from_pdb` + `open_panel` (Jupyter ipywidgets). Polymer γ-path (`declare_same_pattern`) lets users explicitly map a master pattern (longest chain) to shorter copies via main-chain index alignment. See [`fragmenter.md`](fragmenter.md). New in v1.21.0.
- **`abmptools.genesis.mmgbsa`** — GENESIS atdyn-based MM/GBSA single-point ΔG_bind for protein-ligand binding free-energy estimation. 4-stage pipeline: split PDB into receptor + ligand (Biopython) → parameterize 3 systems (acpype with GAFF/GAFF2 + AM1-BCC for the ligand, then `tleap` for complex/ligand/receptor with AMBER ff14SB + DNA.OL15 + RNA.OL3 + TIP3P) → run `mpirun -np 1 atdyn` with `[ENERGY] implicit_solvent=GBSA` + NOBC + 1-step minimize → parse `[STEP4] Compute Single Point Energy` log + compute `ΔG_bind = E_complex - E_ligand - E_receptor` (per GENESIS doc 05_Energy.rst:564, the `ENERGY` column already contains `U_FF + ΔG_solv`) + emit CSV + matplotlib bar plot. Acpype (GPL-3.0) and GENESIS atdyn (LGPL-3.0+) are subprocess-only. CLI: `python -m abmptools.genesis.mmgbsa {example,validate,divide,parameterize,run,analyze,pipeline}` with both JSON config and folder-mode shortcut (`-i / -r / -c` for POC compatibility). See [`mmgbsa.md`](mmgbsa.md) and [`tutorial_mmgbsa.md`](tutorial_mmgbsa.md). Second module under `abmptools.genesis/`, new in v1.22.0.
- **`abmptools.hbond`** — Hydrogen-bond analyzer for COGNAC `.udf` / `.bdf` trajectories with two analysis modes. **imc mode** classifies each COOH into 4 species (cyclic dimer / chain / single COOH→amide / free) matching Yuan 2015 (Mol. Pharm. 12, 4518) NMR deconvolution Table 1; **generic mode** (v1.28+) reports donor-type × acceptor-type pair statistics for arbitrary systems (PVA / peptide / alcohols / mixtures). 12 modules covering UDFManager wrapping (`bdf_reader`), force-field-agnostic functional tagging with 4 built-in FF mappings + element + bond-graph fallback (`func_tags`, `functional_groups`, fallback makes OpenFF SMIRNOFF UDFs work without an antechamber GAFF patch), Luzar-Chandler / strict / custom geometric H-bond detection with orthogonal-PBC minimum image (`hbond_detector`), the IMC classifier and generic pair-type statistics (`classifier`, `pair_type_stats`), continuous / intermittent lifetime + unbiased Luzar-Chandler autocorrelation `C(t)` + τ_HB (`lifetime`), and four visualisation routes via `colorizer` (Mol_Name rename + Draw_Attributes for gourmet, autorun `.act` for gourmet, plain `.py` Python panel script for J-OCTA, per-atom Attributes tagging for J-OCTA filtering). CLI (`python -m abmptools.hbond ... --classify-mode {imc,generic} --donor-groups ... --acceptor-groups ...`) + Python API (`Analyzer`, `AnalyzerConfig`) + Jupyter ipywidgets UI (`open_panel(bdf_path)` with mode dropdown + functional-group checkboxes + RDKit 2D structure preview). Samples: `sample/hbond/imc_amorphous/` (IMC + Yuan 2015 NMR comparison plot) and `sample/amorphous/pva_amorphous/` (PVA 10-mer × 30 OpenFF Sage + AM1-BCC + 5-stage MD + generic-mode demo). See [`hbond.md`](hbond.md). New in v1.25.0; 4-species classifier + Yuan NMR plot in v1.27.0 candidate; generic mode + element fallback + PVA sample in v1.28.0 candidate.
- **`abmptools.crystal`** — organic-crystal FMO pipeline that produces ABINIT-MP AJF inputs and HPC jobscripts from one or more CIF files. 5-stage flow: (1) CIF → supercell PDB via legacy parser (`readcif.py`) or ASE backend (`cif_engine_ase.py` with `ase.io.read` + `Atoms.repeat` + bond-graph `unwrap_molecules` for PBC boundary fix), (2) PDB → `for_abmp/*.{ajf,pdb}` with fragment cut around solute (re-using `setfmo` + `pdb2fmo` with `is_xyz=True` for full-precision `&XYZ` block), (3) HPC jobscripts (PJM / SLURM / PBS / `local`, `string.Template` rendered + `runbatch.sh`), (4) optional `--run-local` invocation (`mpirun -np N abinitmp` for MPI-flat builds), (5) optional `getifiepieda` postprocessing (IFIE/PIEDA CSV + nearest-atom annotation). Driven by a single `CrystalBuildConfig` (one top-level dataclass + 6 leaves: `CIFInputSpec` / `CIFEngineConfig` / `FragmentTemplate` / `FMOMethod` / `HPCJobSpec` / `PostProcessSpec`) with YAML/JSON round-trip. CLI: `abmp-crystal {expand,fragment,jobs,pipeline,postproc,nearest,validate,example}` (8 subcommands, also `python -m abmptools.crystal`). ASE (LGPL-2.1+) and ABINIT-MP (separate distribution) are subprocess-only. Includes `abmptools/crystal/legacy/` namespace re-exporting the flat-layout 5 modules (`readcif`, `pdb2fmo`, `ajf2config`, `pdbmodify`, `getifiepieda`) for backward compatibility. See [`crystal.md`](crystal.md) + [`tutorial_crystal_fmo.md`](tutorial_crystal_fmo.md) + [`crystal_verification.md`](crystal_verification.md) + [`crystal_public_molecule_references.md`](crystal_public_molecule_references.md) (4-molecule MP2/6-31G(d) reference set). New in v1.23.0.
- **`abmptools.core`** — shared dataclasses used across the subpackages above:
  `SystemModel`, `AnnealProtocol`, `ClusterData`, `FixedLabel`, plus the
  `ensemble_family` flag and `classify_ensemble()` helper that distinguish
  COGNAC-only ensembles (e.g. `NPT_Andersen_Nose_Hoover`) from
  GROMACS-representable ones. The `amorphous` and `membrane` subpackages
  build a `SystemModel` from their respective inputs and then use
  `SimulationParams` (semiisotropic Pcoupl, multi-group thermostats, etc.)
  to drive the writer layer.

All subpackages use lazy imports for heavy scientific dependencies (OpenMM, RDKit, PySCF, etc.), so `import abmptools` remains cheap even if those optional stacks are not installed.

## Fortran Extension Integration

`abmptools/f90/src/readifiepiedalib.f90` (219 lines) provides a fast log file parser for IFIE/PIEDA data. It is:

1. Compiled to `readifiepiedalib.so` via `Makefile` (using gfortran with `-shared -fPIC`).
2. Loaded at runtime by `getifiepieda.py` using `ctypes`.
3. Optional — the `-nof90` flag forces pure-Python fallback.

**Note**: MP3/MP4 extraction requires the Fortran module. PB-IFIE, BSSE-IFIE, and monomer/dimer energies require the `-nof90` (pure Python) path.

## CLI Entry Points

All CLIs are invoked as `python -m abmptools.<module>`:

| Module | Purpose | Key Arguments |
|--------|---------|---------------|
| `generateajf` | Generate AJF templates | `-i` (PDB), `-basis`, `--method`, `-pb`, `-np` |
| `getifiepieda` | Extract IFIE/PIEDA | `--frag`, `--mol`, `--ffmatrix`, `--tfmatrix`, `-t`, `-np` |
| `log2cpf` | Convert log to CPF | `-i` (input), `-o` (output) |
| `log2config` | Convert log to config | `-i` (input), `-o` (output), `-np` |
| `pdb2fmo` | PDB to FMO setup | `-i` (PDB), `-p` (parameter file) |
| `generate_difie` | Create DIFIE CPF | `-i` (template), `-t` (range), `-z` (padding), `-np` |
| `convertcpf` | Convert CPF versions | `-i` (input), `-v` (version), `-f` (fragments) |
| `udf2fmo` | UDF to FMO setup | `-i` (UDF), `-p` (param), `-s` (solutes), `-r` (record) |
| `pdbmodify` | Modify PDB files | `-i`, `-move`, `-mode`, `-str` |
| `getcharge` | Extract charges | `-i` (log), `-t` (type), `-f` (fragments) |
| `addsolvfrag` | Add solvation | `-i` (PDB), `-temp` (template), `-solv` |
| `ajfserial` | Numbered AJF files | `-i` (template), `-t` (range), `-str` |
| `cpf2ifielist` | CPF to IFIE list | `-i` (CPF), `-f` (fragment range) |
| `crystal` (also `abmp-crystal`) | Organic-crystal FMO pipeline | subcommands: `expand` / `fragment` / `jobs` / `pipeline` / `postproc` / `nearest` / `validate` / `example`; `--config` (YAML/JSON), `--run-local` |

## Where to Start Reading

1. **`mol_io.py`** → `abinit_io.py` → `pdb_io.py` — Follow the inheritance chain to understand the I/O foundation.
2. **`cpfmanager.py`** — Understand how CPF data is parsed into DataFrames; this is the central data structure.
3. **`getifiepieda.py`** — See how analysis CLIs compose the base classes and data managers.
4. **`generateajf.py`** — Understand the FMO input generation pipeline.
5. **`sample/*/run.sh`** — Run working examples to see inputs and outputs.
