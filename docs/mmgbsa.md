# `abmptools.genesis.mmgbsa` — GENESIS MM/GBSA ΔG_bind Builder

This document is the subsystem reference for the
`abmptools.genesis.mmgbsa` package, introduced in v1.22.0.
It provides a fully-typed pipeline around GENESIS' `atdyn` for
single-point MM/GBSA binding free-energy calculations.

For step-by-step ops see [`tutorial_mmgbsa.md`](./tutorial_mmgbsa.md).

---

## 1. Scope and license

**What it does**: takes a protein-ligand complex PDB + a ligand
residue number, splits it into receptor + ligand PDBs (Biopython),
parameterizes each system with AmberTools (`tleap`) + acpype
(GAFF/GAFF2 + AM1-BCC), runs `atdyn` with `[ENERGY]
implicit_solvent=GBSA` for single-point energies on `complex` /
`ligand` / `receptor`, and computes

```
ΔG_bind = E_complex - E_ligand - E_receptor       [kcal/mol]
```

where `E` is the GENESIS `ENERGY` column from the `[STEP4] Compute
Single Point Energy` row. Per GENESIS doc 05_Energy.rst:564, the
`ENERGY` column is the **total potential** `U = U_FF + ΔG_solv`, so the
`SOLVATION` column is already a sub-component and **must not be added
on top** of `E` (would double-count solvation).

The result is a CSV + bar plot.

**License posture**: GENESIS is **LGPL-3.0-or-later**, acpype is
**GPL-3.0**. abmptools shells out to both via `subprocess` only — no
source modification, dynamic linking, or bundling. This is *mere
aggregation* per the FSF GPL FAQ and LGPL §5/§6, compatible with
abmptools' **Apache-2.0** (v1.23.0+; ≤ v1.22.0 was MIT).

AmberTools (`tleap`) is free for academic + commercial use.
Biopython is BSD-3-Clause + Biopython License (dual).

---

## 2. Quick start

```bash
cd abmptools/
python -m pip install -e .[mmgbsa]
# Dependencies installed: biopython, matplotlib.
# External tools (acpype / tleap / atdyn / mpirun): conda or build manually.

# Folder-mode shortcut (POC compatibility):
python -m abmptools.genesis.mmgbsa pipeline \
    -i ./sample/gbsa/input -r 201 -o ./out

# Config-mode (recommended):
python -m abmptools.genesis.mmgbsa pipeline \
    --config ./sample/gbsa/example_config.json -o ./out

# Re-aggregate ΔG_bind from existing logs (no recomputation):
python -m abmptools.genesis.mmgbsa analyze \
    --target-dirs ./out/3772L-rename ./out/9MM52-rename \
    --out-dir ./out/analysis_v2
```

---

## 3. Pipeline

### 3.1 Build (4 stages)

```
1. divide        Biopython          -> <basename>_{receptor,ligand}_<resno>.pdb
2. parameterize  acpype + tleap     -> {complex,ligand,receptor}.{prmtop,inpcrd,pdb}
3. run_gbsa      mpirun -np 1 atdyn -> {complex,ligand,receptor}.{inp,log,dcd,rst}
4. analyze       matplotlib + csv    -> analysis_results.csv + dg_bind_plot.png
```

Each stage is independently re-runnable via the matching CLI
sub-command, so partial failures can be resumed without redoing the
slow steps.

### 3.2 Output layout

```
output_dir/
├── input/
│   └── config.json                                # frozen config dump
├── dirnames.in                                    # POC compatibility
├── <target>/                                      # one per TargetSpec
│   ├── <target>.pdb                               # original copy
│   ├── <target>_receptor_<resno>.pdb              # Stage 1
│   ├── <target>_ligand_<resno>.pdb                # Stage 1
│   ├── <target>_ligand_<resno>.acpype/            # Stage 2 (acpype)
│   │   ├── <target>_ligand_<resno>_AC.frcmod
│   │   └── <target>_ligand_<resno>_bcc_gaff2.mol2
│   ├── leaprc_{complex,ligand,receptor}           # Stage 2 (tleap inputs)
│   ├── leap_{complex,ligand,receptor}.log         # Stage 2 (tleap logs)
│   ├── {complex,ligand,receptor}.{prmtop,inpcrd,pdb}  # Stage 2 outputs
│   ├── {complex,ligand,receptor}.inp              # Stage 3 (atdyn inputs)
│   └── {complex,ligand,receptor}.{log,dcd,rst}    # Stage 3 outputs
└── analysis/
    ├── analysis_results.csv                       # 1 row per target
    └── dg_bind_plot.png                           # bar plot
```

---

## 4. Public API

```python
from abmptools.genesis.mmgbsa.models import (
    TargetSpec,
    ForceFieldSet,
    LigandParameterization,
    EnergyProtocol,
    MinimizationProtocol,
    MMGBSABuildConfig,
)
from abmptools.genesis.mmgbsa.builder import MMGBSAOrchestrator

cfg = MMGBSABuildConfig(
    targets=[TargetSpec(pdb="3772L.pdb", ligand_resno=201)],
    input_dir="./input",
)
orch = MMGBSAOrchestrator(cfg)

# All 4 stages.
result = orch.run()

# Or stage-by-stage.
splits  = orch.divide()
params  = orch.parameterize()
runs    = orch.run_gbsa()
analyze = orch.analyze()
```

`run()` returns a dict with keys `output_dir`, `config_json`,
`n_targets`, `n_succeeded`, `n_failed`, `failures`, `splits`,
`csv_path`, `png_path`, and `targets` (list of per-target ΔG_bind).

---

## 5. CLI

```
python -m abmptools.genesis.mmgbsa {example, validate, divide, parameterize, run, analyze, pipeline}
```

| Sub-command | Required flags | Behaviour |
|---|---|---|
| `example` | (none) | print default JSON (4 POC targets at ligand_resno=201) |
| `validate` | `--config` | tools + Python deps + targets preview, exits 0/1 |
| `divide` | `--config` (`-o`) | Stage 1 only |
| `parameterize` | `--config` (`-o`) | Stages 1+2 (idempotent) |
| `run` | `--config` (`-o`) | Stages 1+2+3 |
| `analyze` | `--config`+`--run-dir` OR `--target-dirs` | Stage 4 only (replay) |
| `pipeline` | `--config` OR (`-i`+`-r` [`-c`]) | All 4 stages |

The `pipeline` folder-mode shortcut (`-i input_dir -r ligand_resno
[-c chain]`) preserves POC compatibility: it enumerates `*.pdb` in
the input dir and synthesises `TargetSpec` instances on the fly.

---

## 6. Design rationale

### 6.1 Why subprocess (no Python binding) for atdyn / tleap / acpype?

All three are licensed under copyleft (LGPL-3.0+ for GENESIS, GPL-3.0
for acpype) or have non-redistributable installers. Subprocess
invocation is the simplest "mere aggregation" path that keeps
abmptools **Apache-2.0** (v1.23.0+) compatible.

### 6.2 Why `ff14SB`, not `ff19SB`?

POC and tutorial 12.3 use `ff14SB` for MM/GBSA. Existing literature
on GBSA single-point binding scoring is overwhelmingly ff14SB-based
(GBSA force-field re-parameterisations are tied to it). `genesis.grest`
uses `ff19SB` for newer protein REMD; we keep them deliberately
separate.

### 6.3 ΔG_bind formula (GENESIS ENERGY column convention)

Per GENESIS doc `05_Energy.rst:564`:

> The solvation free energy is incorporated into the molecular
> mechanics potential energy function as an effective energy term,
> namely, **U = U_FF + ΔG_elec + ΔG_np**

That is, the `ENERGY` column already contains the SOLVATION
contribution. So the binding free energy is simply the difference
of `ENERGY` columns:

```
ΔG_bind = E_complex - E_ligand - E_receptor
```

POC `4_analyse.py` uses the equivalent decomposed form:

```
ΔG_bind = (egas + S)_complex - (egas + S)_ligand - (egas + S)_receptor
        where egas = ENERGY - SOLVATION
```

Since `egas + S = ENERGY`, both forms are **algebraically identical**.
v1.22.0 ships `compute_dg_bind` (total form) plus `compute_dg_components`
which returns `{"dg_mm": ΔE_MM, "dg_solv": ΔG_solv, "dg_bind": ΔG}` for
users who want the gas-phase / solvation breakdown.

3772L gold value: **−38.27 kcal/mol** (matches the POC `analysis_results.csv`
column-difference).

### 6.4 Why 4 sub-commands + `pipeline`?

acpype runs are slow (AM1-BCC charge calculation can take minutes per
ligand). Splitting `divide` / `parameterize` / `run` / `analyze`
allows re-running a single failed stage without redoing the previous
ones. `parameterize` re-uses cached `<basename>.acpype/` outputs by
default (`LigandParameterization.skip_if_cached=True`).

### 6.5 Why `mpirun -np 1`?

The single-point GBSA evaluation on a small protein-ligand complex
doesn't benefit from spatial decomposition; atdyn MPI rank > 1 can
even produce decomposition artefacts on tiny systems. Override via
`mpi_processes` for very large complexes.

### 6.6 Why store `dirnames.in`?

POC compatibility. The legacy scripts in `sample/gbsa/legacy_scripts/`
read `dirnames.in` to enumerate per-target dirs. Emitting it lets
users mix abmptools.genesis.mmgbsa output with POC scripts during
migration.

---

## 7. Defaults (selected)

| Field | Default | Rationale |
|---|---|---|
| `force_field.leaprc_protein` | `leaprc.protein.ff14SB` | POC |
| `force_field.leaprc_water` | `leaprc.water.tip3p` | POC |
| `ligand.charge_method` | `bcc` | POC (AM1-BCC, gold standard for GBSA) |
| `ligand.atom_type` | `gaff2` | POC |
| `ligand.extra_keys` | `maxcyc=0` | POC (skip antechamber pre-opt) |
| `energy.electrostatic` | `CUTOFF` | POC |
| `energy.cutoffdist_A` | 99.9 | POC (NOBC = effectively no cutoff) |
| `energy.implicit_solvent` | `GBSA` | POC |
| `energy.gbsa_salt_cons` | 0.15 | POC (≈ physiological NaCl) |
| `energy.gbsa_surf_tens` | 0.0072 kcal/(mol·Å²) | POC (Onufriev) |
| `energy.gbsa_vdw_offset` | 0.09 Å | POC |
| `minimize.method` | `SD` | POC (single-point) |
| `minimize.nsteps` | 1 | POC (single-point) |
| `mpi_processes` | 1 | POC |

---

## 8. Caveats and known issues

- **acpype on heteroatoms with H-deficient PDBs**: input PDB must have
  reasonable hydrogen placement. acpype + antechamber will hallucinate
  H positions if missing, sometimes producing wrong tautomers. Pre-clean
  with `pdb4amber` or RDKit before splitting.
- **Net charge auto-detection**: `LigandParameterization.net_charge=None`
  lets acpype infer charge from the structure. For unusual ligands
  (zwitterions, salts), set `net_charge` explicitly.
- **Salt concentration**: only enters via `gbsa_salt_cons` (M); to
  introduce explicit ions for an MM-only run, set
  `implicit_solvent=NONE` and add ions to the PDB pre-split.
- **`forcefield = CHARMM`** would require the CHARMM MM/GBSA
  implementation (different `[ENERGY]` block); v1.22.0 ships AMBER only.
- **POC `4_analyse.py` egas formula**: algebraically identical to
  v1.22.0's total-energy form (since `egas + SOLVATION = ENERGY`).
  Both yield the same ΔG_bind. The legacy script computes the
  decomposed form and reports `egas` separately, which can be useful
  for inspecting the MM vs solvation contributions; v1.22.0 exposes
  the same decomposition via :func:`compute_dg_components`.

---

## 9. References

- POC scripts: `sample/gbsa/legacy_scripts/`
- POC outputs: `~/repos/abmptools-sample/sample/gbsa-genesis/`
- GENESIS: <https://github.com/genesis-release-r-ccs/genesis>
- acpype: <https://github.com/alanwilter/acpype>
- AmberTools: <https://ambermd.org/AmberTools.php>
- MM/GBSA reference: Genheden & Ryde, *Expert Opin. Drug Discov.* 10(5):449
  (2015), <https://doi.org/10.1517/17460441.2015.1032936>
