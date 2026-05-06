# Tutorial: GENESIS MM/GBSA via `abmptools.genesis.mmgbsa`

This tutorial walks through a full MM/GBSA run end-to-end:
environment setup → smoke pipeline → POC reproduction → analysis
review → failure modes. Targets abmptools v1.22.0 + GENESIS 2.1+.

For API + design notes, see [`mmgbsa.md`](./mmgbsa.md).

---

## 1. Environment

### 1.1 Conda env

```bash
mamba create -n mmgbsa -c conda-forge \
    python=3.11 numpy pandas matplotlib biopython acpype \
    ambertools=23 openmpi
mamba activate mmgbsa
```

`acpype` (GPL-3.0) and `ambertools` (academic+commercial) come from
conda-forge. `biopython` + `matplotlib` are pure-Python deps.

### 1.2 GENESIS

Build atdyn from source (LGPL-3.0+):

```bash
git clone --depth=1 https://github.com/genesis-release-r-ccs/genesis.git
cd genesis
./configure CC=gcc FC=gfortran
make -j8
make install
```

Verify:

```bash
which atdyn tleap acpype mpirun
atdyn --version 2>&1 | head -2
```

### 1.3 abmptools

```bash
cd /path/to/abmptools
python -m pip install -e .[mmgbsa]
```

---

## 2. Smoke pipeline (1 target, ~2-5 min)

### 2.1 Generate config

```bash
mkdir -p ~/mmgbsa_smoke && cd ~/mmgbsa_smoke
cp /path/to/abmptools/sample/gbsa/smoke_config.json ./
mkdir input
# Place a small protein-ligand complex into input/smoke.pdb
# (see tests/test_mmgbsa_integration.py for a synthetic example).
```

### 2.2 Validate

```bash
python -m abmptools.genesis.mmgbsa validate --config smoke_config.json
```

Expected:

```
Configuration: OK (project=mmgbsa_smoke, AMBER ff14SB + tip3p + GAFF2,
                   1 targets, implicit_solvent=GBSA)

External tools:
  [OK]    atdyn   (.../atdyn)
  [OK]    tleap   (.../tleap)
  [OK]    acpype  (.../acpype)
  [OK]    mpirun  (.../mpirun)

Python modules ([mmgbsa] extras):
  [OK]    Biopython   -- PDB splitter (Stage 1)
  [OK]    matplotlib  -- ΔG_bind bar plot (Stage 4)

Targets preview:
  T01: smoke.pdb  ligand_resno=201  (any chain)

Force field:
  protein  = leaprc.protein.ff14SB
  dna      = leaprc.DNA.OL15
  rna      = leaprc.RNA.OL3
  water    = leaprc.water.tip3p
  ligand   = leaprc.gaff2 + leaprc.gaff (acpype: BCC charges, atom_type=gaff2)

Protocol:
  energy   = GBSA (CUTOFF, cutoff=99.9 Å)
  minimize = SD × 1 step(s)
  mpi      = 1 rank(s) × atdyn (per system × 1 targets × 3 systems)
```

### 2.3 Run pipeline

```bash
python -m abmptools.genesis.mmgbsa pipeline \
    --config smoke_config.json -o ./out
```

Expected wall-time: 2–5 minutes (most spent in acpype's AM1-BCC).

Outputs:

```
out/
├── input/config.json
├── dirnames.in
├── smoke/
│   ├── smoke.pdb
│   ├── smoke_receptor_201.pdb
│   ├── smoke_ligand_201.pdb
│   ├── smoke_ligand_201.acpype/
│   ├── leaprc_{complex,ligand,receptor}
│   ├── {complex,ligand,receptor}.{prmtop,inpcrd,pdb,inp,log,dcd,rst}
│   └── leap_{complex,ligand,receptor}.log
└── analysis/
    ├── analysis_results.csv
    └── dg_bind_plot.png
```

---

## 3. Stage-by-stage workflow

For debugging or partial re-runs, invoke each stage individually:

```bash
# Stage 1 only
python -m abmptools.genesis.mmgbsa divide --config smoke_config.json -o ./out

# Stage 2 only (re-uses Stage 1 outputs)
python -m abmptools.genesis.mmgbsa parameterize --config smoke_config.json -o ./out

# Stage 3 only (re-uses Stages 1+2)
python -m abmptools.genesis.mmgbsa run --config smoke_config.json -o ./out

# Stage 4 only (any time after Stage 3)
python -m abmptools.genesis.mmgbsa analyze --config smoke_config.json --run-dir ./out
```

`parameterize` is idempotent: if `<target>/<ligand>.acpype/` already
contains the expected outputs, the slow AM1-BCC step is skipped. Set
`ligand.skip_if_cached=False` in the config to force re-parameterisation.

---

## 4. POC reproduction (4 targets, ~30-60 min)

The original POC ran 4 targets (`3772L` / `9MM52` / `L7759` / `M99KZ`)
at `ligand_resno=201`. Reproduce with:

```bash
mkdir -p ~/mmgbsa_poc && cd ~/mmgbsa_poc
cp /path/to/abmptools/sample/gbsa/example_config.json ./
mkdir input
# Stage the 4 PDBs into input/ (research data, not committed).
```

Then either:

```bash
# Folder-mode shortcut (POC compatibility):
python -m abmptools.genesis.mmgbsa pipeline -i ./input -r 201 -o ./out

# Or config-mode (re-runnable):
python -m abmptools.genesis.mmgbsa pipeline --config example_config.json -o ./out
```

The 3772L target reproduces the gold value:

```
ΔG_bind = -38.27 kcal/mol  (matches POC analysis_results.csv,
                            single-frame MM/GBSA single-point)
```

Per GENESIS doc 05_Energy.rst:564, the `ENERGY` column is the total
potential `U = U_FF + ΔG_solv`, so `ΔG_bind = E_complex - E_ligand -
E_receptor`. POC `4_analyse.py` writes the algebraically identical
decomposed form `(egas+S)_c - (egas+S)_l - (egas+S)_r`; both yield
the same number. See `mmgbsa.md` §6.3.

---

## 5. Interpreting the output

### 5.1 `analysis_results.csv`

Columns:
```
Target, Complex_Energy, Complex_Solvation,
        Ligand_Energy,  Ligand_Solvation,
        Receptor_Energy, Receptor_Solvation,
        dG_bind_kcal_per_mol
```

Energies are GENESIS-native kcal/mol (matches AMBER convention).
`dG_bind_kcal_per_mol` is the standard MM/GBSA estimator; **single
point**, no entropy correction, no ensemble averaging.

### 5.2 `dg_bind_plot.png`

A simple bar plot for visual comparison across targets. For
publication-quality figures, regenerate with custom styling using
the CSV.

### 5.3 Sanity checks

A successful run typically gives:

- Each target's `complex.log` shows `[STEP4] Compute Single Point Energy`
  with finite numbers (no NaN, no extreme values).
- `analysis_results.csv` ΔG_bind is in the kcal/mol range
  (expect ~±200 for plausible binding poses; very large or
  very negative may indicate clashes or mis-parameterisation).
- The bar plot shows comparable ΔG values across closely-related
  targets (sanity check for systematic errors in setup).

---

## 6. Failure modes

| Symptom | Cause | Fix |
|---|---|---|
| `acpype: failed (rc=1)` | Bad ligand H placement, missing CONECT, exotic atoms | Pre-clean with `pdb4amber`; check ligand element column |
| `tleap: failed (rc=1)` | Missing residue templates, mismatched chain IDs | Inspect `leap_{complex,ligand,receptor}.log` for the offending residue |
| `atdyn: failed` (NaN energy) | Bad starting structure; ligand-receptor clash | Pre-minimise with OpenMM or AMBER `pmemd` before passing to GENESIS |
| `[STEP4]` not found | atdyn quit early (segfault / OOM) | Re-run with `-v` to capture stderr; check log size |
| Empty `param1.dcd` | atdyn crashed silently | Inspect `<name>.log` tail |
| Same `ΔG_bind` for many targets | Re-used same fixture / cached outputs | Set `ligand.skip_if_cached=False` |
| ΔG values significantly differ from POC `analysis_results.csv` column-diff | Possible setup divergence (different ligand_resno, chain, ff) | Compare per-system ENERGY values to POC logs; v1.22.0 should match POC bit-for-bit when given identical inputs |

---

## 7. Where to look next

- [`docs/mmgbsa.md`](./mmgbsa.md) — package reference (API, defaults, design)
- [`tests/test_mmgbsa_*.py`](../tests/) — 116 unit tests + 1 slow integration
- [`sample/gbsa/`](../sample/gbsa/) — example_config.json + smoke_config.json
- POC scripts at `sample/gbsa/legacy_scripts/`
- POC outputs at `~/repos/abmptools-sample/sample/gbsa-genesis/`
- Genheden & Ryde (2015) MM/GBSA review:
  <https://doi.org/10.1517/17460441.2015.1032936>
