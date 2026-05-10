# Tutorial: Organic-Crystal FMO with `abmptools.crystal`

End-to-end recipe for taking an organic-crystal CIF and producing
ABINIT-MP FMO calculation inputs, jobscripts, and (optionally)
running a smoke calculation locally.

For the API/design reference see [`crystal.md`](./crystal.md).

---

## 1. Environment

```bash
mamba activate abmptoolsenv      # or whatever env hosts your editable install
python -c "import abmptools.crystal; print('OK')"
```

Required Python packages (installed via `pip install abmptools[crystal]`
or `mamba install`):

- `ase >= 3.22` (LGPL-2.1+) — needed for `engine='ase'`
- `pyyaml >= 6.0` — needed for `--config *.yaml` (JSON also works)

Required external tools:

- ABINIT-MP `abinitmp` binary — only when running `--run-local` or
  the rendered HPC jobscripts. Not needed for the input-generation
  stages.

The csp7 sample is bundled at `sample/crystal/csp7_smoke/`:

```
sample/crystal/csp7_smoke/
├── XXXI-MMFF-R00001.cif      ← input crystal (P21/n, Z=4)
├── UNK.ajf                   ← AJF template (any minimal FMO ajf works)
├── crystal.yaml              ← Phase C+ recommended config
├── input_param / segment_data.dat   ← legacy CLI drivers (Phase A/B form)
└── run.sh                    ← legacy 4-step shell recipe (kept for documentation)
```

Copy this directory anywhere writable and use it as the working
directory for the rest of the tutorial.

---

## 2. Smoke run: 1 CIF, 1 command, 5 stages

```bash
cp -r sample/crystal/csp7_smoke /tmp/csp7_demo
cd /tmp/csp7_demo
python -m abmptools.crystal pipeline --config crystal.yaml
```

Expected runtime: ~2 seconds on a modest workstation.

Output (under `./out/XXXI-MMFF-R00001/cifout/layer5/pdb/for_abmp/`):

- `XXXI-MMFF-R00001layer5Zp1-around_ar6.0.ajf` — FMO calculation input
  with `Natom=832` and a full-precision `&XYZ` block
- `XXXI-MMFF-R00001layer5Zp1-around_ar6.0.pdb` — supercell PDB
  (832 atoms, 26 fragments)
- `XXXI-MMFF-R00001layer5Zp1-around_ar6.0_12n-2p-24t.sh` — rendered
  PJM jobscript (12 nodes × 2 procs, OMP=24, queue=small)
- `runbatch.sh` — submitter that runs `pjsub *.sh` for the whole AJF set

The output is byte-equivalent to the v1.22.0 baseline frozen in
`tests/regression/reference/main/crystal_csp7/R00001/`.

---

## 3. csp7 reproduction (3 structures)

Reproduce the Phase B regression fixture by replacing the single-input
config with three CIFs:

```yaml
inputs:
  - cif: XXXI-MMFF-R00001.cif    # P21/n, Z=4 → 832 atoms / for_abmp
    layer: 5
    atoms_in_mol: [32]
  - cif: XXXI-MMFF-R00002.cif    # Z=8 → 520 atoms / for_abmp
    layer: 5
    atoms_in_mol: [32]
  - cif: XXXI-MMFF-R00004.cif    # P-1, Z=2 → 520 atoms / for_abmp
    layer: 5
    atoms_in_mol: [32]
```

Place the three `.cif` files alongside `crystal.yaml` and run
`python -m abmptools.crystal pipeline --config crystal.yaml`. Each
input gets its own `out/<input_name>/` subdirectory; the for_abmp
files match the fixture under
`tests/regression/reference/main/crystal_csp7/{R00001,R00002,R00004}/for_abmp/`.

---

## 4. ASE engine: arbitrary space groups

The default `engine='legacy'` covers the symmetry operators and layer
counts validated for the csp7 P21/n / P-1 systems. Switch to
`engine='ase'` when:

- the CIF uses a space group whose `getsymcoord` lookup fails
  silently (e.g. unusual trigonal denominators)
- you want layer ≥ 6 (the legacy parser only knows 1/4/5)

```bash
python -m abmptools.crystal pipeline --config crystal.yaml --engine ase
```

ASE reads `_symmetry_equiv_pos_as_xyz` directly so all standard
space groups work. Coordinates emerge with cell-origin layout (vs.
the legacy centred layout); molecule order also differs. The AJF
remains valid (`Natom > 0`, `&XYZ` block size matches `Natom`), but
**byte-equivalence with the legacy fixture is not preserved** — this
is by design.

---

## 5. HPC jobs

The `hpc` block of the config drives jobscript rendering:

```yaml
hpc:
  scheduler: PJM      # or SLURM / PBS / local
  queue: small
  group: hp190133
  nodes: 12
  proc_per_node: 2
  omp_threads: 24
  elapse: '24:00:00'
  abinit_dir: /data/hp190133/programs/ABINIT-MP
  binary_name: abinitmp_smp
```

After running `pipeline` (or `jobs`), submit the entire AJF set with
the bundled `runbatch.sh`:

```bash
cd out/XXXI-MMFF-R00001/cifout/layer5/pdb/for_abmp
bash runbatch.sh   # invokes pjsub on each rendered .sh
```

For SLURM / PBS sites, change `scheduler:` and `binary_name:` and
`runbatch.sh` will use `sbatch` / `qsub` automatically.

---

## 6. Local execution (`--run-local`)

For small smoke runs (one or a handful of structures) where you have
`abinitmp` on PATH:

```bash
python -m abmptools.crystal pipeline --config crystal.yaml --run-local
```

The orchestrator pipes each AJF into `abinitmp` directly:

```
abinitmp < <ajf>     >  <log>     2>  <err>
```

Logs land in the same `for_abmp/` directory as `<base>.log` /
`<base>.err`. `--run-local` is strictly for smoke / debug; for
production-scale (csp7 1500 structures, MP2/6-31G*), use HPC submission.

---

## 7. Postprocessing (IFIE/PIEDA + nearest atoms)

Once the FMO logs are available, run the standalone postproc
subcommand:

```bash
python -m abmptools.crystal postproc \
  --config crystal.yaml \
  --logs out/XXXI-MMFF-R00001/cifout/layer5/pdb/for_abmp/*.log \
  --pdb  out/XXXI-MMFF-R00001/cifout/layer5/pdb/for_abmp/XXXI-MMFF-R00001layer5Zp1-around_ar6.0.pdb \
  --z-prime 1
```

Output (under `out/postproc/`):

- `csv/*.csv` — IFIE/PIEDA tables emitted by `getifiepieda`
- `csv/*.csv` (annotated, in place) — appended columns
  `nearest_element_1`, `nearest_distance_1`, ... when
  `postproc.annotate_nearest_atoms=true` and a `--pdb` was passed

---

## 8. Establishing a numeric reference for your crystal

Once `--run-local` works for a smoke calculation (Section 6), the
natural next step is to capture the actual FMO numbers (Total
energy, monomer energies, IFIE table) as a **numeric reference**
that future refactors can be regression-tested against. This
section walks through the recipe that produced the csp7 R00001
layer3 HF/6-31G snapshot under
`tests/regression/reference/main/crystal_csp7/R00001/`.

### 8-1. Estimate the system size before committing wall time

Run the input-only stages first (`pipeline` without `--run-local`)
and read off `Natom` and the dimer pair counts from the rendered
AJF — this lets you predict cost without paying for it.

```bash
python -m abmptools.crystal pipeline --config crystal.yaml
AJF=$(ls out/<input_name>/cifout/layer<L>/pdb/for_abmp/*.ajf | head -1)
grep -E "^Natom|^Method|^BasisSet" "$AJF"
grep -A 4 "FRAGMENT PAIR DISTANCE" out/.../*.log 2>/dev/null \
    || awk '/Natom=/ {n=$1} /XYZ/ {p=1} END{print "atoms:", n}' "$AJF"
```

The orchestrator also prints `nfrag` indirectly via the for_abmp
PDB residue count: `awk '$1=="HETATM"{print $5}' <pdb> | sort -u | wc -l`.

### 8-2. Predict wall time

abinitmp v2r8 wall time scales roughly as `nfrag` (monomer SCF) +
`nfrag^2 / 2 × t_dimer` (dimer SCF), and `t_dimer` scales as
`basis_per_atom^4`. Empirical anchors on a single core
(WSL2 / `/home/okuwaki/bin/abinitmp/v2r8/abinitmp`):

| Layer | Method | Fragments | Pairs | Wall time (1 core) |
|---|---|---:|---:|---:|
| 2 | HF / STO-3G | 32 | 496 | ~30 s |
| 3 | HF / 6-31G | 24 | 276 | ~9 h 14 min |

**Choosing the binary** (see memory `reference_abinitmp_parallelism`):

| Binary name | Parallelism | How `--run-local` invokes it |
|---|---|---|
| `abinitmp` | MPI flat (1 thread per rank) | `mpirun -np <proc_per_node> abinitmp` (or direct, with `mpi_launcher: ""`) |
| `abinitmp_omp` | MPI + OpenMP hybrid | `mpirun -np <proc_per_node> abinitmp_omp` + `OMP_NUM_THREADS=<omp_threads>` |

`OMP_NUM_THREADS` only takes effect with `abinitmp_omp`. Setting it on
plain `abinitmp` is silently ignored — the layer3 / 6-31G run that
took ~9 h on 1 core in the table above was the result of using flat
`abinitmp` with `proc_per_node=1`. For local reference calculations
with multiple cores, switch to `abinitmp_omp` or raise
`proc_per_node` (and have an MPI runtime on PATH). The AJF
`Memory` field is per-MPI-process megabytes.

### 8-3. Run with a persistent work dir

The bundled `tests/integration/run_crystal_smoke.sh` writes to
`mktemp -d`, so its log disappears after the run. To keep a
reference run, write to a fixed `.smoke_artifacts/` (which the
repo `.gitignore` already excludes):

```bash
WORK=.smoke_artifacts/<system>_l<L>_<method>_<basis>
mkdir -p "$WORK"
cp sample/crystal/csp7_smoke/UNK.ajf "$WORK/"      # or your own template
cp <your_input>.cif "$WORK/"
cat > "$WORK/crystal.yaml" <<YAML
project_name: <system>_l<L>_<method>_<basis>
output_dir: ./out
inputs:
  - {cif: <your_input>.cif, layer: <L>, atoms_in_mol: [<N>]}
cif_engine: {engine: legacy}
fragment:
  cutmode: around
  solutes: [0]
  criteria: 6.0
  molname: [UNK]
  template_ajf: ./UNK.ajf
fmo: {method: <HF|MP2>, basis_set: <STO-3G|6-31G|...>, memory: 4000, abinit_ver: rev23, npro: 1, is_xyz: true}
hpc:
  scheduler: local
  abinit_dir: /home/okuwaki/bin/abinitmp/v2r8
  binary_name: abinitmp           # or abinitmp_omp for hybrid
  nodes: 1
  proc_per_node: 1                # MPI rank count for --run-local
  omp_threads: 1                  # OMP only effective with abinitmp_omp
  mpi_launcher: ""                # "" = direct, "mpirun" / "mpiexec" / "srun" etc.
  elapse: '12:00:00'
postproc: {enable: false}
YAML
cd "$WORK"
python -m abmptools.crystal pipeline --config crystal.yaml --run-local
```

The log lands at `out/<system>/cifout/layer<L>/pdb/for_abmp/*.log`
and persists for extraction.

### 8-4. Extract numbers via `getifiepieda` (CSV reference)

Phase D-3 (revised, 2026-05-09) replaces the earlier ad-hoc Python
regex extractor with a thin wrapper around the in-tree
`abmptools.getifiepieda` post-processor — single source of truth,
shared with production csp7 1500-structure runs. The shipped
extract script is therefore a `subprocess.run` shell:

```bash
LOG=$(ls .smoke_artifacts/<your_run>/out/.../*.log | head -1)
python tests/regression/reference/main/crystal_csp7/R00001/extract_layer3_hf_631g.py \
    "$LOG" \
    --out-dir tests/regression/reference/main/crystal_<your_system>/<run_id>
```

It runs (under the hood)

```bash
python -m abmptools.getifiepieda \
    --multi 1 -dimeres -imd -zp 5 -t <id> <id> 1 -nof90 \
    -i '["<dir>/<prefix-without-id>", "<suffix-with-extension>"]'
```

and copies the canonical outputs to your `--out-dir`:

| File | Source CSV | Contents |
|---|---|---|
| `expected_<run_id>_ifiesum.csv` | `csv/frag1-dimer-es-false-ifiesum.csv` | 1 row: target-frag-1 sums (HF-IFIE / ES / EX / CT-mix / **`MonomerEnergy(1)`** / etc) |
| `expected_<run_id>_ifiedt.csv`  | `csv/frag1-dimer-es-false-ifiedt.csv`  | N rows: per-pair detail with `UNK<i>(<i>)` resname annotation |

The `MonomerEnergy(<frag>)` column is the **in-crystal monomer
energy** of the target fragment: the polarized monomer SCF/MP2
total (HF + MP2 correlation, hartree) computed inside the
supercell's electrostatic field. Combine it with the
gas-phase isolated-monomer total (Section 8-5, same fragment
coordinates) to extract a polarization-plus-correlation
contribution. Note that this difference is **not** a structural
deformation energy (the geometry is unchanged between the two
calculations); a clean structural-deformation analysis needs a
separate gas-phase relaxed-geometry calculation, e.g. a
polymorph-vs-polymorph comparison of `MonomerEnergy(<frag>)`.

**Adapting to a different system**: copy the extract script as-is,
then edit only the destination filename hardcoded in `main()` (and
`split_log_for_getifiepieda` if your structure-id zero-padding
differs from 5 digits). The CLI flags passed to getifiepieda are
fixed for the layer3-style run; for a different cutmode / target
fragment, change `--multi` and `-dimeres`/`-imd` accordingly.

### 8-5. Isolated monomer reference (gas-phase)

For the gas-phase monomer reference you need a **separate ajf** that
contains only fragment 0 (the central molecule) with no surrounding
neighbours. The shipped public-molecule samples
(`sample/crystal/{urea,glycine,benzene,naphthalene}/monomer.yaml`)
do this with `cutmode='around', criteria=0.0, layer=1`: the
fragment-0 atoms are kept and every other molecule is screened
out by the 0 Å cutoff, leaving a 1-fragment AJF whose atom
coordinates match those of frag 1 in the supercell exactly.

Compare the resulting log's `Total energy` to the
`MonomerEnergy(<frag>)` column of `ifiesum.csv`. The difference is
a polarization-plus-correlation shift between the in-crystal and
isolated monomer at the **same geometry** — it is **not** a
structural deformation energy. A structural deformation analysis
needs separate gas-phase relaxation, e.g. a polymorph-vs-polymorph
comparison of relaxed monomer energies.

### 8-6. Add the regression test

`tests/test_crystal_numeric_regression.py` is the template. The cheap
test (`test_expected_csvs_have_committed_shape`) sanity-checks the
`test_expected_json_shape`) generalise by parametrising on the
`(extract_script, excerpt_log, expected_json)` triple — extend
`pytest.mark.parametrize` with the new system instead of writing
a new test file. Copy the live regression test
(`test_live_*_against_reference`) only if the wall time is
manageable and the abinitmp build is reproducible.

committed CSVs (column presence, MP2-padding-with-zero, monomer
energy range, dimer-es=False filter). The slow live test
(`test_live_layer3_hf_631g_against_reference`) is gated by
`ENABLE_FMO_LIVE_REGRESSION=1` and reproduces the full chain
(`pipeline --run-local` → 9 h abinitmp → extract script → byte-
compare CSVs) — useful for verifying that an `anlfmo.py` change
hasn't drifted the numbers.

```bash
pytest tests/test_crystal_numeric_regression.py -v
# expect: test_expected_csvs_have_committed_shape PASSED
#         test_live_layer3_hf_631g_against_reference SKIPPED (gated)

# Manual full reproduction (~9 h):
ENABLE_FMO_LIVE_REGRESSION=1 \
  pytest tests/test_crystal_numeric_regression.py::test_live_layer3_hf_631g_against_reference -v
```

### 8-7. Stash the full artifacts off-repo

Per the abmptools artifact-placement convention, the full
`.ajf / .pdb / .log / .err / crystal.yaml / run_local.sh` set
goes to OneDrive (or wherever your team's `abmptools-dump/`
lives), not the repository:

```bash
DUMP=/mnt/c/Users/<you>/OneDrive/abmptools-dump/crystal/<run_id>
mkdir -p "$DUMP"
cp .smoke_artifacts/<your_run>/out/.../*.{ajf,pdb,log,err} "$DUMP/"
cp .smoke_artifacts/<your_run>/{crystal.yaml,run_local.sh,run.stdout} "$DUMP/"
```

The repository keeps only the two reference CSVs (~2 KB) and the
extract script wrapper (~3 KB). `.smoke_artifacts/` itself is
`.gitignore`'d at the repo root, so the multi-MB log stays local.

The legacy `expected_layer3_hf_631g.json` and
`excerpt_layer3_hf_631g.log` from the original Phase D-3 (regex
extractor) remain in the reference directory for archival but are
not read by the test suite anymore.

---

## 9. Failure modes

| Symptom | Likely cause | Fix |
|---|---|---|
| `pdb2fmo` raises with "no atoms in mol" | Wrong `atoms_in_mol` for the CIF | Inspect the CIF; set `atoms_in_mol=[N1]` (or `[N1, N2]` for Z'≥2) |
| `detect_molecules` raises with "do not match atoms_in_mol" (`engine='ase'`) | `bond_tolerance` too tight or wrong `atoms_in_mol` | Increase `cif_engine.bond_tolerance` (default 0.4 Å) or fix `atoms_in_mol` |
| `Natom=0` and `&XYZ` block missing | `is_xyz=true` in config but matching `<base>.xyz` not staged | The orchestrator stages it automatically; if calling `pdb2fmo` by hand, copy `cifout/layer<L>/xyz/<base>.xyz` next to the PDB |
| `--run-local` fails with `FileNotFoundError: abinitmp binary not found` | `abinitmp` not on PATH and `abinit_dir`/`binary_name` not set or wrong | Set both fields in the config or `export PATH=...:$PATH` |
| `--run-local` runs single-core despite `omp_threads: 4` | Plain `abinitmp` is the MPI flat build — it ignores `OMP_NUM_THREADS`. | Switch `binary_name: abinitmp_omp` (the hybrid build) or raise `proc_per_node` so multiple MPI ranks are launched via `mpi_launcher`. |
| `--run-local` fails with `mpirun: command not found` (or similar) | `mpi_launcher: "mpirun"` (default) but no MPI runtime on PATH | Either set `mpi_launcher: ""` for direct flat invocation (only valid with `proc_per_node=1` and a non-OMP binary), or install/load MPI (`module load openmpi` / activate a conda env that ships `mpirun`). |
| `_build_run_command` raises `ValueError: mpi_launcher must be non-empty when ...` | Tried to launch a multi-rank or `_omp` build without a launcher | Set `mpi_launcher: "mpirun"` (or your site's launcher). Direct invocation is only OK for single-rank flat. |
| Live FMO calc takes >> predicted wall time | Almost certainly running on 1 core | See Section 8-2; switch to MPI or downgrade method/basis |
| Legacy fixture diff after refactor | Phase B regression test catches this | `pytest tests/test_crystal_regression.py -m slow` shows the offending file pair |
| `extract_*.py` regex misses the FMO TOTAL ENERGY block | (legacy regex extractor) Anchored at `## TIME PROFILE` and matched the early MONOMER SCC profile instead | Use the new `extract_*.py` shipped in Phase D-3 revised (calls `getifiepieda` directly) — no regex maintenance |
| `getifiepieda` raises `IndexError: list index out of range` at `anlfmo.py:649` (`getlogorpdbfrag`) | Old anlfmo: `ReadGeom = ` was empty (`is_xyz=True` route) and the parser blindly indexed `items[2]` | Fixed in 2026-05-09 (`anlfmo.py:648-660` falls back to sibling PDB basename). Update `abmptools` to ≥ that commit |
| `getifiepieda` raises `ValueError: 10 columns passed, passed data had 7 columns` | Old anlfmo: `readmultiifie` didn't set `self.logMethod` / `self.icolumn` for HF logs | Fixed 2026-05-09 — `readmultiifie` now mirrors `readsingleifie`'s method-aware setup |
| `getifiepieda` raises `KeyError: 'MP2-IFIE'` (or `'PR-TYPE1'`/`'GRIMME'`/...) | Old anlfmo: `getfiltifpifd` indexed MP2 columns that HF logs don't have | Fixed 2026-05-09 — `getfiltifpifd` now zero-pads MP2 sums when `self.logMethod == 'HF'` |
| `getifiepieda --imd` raises `ValueError: 3 columns passed, passed data had 2 columns` | Old anlfmo: `getmomenedf`/`getdimenedf` assumed 3-column / 4-column rows (MP2 layout) | Fixed 2026-05-09 — both DataFrame builders branch on `self.logMethod` for HF (2 / 3 cols, MP2 padded with 0) |
| `--engine ase` Monomer SCC oscillates by 100s of hartree, never converges | Pre-2026-05-10 ASE backend: `Atoms.repeat` carried PBC forward but didn't unwrap molecules — atoms ended up scattered to opposite faces of the supercell, breaking the FMO initial guess | Fixed 2026-05-10 — `cif_engine_ase.run_ase()` now calls `unwrap_molecules` (BFS along the bond graph + minimum-image translation) so each detected molecule sits contiguously in real space before fragment cutting |
