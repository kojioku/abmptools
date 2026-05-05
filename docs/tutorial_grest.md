# Tutorial: GENESIS gREST_SSCR via `abmptools.genesis.grest`

This tutorial walks through a full gREST_SSCR run end-to-end:
environment setup → smoke build → real MD → analysis (replica
transition, acceptance ratio, 1D distance PMF). It targets
abmptools v1.20.0 + GENESIS 2.1+.

For the package's API + design notes, see [`grest.md`](./grest.md).

---

## 1. Environment

### 1.1 Conda env (recommended)

```bash
mamba create -n grest -c conda-forge \
    python=3.11 numpy pandas matplotlib \
    ambertools=23 openmpi
mamba activate grest
```

### 1.2 abmptools

```bash
cd /path/to/abmptools
python -m pip install -e .[grest]
```

`pip install abmptools[grest]` pulls in `matplotlib` for plotting.
`numpy` is already a core dependency.

### 1.3 GENESIS build (LGPL-3.0-or-later)

```bash
git clone --depth=1 https://github.com/genesis-release-r-ccs/genesis.git
cd genesis
./configure CC=gcc FC=gfortran   # or CC=icx FC=ifx for Intel
make -j8
make install                     # binaries land in $PREFIX/bin/
```

If `make` fails on `analysis/src/lib/fileio_data_.c` under icx with
`ftello64 / fseeko64` undeclared, apply this patch:

```bash
sed -i 's/ftello64/ftello/g; s/fseeko64/fseeko/g' \
    analysis/src/lib/fileio_data_.c
make
```

Sanity check:

```bash
which spdyn atdyn remd_convert tleap mpirun
spdyn --version 2>&1 | head -3
```

---

## 2. Smoke build (KGG, ~2-5 min on 16 cores)

### 2.1 Generate config

```bash
mkdir -p ~/grest_smoke && cd ~/grest_smoke
cp /path/to/abmptools/sample/grest/kgg_smoke.json ./

# Or generate from scratch:
python -m abmptools.genesis.grest example > cfg.json
```

The shipped `kgg_smoke.json` reduces nsteps to 4000 in step3 and uses
the POC ladder (300, 318.11, 337.11, 357.10 K) so the full pipeline
runs in <5 minutes on a small workstation.

### 2.2 Prepare input PDB

Create a minimal KGG tripeptide PDB (or use any 3-residue protein).
For testing, the package has a synthesised KGG PDB used in the
integration test:

```bash
cat > kgg.pdb <<'EOF'
ATOM      1  N   LYS A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  LYS A   1       1.500   0.000   0.000  1.00  0.00           C
...                                # see tests/test_grest_integration.py
TER
END
EOF
```

### 2.3 Validate config

```bash
python -m abmptools.genesis.grest validate --config kgg_smoke.json
```

Expected output:

```
Configuration: OK (project=kgg_smoke, AMBER ff19SB + tip3p, 4 replicas, REST=explicit (3 residues), T=[300.00, ..., 357.10] K)

External tools:
  [OK]    tleap         (.../tleap)
  [OK]    spdyn         (.../spdyn)
  [OK]    atdyn         (.../atdyn)
  [OK]    remd_convert  (.../remd_convert)
  [OK]    mpirun        (.../mpirun)
  [WARN]  cpptraj       -- AmberTools -- around-mode REST residue resolution (optional)

REST selection preview:
  mode       = explicit
  selection  = rno:1-3
  n_residues = 3
  param_type = C L  (GENESIS Setup_Remd_Solute_Tempering tokens)

Replica temperature ladder (manual/geometric):
  rep01 ->  300.000 K     ratio T_i/T_0 = 1.000
  rep02 ->  318.110 K     ratio T_i/T_0 = 1.0604  adj T_i/T_0 = 1.0604
  rep03 ->  337.110 K     ratio T_i/T_0 = 1.1237  adj T_i/T_1 = 1.0597
  rep04 ->  357.100 K     ratio T_i/T_0 = 1.1903  adj T_i/T_2 = 1.0593

MPI mapping preview:
  total processes = 4 replicas x 2 MPI/replica = 8
  threads/proc    = 2 (OMP_NUM_THREADS)
  total cores     = 16
```

### 2.4 Build the input set

```bash
python -m abmptools.genesis.grest build --config kgg_smoke.json -o ./out
```

Output:

```
out/
├── input/config.json
├── build/system.{tleap,prmtop,coor}, system_ref.pdb, system_tleap.log
├── inp/step{1,2,3,5}_*.inp
├── logs/                  # populated by run.sh
├── run.sh                 # +x
├── run_pjsub.sh
└── run_sbatch.sh
```

### 2.5 Run MD

```bash
cd out
bash run.sh 2>&1 | tee logs/run.log
```

Expected wall-clock: 2–5 minutes for the smoke profile (4000 nsteps in
step3). Look for the `Setup_Remd>` block in `logs/step3.log` followed
by ~40 `REMD>` exchange entries.

---

## 3. Smoke vs production

The smoke run is a *connectivity test*: it verifies that tleap accepts
the PDB, atdyn minimises without nonbonded clashes, spdyn equilibrates,
and gREST_SSCR exchanges replicas. It does **not** sample meaningful
PMFs (only ~14 ps of total dynamics across 4 replicas).

For research-quality results:

- Use ~12 replicas spanning T_min=300 K, T_max=420 K (geometric).
- Use POC's nsteps = 60 M (~210 ns at dt=3.5 fs) per replica.
- Plan for ~12 nodes × ~12 hours on Fugaku (POC reference).

The shipped `sample/grest/poc_reproduction.json` matches the POC scale.

---

## 4. Production build (optional)

```bash
cp /path/to/abmptools/sample/grest/poc_reproduction.json ./
# Edit input_pdb to point at your protein.
python -m abmptools.genesis.grest validate --config poc_reproduction.json
python -m abmptools.genesis.grest build --config poc_reproduction.json -o ./prod
```

The `validate` step prints the geometric ladder ratios so you can
inspect that adjacent T_i/T_{i-1} sit in the 1.05–1.07 sweet-spot for
healthy acceptance.

---

## 5. Analysis

After `bash run.sh` completes (or even after a partial run), use the
`analyze` subcommand. Tasks can be requested individually or as `all`.

### 5.1 Parameter sort + replica transition + acceptance

```bash
python -m abmptools.genesis.grest analyze \
    --config kgg_smoke.json \
    --run-dir ./out \
    --what sort transition acceptance
```

Outputs (under `out/analysis/`):

- `replica_transition.png` — random walks of each replica through
  parameter space (cf. POC images 3 / 19).
- `replica_transition.csv` — raw `(step, replica, parameter)` triples.
- `acceptance_ratio.png` — cumulative acceptance ratio per replica
  pair, with a 100-exchange burn-in skip.
- `acceptance_ratio.csv` — raw `(step, rep_a, rep_b, ratio)` rows.

`build/param1.dcd` is the lowest-T parameter-sorted trajectory.

### 5.2 1D distance PMF

The POC validated a contact PMF between Lys96-NZ (antigen) and
Asn274-OD1 (antibody light-chain residue 31). Reproduce:

```bash
python -m abmptools.genesis.grest analyze \
    --config kgg_smoke.json \
    --run-dir ./out \
    --what pmf \
    --mask1 "rno:96 and an:NZ" \
    --mask2 "rno:274 and an:OD1" \
    --temperature 300.0
```

Outputs (under `out/analysis/`):

- `dist_pmf.png` — `-kT log P(r)` curve.
- `dist_pmf.xvg` — two-column `(r_A, pmf_kJ_per_mol)`.
- `distance.dat` — raw cpptraj distance time-series (intermediate).

The expected POC qualitative feature is a sharp peak around
`r = 3.0 Å` (the contact distance for the Lys-NZ ... Asn-OD1
H-bond network).

### 5.3 Mask syntax

Masks use GENESIS-style `rno:` (residue number) + `an:` (atom name)
selectors:

| Mask | Meaning |
|---|---|
| `rno:96` | all atoms of residue 96 |
| `rno:96 and an:NZ` | NZ of residue 96 |
| `rno:96 or rno:274` | both residues |

The package translates `rno:N and an:X` to cpptraj `:N@X` internally.

---

## 6. Around-mode REST selection

Instead of listing residues explicitly, ask cpptraj to find them:

```json
"rest_selection": {
  "mode": "around",
  "center": "rno:96",
  "radius_A": 5.0
}
```

At build time, abmptools renders a cpptraj input script that selects
residues with at least one atom within 5 Å of any atom in residue 96
(excluding water / common counterions), then parses the residue
list from cpptraj's `resinfo` output. This is helpful for rapid
iteration when you don't yet know which residues form the binding
interface.

`validate` shows the resolved list:

```
REST selection preview:
  mode       = around
  centre     = rno:96, radius     = 5.00 Å
  (resolved via cpptraj at build time)
```

`build` populates `output_dir/build/rest_residues.txt` with the
expanded list.

---

## 7. Failure modes

| Symptom | Cause | Fix |
|---|---|---|
| `ValueError: GrestStage.nsteps must be a multiple of 2*exchange_period` | config violates GENESIS Setup_Remd | Adjust `nsteps` to a multiple of `2 * exchange_period` |
| `tleap: failed (rc=1)` in `system_tleap.log` | malformed PDB, missing residue templates, or TER lines absent | Inspect `system_tleap.log` for the offending residue; pre-clean with `pdb4amber` |
| `Setup_Remd> mod(nsteps, rstout_period) is not ZERO` from spdyn | restart-period mismatch | dataclass enforces this at config-construct time; if you bypass via direct `.inp` edit, fix manually |
| `acceptance_ratio.png` shows ratios near 0.0 | temperature ladder spread too wide | Reduce `T_max` or increase `n_replicas`; aim for adjacent ratios 1.05–1.07 |
| `Setup_Mpi_Remd> total MPI != n_replicas * mpi/replica` | mpirun mismatch | Ensure `mpirun -np` matches `n_replicas * mpi_processes_per_replica` |
| `remd_convert` produces empty `param1.dcd` | one or more `step3_repN.rem` files truncated | Wait for the run to complete cleanly; `.rem` files are flushed only at the end |

---

## 8. Where to look next

- `docs/grest.md` — package reference, API, defaults.
- `tests/test_grest_*.py` — 159 unit tests + 1 slow integration smoke.
- `sample/grest/` — three configs (smoke / POC reproduction / around).
- POC notes: `/home/okuwaki/llm-project/SI/grest-POC.md`.
- gREST_SSCR original paper:
  Re, Kim & Sugita, *Int. J. Mol. Sci.* 22:270 (2021).
- Tutorial 12.3: <https://mdgenesis.org/tutorials/genesis_tutorial_12.3_2022/>
