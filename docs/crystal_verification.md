# `abmptools.crystal` Verification Record

This document records what is verified for the `abmptools.crystal`
subpackage and how to reproduce it. For API/design see
[`crystal.md`](./crystal.md); for the user tutorial see
[`tutorial_crystal_fmo.md`](./tutorial_crystal_fmo.md).

Last updated: 2026-05-08 (commit on `feature/crystal`).

---

## 1. Test matrix

63 collected unit + regression tests (62 passed + 1 gated live) cover the subpackage. Files under `tests/`:

| File | Tests | Verifies |
|---|---:|---|
| `test_crystal_models.py` | 25 | Dataclass round-trip (YAML/JSON), default values, constructor validation for `CrystalBuildConfig` and 7 leaf dataclasses |
| `test_crystal_cif_legacy.py` | 2 | `cif_engine_legacy.run_legacy` produces the same supercell PDB as the legacy `python -m abmptools.readcif` CLI |
| `test_crystal_cif_ase.py` | 13 | ASE backend: `read_cif_to_atoms` / `expand_supercell` / `detect_molecules` / PDB+XYZ writers / `run_ase_pipeline` end-to-end |
| `test_crystal_job_templates.py` | 8 | PJM / SLURM / PBS / local templates render with placeholder substitution and `template_override` |
| `test_crystal_atom_distance.py` | 6 | `find_nearest_atoms` against a known PDB, edge cases (single residue, multi-residue) |
| `test_crystal_builder.py` | 10 | `CrystalOrchestrator` 5-stage pipeline, `_emit_input_param` / `_emit_segment_data` driver synthesis, `_resolve_abinit_binary` 3-step lookup, `run_abinit` fake-binary execution and non-zero exit propagation |
| `test_crystal_cli.py` | 5 | 8 subcommands surfaced under `abmp-crystal`, `--help`, YAML-driven `pipeline` byte-equivalence |
| `test_crystal_regression.py` | 4 | **Phase B regression**: byte-equivalence with the v1.22.0 csp7 baseline (R00001 / R00002 / R00004 layer5) + legacy namespace identity |
| `test_crystal_numeric_regression.py` | 3 | **Phase D-3 numeric reference** (csp7 R00001 layer3 HF/6-31G, abinitmp v2r8): extract-script roundtrip, JSON shape validation, live regression (gated, ~9 h) |

Aggregate runtime: ~17 s for the regular layer (ase / abinitmp not exercised).

---

## 2. How to run

### Unit + light regression (default)

```bash
mamba activate abmptoolsenv
cd /home/okuwaki/llm-project/fcews-workspace/abmptools
pytest tests/test_crystal_*.py -v
```

Expected: 62 passed, 0 failed, 1 skipped (live numeric regression is
opt-in). Wall time ≤ 30 s. No abinitmp binary needed.

### Phase B byte-equivalence regression (slow)

```bash
pytest tests/test_crystal_regression.py -v -m slow
```

This re-runs the full pipeline through the legacy engine for csp7
R00001 / R00002 / R00004 layer5 and byte-compares for_abmp/*.{ajf,pdb}
against the fixture under
`tests/regression/reference/main/crystal_csp7/`. Wall time ~7 s.

### Integration smoke (real abinitmp, layer2 + HF/STO-3G)

```bash
bash tests/integration/run_crystal_smoke.sh
```

Runs `abmp-crystal pipeline --config crystal.yaml --run-local` on
csp7 R00001 layer2 (32 mol, 4 frag/cell) + HF/STO-3G + 1 thread.
Skips when `abinitmp` is not on PATH. Wall time ~30 s on a workstation.
Verifies that the rendered AJF is accepted by `abinitmp` and a log
banner is emitted.

### Live numeric regression (real abinitmp, layer3 + HF/6-31G)

```bash
ENABLE_FMO_LIVE_REGRESSION=1 \
ABINIT_DIR=/home/okuwaki/bin/abinitmp/v2r8 \
ABINIT_BIN=abinitmp \
pytest tests/test_crystal_numeric_regression.py::test_live_layer3_hf_631g_against_reference -v -m slow
```

Wall time ~9 hours on 1 core (current abinitmp v2r8 build does not
respond to `OMP_NUM_THREADS`; MPI parallelism is required for
production-scale runs). Compares Total energy, 24 monomer energies,
and 95 ESP-AOC IFIE values against the committed reference within the
tolerances stored in the JSON itself.

---

## 3. Verified environment (2026-05-08 snapshot)

| Component | Version |
|---|---|
| Python | 3.11.11 |
| pytest | 9.0.2 |
| ase | 3.28.0 (LGPL-2.1+) |
| pyyaml | 6.0.3 |
| numpy / scipy | per `abmptoolsenv` mamba env |
| abinitmp | v2r8 at `/home/okuwaki/bin/abinitmp/v2r8/abinitmp` |
| Host | WSL2 Ubuntu 22.04 |

The Phase A/B/C/D commits land on `feature/crystal` (`108722a` for
Phase A-C, `342e7c8` for Phase D, plus the Phase D-3 numeric reference
commit added 2026-05-08).

---

## 4. Numeric reference: csp7 R00001 layer3 HF/6-31G

Frozen in `tests/regression/reference/main/crystal_csp7/R00001/expected_layer3_hf_631g.json`
(21 KB). Companion excerpt log
`excerpt_layer3_hf_631g.log` (10 KB) lets the extract script roundtrip
without staging the full 322 KB log.

| Quantity | Value |
|---|---|
| Method | FMO2-HF / 6-31G |
| Fragments | 24 (cutmode `around` criteria 6.0 Å) |
| Atoms | 768 |
| Dimer pairs | 276 (95 ESP-AOC + 181 ESP-PTC + 0 non-approx) |
| Nuclear repulsion | +340271.8399768868 hartree |
| Electronic energy | -374722.6033745366 hartree |
| **Total energy** | **-34450.7633976498 hartree** |
| Monomer SCF time | 10280.1 s |
| Dimer SCF time | 19645.5 s |
| Total elapsed | 29925.7 s (8h 18min, 1 core abinitmp v2r8) |

Tolerances (live regression): ±1e-5 hartree on total / monomer
energies; ±1e-3 kcal/mol on IFIE values.

The full log + IFIE block is archived under
OneDrive `abmptools-dump/crystal/csp7_r00001_l3_hf_631g/` (not
checked into the repository — see memory `feedback_artifact_placement`).

---

## 5. Failure modes → test reverse index

| Symptom | Caught by |
|---|---|
| ASE backend produces wrong supercell origin or atom ordering | `test_crystal_cif_ase::test_run_ase_pipeline_*` |
| Supercell PDB drifts from v1.22.0 baseline | `test_crystal_regression::test_orchestrator_engine_legacy_*` |
| AJF format regression (Natom / `&XYZ` block) | `test_crystal_regression::test_orchestrator_engine_legacy_*` (byte-compare) |
| `_resolve_abinit_binary` resolution rule changes | `test_crystal_builder::test_resolve_abinit_*` |
| `--run-local` fails to propagate non-zero exits | `test_crystal_builder::test_run_abinit_fail_propagation` |
| Reference JSON / extractor drift | `test_crystal_numeric_regression::test_extract_script_roundtrip` |
| HPC template token escape regressions | `test_crystal_job_templates::test_render_*` |
| CLI subcommand removed or renamed | `test_crystal_cli::test_help_lists_subcommands` |
| YAML `cif_engine.engine` validation | `test_crystal_models::test_cif_engine_*` |

---

## 6. Known limitations of current verification

1. **abinitmp parallelism** — the v2r8 binary at the recorded path
   ignored `OMP_NUM_THREADS=4` and ran on a single core. The live
   regression therefore takes ~9 hours and is not suitable for CI.
   Production use requires the MPI build (Fugaku / cluster); reference
   data was collected on a single core for reproducibility on
   developer workstations.
2. **Single CIF for numeric reference** — only csp7 R00001 has a
   numeric snapshot. R00002 (Z=8) and R00004 (P-1, Z=2) have only
   byte-equivalence input fixtures, no FMO outputs.
3. **No MP2 reference yet** — HF/6-31G only. Production calculations
   for csp7 use MP2/6-31G(d) on Fugaku; that reference is on the
   roadmap but not yet checked in.
4. **Single space group group covered numerically** — only P21/n. Other
   space groups (P-1, P21/c, P212121) have input-generation tests but
   no abinitmp-verified outputs.

---

## 7. Reproducing the numeric reference from scratch

```bash
cd /home/okuwaki/llm-project/fcews-workspace/abmptools
mkdir -p .smoke_artifacts/csp7_r00001_l3_hf_631g
cp sample/crystal/csp7_smoke/{XXXI-MMFF-R00001.cif,UNK.ajf} \
   .smoke_artifacts/csp7_r00001_l3_hf_631g/

cat > .smoke_artifacts/csp7_r00001_l3_hf_631g/crystal.yaml <<'YAML'
project_name: csp7_r00001_l3_hf_631g
output_dir: ./out
inputs:
  - cif: XXXI-MMFF-R00001.cif
    layer: 3
    atoms_in_mol: [32]
cif_engine: {engine: legacy}
fragment:
  cutmode: around
  solutes: [0]
  criteria: 6.0
  molname: [UNK]
  template_ajf: ./UNK.ajf
fmo: {method: HF, basis_set: 6-31G, memory: 4000, abinit_ver: rev23, npro: 1, is_xyz: true}
hpc:
  scheduler: local
  abinit_dir: /home/okuwaki/bin/abinitmp/v2r8
  binary_name: abinitmp
  nodes: 1
  proc_per_node: 1
  omp_threads: 1
  elapse: '12:00:00'
postproc: {enable: false}
YAML

cd .smoke_artifacts/csp7_r00001_l3_hf_631g
python -m abmptools.crystal pipeline --config crystal.yaml --run-local

# ~9 h later, regenerate the JSON reference:
LOG=$(ls out/XXXI-MMFF-R00001/cifout/layer3/pdb/for_abmp/*.log)
python ../../tests/regression/reference/main/crystal_csp7/R00001/extract_layer3_hf_631g.py "$LOG" \
    --out ../../tests/regression/reference/main/crystal_csp7/R00001/expected_layer3_hf_631g.json
```
