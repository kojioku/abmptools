# csp7 single-structure smoke pipeline (private inputs)

The csp7 series cif and the matching `UNK.ajf` template are **internal
data** and are no longer bundled with the public abmptools repo. They
live in the abmptools-sample sibling repo at:

```
$ABMPTOOLS_SAMPLE_DIR/sample/csp7_ciftest/crystal_reference/R00001_layer3_hf_631g/
    XXXI-MMFF-R00001.cif
    UNK.ajf
    crystal.yaml
    run_local.sh
```

This directory now keeps only the **driver / config files** that
remain reusable across systems (and `run.sh` for the Phase A/B legacy
CLI). To run the smoke pipeline, first stage the cif and AJF template
from abmptools-sample:

```bash
export ABMPTOOLS_SAMPLE_DIR=~/repos/abmptools-sample
SRC="${ABMPTOOLS_SAMPLE_DIR}/sample/csp7_ciftest/crystal_reference/R00001_layer3_hf_631g"
cp "$SRC/XXXI-MMFF-R00001.cif" .
cp "$SRC/UNK.ajf" .

# Then run either Phase A/B legacy CLI…
bash run.sh
# …or the Phase C+ YAML-driven pipeline:
python -m abmptools.crystal pipeline --config crystal.yaml
```

The integration smoke at `tests/integration/run_crystal_smoke.sh` does
the staging automatically (skips with a friendly message if
`ABMPTOOLS_SAMPLE_DIR` is unset / the cif is missing).

## What's in this directory

- `crystal.yaml` — recommended Phase C+ pipeline config (declarative,
  folds `input_param` / `segment_data.dat` / `UNK.ajf` into one schema).
- `input_param` / `segment_data.dat` — Phase A/B legacy CLI drivers
  (Python-dict literal + per-fragment graph). Reusable across csp7
  systems since the cocrystal molecule layout is identical.
- `run.sh` — 4-step Phase A/B legacy recipe, kept for documentation.

## What it does

The legacy `run.sh` runs (after the cif/UNK.ajf are staged):

```bash
python -m abmptools.readcif -i XXXI-MMFF-R00001.cif -an 32 -l 5
cp input_param segment_data.dat UNK.ajf cifout/layer5/pdb/
cp cifout/layer5/xyz/XXXI-MMFF-R00001layer5Zp1.xyz cifout/layer5/pdb/
cd cifout/layer5/pdb
python -m abmptools.pdb2fmo -i XXXI-MMFF-R00001layer5Zp1.pdb -p input_param -xyz
```

The `-xyz` flag instructs `pdb2fmo` to emit the AJF in
**direct-coordinate mode**: the `&XYZ` block carries full
floating-point precision (e.g. `2.891228423347176`) instead of the
PDB-truncated `%8.3f` (3-digit decimal) format.

## Phase C+: YAML-driven pipeline (recommended)

```bash
python -m abmptools.crystal pipeline --config crystal.yaml
```

Output ends up in `./out/XXXI-MMFF-R00001/cifout/layer5/pdb/for_abmp/`
with the `.ajf` carrying `Natom=832` + a full-precision `&XYZ` block,
byte-identical to the abmptools-sample fixture under
`crystal_reference/phase_b_layer5/R00001/`.

The PJM jobscripts (12 nodes × 2 procs, OMP=24) are emitted alongside
the AJF; switch `hpc.scheduler` to `SLURM` / `PBS` / `local` to render
the corresponding scheduler templates instead.
