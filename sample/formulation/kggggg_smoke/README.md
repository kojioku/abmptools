# `kggggg_smoke` — Smallest formulation smoke build

Tiny **AA mixed-solution** build modeled after Hossain et al. 2023
(*Nanoscale* 15, 19180) but using only commercial-permissive force
fields (AMBER ff14SB + GAFF2 + TIP3P + Joung-Cheatham).

## Composition

| Species | Count | Notes |
|---|---|---|
| KGGGG peptide (Ac-KGGGG-NMe) | 2 | natural one-letter sequence, ff14SB |
| Na-caprate (decanoate, neutral) | 8 | GAFF2 via acpype/AM1-BCC |
| Na-caprate (decanoate, charged COO⁻) | 8 | GAFF2 via acpype/AM1-BCC |
| Taurocholate | 2 | GAFF2 via acpype/AM1-BCC, charge -1 |
| Box | 6 nm cubic | smaller than Hossain (10 nm) for smoke runtime |
| Salt | 0.15 M NaCl | Joung-Cheatham |

Expected: ~8,000 atoms.

## Run

```bash
# Validate config + tool availability
micromamba run -n abmptoolsenv python -m abmptools.formulation validate \
    --config sample/formulation/kggggg_smoke/config.json

# Build (~1-3 min depending on acpype + RDKit speed)
micromamba run -n abmptoolsenv python -m abmptools.formulation build \
    --config sample/formulation/kggggg_smoke/config.json \
    --output-dir /tmp/kggggg_smoke

# Run MD (em + nvt + npt + 200 ps prod), ~10 min on 4-core CPU
cd /tmp/kggggg_smoke
NT=4 bash run.sh
```

## Notes

- This sample is for **smoke testing the builder pipeline end-to-end**.
  For statistically meaningful aggregation analysis use the 10 nm /
  500 ns scale (`insulin_smoke/config.json` as starting point, then
  bump `production.nsteps` to 250_000_000).
- KGGGG is not from the Hossain paper; substitute the desired
  peptide via `peptides[0].sequence` (single-letter, natural AAs
  only) or `peptides[0].pdb_path` (pre-built PDB).
