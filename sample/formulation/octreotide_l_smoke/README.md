# `octreotide_l_smoke` — Octreotide-L analog smoke build

L-amino approximation of Hossain 2023 octreotide system (2 × peptide
+ 32 × Na-caprate + 2 × taurocholate, 10 nm cubic). All-L for tleap
+ ff14SB compatibility (the original octreotide has D-Phe1 + D-Trp4
+ Thr-ol C-term; this approximation uses Phe1 + Trp4 + Thr-ol → Thr).

## Sequence

`Ac-F-C-F-W-K-T-C-T-NMe`, **Cys2-Cys7 disulfide bridge**.

| Field | Value |
|---|---|
| Original octreotide | D-Phe-Cys-Phe-D-Trp-Lys-Thr-Cys-Thr-ol |
| This L-approximation | L-Phe-Cys-L-Phe-L-Trp-L-Lys-L-Thr-Cys-L-Thr |
| Disulfide | Cys2-Cys7 (declared in `disulfide_bonds`) |
| Caps | ACE (N-term) + NMe (C-term) |
| Net charge | +1 (Lys side chain) |

## Composition

| Species | Count | resname (after mangling) |
|---|---|---|
| octreotide_l | 2 (×8 residues each) | — |
| Na-caprate neutral | 16 | CPN |
| Na-caprate charged (COO⁻) | 16 | CPC |
| Taurocholate | 2 | TCH |
| Water (TIP3P) | filled to 10 nm cubic | WAT |
| Na+ / Cl- | neutralize + 0.15 M NaCl | NA / CL |

Expected: ~31,000 atoms.

## Run

```bash
# Validate
micromamba run -n abmptoolsenv python -m abmptools.formulation validate \
    --config sample/formulation/octreotide_l_smoke/config.json

# Build (~30 s — acpype on 3 small molecules + tleap + packmol)
micromamba run -n abmptoolsenv python -m abmptools.formulation build \
    --config sample/formulation/octreotide_l_smoke/config.json \
    --output-dir /tmp/octreotide_l_smoke

# Run MD (em + nvt + npt + 1 ns prod). ~30 min GPU / ~3 h CPU 4-core.
cd /tmp/octreotide_l_smoke
NT=4 bash run.sh
```

## Force-field caveats vs. Hossain 2023

- D-Phe1, D-Trp4 → L-Phe, L-Trp (β-turn around D-Trp not stabilized)
- Thr-ol C-term → Thr (C-term carboxyl, not alcohol)
- CGenFF 1.0.0 → GAFF2/AM1-BCC (taurocholate, caprate)
- CHARMM TIP3P → AMBER TIP3P

Qualitative aggregation / release direction is expected to reproduce
within ~5-15 kJ/mol PMF drift from the published values. For
rank-order or absolute-value comparison vs. Hossain, use the
D-isomer build via `octreotide_d_gaff_smoke/` (Phase 2; D-amino +
non-standard residues parameterised via acpype/GAFF2 on the whole
peptide).
