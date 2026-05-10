# Urea MP2/6-31G(d) reference (abinitmp v2r8, 2026-05-10)

Numeric reference for csp7-style organic-crystal FMO on urea built
**entirely from public-domain literature values** (Worsham, Levy &
Peterson, *Acta Cryst.* 10, 319, 1957). All inputs and outputs are
shippable in the public abmptools repo.

## Pipeline

| Step | Command |
|---|---|
| 1. Build `urea.cif` | `python build_urea_v2.py` (P-421m, Z=2, hand-crafted) |
| 2. Stage `UNK.ajf` template | `cp <csp7_smoke>/UNK.ajf .` (1-fragment FMO template, generic) |
| 3. Crystal layer3 (in-crystal IFIE) | `abmp-crystal pipeline --config crystal_layer3.yaml --run-local` |
| 4. Isolated monomer | `abmp-crystal pipeline --config monomer.yaml --run-local` |
| 5. getifiepieda extraction | `extract_layer3_hf_631g.py … --out-dir reference/` |

`mpirun -np 4 abinitmp` (flat MPI). Layer3 finishes in ~2.5 min on
4 cores; the isolated monomer finishes in seconds.

## Files

| File | Contents |
|---|---|
| `expected_layer3_mp2_631gd_ifiesum.csv` | Target frag 1 sum: HF-IFIE / MP2-IFIE / ES / EX / CT-mix / **`MonomerEnergy(1)`** (in-crystal HF + MP2 total, hartree) |
| `expected_layer3_mp2_631gd_ifiedt.csv`  | 1:1 IFIE detail for dimer-es=False pairs around frag 1, with `UNK<i>(<i>)` resname annotation |
| `expected_isolated_monomer_mp2_631gd.txt` | Isolated 1-molecule HF / MP2 total energies (hartree) — gas-phase reference at the same fragment geometry |

## Monomer energy: in crystal vs isolated

The fragment-1 atomic coordinates are identical in both calculations.
The "in-crystal" value is the polarized monomer SCF/MP2 of frag 1
inside the supercell; the "isolated" value is the same monomer
calculated alone in vacuum. The difference is therefore **not** a
structural deformation energy (which would require a separate
gas-phase relaxed geometry); it is the polarization-plus-correlation
contribution from sitting in the crystal field.

| Quantity | In crystal (frag 1) | Isolated | Δ (hartree) | Δ (kcal/mol) |
|---|---|---|---|---|
| HF total  | -223.9659 | -223.9704 | +0.0045 | +2.8 |
| MP2 total | -224.5950 | -224.6033 | +0.0083 | +5.2 |

For urea this Δ is in the ~few kcal/mol range, consistent with a
modest polarization shift. The crystal-environment IFIEs around
frag 1 sum to **-28.9 kcal/mol (HF)** and **-11.3 kcal/mol (MP2
corr)**, so the per-molecule cohesive contribution
(≈ ΔE_monomer + ∑IFIE/2) is in the −10 to −20 kcal/mol band, in
the right ballpark for urea's experimental cohesive energy
~16 kcal/mol.

## Reproducing

```bash
mamba activate abmptoolsenv
cd sample/crystal/urea
abmp-crystal pipeline --config crystal_layer3.yaml --run-local   # ~2.5 min
abmp-crystal pipeline --config monomer.yaml --run-local            # ~5 s

# Symlink the layer3 log to the 5-digit zero-padded name expected by
# the extract script:
DIR=out_layer3_v2/urea/cifout/layer3/pdb/for_abmp
for ext in log ajf pdb xyz; do
    ln -sf urealayer3Zp1-around_ar6.0.${ext} \
           ${DIR}/urea00001layer3Zp1-around_ar6.0.${ext}
done

python ../../../tests/regression/reference/main/crystal_csp7/R00001/extract_layer3_hf_631g.py \
    ${DIR}/urea00001layer3Zp1-around_ar6.0.log \
    --out-dir reference/
```

The `extract` step renames the canonical
`csv/frag1-dimer-es-false-{ifiesum,ifiedt}.csv` outputs of
`python -m abmptools.getifiepieda --multi 1 -dimeres -imd -zp 5
-t 1 1 1 -nof90` to the `expected_layer3_mp2_631gd_*.csv` names
shipped in this directory. They should be byte-equivalent to the
committed reference.
