# Naphthalene MP2/6-31G(d) reference (abinitmp v2r8, 2026-05-10)

Numeric reference for organic-crystal FMO on naphthalene, using
the room-temperature X-ray structure from the **Crystallography
Open Database** (CC0):

> Hoser, A. A. and Madsen, A. Ø.
> *Acta Crystallographica Section A: Foundations and Advances* 73,
> 102-114 (2017).
> COD entry: [2311088](https://www.crystallography.net/cod/2311088.html)

The cif is shipped verbatim under `sample/crystal/naphthalene/naphthalene.cif`,
preserving the COD header attribution as required by the database.

## Files

| File | Contents |
|---|---|
| `expected_layer3_mp2_631gd_ifiesum.csv` | Target frag 1 sum: HF-IFIE / MP2-IFIE / ES / EX / CT-mix / `MonomerEnergy(1)` |
| `expected_layer3_mp2_631gd_ifiedt.csv`  | 1:1 IFIE detail for dimer-es=False pairs around frag 1 |
| `expected_isolated_monomer_mp2_631gd.txt` | Isolated naphthalene HF/MP2 totals |

## Layer3 setup

| Quantity | Value |
|---|---|
| Space group | P21/c (Z=2) |
| Cell | a=7.825, b=5.935, c=8.100 Å, β=114.44° |
| Temperature | 293 K |
| atoms_in_mol | 18 |
| `bond_tolerance` | **0.2 Å** (lower than default 0.4 to avoid aromatic stack false-bonds) |
| MPI ranks | 4 (`mpirun -np 4 abinitmp`) |

| Layer 3 cutmode 6.0 Å | 16 fragments / 288 atoms |
| Wall time | ~76 minutes |

## Sum-IFIE around frag 1 (kcal/mol)

| Quantity | Value |
|---|---|
| HF-IFIE | -1.45 |
| **MP2-IFIE (dispersion correction)** | **-26.34** |
| ES | -6.82 |
| EX | +12.35 |
| CT-mix | -6.98 |

Naphthalene's HF-IFIE is small (-1.5 kcal/mol) — much like
benzene, the electrostatic + exchange terms nearly cancel for the
apolar aromatic core. **MP2 contributes the dominant binding via
π-stacking dispersion (-26.3 kcal/mol)**, which is roughly 2× the
benzene MP2-IFIE — consistent with the larger polarizable
π-system. Together with benzene this is the textbook hydrocarbon
"dispersion is essential" demonstration.

## Deformation energy

Per-fragment deformation = `E_in_crystal − E_isolated`:

| Quantity | In-crystal (frag 1) | Isolated | Δ (hartree) | Δ (kcal/mol) |
|---|---|---|---|---|
| MP2 total | -384.5579 | -384.5600 | **+0.0021** | **+1.3** |

A small positive shift, in line with weak polarization of an
apolar molecule in the crystal field, similar to benzene.

## Reproducing

```bash
mamba activate abmptoolsenv
cd sample/crystal/naphthalene
abmp-crystal pipeline --config crystal_layer3.yaml --run-local
abmp-crystal pipeline --config monomer.yaml --run-local

DIR=out_layer3_v2/naphthalene/cifout/layer3/pdb/for_abmp
for ext in log ajf pdb xyz; do
    ln -sf naphthalenelayer3Zp1-around_ar6.0.${ext} \
           ${DIR}/naphthalene00001layer3Zp1-around_ar6.0.${ext}
done

python ../../../tests/regression/reference/main/crystal_csp7/R00001/extract_layer3_hf_631g.py \
    ${DIR}/naphthalene00001layer3Zp1-around_ar6.0.log \
    --out-dir reference/
```
