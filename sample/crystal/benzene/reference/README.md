# Benzene MP2/6-31G(d) reference (abinitmp v2r8, 2026-05-10)

Numeric reference for organic-crystal FMO on benzene (the textbook
apolar molecular crystal), built from public-domain literature
(Cox, Cruickshank & Smith, *Proc. R. Soc. Lond. A* 247, 1, 1958).

## Files

| File | Contents |
|---|---|
| `expected_layer3_mp2_631gd_ifiesum.csv` | Target frag 1 sum: HF-IFIE / MP2-IFIE / ES / EX / CT-mix / `MonomerEnergy(1)` |
| `expected_layer3_mp2_631gd_ifiedt.csv`  | 1:1 IFIE detail for dimer-es=False pairs around frag 1 |
| `expected_isolated_monomer_mp2_631gd.txt` | Isolated benzene HF/MP2 totals |

## Layer3 setup

| Quantity | Value |
|---|---|
| Space group | Pbca (Z=4) |
| Cell | a=7.46, b=9.66, c=7.03 Å |
| atoms_in_mol | 12 |
| `bond_tolerance` | 0.4 (default) |
| Layer 3 cutmode 6.0 Å | 26 fragments / 312 atoms |
| MPI ranks | 4 (`mpirun -np 4 abinitmp`) |
| Wall time | ~7 minutes |

## Sum-IFIE around frag 1 (kcal/mol)

| Quantity | Value |
|---|---|
| HF-IFIE | -2.50 |
| MP2-IFIE (dispersion correction) | **-12.14** |
| ES | -4.54 |
| EX | +6.21 |
| CT-mix | -4.17 |

Benzene's Hartree-Fock IFIE is small (-2.5 kcal/mol) — the
electrostatic + exchange + CT terms nearly cancel for an apolar
aromatic. The **MP2 correction (-12.1 kcal/mol)** is the
dispersion contribution that actually binds the crystal; this is
the pedagogical example of why dispersion is essential for
hydrocarbon crystals.

## Deformation energy

Per-fragment deformation = `E_in_crystal − E_isolated`:

| Quantity | In-crystal (frag 1) | Isolated | Δ (hartree) | Δ (kcal/mol) |
|---|---|---|---|---|
| MP2 total | -231.4560 | -231.4575 | **+0.0015** | **+0.94** |

A small positive shift, consistent with weak polarization of an
apolar molecule in the crystal field.

## Reproducing

```bash
mamba activate abmptoolsenv
cd sample/crystal/benzene
abmp-crystal pipeline --config crystal_layer3.yaml --run-local   # ~7 min
abmp-crystal pipeline --config monomer.yaml --run-local            # ~5 s

DIR=out_layer3_v2/benzene/cifout/layer3/pdb/for_abmp
for ext in log ajf pdb xyz; do
    ln -sf benzenelayer3Zp1-around_ar6.0.${ext} \
           ${DIR}/benzene00001layer3Zp1-around_ar6.0.${ext}
done

python ../../../tests/regression/reference/main/crystal_csp7/R00001/extract_layer3_hf_631g.py \
    ${DIR}/benzene00001layer3Zp1-around_ar6.0.log \
    --out-dir reference/
```
