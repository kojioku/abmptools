# Glycine α-form MP2/6-31G(d) reference (abinitmp v2r8, 2026-05-10)

Numeric reference for organic-crystal FMO on the α polymorph of
glycine (zwitterion crystal), built from public-domain literature
(Marsh, *Acta Cryst.* 11, 654, 1958). All inputs are shippable.

## Files

| File | Contents |
|---|---|
| `expected_layer3_mp2_631gd_ifiesum.csv` | Target frag 1 sum: HF-IFIE / MP2-IFIE / ES / EX / CT-mix / `MonomerEnergy(1)` |
| `expected_layer3_mp2_631gd_ifiedt.csv`  | 1:1 IFIE detail for dimer-es=False pairs around frag 1 |
| `expected_isolated_monomer_mp2_631gd.txt` | Isolated 1-zwitterion HF/MP2 totals |

## Layer3 setup

| Quantity | Value |
|---|---|
| Space group | P21/n (Z=4) |
| Cell | a=5.105, b=11.969, c=5.464, β=111.65° |
| atoms_in_mol | 10 (full zwitterion: NH3+CH2COO−) |
| `bond_tolerance` | **0.2 Å** (lower than default 0.4 to avoid hydrogen-bond pair-off detection) |
| Layer 3 cutmode 6.0 Å | 21 fragments / 210 atoms |
| MPI ranks | 4 (`mpirun -np 4 abinitmp`) |
| Wall time | ~22.5 minutes |

## Sum-IFIE around frag 1 (kcal/mol)

| Quantity | Value |
|---|---|
| HF-IFIE | -711.27 |
| MP2-IFIE (dispersion correction) | +25.41 |
| ES | -1110.54 |
| EX | +1357.89 |
| CT-mix | -958.62 |

Glycine α-crystal is a **zwitterion lattice**. Each NH3+ / COO−
pair generates very large electrostatic stabilisation (ES) and
exchange repulsion (EX) that nearly cancel; the residual ~−700
kcal/mol HF-IFIE reflects the remaining net cohesive contribution.
This is normal for ionic-like organic crystals.

## Monomer energy: in crystal vs isolated (caveat!)

| Quantity | E' (in crystal, hartree) | MP2 total (isolated, hartree) | Δ (hartree) | Δ (kcal/mol) |
|---|---:|---:|---:|---:|
| monomer energy | -282.0075 | -282.5508 | +0.5433 | +341 ※ |

The fragment-1 atomic coordinates are identical in both runs.
The Δ above is **not** a structural deformation energy and is
**not** a clean polarization shift either. It is dominated by the
**zwitterion penalty in vacuum**: isolating a single NH3+/COO−
zwitterion at the crystal geometry, without the surrounding
electrostatic field, puts the SCF in a high-energy
charge-separated state that proton-transfer back to the neutral
NH2/COOH form would relieve. A clean polarization estimate would
require a separate gas-phase optimization on the neutral form, or
a constrained QM/QM cluster calculation; both are outside this
reference.

`expected_layer3_mp2_631gd_ifiesum.csv` and the IFIE detail are
still meaningful as raw numerics — they just shouldn't be combined
naïvely with the isolated-monomer total to give a single "lattice
energy".

## Reproducing

```bash
mamba activate abmptoolsenv
cd sample/crystal/glycine
abmp-crystal pipeline --config crystal_layer3.yaml --run-local   # ~22 min
abmp-crystal pipeline --config monomer.yaml --run-local            # ~5 s

DIR=out_layer3_v2/glycine/cifout/layer3/pdb/for_abmp
for ext in log ajf pdb xyz; do
    ln -sf glycinelayer3Zp1-around_ar6.0.${ext} \
           ${DIR}/glycine00001layer3Zp1-around_ar6.0.${ext}
done

python ../../../tests/regression/reference/main/crystal_csp7/R00001/extract_layer3_hf_631g.py \
    ${DIR}/glycine00001layer3Zp1-around_ar6.0.log \
    --out-dir reference/
```
