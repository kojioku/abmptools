# Public-molecule MP2/6-31G(d) reference summary

Cross-molecule summary of the four open-data references shipped under
`sample/crystal/{urea,glycine,benzene,naphthalene}/`. All values from
`abinitmp v2r8` runs with `mpirun -np 4` (MPI flat).

The same fragment-1 atomic coordinates are used for both columns of
each row — the in-crystal value is the polarized monomer of fragment 1
inside the layer-3 supercell (ABINIT-MP `E'(I)`-style `MonomerEnergy(1)`
column from `frag1-dimer-es-false-ifiesum.csv`), and the isolated value
is the MP2 total of the same monomer calculated alone in vacuum
(`expected_isolated_monomer_mp2_631gd.txt`). The Δ is the polarization-
plus-correlation contribution at fixed geometry; **it is not a
structural deformation energy** (the geometry is unchanged between
the two calculations).

| Molecule | E' (in crystal, hartree) | MP2 total (isolated, hartree) | Δ (hartree) | Δ (kcal/mol) |
|---|---:|---:|---:|---:|
| Urea           | -224.5950 | -224.6033 | +0.0083 | +5.2 |
| Benzene        | -231.4560 | -231.4575 | +0.0015 | +0.94 |
| Naphthalene    | -384.5579 | -384.5600 | +0.0021 | +1.3 |
| Glycine α      | -282.0075 | -282.5508 | +0.5433 | +341 ※ |

※ The glycine entry is dominated by the **vacuum zwitterion penalty**
(NH3+ / COO− at the crystal geometry without the surrounding
electrostatic field is unstable in vacuum), not by polarization. See
`sample/crystal/glycine/reference/README.md` for the full caveat.

## Sum-IFIE around frag 1 (kcal/mol, layer3 cutmode 6.0 Å)

| Molecule | HF-IFIE | MP2-IFIE | ES | EX | CT-mix |
|---|---:|---:|---:|---:|---:|
| Urea           |  -28.9 | -11.3 |   -37.7 |   +18.6 |   -9.8 |
| Benzene        |   -2.5 | -12.1 |    -4.5 |    +6.2 |   -4.2 |
| Naphthalene    |   -1.5 | -26.3 |    -6.8 |   +12.4 |   -7.0 |
| Glycine α      | -711.3 | +25.4 | -1110.5 | +1357.9 | -958.6 |

## Source files

For each molecule, `sample/crystal/<molecule>/reference/` contains:

- `expected_layer3_mp2_631gd_ifiesum.csv` — target frag-1 sum row (the
  source of the `E' (in crystal)` and HF-IFIE / MP2-IFIE / ES / EX /
  CT-mix columns above)
- `expected_layer3_mp2_631gd_ifiedt.csv`  — 1:1 IFIE detail
- `expected_isolated_monomer_mp2_631gd.txt` — `MP2 total (isolated)`
- `README.md` — per-molecule notes, reproduction recipe

## Phase / wall time

| Molecule | space group | Z | atoms_in_mol | bond_tolerance | layer3 frag/atom | wall time (mpirun -np 4) |
|---|---|---:|---:|---:|---|---:|
| Urea           | P-421m | 2 |  8 | 0.4 | 14 / 112 |  ~2.5 min |
| Glycine α      | P21/n  | 4 | 10 | 0.2 | 21 / 210 |  ~22 min |
| Benzene        | Pbca   | 4 | 12 | 0.4 | 26 / 312 |   ~7 min |
| Naphthalene    | P21/c  | 2 | 18 | 0.2 | 16 / 288 |  ~76 min |

## Source attribution

| Molecule | cif source |
|---|---|
| Urea           | Worsham, Levy & Peterson, *Acta Cryst.* 10, 319 (1957) — public-domain literature, hand-built via ASE `crystal()` |
| Glycine α      | Marsh, *Acta Cryst.* 11, 654 (1958) — public-domain literature, hand-built via ASE `crystal()` |
| Benzene        | Cox, Cruickshank & Smith, *Proc. R. Soc. Lond. A* 247, 1 (1958) — public-domain literature, hand-built via ASE `crystal()` |
| Naphthalene    | Hoser & Madsen, *Acta Cryst. A* 73, 102 (2017) — Crystallography Open Database entry [2311088](https://www.crystallography.net/cod/2311088.html), CC0 license, cif retained verbatim |
