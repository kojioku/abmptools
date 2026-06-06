# `abmptools.formulation` — AA peptide-formulation mixed-solution builder (v1.30.0)

A builder + analysis stack for atomistic MD simulations of *peptide
drug × permeation enhancer × bile salt* mixed solutions in a cubic
water box, modeled after Hossain et al. 2023 (*Nanoscale* 15,
19180-19195).

> **Platform**: 現行 (Amber route) は Linux / macOS / WSL2 のみ
> (AmberTools が Windows native install 不可のため)。 Windows native
> 動作のための **OpenFF route (Phase 1 開発中)** の詳細は
> [`platform_support.md`](platform_support.md) を参照。

## Why this sub-package

- `amorphous` packs small organics (single composition class) for
  OpenFF — no peptide support, no AMBER backend.
- `membrane` builds peptide × lipid bilayer with AMBER/CHARMM — no
  small molecule mixing.
- `cg.membrane` is Martini coarse-grained.

Hossain 2023 — and a wide class of formulation studies — needs all
of those at once: **peptide + small-molecule enhancer + bile salt +
water + ions in a single cubic box**. `formulation` fills the gap.

## License-clean force-field choice

| Force field | License | Used for |
|---|---|---|
| AMBER ff14SB (`leaprc.protein.ff14SB`) | public domain | peptide backbone + sidechain |
| GAFF2 (via acpype) | public domain (AMBER); acpype = GPL-3.0 subprocess only | small molecules |
| AM1-BCC charges | acpype/sqm | partial charges |
| TIP3P (`leaprc.water.tip3p`) | public domain | water |
| Joung-Cheatham Na/Cl | public domain | ions |

**No CHARMM36 / CGenFF**: Hossain 2023 used those, but they require
academic licensing. Our reproduction uses commercial-permissive OSS
force fields. Qualitative aggregation / release-order trends are
expected to reproduce; absolute PMF values may differ by ~5-15 kJ/mol.

## Pipeline

```
FormulationBuildConfig (JSON)
        │
        ▼
┌────────────────────────────────────────────────────────────────┐
│  FormulationBuilder.build() — 7 stages                          │
├────────────────────────────────────────────────────────────────┤
│  0  ff_staging              [output_dir/build, md, analysis]    │
│  1  peptide_atomistic       tleap (ff14SB) → <name>.pdb         │
│  2  small_mol_parameterize  acpype GAFF2/AM1-BCC × N species    │
│  3  packmol                 multi-component cubic box           │
│  4  solvate_ions            tleap solvatebox + parmed → top/gro │
│  5  index                   resname-table → system.ndx          │
│  6  mdp_render              em/nvt/npt/prod (Hossain protocol)  │
│  7  run_script              em → nvt → npt → prod bash wrapper  │
└────────────────────────────────────────────────────────────────┘
```

Output layout:

```
<output_dir>/
├── build/
│   ├── input/
│   │   ├── <peptide>.pdb
│   │   ├── <resname>.pdb           per small molecule monomer
│   │   └── <resname>.acpype/       per small molecule (itp/mol2/frcmod/gro)
│   ├── mixture.pdb                 packmol output
│   ├── system.tleap / system_tleap.log / system.prmtop / system.inpcrd
│   └── system.top / system.gro     parmed gromacs output (copies)
├── system.top / system.gro / system.ndx
├── md/em.mdp / nvt.mdp / npt.mdp / prod.mdp
├── analysis/                       populated by `formulation analyze`
├── run.sh
└── config.json                     snapshot of input config
```

## Data classes

| Class | Field highlights |
|---|---|
| `PeptideSpec` | `sequence` OR `pdb_path`, `n_copies`, `disulfide_bonds` (tleap `bond` directives), `cap_n`/`cap_c` |
| `EnhancerSpec` | `resname`, neutral + charged forms with independent `smiles_*` and counts (Hossain 50:50 fatty-acid convention) |
| `BileSaltSpec` | `resname`, `smiles`, `net_charge` (passed to acpype `-n`) |
| `SystemSpec` | `peptides`, `enhancers`, `bile_salts`, `box_size_nm`, `salt_concentration_M`, `neutralize` |
| `EquilibrationProtocol` | em/nvt/npt configurable (defaults Hossain) |
| `ProductionProtocol` | Parrinello-Rahman isotropic NPT |
| `USProtocol` | n_windows, window spacing, k, pull rate (for `release_us`) |
| `FormulationBuildConfig` | wraps SystemSpec + protocols + tool paths + force-field selection |

JSON round-trip via `cfg.to_json(path)` and
`FormulationBuildConfig.from_json(path)`.

## CLI

```bash
python -m abmptools.formulation example [--output out.json]
python -m abmptools.formulation validate --config cfg.json
python -m abmptools.formulation build    --config cfg.json [--output-dir DIR]
python -m abmptools.formulation analyze  --traj prod.xtc --top system.top --out analysis/ \
                                          [--enhancer-resnames CPRN,CPRC] \
                                          [--bile-salt-resnames TCH]
python -m abmptools.formulation release_us --config cfg.json --target-peptide-idx 0
```

## Force-field trade-off

| Axis | Hossain 2023 (CHARMM36 + CGenFF) | This module (ff14SB + GAFF2) | Impact |
|---|---|---|---|
| Peptide backbone | CHARMM36m | AMBER ff14SB | minor; ff14SB slightly under-stabilises α-helix |
| Fatty acid neutral/charged | CGenFF (LJ tuned) | GAFF2 + AM1-BCC | CMC may shift by ~10-20% |
| Bile salt (steroid scaffold) | CGenFF | GAFF2 generic | macrocycle stiffness less accurate |
| Water | CHARMM-modified TIP3P | TIP3P plain | negligible |
| Ions | CHARMM Na/Cl | Joung-Cheatham | JC more accurate for ion pairing |
| License | academic | full commercial | enables industrial use |

Use this module when **qualitative** aggregation / release order is
the question, or for industrial pre-screening. For absolute binding
free energies that must match the published values, fall back to
the academic CHARMM36 + CGenFF stack (out of scope here).

## License posture

- acpype (GPL-3.0) is invoked via subprocess only (mere aggregation
  per GPL FAQ; precedent in `genesis/mmgbsa/ligand_parameterize.py`).
- MDAnalysis (GPL-2.0) is imported lazily inside
  `analysis/aggregate_transition.py` and `analysis/contact_map.py`,
  and is shipped only as an opt-in `[formulation-analysis]` extra. The
  abmptools wheel does **not** install or link MDAnalysis by default.
- The default install (`pip install abmptools` + `[formulation]`
  extras: parmed, rdkit-pypi, matplotlib) is fully Apache-2.0 compatible.

## See also

- `sample/formulation/kggggg_smoke/` — smallest end-to-end build
- `sample/formulation/insulin_smoke/` — Hossain-style 99k-atom build
  (PDB 2G4M downloaded on demand)
- `docs/tutorial_formulation_smoke.md` — step-by-step setup + run
