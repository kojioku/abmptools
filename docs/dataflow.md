# Data Flow

## Overview

ABMPTools processes molecular data through three main pipelines:

1. **FMO Setup** — Structural data → AJF input files for ABINIT-MP.
2. **Analysis** — ABINIT-MP outputs (LOG/CPF) → CSV results with IFIE/PIEDA data.
3. **DIFIE Averaging** — Multiple CPF snapshots → Single averaged CPF.

## Main Data Flow Diagram

```
                        ┌─────────────┐
                        │  PDB file   │
                        └──────┬──────┘
                               │
              ┌────────────────┼────────────────┐
              │                │                │
              ▼                ▼                ▼
        ┌──────────┐    ┌──────────┐    ┌──────────────┐
        │ pdb2fmo  │    │generateajf│   │  pdbmodify   │
        │ (config) │    │          │    │  (editing)   │
        └────┬─────┘    └────┬─────┘    └──────┬───────┘
             │               │                 │
             ▼               ▼                 ▼
        ┌──────────────────────┐          ┌─────────┐
        │    AJF input file    │          │ PDB out │
        └──────────┬───────────┘          └─────────┘
                   │
                   ▼
          ┌────────────────┐
          │   ABINIT-MP    │  (external calculation)
          │   execution    │
          └───────┬────────┘
                  │
         ┌────────┴────────┐
         ▼                 ▼
   ┌──────────┐      ┌──────────┐
   │ LOG file │      │ CPF file │
   └────┬─────┘      └────┬─────┘
        │                  │
   ┌────┴────┐        ┌────┴──────────┐
   ▼         ▼        ▼               ▼
┌───────┐ ┌───────┐ ┌───────────┐ ┌────────────┐
│log2cpf│ │getifie│ │convertcpf │ │generate    │
│       │ │pieda  │ │           │ │  _difie    │
└───┬───┘ └───┬───┘ └─────┬─────┘ └─────┬──────┘
    ▼         ▼           ▼              ▼
  CPF      CSV files    CPF (filtered)  DIFIE CPF
           (IFIE/PIEDA   (version       (averaged)
            matrices)     converted)
```

## FMO Setup Pipeline (PDB → AJF)

### Protein Systems

```
PDB ──→ pdb_io.readpdb() ──→ internal atom/residue data
                                      │
                                      ▼
                              generateajf.py
                              (method, basis, solvation, etc.)
                                      │
                                      ▼
                                  AJF output
```

Modules involved: `pdb_io.py` → `generateajf.py`

### Non-Protein Systems (Molecular Aggregates)

```
                ┌───────────────┐
                │ 1-molecule AJF│  (manual fragment definition)
                └───────┬───────┘
                        │
                        ▼
                   ajf2config.py ──→ segment_data.dat (config dict)
                                           │
PDB (full system) ──→ pdb2fmo.py ◀─────────┘
                          │
                          ▼
                   setfmo.setrfmoparam()
                   (cut mode: sphere/cube/around/neutral/none)
                          │
                          ▼
                    AJF + PDB output
```

Modules involved: `ajf2config.py` → `pdb2fmo.py` → `setfmo.py` → `generateajf.py`

### UDF (OCTA COGNAC MD) Systems

```
UDF file ──→ udfrm_io.convert_udf_pdb() ──→ PDB (per-frame)
                                                   │
                                                   ▼
                                          udf2fmo.py / pdb2fmo.py
                                                   │
                                                   ▼
                                              AJF output
```

Modules involved: `udf_io.py` → `udfrm_io.py` → `udf2fmo.py`

## Analysis Pipeline (LOG/CPF → CSV)

### Single-Point IFIE Analysis

```
LOG file ──→ anlfmo.readlog()
                 │
                 ├──→ version detection (getversion)
                 ├──→ fragment info (getfraginfo)
                 ├──→ IFIE/PIEDA data extraction
                 │         │
                 │    [optional: Fortran SO for speed]
                 │         │
                 ▼         ▼
           getifiepieda.py (mode selection)
                 │
         ┌───────┼────────┬──────────┐
         ▼       ▼        ▼          ▼
      frag     mol      ffmatrix   fragids
      mode     mode     mode       mode
         │       │        │          │
         ▼       ▼        ▼          ▼
       CSV     CSV    CSV files    CSV
    (per-frag) (per-mol) (per-PIEDA  (per-frag)
                         component)
```

### Time-Series / Multi-Sample Analysis

```
LOG files (numbered) ──→ getifiepieda --multi / --tfmatrix
                              │
                         multiprocessing.Pool
                         (parallel log reading)
                              │
                              ▼
                         Aggregated CSV
                    (IFIE over time, sum, individual)
```

Key output files per mode:
- `--frag`: `log-fragi-fragj-ifie.csv`
- `--ffmatrix`: `{ES,EX,CT,DI,HF,MP2,...}-ffmatrix.csv`, `Distance-ffmatrix.csv`
- `--tfmatrix`: Time × fragment matrix CSV (one file per reference fragment)
- `--multi`: `fragi-distr-ifiedt.csv` (per-pair), `fragi-distr-ifiesum.csv` (totals)

### CPF-Based Analysis

```
CPF file ──→ CPFManager.parse()
                  │
                  ├──→ atominfo DataFrame
                  ├──→ fraginfo DataFrame
                  ├──→ diminfo DataFrame (IFIE/PIEDA components)
                  ├──→ mominfo DataFrame (monomer energies)
                  └──→ static_data dict
                          │
                          ▼
                    cpf2ifielist.py
                    (filter by fragment range)
                          │
                          ▼
                   Formatted CSV (kcal/mol, Å)
```

## DIFIE Pipeline (Multiple CPFs → Averaged CPF)

```
CPF files (trajectory snapshots)
    file-001.cpf, file-002.cpf, ..., file-N.cpf
                        │
                        ▼
              generate_difie.py
              -t start end interval
              -np (parallel reading)
                        │
              ┌─────────┴─────────┐
              ▼                   ▼
    getcpfobj() per timestep   setoutcpfstatic()
    (parallel via Pool)        (prepare output structure)
              │                   │
              ▼                   ▼
         getavestddf()      representative
         (mean, std)        structure coords
              │                   │
              └─────────┬─────────┘
                        ▼
                  DIFIE CPF output
            (M-/S- prefixed columns for
             mean and std of each property)
```

## File Conversion Paths

```
LOG ──→ log2cpf.py    ──→ CPF
LOG ──→ log2config.py ──→ Python config dict
AJF ──→ ajf2config.py ──→ segment_data.dat
CIF ──→ readcif.py    ──→ Cartesian coords (XYZ/PDB)
CPF ──→ convertcpf.py ──→ CPF (different version / fragment subset)
UDF ──→ udfrm_io.py   ──→ PDB / XYZ
PDB ──→ pdbmodify.py  ──→ PDB (edited: renamed, moved, re-numbered)
AJF ──→ ajfserial.py  ──→ Numbered AJF series (for trajectory FMO)
```

## Amorphous Build Pipeline (SMILES/SDF → GROMACS MD)

```
PubChem CID / name ──→ amorphous.pubchem ──→ 3D SDF (local cache)
                                  │
                                  ↓   (or direct SMILES / SDF input)
SMILES / SDF  ──→ molecule_prep.py     ──→ OpenFF Molecule + single-mol PDB
                      │
                      ↓
                 density.py             ──→ count & cubic-box size from target density
                      │
                      ↓
                 packing.py (packmol)   ──→ mixture.pdb (multi-component packed)
                      │
                      ↓
                 parameterizer.py       ──→ system.gro / system.top (OpenFF + AM1-BCC)
                      │                    via openff-interchange
                      ↓
                 ndx_writer.py          ──→ system.ndx  (per-component index groups)
                      │
                      ↓
                 mdp_protocol.py        ──→ 01_em.mdp / 02_nvt / 03_npt / 04_anneal /
                      │                    05_npt_final.mdp  + run_all.sh + wrap_pbc.sh
                      ↓
                 (run_all.sh)           ──→ 5-stage GROMACS MD
                      ↓
                 (wrap_pbc.sh)          ──→ *_pbc.xtc / 05_npt_final_pbc.gro (VMD ready)
```

Orchestrated end-to-end by `abmptools.amorphous.builder.AmorphousBuilder.build()`.

## CG Builder Pipeline (Martini 3, 1.18.0+)

Two CG sub-packages live under `abmptools.cg/`. Both treat external
tools (`vermouth`, `insane`, `gmx`, `tleap`) as subprocess-only black
boxes — no source bundled or modified.

### `abmptools.cg.peptide` — Martini 3 peptide CG box (1.18.0)

```
PeptideBuildConfig (JSON)
    │
    ↓
peptide_atomistic.py     ──→ atomistic peptide PDB
   (tleap or extended-backbone fallback)
    │
    ↓
martinize_runner.py      ──→ <name>_cg.pdb + <name>.itp
   (martinize2 -ff martini3001, vermouth Apache-2.0 subprocess)
    │
    ↓
system_packer.py         ──→ packed.gro
   (gmx insert-molecules)
    │
    ↓
top_writer.py            ──→ topol.top  (M3 ITPs + peptide ITP)
    │
    ↓
water_box.py             ──→ martini_v3.0.0_water.gro (auto-generated
   if not in ff/, gmx insert-molecules with W density 8.36/nm³)
    │
    ↓
system_packer.solvate    ──→ system_solv.gro  (gmx solvate -cs)
    │
    ↓
system_packer.add_ions   ──→ system_ions.gro  (gmx grompp + genion,
                              NaCl 0.15 M, neutralize)
    │
    ↓
mdp_templates.py + run script writer
                         ──→ mdp/{em,nvt,npt,md}.mdp + index.ndx + run.sh
```

Orchestrated end-to-end by `abmptools.cg.peptide.builder.PeptideCGBuilder.build()`.

### `abmptools.cg.membrane` — Martini 3 peptide-membrane PMF (1.19.0)

```
MembraneCGBuildConfig (JSON)
    │
    ↓
_copy_ff_files            ──→ output_dir/martini_v3.0.0*.itp (4 files)
    │
    ↓
PeptideCGBuilder (sub-call,    ──→ molecules/<name>/{name}_cg.pdb + .itp
solvent_enabled=False,
mdp_*=False)
    │
    ↓
insane_runner.run_insane  ──→ bilayer.gro + insane_topol.top
   (insane GPL-2.0 subprocess; peptide -dm <z_offset> in POPC + W + NaCl)
    │
    ↓
topology_composer.compose_topology
                          ──→ topol.top  (4 M3 ITPs + peptide ITP +
                              "Protein" → moleculetype name +
                              NA+/CL- → NA/CL normalisation)
    │
    ↓
topology_composer.normalize_ion_atom_names_gro
                          ──→ system_ions.gro  (NA+/CL- → NA/CL in .gro)
    │
    ↓
system_packer.write_ndx_from_gro_cg
                          ──→ index.ndx  (Bilayer / Peptide / W / NA / CL /
                                         Non_Bilayer; no gmx make_ndx)
    │
    ↓
pulling.find_pbc_center_atom (imported from abmptools.membrane)
                          ──→ bilayer pbc-atom index
    │
    ↓
mdp_templates + pulling.write_pulling_mdp_cg + umbrella.write_window_mdps
                          ──→ mdp/{em,nvt,npt}.mdp +
                              pull/pull.mdp (NVT-chassis + direction-periodic) +
                              windows/win_NNN/window.mdp × N
                              (NPT-semiisotropic + direction + pbcatom)
    │
    ↓
umbrella.write_run_script ──→ run.sh
    │
    ↓
(bash run.sh)   em → nvt → npt → pull →
                python -m abmptools.cg.membrane make-windows →
                per-window grompp + mdrun (× N) →
                python -m abmptools.cg.membrane wham
                          ──→ analysis/pmf.xvg + histo.xvg
```

Orchestrated end-to-end by `abmptools.cg.membrane.builder.MembraneCGBuilder.build()` (build phase only) + `run.sh` (MD execution + post-analysis).

`pulling.parse_pullx_xvg / extract_window_frames / find_pbc_center_atom`,
`mdp_us_protocol.render_pull_block`, and `pmf.run_wham` are imported from
`abmptools.membrane.*` (duck-typed via `UmbrellaCGProtocol` field-name
match with `USProtocol`) — no helper duplication between CG and AA.

## GENESIS Builder Pipelines (1.20.0+)

GENESIS-based sub-packages live under `abmptools.genesis/`. Both treat
external binaries (`atdyn` / `spdyn` / `remd_convert`, LGPL-3.0+) as
subprocess invocations only — no source bundling, no linking.

### `abmptools.genesis.grest` — gREST_SSCR replica exchange (1.20.0)

```
input.pdb (protein)
        │
        ▼
[Stage 1] tleap                                ──→ system.{prmtop,coor,ref.pdb}
   AMBER ff19SB + TIP3P + solvateBox + addions
        │
        ▼
[Stage 2] REST residue resolution              ──→ rest_residues.txt
   ┌─ explicit: parse "1-138" / "21,96,274-275" → List[int]
   └─ around:   cpptraj `(:96)<@5.0 & !:WAT` mask → List[int]
        │
        ▼
[Stage 3] temperature ladder                   ──→ temperature_ladder.txt
   ┌─ auto:   T_i = T_min × (T_max/T_min)^(i/(N-1))  (geometric)
   └─ manual: passthrough
        │
        ▼
[Stage 4] inp_writer (4 GENESIS .inp files)    ──→ inp/{step1,step2,step3,step5}.inp
   step1_minimize        atdyn SD
   step2_equilibrate     spdyn VVER NPT
   step3_grest           spdyn VRES NVT + [REMD] type1=REST + ladder
   step5_remd_convert    parameter sort
        │
        ▼
[Stage 5] grest_runner (run.sh + HPC scaffolds)
   bash run.sh: atdyn → mpirun spdyn (equil) → mpirun spdyn (gREST) → remd_convert
        │
        ▼
analyze ─→ analysis/{replica_transition, acceptance_ratio, dist_pmf}.{png,csv}
   ・remd_convert で param sort (param1.dcd = 最低 T)
   ・replica transition plot (.rem files → matplotlib)
   ・acceptance ratio plot (REMD log 集計、burn_in=100)
   ・1D 距離 PMF (-kT log P(r)、cpptraj 距離 → numpy histogram)
```

Orchestrated by `abmptools.genesis.grest.builder.GrestBuilder.build()` (build phase) + `run.sh` (MD) + `analyze` subcommand (post-MD).

### `abmptools.genesis.mmgbsa` — MM/GBSA single-point ΔG_bind (1.22.0)

```
input/<target>.pdb (protein-ligand complex)
        │
        ▼
[Stage 1] pdb_splitter (Biopython)             ──→ <target>/<target>_{receptor,ligand}_<resno>.pdb
   - residue resno を ligand として抽出、それ以外を receptor
   - 任意の chain filter
        │
        ▼
[Stage 2] parameterize (acpype + tleap × 3)    ──→ <target>/{complex,ligand,receptor}.{prmtop,inpcrd,pdb}
   ┌─ acpype -i ligand.pdb -c bcc -k maxcyc=0
   │     → ligand.acpype/{*_AC.frcmod, *_bcc_gaff2.mol2}
   ├─ tleap -f leaprc_complex   (loadAmberParams + loadmol2 + loadpdb + combine)
   ├─ tleap -f leaprc_ligand    (mol2 only)
   └─ tleap -f leaprc_receptor  (pdb only)
   FF: ff14SB + DNA.OL15 + RNA.OL3 + TIP3P + GAFF/GAFF2
        │
        ▼
[Stage 3] gbsa_runner (mpirun atdyn × 3)       ──→ <target>/{complex,ligand,receptor}.{inp,log,dcd,rst}
   - [ENERGY] implicit_solvent=GBSA + NOBC + cutoffdist=99.9 Å
   - [MINIMIZE] method=SD nsteps=1 (single-point energy)
        │
        ▼
[Stage 4] analysis (log parser + plot)         ──→ analysis/{analysis_results.csv, dg_bind_plot.png}
   - [STEP4] Compute Single Point Energy 行から
     ENERGY 列 (= U_FF + ΔG_solv) と SOLVATION 列を抽出
   - ΔG_bind = E_complex - E_ligand - E_receptor (kcal/mol)
   - matplotlib bar plot
```

Orchestrated by `abmptools.genesis.mmgbsa.builder.MMGBSAOrchestrator.run()` (`fail_fast` 制御、N targets 逐次)。

### `abmptools.crystal` — Organic-crystal FMO pipeline (1.23.0)

```
inputs[*].cif (one or more CIF files)
        │
        ▼
[Stage 1] CIF → supercell PDB                    ──→ <project>/<name>/cifout/layer<N>/<name>layer<N>.pdb
   ┌─ cif_engine.engine = "legacy" : readcif.py ハードコード対称展開
   │     getsymcoord (1/3, 2/3, 1/4, 3/4, 1/6, 5/6 のみ) で symmetry 展開
   └─ cif_engine.engine = "ase"    : ase.io.read(cif, format="cif")
                                     + Atoms.repeat((2*layer+1,)*3)
                                     + ase.neighborlist で connected components
                                     + unwrap_molecules (BFS + find_mic で PBC 境界またぎ修正)
        │
        ▼
[Stage 2] PDB → for_abmp/                       ──→ <project>/<name>/cifout/layer<N>/pdb/for_abmp/<name>layer<N><Zp>-<cutmode>_<solute>r<criteria>.{ajf,pdb}
   - setfmo (multiple inheritance) で fragment cut around solute
     (cutmode=around / sphere / cube / neutral / none)
   - pdb2fmo + setfmo で `&BASIS BasisSet=<name>` (例: 6-31Gdag),
     `&SCF Method=<MP2/HF>` を render
   - is_xyz=True (default) で `&XYZ` block にフル精度浮動小数点座標を埋め込み
     (PDB 経由の %8.3f 切捨を回避、`&FMOXYZ ReadGeom=` を空にして LOG 解析互換維持)
        │
        ▼
[Stage 3] HPC jobscripts                         ──→ <project>/<name>/.../runab23q2_*.sh + runbatch.sh
   - job_templates.render_jobscript(spec, ajf_path) で string.Template
     系統別 4 種 (PJM / SLURM / PBS / local)
   - bash の ${var} と Python Template の ${var} 衝突は $$var エスケープで回避
        │
        ▼
[Stage 4] (optional) --run-local                 ──→ <project>/<name>/.../<ajf>.log
   - mpirun -np <npro> abinitmp <ajf> > <log>
     (abinitmp は MPI flat 想定; abinitmp_omp は MPI/OMP 混成、別バイナリ)
   - HPCJobSpec.mpi_launcher で mpirun コマンド override 可
        │
        ▼
[Stage 5] (optional) postprocess                 ──→ <project>/<name>/postproc/{ifie.csv, ifie_annotated.csv}
   - getifiepieda --frag <target> --multi --csv で IFIE/PIEDA 抽出
   - atom_distance.find_nearest_atoms で近接原子情報を CSV に annotate
```

Orchestrated end-to-end by `abmptools.crystal.builder.CrystalOrchestrator.run()` (Stage 1〜5 を逐次、`fail_fast` 制御、N inputs 並列)。Default は `legacy` engine (既存実績尊重)、`ase` は opt-in (`cif_engine.engine: ase`)。座標は `is_xyz=True` で `&XYZ` block 埋め込み。

## FMO Fragment Auto-splitter Pipeline (1.21.0+)

`abmptools.fragmenter` は MD pipeline ではなく、PDB → fragment
definition (`segment_data.dat`) のワンショット静的解析パイプライン。

```
input.pdb
        │
        ▼
pdb_loader (RDKit、4-strategy bond perception) ──→ List[Mol] (連結成分)
   1. proximityBonding=True,  sanitize=True   (default、CONECT 優先 + 近接補完)
   2. proximityBonding=False, sanitize=True   (CONECT のみ信用、誤結合回避)
   3. proximityBonding=True,  sanitize=False  (best effort、valence エラー時)
   4. obabel 前処理 → 再ロード                (subprocess fallback)
        │
        ▼
grouping (heavy-atom-only canonical SMILES)    ──→ List[MoleculeGroup]
   - MoleculeGroup: pattern, members, atom_to_member
        │
        ▼
auto_split (graph diameter MW walk)            ──→ List[CutSite]
   ・graph diameter (heavy-atom-only 2-pass BFS) で主鎖検出
   ・主鎖 walk + 側鎖 MW を累積
   ・累積 ≥ target_mw (default 200) で C-C 切断 candidate
   ・filter: 環内 / 多重結合 / ヘテロ隣接 (N/O/S/P/F/Cl/Br/I) 除外
        │
        ▼
[A 経路 = Jupyter UI]   open_panel: ipywidgets dropdown / SVG / checkbox
[C 経路 = ヘッドレス]   suggest → group_NNN.svg + group_NNN.json
                        (cut_sites.enabled 編集 → apply)
        │
        ▼
cut_apply (RWMol で破壊的に bond 削除)         ──→ List[FragmentResult]
        │
        ▼
expand_to_system (n_copies 全コピーへ展開)     ──→ segment_data.dat
                                                  (log2config 互換、pdb2fmo がそのまま読む)
```

ポリマー γ 経路 (`declare_same_pattern`) はオプションで、異なる SMILES (PE N=10 / N=11 等) を明示的に同一視し、master (最も cut の多い group) のパターンを atom-path-index 対応で短い chain にも転送する。

## Internal Data Structures

All modules converge on **pandas DataFrames** for structured data:

| DataFrame | Source Module | Key Columns |
|-----------|--------------|-------------|
| `atominfo` | `cpfmanager` | alabels, elems, elemtypes, resnames, resids, fragids, xcoords, ycoords, zcoords, MUL-HF |
| `fraginfo` | `cpfmanager` | Fragment definitions, atom ranges |
| `diminfo` | `cpfmanager` | fragi, fragj, min-dist, NR, HF, ES, EX, CT, DQ |
| `mominfo` | `cpfmanager` | fragi, DPM-HF-{X,Y,Z}, NR, HF |
| `static_data` | `cpfmanager`, `logmanager` | nuclear_repulsion_energy, total_energy, natom, nfrag, ndimer |

Unit conventions:
- Energies in IFIE/PIEDA tables: **kcal/mol** (conversion factor: 627.5095 from Hartree).
- Distances: **Angstrom** (conversion factor: 0.529177 from Bohr).
