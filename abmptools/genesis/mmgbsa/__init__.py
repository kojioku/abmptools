# -*- coding: utf-8 -*-
"""
abmptools.genesis.mmgbsa
------------------------
GENESIS-based MM/GBSA single-point ΔG_bind calculator for protein-ligand
binding free energy estimation.

The 4-stage pipeline (`MMGBSAOrchestrator.run()`):

1. ``_stage1_split_pdb``     -- Biopython で receptor (指定残基除外) +
                                ligand (指定残基のみ) に分割
2. ``_stage2_parameterize``  -- ``acpype`` (GAFF/GAFF2 + AM1-BCC) +
                                ``tleap`` で complex / ligand / receptor
                                の 3 系 prmtop+inpcrd 生成
3. ``_stage3_run_gbsa``      -- ``mpirun -np 1 atdyn`` で
                                ``[ENERGY] implicit_solvent=GBSA`` +
                                NOBC + 1-step minimize の single-point
                                energy 計算 (3 系)
4. ``_stage4_analyze``       -- ``[STEP4] Compute Single Point Energy``
                                log パース → Energy + Solvation 抽出 →
                                ``ΔG_bind = (E+S)_complex − (E+S)_ligand
                                − (E+S)_receptor`` → CSV + matplotlib
                                棒グラフ

Force field default: AMBER **ff14SB** + DNA.OL15 + RNA.OL3 + TIP3P +
GAFF/GAFF2 (POC 通り、MM/GBSA で広く使われる安定 default。
``abmptools.genesis.grest`` の ff19SB とは意図的に別)。

External dependencies (subprocess only, never bundled):

- GENESIS atdyn >= 2.1 (LGPL-3.0-or-later) -- subprocess only
- AmberTools tleap (free academic + commercial) -- conda 推奨
- acpype (GPL-3.0) -- subprocess only (mere aggregation per GPL FAQ)
- mpirun (OpenMPI / MPICH) -- ``mpirun -np 1 atdyn``
- Biopython (Biopython License + BSD-3-Clause) -- PDB splitter
- matplotlib (PSF/BSD) -- ΔG_bind 棒グラフ

CLI: ``python -m abmptools.genesis.mmgbsa
{example,validate,divide,parameterize,run,analyze,pipeline}``

Sibling: ``abmptools.genesis.grest`` (1.20.0、REMD/REST sampling).
"""
