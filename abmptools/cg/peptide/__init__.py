# -*- coding: utf-8 -*-
"""
abmptools.cg.peptide
---------------------
Martini 3 peptide CG system builder.

vermouth-martinize (Apache-2.0) を ``subprocess`` で呼び、AmberTools tleap
(optional) で生成した atomistic PDB を CG 化する。GROMACS で solvate /
genion / MDP / run script まで生成する end-to-end builder。

主要 API:
    PeptideSpec
        1 ペプチドの仕様 (name / sequence / count)
    PeptideBuildConfig
        全システム設定 (JSON 往復可能、YAML optional)

CLI:
    python -m abmptools.cg.peptide {build,validate,example}
"""
from .models import PeptideBuildConfig, PeptideSpec

__all__ = ["PeptideBuildConfig", "PeptideSpec"]
