# -*- coding: utf-8 -*-
"""
abmptools.crystal
------------------
Crystal-FMO workflow: CIF を入力に有機物結晶のスーパーセル展開、FMO
入力 (ajf) 生成、HPC ジョブスクリプト出力、IFIE/PIEDA 解析までを
end-to-end に組める。

Phase A では legacy ファイル (abmptools.{readcif,pdb2fmo,ajf2config,
pdbmodify,getifiepieda}) を ``abmptools.crystal.legacy`` 経由で再 export
するのみ。本体実装 (CrystalOrchestrator / ase バックエンド / YAML config /
HPC テンプレ) は Phase C で追加。

CLI: ``python -m abmptools.crystal {example,validate}`` (Phase A の最小
セット。Phase C で expand/fragment/jobs/pipeline/postproc/nearest を追加)。

External dependencies (Phase C 以降、optional [crystal] extras):

- ``ase>=3.22`` (LGPL-2.1+) — CIF 読み込み、対称展開、スーパーセル
- ``pyyaml>=6.0`` — YAML config 読み込み (JSON でも動く)
- ``abinitmp`` — ローカル実行する場合のみ (HPC 提出は subprocess なし)

Sibling subpackages: ``cg.peptide``, ``cg.membrane``, ``genesis.grest``,
``genesis.mmgbsa``, ``fragmenter``。
"""
from . import legacy
from .models import (
    CIFEngineConfig,
    CIFInputSpec,
    CrystalBuildConfig,
    FMOMethod,
    FragmentTemplate,
    HPCJobSpec,
    PostProcessSpec,
)

# Heavy modules (depend on optional ase / yaml) are imported lazily by
# the CLI; the orchestrator is exposed eagerly because it only depends
# on the legacy engine when ``cif_engine.engine='legacy'`` (the default).
from .builder import CrystalOrchestrator

__all__ = [
    "legacy",
    "CIFEngineConfig",
    "CIFInputSpec",
    "CrystalBuildConfig",
    "CrystalOrchestrator",
    "FMOMethod",
    "FragmentTemplate",
    "HPCJobSpec",
    "PostProcessSpec",
]
