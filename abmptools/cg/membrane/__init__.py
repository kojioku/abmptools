# -*- coding: utf-8 -*-
"""
abmptools.cg.membrane
----------------------
Martini 3 peptide-membrane PMF (umbrella sampling) builder.

`insane <https://github.com/Tsjerk/Insane>`_ (GPL-2.0) を ``subprocess`` で
呼び、Martini 3 lipid bilayer にペプチドを埋め込み、umbrella sampling 用
13 windows + WHAM 解析まで end-to-end に組む CG 系 2 番目のサブパッケージ。

ペプチド CG 化は :mod:`abmptools.cg.peptide` を内部 sub-call、umbrella +
WHAM の generic helper は :mod:`abmptools.membrane.{pulling,pmf,mdp_us_protocol}`
を import 経由で再利用する (コード重複なし)。

主要 API:
    LipidMix
        単一 lipid 種の指定 (POPC v1)
    PeptideMembraneSpec
        permeant ペプチド (sequence または既存 CG PDB)
    EquilibrationCGProtocol / PullingCGProtocol / UmbrellaCGProtocol
        プロトコルパラメータ群
    MembraneCGBuildConfig
        全システム設定 (JSON 往復可)

CLI:
    python -m abmptools.cg.membrane {build,validate,example,make-windows,wham}
"""
from .models import (
    EquilibrationCGProtocol,
    LipidMix,
    MembraneCGBuildConfig,
    PeptideMembraneSpec,
    PullingCGProtocol,
    UmbrellaCGProtocol,
)

__all__ = [
    "LipidMix",
    "PeptideMembraneSpec",
    "EquilibrationCGProtocol",
    "PullingCGProtocol",
    "UmbrellaCGProtocol",
    "MembraneCGBuildConfig",
]
