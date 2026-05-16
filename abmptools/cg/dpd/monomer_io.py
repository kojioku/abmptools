# -*- coding: utf-8 -*-
"""
abmptools.cg.dpd.monomer_io
---------------------------
cg_segmenter dpdgen_exporter が生成する `{name}_monomer` ファイルの読み込み。

`{name}_monomer` は Python script で以下を定義する:

    bond12      = [[0, 1], [1, 2], ...]
    bond12h     = [[1, 0, 1, 0.86, 50, 0], ...]
    bond13_150  = [[0, 2], ...]
    bond13_150h = [...]
    bond14_150  = [[0, 3], ...]
    bond14_150h = [...]
    angle13     = [[0, 1, 2], ...]
    angle13data = [[0, 1, 2, 30, 5.0], ...]

dpdgen 互換 monomer file (= ``monomers = [{'name': ..., 'particle': [...], 'bond': [...], 'angle': [...]}]``)
も読み込めるよう、 :func:`read_monomers_dict` 経由でサポートする。
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

from .models import MonomerSpec

logger = logging.getLogger(__name__)


def read_monomer(path: Union[str, Path], name: Optional[str] = None) -> MonomerSpec:
    """cg_segmenter dpdgen_exporter 形式の `{name}_monomer` を読む。

    Parameters
    ----------
    path : str | Path
        monomer ファイルパス (Python script 形式)。
    name : str | None
        monomer 名 (省略時は path stem から `_monomer` を除いた部分を採用)。

    Returns
    -------
    MonomerSpec

    Notes
    -----
    cg_segmenter の出力にはまだ ``particle_names`` (segment ラベル) が含まれて
    いない (= 単に 0..n の index のみ)。 ここでは bond12 / angle13 から登場する
    粒子 index の最大値 + 1 を n_particles とし、 ``particle_names`` は
    ``["P0", "P1", ...]`` で仮置きする。
    実 segment 名 (例: "segA", "WAT") との対応は、 :class:`DpdSystemSpec` を組む
    側で `MonomerSpec.particle_names` を上書きする責務。
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"monomer file not found: {path}")
    if name is None:
        stem = path.stem
        name = stem[: -len("_monomer")] if stem.endswith("_monomer") else stem

    ns: Dict[str, Any] = {}
    exec(path.read_text(encoding="utf-8"), ns)  # noqa: S102

    bond12 = [tuple(b) for b in ns.get("bond12", [])]
    bond13_150 = [tuple(b) for b in ns.get("bond13_150", [])]
    bond14_150 = [tuple(b) for b in ns.get("bond14_150", [])]
    angle13 = [tuple(a) for a in ns.get("angle13", [])]

    bond12h = list(ns.get("bond12h", []))
    bond13_150h = list(ns.get("bond13_150h", []))
    bond14_150h = list(ns.get("bond14_150h", []))
    angle13data = list(ns.get("angle13data", []))

    # 粒子数推定 (bond12 / angle13 の最大 index + 1)
    max_idx = -1
    for src in (bond12, bond13_150, bond14_150):
        for (i, j) in src:
            max_idx = max(max_idx, i, j)
    for (a, b, c) in angle13:
        max_idx = max(max_idx, a, b, c)
    n = max_idx + 1 if max_idx >= 0 else 0

    spec = MonomerSpec(
        name=name,
        particle_names=[f"P{i}" for i in range(n)],
        bond12=bond12,
        bond12h=bond12h,
        bond13_150=bond13_150,
        bond13_150h=bond13_150h,
        bond14_150=bond14_150,
        bond14_150h=bond14_150h,
        angle13=angle13,
        angle13data=angle13data,
    )
    logger.info(
        "read_monomer: %s (n_particles=%d, bond12=%d, bond13_150=%d, bond14_150=%d, angle13=%d)",
        name, n, len(bond12), len(bond13_150), len(bond14_150), len(angle13),
    )
    return spec


def read_monomers_dict(path: Union[str, Path]) -> List[Dict[str, Any]]:
    """dpdgen 互換の ``monomers = [{...}]`` を読み込む。

    cg_segmenter は ``monomers`` 変数を出力しないが、 dpdgen 既存資産との
    互換のため、 ``monomers`` 変数を含む monomer file も読み込めるようにする。

    Returns
    -------
    List of monomer dict (元の Python dict そのまま)
    """
    path = Path(path)
    ns: Dict[str, Any] = {}
    exec(path.read_text(encoding="utf-8"), ns)
    monomers = ns.get("monomers", [])
    if not monomers:
        raise ValueError(
            f"'monomers' variable not found in {path}; "
            f"use read_monomer() for cg_segmenter-style files"
        )
    return list(monomers)


def assign_particle_names(spec: MonomerSpec, names: List[str]) -> MonomerSpec:
    """``MonomerSpec.particle_names`` を差し替える (新規 spec を返す)。

    Parameters
    ----------
    spec : MonomerSpec
        対象 spec (元 spec は変更しない)。
    names : List[str]
        n_particles と同じ長さの新しい segment label list。
    """
    if len(names) != spec.n_particles:
        raise ValueError(
            f"name list length {len(names)} != n_particles {spec.n_particles}"
        )
    import dataclasses
    return dataclasses.replace(spec, particle_names=list(names))
