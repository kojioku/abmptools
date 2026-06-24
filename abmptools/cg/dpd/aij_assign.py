# -*- coding: utf-8 -*-
"""
abmptools.cg.dpd.aij_assign
---------------------------
**既存の Cognac DPD 入力 UDF** の相互作用パラメータ a を、aij.dat の内容で
割り当て直す (build-udf のように新規生成するのではなく、既存 UDF をパッチする)。

割り当て規則:
    - aij.dat が `aij` モード … 値をそのまま a に入れる
    - aij.dat が `chi` モード … `a = aii + chi / 0.306` (= aii + 3.268·χ、
      Groot-Warren。`AijMatrix.to_a_values()` が実施)。`aii` は基準同種 a (引数)

照合:
    UDF の `Interactions.Pair_Interaction[]` 各エントリの `Site1_Name` /
    `Site2_Name` (= 粒子名ペア) を読み、aij.dat のペアと **順不同で照合**して
    `Interactions.Pair_Interaction[].DPD.a` を書き換える。

挿入パス (Cognac 標準。dpdgen `udfdpd_io` と同じ):
    Interactions.Pair_Interaction[N].Site1_Name / Site2_Name / Potential_Type
    Interactions.Pair_Interaction[N].DPD.a          ← ここに a を put

I/O は `UDFManager` round-trip (OCTA 同梱、`abmptools.udfcharge` と同方式)。
`UDFManager.put` は numpy 値をサイレントに 0 化するため **`float()` cast 必須**
(`reference_udfmanager_put_numpy_silent_zero`)。

照合・変換ロジック (`build_a_lookup` / `match_aij_to_pairs`) は UDFManager 非依存で
単体テスト可能。`assign_aij_to_udf` のみ UDFManager を lazy import する。
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

from .models import AijMatrix

logger = logging.getLogger(__name__)

_PI = "Interactions.Pair_Interaction[]"


def _pair_key(a: str, b: str) -> frozenset:
    """順不同のペアキー (同種 a==b は 1 要素 frozenset)。"""
    return frozenset((a, b))


def build_a_lookup(aij: AijMatrix) -> Dict[frozenset, float]:
    """AijMatrix → 順不同ペアキー dict ``{frozenset({i, j}): a}``。

    chi モードは ``to_a_values()`` で a に変換済みの値が入る。
    """
    lut: Dict[frozenset, float] = {}
    for i, j, a in aij.to_a_values():
        lut[_pair_key(i, j)] = float(a)
    return lut


def match_aij_to_pairs(
    udf_pairs: List[Tuple[int, str, str]],
    aij: AijMatrix,
) -> Tuple[List[Tuple[int, float]], List[Tuple[int, str, str]]]:
    """UDF の Pair_Interaction ペアと aij を照合し、割当値と未照合を返す。

    Parameters
    ----------
    udf_pairs
        ``(pair_index, site1_name, site2_name)`` のリスト (UDF から読んだ順)。
    aij
        aij.dat 由来の AijMatrix。

    Returns
    -------
    (assignments, unmatched)
        assignments : ``(pair_index, a_value)`` のリスト (aij.dat に該当あり)
        unmatched   : ``(pair_index, site1, site2)`` のリスト (aij.dat に無し)
    """
    lut = build_a_lookup(aij)
    assignments: List[Tuple[int, float]] = []
    unmatched: List[Tuple[int, str, str]] = []
    for idx, s1, s2 in udf_pairs:
        key = _pair_key(s1, s2)
        if key in lut:
            assignments.append((idx, lut[key]))
        else:
            unmatched.append((idx, s1, s2))
    return assignments, unmatched


def assign_aij_to_udf(
    udf_path: Union[str, Path],
    aij: Union[AijMatrix, str, Path],
    out_path: Optional[Union[str, Path]] = None,
    aii: float = 25.0,
    only_dpd: bool = True,
    record: int = -1,
) -> Dict[str, Any]:
    """既存 Cognac DPD UDF の Pair_Interaction の a を aij で割り当てる。

    Parameters
    ----------
    udf_path
        既存 DPD 入力 UDF (Cognac、`Interactions.Pair_Interaction[]` を持つ)。
    aij
        `AijMatrix` か、aij.dat のパス (パスなら `read_aij(path, aii)` で読む)。
    out_path
        出力先。None なら ``udf_path`` を上書き。
    aii
        chi モード変換の基準同種 a (aij がパスのときのみ使用)。
    only_dpd
        True なら `Potential_Type == "DPD"` のペアのみ対象。
    record
        UDFManager の jump 先レコード (default -1 = static/最終)。

    Returns
    -------
    dict
        {assigned, unmatched, total_pairs, out_path, unused_aij_pairs}
    """
    from UDFManager import UDFManager  # OCTA 同梱 (PyPI 非配布)

    if not isinstance(aij, AijMatrix):
        from .aij_io import read_aij
        aij = read_aij(aij, aii=aii)

    u = UDFManager(str(udf_path))
    u.jump(record)

    n = int(u.size(_PI) or 0)
    udf_pairs: List[Tuple[int, str, str]] = []
    skipped_non_dpd = 0
    for k in range(n):
        ptype = u.get(f"{_PI}.Potential_Type", [k])
        if only_dpd and ptype != "DPD":
            skipped_non_dpd += 1
            continue
        s1 = u.get(f"{_PI}.Site1_Name", [k])
        s2 = u.get(f"{_PI}.Site2_Name", [k])
        udf_pairs.append((k, str(s1), str(s2)))

    assignments, unmatched = match_aij_to_pairs(udf_pairs, aij)

    for idx, a in assignments:
        u.put(float(a), f"{_PI}.DPD.a", [idx])  # float() cast 必須 (numpy silent-0)

    out = Path(out_path) if out_path else Path(udf_path)
    u.write(str(out))

    # aij.dat 側で UDF に使われなかったペア (情報提供用)
    udf_keys = {_pair_key(s1, s2) for _, s1, s2 in udf_pairs}
    unused = [
        (i, j) for (i, j, _a) in aij.to_a_values()
        if _pair_key(i, j) not in udf_keys
    ]

    logger.info(
        "assign_aij_to_udf: %d/%d pair(s) assigned, %d unmatched, "
        "%d non-DPD skipped -> %s",
        len(assignments), len(udf_pairs), len(unmatched), skipped_non_dpd, out,
    )
    if unmatched:
        logger.warning(
            "aij.dat に無く未割当のペア (UDF 側): %s",
            [(s1, s2) for _, s1, s2 in unmatched],
        )

    return {
        "assigned": assignments,
        "unmatched": unmatched,
        "total_pairs": len(udf_pairs),
        "skipped_non_dpd": skipped_non_dpd,
        "unused_aij_pairs": unused,
        "out_path": str(out),
    }
