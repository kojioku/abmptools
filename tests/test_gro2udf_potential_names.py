# -*- coding: utf-8 -*-
"""ポテンシャル名は正規形 ``<型>-<型>`` にする。

以前は index を無条件に付けて ``c3-hc-0`` / ``c3-c3-1`` のようにしていた。
UDF の中では表と参照が揃うので解決はでき、`Export_GROMACS.py` も分子側の
``Potential_Name`` で引くだけなので通る。 ところが **J-OCTA で NPT を流すと
落ち、「力場を取得しなおす」と通るようになる** (2026-09-12 に実機で確認)。
その操作の後のファイルを見ると、名前が連番の無い正規形に書き直されていた。

型の組み合わせが重なるポテンシャルが複数あるときだけ連番で区別する。
名前は UDF のキーなので、重複させるわけにはいかないため。
"""
from __future__ import annotations

from abmptools.gro2udf.top_adapter import TopAdapter, _unique


def test_a_plain_name_keeps_its_shape():
    assert _unique("c3-hc", {}) == "c3-hc"


def test_a_repeat_gets_a_suffix():
    used = {}
    assert [_unique("c3-hc", used) for _ in range(3)] == [
        "c3-hc", "c3-hc-1", "c3-hc-2"]


def test_different_names_do_not_interfere():
    used = {}
    _unique("c3-hc", used)
    assert _unique("c3-c3", used) == "c3-c3"


def test_bond_names_have_no_index():
    """PE の 2 種類の結合が c3-hc / c3-c3 になる (以前は c3-hc-0 / c3-c3-1)。"""
    specs = TopAdapter._build_bond_type_specs([
        ("c3", "hc", 1, [0.1097, 661.2]),
        ("c3", "c3", 1, [0.1538, 601.8]),
    ])
    assert [s.name for s in specs] == ["c3-hc", "c3-c3"]


def test_angle_names_have_no_index():
    specs = TopAdapter._build_angle_type_specs([
        ("c3", "c3", "hc", 1, [110.0, 388.0]),
        ("hc", "c3", "hc", 1, [109.5, 329.0]),
    ])
    assert [s.name for s in specs] == ["c3-c3-hc", "hc-c3-hc"]


def test_a_genuine_duplicate_still_gets_distinguished():
    """同じ型の組で別パラメータが 2 つあっても、名前は衝突させない。"""
    specs = TopAdapter._build_bond_type_specs([
        ("c3", "hc", 1, [0.1097, 661.2]),
        ("c3", "hc", 1, [0.1100, 700.0]),
    ])
    names = [s.name for s in specs]
    assert names == ["c3-hc", "c3-hc-1"]
    assert len(set(names)) == 2


def test_torsion_multiplicity_terms_keep_the_colon_form():
    """多重度は ``:0`` / ``:1`` で分ける。 ここは J-OCTA も同じ形だった。"""
    specs = TopAdapter._build_torsion_type_specs([
        ("c3", "c3", "c3", "c3", 9, False,
         [0.0, 0.6, 3.0, 180.0, 0.25, 2.0, 0.0, 0.18, 1.0]),
    ])
    assert [s.name for s in specs] == [
        "c3-c3-c3-c3:0", "c3-c3-c3-c3:1", "c3-c3-c3-c3:2"]
