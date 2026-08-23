# -*- coding: utf-8 -*-
"""udf2gro の単位換算ガードの検証。

`udf2gro` は値を ``udf.get(..., "[nm]")`` のように**単位を指定して読み**、換算は
UDFManager が UDF の ``Unit_Parameter`` を基準に行う。**``Unit_Parameter`` が無いと
換算は黙って素通りし**、Å の座標が nm として書き出される。`gmx grompp` は形式が
正しいので通してしまい、箱が 10 倍 (= 密度 1/1000) の系が警告なしに走る。

そのため宣言が無ければエラーにする。ここではそのガードを検証する。
"""
from __future__ import annotations

import pytest

from abmptools.udf2gro.udf_adapter import UdfAdapter


class _FakeUdf:
    """Unit_Parameter の有無だけを再現する最小のスタブ。"""

    def __init__(self, length=None):
        self.values = {}
        if length is not None:
            self.values = {
                "Unit_Parameter.Mass": 1.0,
                "Unit_Parameter.Energy": 4.184,
                "Unit_Parameter.Length": length,
            }

    def get(self, path, *args):
        return self.values.get(path)

    def put(self, value, path):
        self.values[path] = value


class TestUnitParameterGuard:
    def test_raises_when_undeclared(self):
        """宣言が無く指定も無ければ、黙って変換せずエラーにする。"""
        with pytest.raises(RuntimeError, match="Unit_Parameter"):
            UdfAdapter(_FakeUdf())._ensure_unit_parameter()

    def test_error_message_is_actionable(self):
        """対処法 (--unit all_atom) がメッセージに出ること。"""
        with pytest.raises(RuntimeError) as exc:
            UdfAdapter(_FakeUdf())._ensure_unit_parameter()
        assert "all_atom" in str(exc.value)

    def test_accepts_declared_unit_parameter(self):
        """UDF が宣言していればそのまま使い、上書きしない。"""
        udf = _FakeUdf(length=0.1)
        UdfAdapter(udf)._ensure_unit_parameter()
        assert udf.values["Unit_Parameter.Length"] == 0.1

    def test_all_atom_shortcut_injects_angstrom_kcal(self):
        """'all_atom' = 長さ Å (0.1 nm) / エネルギー kcal/mol (4.184 kJ/mol)。"""
        udf = _FakeUdf()
        UdfAdapter(udf, unit_parameter="all_atom")._ensure_unit_parameter()
        assert udf.values["Unit_Parameter.Length"] == pytest.approx(0.1)
        assert udf.values["Unit_Parameter.Energy"] == pytest.approx(4.184)
        assert udf.values["Unit_Parameter.Mass"] == pytest.approx(1.0)

    def test_explicit_tuple(self):
        udf = _FakeUdf()
        UdfAdapter(udf, unit_parameter=(2.0, 8.0, 0.5))._ensure_unit_parameter()
        assert udf.values["Unit_Parameter.Mass"] == pytest.approx(2.0)
        assert udf.values["Unit_Parameter.Energy"] == pytest.approx(8.0)
        assert udf.values["Unit_Parameter.Length"] == pytest.approx(0.5)

    def test_unknown_string_rejected(self):
        with pytest.raises(ValueError, match="unit_parameter"):
            UdfAdapter(_FakeUdf(), unit_parameter="reduced")

    def test_all_atom_constant(self):
        assert UdfAdapter.ALL_ATOM_UNIT == (1.0, 4.184, 0.1)
