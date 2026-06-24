# -*- coding: utf-8 -*-
"""tests/cg_dpd/test_aij_assign.py

既存 DPD UDF への aij 割り当て (assign_aij_to_udf) の検証。

- 照合 / chi→a 変換ロジック (build_a_lookup / match_aij_to_pairs): UDFManager 非依存
- assign_aij_to_udf の I/O: モック UDFManager で配線確認 (実 OCTA 不要)
"""
from __future__ import annotations

import sys
import types

import pytest

from abmptools.cg.dpd import AijMatrix
from abmptools.cg.dpd.aij_assign import (
    build_a_lookup, match_aij_to_pairs, _pair_key,
)
from .conftest import requires_cognac


# ---------------------------------------------------------------------------
# build_a_lookup / chi→a 変換
# ---------------------------------------------------------------------------


def test_build_a_lookup_a_mode_direct():
    aij = AijMatrix(segments=["A", "B"],
                    pairs=[("A", "A", 25.0), ("A", "B", 40.0), ("B", "B", 25.0)],
                    mode="a", aii=25.0)
    lut = build_a_lookup(aij)
    assert lut[_pair_key("A", "B")] == 40.0
    assert lut[_pair_key("A", "A")] == 25.0


def test_build_a_lookup_chi_mode_conversion():
    """chi モードは a = aii + chi/0.306 (= aii + 3.268·χ)."""
    aij = AijMatrix(segments=["A", "B"],
                    pairs=[("A", "B", -4.108)], mode="chi", aii=25.0)
    lut = build_a_lookup(aij)
    expected = 25.0 + (-4.108) / 0.306        # = 25 + 3.268·(-4.108)
    assert lut[_pair_key("A", "B")] == pytest.approx(expected)
    # ユーザー指定式 aii + 3.268·χ とも一致 (1/0.306 ≈ 3.268)
    assert lut[_pair_key("A", "B")] == pytest.approx(25.0 + 3.2679738 * -4.108, rel=1e-4)


def test_pair_key_order_insensitive():
    assert _pair_key("A", "B") == _pair_key("B", "A")
    assert _pair_key("A", "A") == frozenset({"A"})


# ---------------------------------------------------------------------------
# match_aij_to_pairs
# ---------------------------------------------------------------------------


def test_match_order_insensitive_and_unmatched():
    aij = AijMatrix(segments=["A", "B", "W"],
                    pairs=[("A", "B", 40.0), ("A", "A", 25.0)],
                    mode="a", aii=25.0)
    # UDF 側は B-A の順 (順不同で照合できること) + aij に無い W-W
    udf_pairs = [(0, "B", "A"), (1, "A", "A"), (2, "W", "W")]
    assignments, unmatched = match_aij_to_pairs(udf_pairs, aij)
    assert (0, 40.0) in assignments          # B-A ↔ A-B 照合
    assert (1, 25.0) in assignments
    assert unmatched == [(2, "W", "W")]      # W-W は aij.dat に無し


# ---------------------------------------------------------------------------
# assign_aij_to_udf (モック UDFManager)
# ---------------------------------------------------------------------------


class _MockUDF:
    """Interactions.Pair_Interaction[] を持つ最小モック UDFManager。"""

    _instances = []

    def __init__(self, path):
        self.path = path
        # (Site1, Site2, Potential_Type, DPD.a)
        self.pairs = [
            ["A", "B", "DPD", 0.0],
            ["A", "A", "DPD", 0.0],
            ["X", "Y", "Coulomb", 0.0],   # 非 DPD (only_dpd でスキップされる)
        ]
        self.written = None
        _MockUDF._instances.append(self)

    def jump(self, rec):
        self.rec = rec

    def size(self, path):
        assert path == "Interactions.Pair_Interaction[]"
        return len(self.pairs)

    def get(self, path, idx):
        k = idx[0]
        if path.endswith(".Potential_Type"):
            return self.pairs[k][2]
        if path.endswith(".Site1_Name"):
            return self.pairs[k][0]
        if path.endswith(".Site2_Name"):
            return self.pairs[k][1]
        raise KeyError(path)

    def put(self, value, path, idx):
        assert path.endswith(".DPD.a")
        # numpy が来ていないこと (float cast 済) を確認
        assert isinstance(value, float)
        self.pairs[idx[0]][3] = value

    def write(self, path):
        self.written = path


@pytest.fixture
def mock_udfmanager(monkeypatch):
    _MockUDF._instances.clear()
    mod = types.ModuleType("UDFManager")
    mod.UDFManager = _MockUDF
    monkeypatch.setitem(sys.modules, "UDFManager", mod)
    return _MockUDF


def test_assign_aij_to_udf_dpd_only(mock_udfmanager, tmp_path):
    from abmptools.cg.dpd.aij_assign import assign_aij_to_udf
    aij = AijMatrix(segments=["A", "B"],
                    pairs=[("A", "B", 40.0), ("A", "A", 25.0)],
                    mode="a", aii=25.0)
    out = tmp_path / "patched.udf"
    result = assign_aij_to_udf("dummy.udf", aij, out_path=str(out))

    inst = mock_udfmanager._instances[-1]
    # DPD ペアのみ (A-B, A-A) に a が入り、Coulomb はスキップ
    assert inst.pairs[0][3] == 40.0     # A-B
    assert inst.pairs[1][3] == 25.0     # A-A
    assert inst.pairs[2][3] == 0.0      # Coulomb (skip)
    assert result["total_pairs"] == 2
    assert result["skipped_non_dpd"] == 1
    assert len(result["assigned"]) == 2
    assert inst.written == str(out)


def test_assign_aij_to_udf_chi_mode(mock_udfmanager, tmp_path):
    from abmptools.cg.dpd.aij_assign import assign_aij_to_udf
    aij = AijMatrix(segments=["A", "B"],
                    pairs=[("A", "B", -4.108)], mode="chi", aii=25.0)
    out = tmp_path / "patched_chi.udf"
    assign_aij_to_udf("dummy.udf", aij, out_path=str(out))
    inst = mock_udfmanager._instances[-1]
    assert inst.pairs[0][3] == pytest.approx(25.0 + (-4.108) / 0.306)  # A-B chi→a
    assert inst.pairs[1][3] == 0.0   # A-A は aij に無し → 未割当


def test_assign_aij_unmatched_reported(mock_udfmanager, tmp_path):
    from abmptools.cg.dpd.aij_assign import assign_aij_to_udf
    aij = AijMatrix(segments=["A"], pairs=[("A", "A", 25.0)], mode="a", aii=25.0)
    result = assign_aij_to_udf("dummy.udf", aij, out_path=str(tmp_path / "o.udf"))
    # UDF の A-B は aij.dat に無い → unmatched
    assert any(s1 == "A" and s2 == "B" for _, s1, s2 in result["unmatched"])


# ---------------------------------------------------------------------------
# 実 cognac UDF への割り当て (OCTA 有り環境のみ)
# ---------------------------------------------------------------------------


@requires_cognac
def test_assign_aij_real_udf_roundtrip(tmp_path, cholesterol_cg, sample_aij_a_mode):
    """build-udf で実 UDF 生成 → assign-aij で特定 pair を書き換え → 反映確認。"""
    from abmptools.cg.dpd import CGDpdBuilder, AijMatrix
    from abmptools.cg.dpd.aij_assign import assign_aij_to_udf
    from UDFManager import UDFManager

    builder = CGDpdBuilder.from_files(
        monomer=cholesterol_cg["monomer"], aij=sample_aij_a_mode,
        calc_sett=cholesterol_cg["calc_sett"],
    )
    udf = builder.build_udf(tmp_path / "real.udf")

    # P0-P1 を 99.0 に上書きする aij
    new_aij = AijMatrix(segments=["P0", "P1"], pairs=[("P0", "P1", 99.0)],
                        mode="a", aii=25.0)
    out = tmp_path / "real_patched.udf"
    result = assign_aij_to_udf(udf, new_aij, out_path=out)
    assert len(result["assigned"]) == 1

    u = UDFManager(str(out))
    n = u.size("Interactions.Pair_Interaction[]")
    got = None
    for k in range(n):
        s1 = u.get(f"Interactions.Pair_Interaction[{k}].Site1_Name")
        s2 = u.get(f"Interactions.Pair_Interaction[{k}].Site2_Name")
        if {s1, s2} == {"P0", "P1"}:
            got = u.get(f"Interactions.Pair_Interaction[{k}].DPD.a")
    assert got == 99.0
