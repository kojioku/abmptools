# -*- coding: utf-8 -*-
"""
Byte / field-level equivalence tests between legacy gen_udf and gen_udf_v2.

Phase 2c-E (D-4-3): legacy `abmptools.udfcreate.udfcreate.gen_udf` と新規
`abmptools.udfcreate_v2.gen_udf_v2` を同一入力で実行し、出力 UDF が
**主要 fields で一致** することを確認する。byte-identical までは要求
しない (空白・順序・浮動小数点表現の差異が許容される)。これが pass
すれば fcewsmb の `use_gen_udf_v2=True` を default 化できる。

Test scope (現状):
  - gen_udf_v2 が cognac-readable な UDF を出すこと (UDFManager で開ける)
  - 主要 fields の round-trip 確認 (Time triplet / Cell_Size / Atom_Type[]
    / Set_of_Molecules / Structure.Position / charges)

**legacy gen_udf との直接比較は今後の拡張**:
gen_udf の入力 (`udf_param`, `som_param`) は ~100 個のフィールドを要求し、
fcewsmb 内部で構築される。そのため "stand-alone fixture" を作ること自体が
困難 (per-atom force field row, ffname の二重リスト, etc.)。

実用的な比較手順 (integration test):
1. fcews-manybody で `use_gen_udf_v2=False` で 2_mkudf 実行 → legacy npt.udf
2. 同じ入力で `use_gen_udf_v2=True` で再実行 → v2 npt.udf
3. 両 UDF を UDFManager で開いて key fields を比較

このパターンは fcews-manybody/tests/integration/ の oligomer / methanol-acetone
fixtures で実装する方が筋が良い。本ファイルは v2 単独の sanity test に専念
する。

UDFManager がない環境では skip。
"""
from __future__ import annotations

import os

import pytest


udfm = pytest.importorskip("UDFManager")


def _build_minimal_acetone_inputs():
    """gen_udf_v2 用の最小入力 (acetone monomer × 1 instance)。

    legacy gen_udf も理論上は通るが、入力形式の細かい規約 (atom[i][1] が
    str expected 等) で TypeError が出やすい。本 fixture は **v2 が正しく
    UDF を書けるか** の自己完結 sanity test 用。

    Returns:
        (udf_param, som_param) tuple.
    """
    cellsize = [2.5, 2.5, 2.5]
    ljparam = [
        ["c", 12.011, 0.0, 0.36, 0.34, 0.0, 0.0],
        ["o", 15.999, 0.0, 0.65, 0.30, 0.0, 0.0],
    ]
    bondparam = [["c", "o", 4.5e5, 0.122]]
    angleparam = []
    torsionparam = []
    atom_list = [["C1", 0, "c"], ["O1", 0, "o"]]
    totalstep = 5000
    outstep = 50
    totalmass = 28.01
    algo = ["NVT_Nose_Hoover"]
    poslist = [
        [
            [[0.5, 0.5, 0.5], [0.6, 0.5, 0.5]],   # 1 mol with 2 atoms
        ],
        [],
    ]
    udf_param = [
        cellsize, ljparam, bondparam, angleparam, torsionparam,
        atom_list, totalstep, outstep, totalmass, algo, poslist,
    ]
    fname    = ["co.xyz", ""]
    atom     = [[["C", 0, "c"], ["O", 0, "o"]], []]
    ffname   = [["c", "o"], []]
    molname  = ["CO", "CO2"]
    bondff   = [["c-o"], []]
    bond     = [[[1, 2]], []]
    angleff  = [[], []]
    angle    = [[], []]
    torsionff = [[], []]
    torsion  = [[], []]
    chg      = [[0.5, -0.5], []]
    som_param = [fname, atom, ffname, molname,
                 bondff, bond, angleff, angle,
                 torsionff, torsion, chg]
    return udf_param, som_param


@pytest.fixture
def v2_output_udf(tmp_path):
    """gen_udf_v2 の出力 UDF を返す fixture。"""
    from abmptools.udfcreate_v2 import gen_udf_v2

    udf_param, som_param = _build_minimal_acetone_inputs()
    out_path = str(tmp_path / "v2.udf")
    gen_udf_v2(udf_param, out_path, som_param)
    assert os.path.isfile(out_path)
    return out_path


# ---------------------------------------------------------------------------
# v2 single-writer sanity tests (legacy 比較は別途 fcewsmb 経由で)
# ---------------------------------------------------------------------------


def test_v2_output_opens_with_udfmanager(v2_output_udf):
    """gen_udf_v2 出力が UDFManager で開ける。"""
    u = udfm.UDFManager(v2_output_udf)
    u.jump(-1)


def test_v2_output_simulation_time(v2_output_udf):
    """Simulation_Conditions.Time triplet が入力値で round-trip。"""
    u = udfm.UDFManager(v2_output_udf); u.jump(-1)
    base = "Simulation_Conditions.Dynamics_Conditions.Time"
    assert u.get(f"{base}.Total_Steps") == 5000
    assert u.get(f"{base}.Output_Interval_Steps") == 50


def test_v2_output_cell_size(v2_output_udf):
    """Initial_Unit_Cell.Cell_Size が round-trip。"""
    u = udfm.UDFManager(v2_output_udf); u.jump(-1)
    base = "Initial_Structure.Initial_Unit_Cell.Cell_Size"
    for axis in ("a", "b", "c"):
        assert u.get(f"{base}.{axis}") == pytest.approx(2.5, rel=1e-9)


def test_v2_output_atom_types(v2_output_udf):
    """Molecular_Attributes.Atom_Type[] が入力 ljparam と一致。"""
    u = udfm.UDFManager(v2_output_udf); u.jump(-1)
    n = u.size("Molecular_Attributes.Atom_Type[]")
    assert n == 2
    pairs = {
        (u.get("Molecular_Attributes.Atom_Type[].Name", [i]),
         u.get("Molecular_Attributes.Atom_Type[].Mass", [i]))
        for i in range(n)
    }
    assert pairs == {("c", 12.011), ("o", 15.999)}


def test_v2_output_set_of_molecules(v2_output_udf):
    """Set_of_Molecules.molecule[] count + 一部 atom name が round-trip。"""
    u = udfm.UDFManager(v2_output_udf); u.jump(-1)
    n_mol = u.size("Set_of_Molecules.molecule[]")
    assert n_mol == 1   # 1 instance of 'CO'
    assert u.get("Set_of_Molecules.molecule[].Mol_Name", [0]) == "CO"
    n_atom = u.size("Set_of_Molecules.molecule[].atom[]", [0])
    assert n_atom == 2


def test_v2_output_positions_round_trip(v2_output_udf):
    """Structure.Position.mol[].atom[].{x,y,z} が入力 poslist と一致 (nm)。"""
    u = udfm.UDFManager(v2_output_udf); u.jump(-1)
    # 入力は [nm]、UDFManager.get("[nm]") で同じ単位で読み戻す
    # mol[0].atom[0]: [0.5, 0.5, 0.5] nm
    for axis, expected in zip(("x", "y", "z"), (0.5, 0.5, 0.5)):
        v = u.get(f"Structure.Position.mol[].atom[].{axis}",
                  [0, 0], "[nm]")
        assert v == pytest.approx(expected, rel=1e-6)
    # mol[0].atom[1]: [0.6, 0.5, 0.5] nm
    for axis, expected in zip(("x", "y", "z"), (0.6, 0.5, 0.5)):
        v = u.get(f"Structure.Position.mol[].atom[].{axis}",
                  [0, 1], "[nm]")
        assert v == pytest.approx(expected, rel=1e-6)


def test_v2_output_charges_via_es_element(v2_output_udf):
    """electrostatic_Site[].ES_Element = charge × e2Q (gen_udf と同係数)。"""
    from abmptools.udfcreate_v2 import _E2Q

    u = udfm.UDFManager(v2_output_udf); u.jump(-1)
    es_0 = u.get(
        "Set_of_Molecules.molecule[].electrostatic_Site[].ES_Element",
        [0, 0],
    )
    es_1 = u.get(
        "Set_of_Molecules.molecule[].electrostatic_Site[].ES_Element",
        [0, 1],
    )
    assert es_0 == pytest.approx(0.5 * _E2Q, rel=1e-6)
    assert es_1 == pytest.approx(-0.5 * _E2Q, rel=1e-6)


# ---------------------------------------------------------------------------
# Legacy gen_udf vs v2 — 入力 shape の制約で skip-ed (将来 expand)
# ---------------------------------------------------------------------------


def test_legacy_vs_v2_via_fcewsmb_integration():
    """legacy gen_udf vs gen_udf_v2 の field-level 比較は fcews-manybody
    integration test 経由で行う。

    rationale: gen_udf の入力 (`udf_param`, `som_param`) は per-atom
    force field row や ffname の二重リスト等 ~100 個の field を持ち、
    standalone fixture を作ると入力 shape の細かい規約違反で TypeError
    が頻発する。fcewsmb の `_setupudf_legacy` 内で実際に組まれる入力を
    使う方が信頼性が高い。

    したがって本テストは pytest.skip マーカーで future-work プレースホルダ
    として残置する。具体的な実装は:
      fcews-manybody/tests/integration/test_use_gen_udf_v2.py で
      use_gen_udf_v2=False と True で 2_mkudf を回し、出力 npt.udf の
      key fields を比較する。
    """
    pytest.skip(
        "moved to fcews-manybody/tests/integration/ — see test docstring"
    )
