# -*- coding: utf-8 -*-
"""abmptools.udfcharge の単体テスト。

フィクスチャは gro2udf の default_template.udf (空 Set_of_Molecules を持つ
COGNAC schema base) から **プログラムで生成**する:
- 単分子 UDF (電荷あり)  = template
- バルク UDF (同名分子 N 個、 電荷なし) = 割り当て先

UDFManager (OCTA 同梱) が無い環境では skip。
"""
from __future__ import annotations

import shutil
from pathlib import Path

import pytest

pytest.importorskip("UDFManager", reason="OCTA UDFManager not installed")

from abmptools.udfcharge import (  # noqa: E402
    CHARGE_UNIT,
    assign_charges_to_bulk,
    read_molecule_charges,
)

# 4-atom methyl-like 分子 (電荷和 0)
MOL_NAME = "MOL"
ATOM_TYPES = ["c3", "h1", "h1", "h1"]
ATOM_ELEMS = ["C", "H", "H", "H"]
CHARGES = [-0.6, 0.2, 0.2, 0.2]


def _template_udf_path() -> Path:
    import abmptools.gro2udf as g
    return Path(g.__file__).parent / "default_template.udf"


def _build_udf(out_path: Path, *, n_copies: int, with_charges: bool,
               mol_name=MOL_NAME, atom_types=ATOM_TYPES, atom_elems=ATOM_ELEMS,
               charges=CHARGES) -> Path:
    """default_template から n_copies 個の同名分子を持つ UDF を作る。"""
    from UDFManager import UDFManager

    shutil.copy(_template_udf_path(), out_path)
    u = UDFManager(str(out_path))
    u.jump(-1)
    n_atoms = len(atom_types)
    for imol in range(n_copies):
        u.put(mol_name, "Set_of_Molecules.molecule[].Mol_Name", [imol])
        for i in range(n_atoms):
            u.put(i, "Set_of_Molecules.molecule[].atom[].Atom_ID", [imol, i])
            u.put(atom_elems[i], "Set_of_Molecules.molecule[].atom[].Atom_Name", [imol, i])
            u.put(atom_types[i], "Set_of_Molecules.molecule[].atom[].Atom_Type_Name", [imol, i])
            if with_charges:
                u.put("POINT_CHARGE",
                      "Set_of_Molecules.molecule[].electrostatic_Site[].Type_Name", [imol, i])
                u.put(float(charges[i] * CHARGE_UNIT),
                      "Set_of_Molecules.molecule[].electrostatic_Site[].ES_Element", [imol, i])
                u.put(i,
                      "Set_of_Molecules.molecule[].electrostatic_Site[].atom[]", [imol, i, 0])
    u.write(str(out_path))
    return out_path


def _read_all_charges(udf_path: Path, mol_name=MOL_NAME):
    """UDF 内 mol_name 全分子の per-atom 電荷 [e] を返す (list[list[float]])。"""
    from UDFManager import UDFManager

    u = UDFManager(str(udf_path))
    u.jump(-1)
    nmol = u.size("Set_of_Molecules.molecule[]")
    out = []
    for i in range(nmol):
        if u.get("Set_of_Molecules.molecule[].Mol_Name", [i]) != mol_name:
            continue
        na = u.size("Set_of_Molecules.molecule[].atom[]", [i])
        ne = u.size("Set_of_Molecules.molecule[].electrostatic_Site[]", [i])
        q = [0.0] * na
        for j in range(ne or 0):
            es = u.get("Set_of_Molecules.molecule[].electrostatic_Site[].ES_Element", [i, j])
            q[j] = es / CHARGE_UNIT
        out.append(q)
    return out


# --------------------------------------------------------------------------
# fixtures
# --------------------------------------------------------------------------

@pytest.fixture
def single_udf(tmp_path):
    return _build_udf(tmp_path / "single.udf", n_copies=1, with_charges=True)


@pytest.fixture
def bulk_udf(tmp_path):
    return _build_udf(tmp_path / "bulk.udf", n_copies=5, with_charges=False)


# --------------------------------------------------------------------------
# tests
# --------------------------------------------------------------------------

def test_read_template(single_udf):
    tmpl = read_molecule_charges(single_udf)
    assert tmpl.mol_name == MOL_NAME
    assert tmpl.n_atoms == 4
    assert tmpl.atom_type_names == ATOM_TYPES
    assert tmpl.atom_names == ATOM_ELEMS
    for got, exp in zip(tmpl.charges, CHARGES):
        assert abs(got - exp) < 1e-6
    assert abs(tmpl.net_charge) < 1e-6


def test_bulk_starts_uncharged(bulk_udf):
    all_q = _read_all_charges(bulk_udf)
    assert len(all_q) == 5
    for q in all_q:
        assert all(abs(c) < 1e-9 for c in q)  # 全原子 0


def test_assign_charges_to_all_bulk_molecules(single_udf, bulk_udf, tmp_path):
    tmpl = read_molecule_charges(single_udf)
    out = tmp_path / "bulk_charged.udf"
    res = assign_charges_to_bulk(bulk_udf, tmpl, out)

    assert res.n_molecules_assigned == 5
    assert res.n_molecules_total == 5
    assert res.n_atoms_per_mol == 4
    assert Path(res.out_path) == out

    all_q = _read_all_charges(out)
    assert len(all_q) == 5
    for q in all_q:  # 全分子が template と一致
        for got, exp in zip(q, CHARGES):
            assert abs(got - exp) < 1e-6


def test_input_bulk_not_modified(single_udf, bulk_udf, tmp_path):
    """出力は別ファイル、 入力 bulk は無改変 (電荷 0 のまま)。"""
    tmpl = read_molecule_charges(single_udf)
    assign_charges_to_bulk(bulk_udf, tmpl, tmp_path / "out.udf")
    for q in _read_all_charges(bulk_udf):
        assert all(abs(c) < 1e-9 for c in q)


def test_default_out_path(single_udf, bulk_udf):
    """out_path 省略時は <bulk>_charged.udf。"""
    tmpl = read_molecule_charges(single_udf)
    res = assign_charges_to_bulk(bulk_udf, tmpl, None)
    assert res.out_path.endswith("bulk_charged.udf")
    assert Path(res.out_path).is_file()


def test_atom_count_mismatch_raises(single_udf, tmp_path):
    """atom 数の違う同名分子へは strict で例外。"""
    tmpl = read_molecule_charges(single_udf)  # 4 atoms
    bad = _build_udf(tmp_path / "bad.udf", n_copies=2, with_charges=False,
                     atom_types=["c3", "h1"], atom_elems=["C", "H"])  # 2 atoms
    with pytest.raises(ValueError, match="atom 数"):
        assign_charges_to_bulk(bad, tmpl, tmp_path / "x.udf")


def test_type_mismatch_raises(single_udf, tmp_path):
    """Atom_Type_Name 列の違う同名分子へは strict で例外。"""
    tmpl = read_molecule_charges(single_udf)
    bad = _build_udf(tmp_path / "bad2.udf", n_copies=2, with_charges=False,
                     atom_types=["XX", "h1", "h1", "h1"])
    with pytest.raises(ValueError, match="Atom_Type_Name"):
        assign_charges_to_bulk(bad, tmpl, tmp_path / "x.udf")


def test_no_matching_molecule_raises(single_udf, tmp_path):
    """Mol_Name 不一致なら strict で例外。"""
    tmpl = read_molecule_charges(single_udf)
    other = _build_udf(tmp_path / "other.udf", n_copies=3, with_charges=False,
                       mol_name="OTHER")
    with pytest.raises(ValueError, match="一致する分子"):
        assign_charges_to_bulk(other, tmpl, tmp_path / "x.udf")


def test_cli(single_udf, bulk_udf, tmp_path):
    from abmptools.udfcharge.__main__ import main

    out = tmp_path / "cli_out.udf"
    rc = main(["--template", str(single_udf), "--bulk", str(bulk_udf),
               "--out", str(out)])
    assert rc == 0
    assert out.is_file()
    for q in _read_all_charges(out):
        for got, exp in zip(q, CHARGES):
            assert abs(got - exp) < 1e-6
