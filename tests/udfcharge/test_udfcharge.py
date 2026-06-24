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
# 分子内座標 [nm] (C + 3H、 適当な四面体配置)
MOL_XYZ = [
    (0.000, 0.000, 0.000),
    (0.063, 0.063, 0.063),
    (-0.063, -0.063, 0.063),
    (0.063, -0.063, -0.063),
]
# 結合トポロジー (atom 0-based): C-H ×3
MOL_BONDS = [(0, 1), (0, 2), (0, 3)]
MOL_ANGLES = [(1, 0, 2), (1, 0, 3), (2, 0, 3)]


def _template_udf_path() -> Path:
    import abmptools.gro2udf as g
    return Path(g.__file__).parent / "default_template.udf"


def _build_udf(out_path: Path, *, n_copies: int, with_charges: bool,
               mol_name=MOL_NAME, atom_types=ATOM_TYPES, atom_elems=ATOM_ELEMS,
               charges=CHARGES, with_coords: bool = False,
               with_bonds: bool = False) -> Path:
    """default_template から n_copies 個の同名分子を持つ UDF を作る。

    with_coords=True で各分子を 1.0 nm 間隔で並べた座標 + cell を dynamic record に、
    with_bonds=True で bond/angle を static record に書き込む (charge 転写が座標・
    トポロジーを壊さないことの検証用)。
    """
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
        if with_bonds:
            for k, (a1, a2) in enumerate(MOL_BONDS):
                u.put("", "Set_of_Molecules.molecule[].bond[].Potential_Name", [imol, k])
                u.put(a1, "Set_of_Molecules.molecule[].bond[].atom1", [imol, k])
                u.put(a2, "Set_of_Molecules.molecule[].bond[].atom2", [imol, k])
                u.put(1.0, "Set_of_Molecules.molecule[].bond[].Order", [imol, k])
            for k, (a1, a2, a3) in enumerate(MOL_ANGLES):
                u.put("", "Set_of_Molecules.molecule[].angle[].Potential_Name", [imol, k])
                u.put(a1, "Set_of_Molecules.molecule[].angle[].atom1", [imol, k])
                u.put(a2, "Set_of_Molecules.molecule[].angle[].atom2", [imol, k])
                u.put(a3, "Set_of_Molecules.molecule[].angle[].atom3", [imol, k])

    if with_coords:
        u.eraseRecord(0, u.totalRecord())
        u.newRecord()
        for imol in range(n_copies):
            for i in range(n_atoms):
                bx, by, bz = MOL_XYZ[i % len(MOL_XYZ)]
                u.put(float(imol) + bx, "Structure.Position.mol[].atom[].x", [imol, i], "[nm]")
                u.put(by, "Structure.Position.mol[].atom[].y", [imol, i], "[nm]")
                u.put(bz, "Structure.Position.mol[].atom[].z", [imol, i], "[nm]")
        for ax in ("a", "b", "c"):
            u.put(float(n_copies + 2), f"Structure.Unit_Cell.Cell_Size.{ax}", "[nm]")
        for ax in ("alpha", "beta", "gamma"):
            u.put(90.0, f"Structure.Unit_Cell.Cell_Size.{ax}")
        u.put(0, "Steps")

    u.write(str(out_path))
    return out_path


def _read_all_positions(udf_path: Path, mol_name=MOL_NAME):
    """record の per-atom 座標 [nm] を返す (list[list[(x,y,z)]])。"""
    from UDFManager import UDFManager

    u = UDFManager(str(udf_path))
    u.jump(u.totalRecord() - 1)
    nmol = u.size("Set_of_Molecules.molecule[]")
    out = []
    for i in range(nmol):
        u.jump(-1)
        if u.get("Set_of_Molecules.molecule[].Mol_Name", [i]) != mol_name:
            continue
        u.jump(u.totalRecord() - 1)
        na = u.size("Structure.Position.mol[].atom[]", [i])
        out.append([
            tuple(u.get(f"Structure.Position.mol[].atom[].{ax}", [i, j], "[nm]") for ax in "xyz")
            for j in range(na or 0)
        ])
    return out


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


def test_coordinates_preserved(single_udf, tmp_path):
    """電荷転写後も Structure.Position の座標が保持される (record を壊さない)。"""
    tmpl = read_molecule_charges(single_udf)
    bulk = _build_udf(tmp_path / "bulk_xyz.udf", n_copies=4, with_charges=False,
                      with_coords=True)
    before = _read_all_positions(bulk)
    assert len(before) == 4
    assert all(len(p) == 4 for p in before)  # 4 atoms each

    out = tmp_path / "bulk_xyz_charged.udf"
    assign_charges_to_bulk(bulk, tmpl, out)

    after = _read_all_positions(out)
    assert after == before  # 座標は完全一致 (無改変)

    # かつ電荷は入っている
    for q in _read_all_charges(out):
        for got, exp in zip(q, CHARGES):
            assert abs(got - exp) < 1e-6


def _read_all_bonds(udf_path: Path, mol_name=MOL_NAME):
    """static record の per-molecule bond (atom1,atom2) を返す。"""
    from UDFManager import UDFManager

    u = UDFManager(str(udf_path))
    u.jump(-1)
    nmol = u.size("Set_of_Molecules.molecule[]")
    out = []
    for i in range(nmol):
        if u.get("Set_of_Molecules.molecule[].Mol_Name", [i]) != mol_name:
            continue
        nb = u.size("Set_of_Molecules.molecule[].bond[]", [i]) or 0
        out.append([
            (u.get("Set_of_Molecules.molecule[].bond[].atom1", [i, k]),
             u.get("Set_of_Molecules.molecule[].bond[].atom2", [i, k]))
            for k in range(nb)
        ])
    return out


def test_topology_preserved(single_udf, tmp_path):
    """電荷転写後も bond/angle トポロジーが保持される (static record を壊さない)。"""
    tmpl = read_molecule_charges(single_udf)
    bulk = _build_udf(tmp_path / "bulk_topo.udf", n_copies=3, with_charges=False,
                      with_coords=True, with_bonds=True)
    before = _read_all_bonds(bulk)
    assert len(before) == 3
    assert all(b == MOL_BONDS for b in before)

    out = tmp_path / "bulk_topo_charged.udf"
    assign_charges_to_bulk(bulk, tmpl, out)

    after = _read_all_bonds(out)
    assert after == before                       # bond 無改変
    assert all(len(p) == 4 for p in _read_all_positions(out))  # 座標も保持
    for q in _read_all_charges(out):             # 電荷は入っている
        for got, exp in zip(q, CHARGES):
            assert abs(got - exp) < 1e-6


def test_cli(single_udf, bulk_udf, tmp_path):
    # 旧フラット呼び出し (transfer 省略) も後方互換で動く
    from abmptools.udfcharge.__main__ import main

    out = tmp_path / "cli_out.udf"
    rc = main(["--template", str(single_udf), "--bulk", str(bulk_udf),
               "--out", str(out)])
    assert rc == 0
    assert out.is_file()
    for q in _read_all_charges(out):
        for got, exp in zip(q, CHARGES):
            assert abs(got - exp) < 1e-6


# --------------------------------------------------------------------------
# restore_formal_charge (中和の逆変換)
# --------------------------------------------------------------------------

# 元電荷 A (Σ = +1)、 5 atoms
RESTORE_A = [0.50, 0.30, 0.40, -0.10, -0.10]
RESTORE_TYPES = ["c3", "c3", "c3", "c3", "c3"]
RESTORE_ELEMS = ["C", "C", "C", "C", "C"]


def _forward_neutralize(A, S):
    """A (Σ=S) を |A| 比例で中和して B (Σ≈0) を返す (機能の forward 側)。"""
    sa = sum(abs(a) for a in A)
    return [a - S * abs(a) / sa for a in A]


def _neutral_udf(tmp_path, A, S, name="neutral"):
    B = _forward_neutralize(A, S)
    assert abs(sum(B)) < 1e-9
    return _build_udf(tmp_path / f"{name}.udf", n_copies=1, with_charges=True,
                      mol_name="ION", atom_types=RESTORE_TYPES,
                      atom_elems=RESTORE_ELEMS, charges=B), B


def test_restore_recovers_original_charges(tmp_path):
    from abmptools.udfcharge import restore_formal_charge

    S = 1
    udf, B = _neutral_udf(tmp_path, RESTORE_A, S)
    out = tmp_path / "restored.udf"
    res = restore_formal_charge(udf, S, out)

    assert res.formal_charge == 1
    assert abs(res.input_total) < 1e-9          # 入力は中和済み
    assert abs(res.output_total - 1.0) < 1e-6   # 出力は形式電荷
    a_hat = read_molecule_charges(out).charges
    for got, exp in zip(a_hat, RESTORE_A):      # 元電荷を復元
        assert abs(got - exp) < 1e-6


def test_restore_negative_formal_charge(tmp_path):
    from abmptools.udfcharge import restore_formal_charge

    A = [-0.6, -0.3, -0.4, 0.15, 0.15]          # Σ = -1
    S = -1
    udf, _ = _neutral_udf(tmp_path, A, S, name="anion")
    res = restore_formal_charge(udf, S, tmp_path / "anion_restored.udf")
    assert abs(res.output_total - (-1.0)) < 1e-6
    a_hat = read_molecule_charges(res.out_path).charges
    for got, exp in zip(a_hat, A):
        assert abs(got - exp) < 1e-6


def test_restore_default_out_path(tmp_path):
    from abmptools.udfcharge import restore_formal_charge

    udf, _ = _neutral_udf(tmp_path, RESTORE_A, 1)
    res = restore_formal_charge(udf, 1, None)
    assert res.out_path.endswith("neutral_q+1.udf")
    assert Path(res.out_path).is_file()


def test_restore_zero_formal_charge_is_noop(tmp_path):
    """形式電荷 0 (既に中和) は電荷を変えない。"""
    from abmptools.udfcharge import restore_formal_charge

    udf, B = _neutral_udf(tmp_path, RESTORE_A, 1)   # B は Σ=0
    res = restore_formal_charge(udf, 0, tmp_path / "z.udf")
    assert abs(res.output_total) < 1e-6
    a_hat = read_molecule_charges(res.out_path).charges
    for got, exp in zip(a_hat, B):
        assert abs(got - exp) < 1e-6


def test_restore_preserves_coords_and_bonds(tmp_path):
    """restore も electrostatic_Site のみ更新し座標・結合を壊さない。"""
    from abmptools.udfcharge import restore_formal_charge

    B = _forward_neutralize(RESTORE_A, 1)
    udf = _build_udf(tmp_path / "full.udf", n_copies=1, with_charges=True,
                     mol_name="ION", atom_types=RESTORE_TYPES,
                     atom_elems=RESTORE_ELEMS, charges=B,
                     with_coords=True, with_bonds=True)
    before_pos = _read_all_positions(udf, mol_name="ION")
    before_bonds = _read_all_bonds(udf, mol_name="ION")
    out = tmp_path / "full_restored.udf"
    restore_formal_charge(udf, 1, out)
    assert _read_all_positions(out, mol_name="ION") == before_pos
    assert _read_all_bonds(out, mol_name="ION") == before_bonds


def test_cli_restore(tmp_path):
    from abmptools.udfcharge.__main__ import main

    udf, _ = _neutral_udf(tmp_path, RESTORE_A, 1)
    out = tmp_path / "cli_restored.udf"
    rc = main(["restore", "--udf", str(udf), "--formal-charge", "1",
               "--out", str(out)])
    assert rc == 0
    a_hat = read_molecule_charges(out).charges
    for got, exp in zip(a_hat, RESTORE_A):
        assert abs(got - exp) < 1e-6


# --- uniform (均等分配) mode ---

def _forward_uniform(A, S):
    """A (Σ=S) を均等に中和して B (Σ≈0) を返す: B_i = A_i - S/N。"""
    shift = S / len(A)
    return [a - shift for a in A]


def _uniform_neutral_udf(tmp_path, A, S, name="uniform"):
    B = _forward_uniform(A, S)
    assert abs(sum(B)) < 1e-9
    types = ["c3"] * len(A)
    elems = ["C"] * len(A)
    return _build_udf(tmp_path / f"{name}.udf", n_copies=1, with_charges=True,
                      mol_name="ION", atom_types=types, atom_elems=elems,
                      charges=B), B


def test_restore_uniform_recovers_exactly(tmp_path):
    """uniform で中和した電荷は mode='uniform' で厳密復元 (A_i = B_i + S/N)。"""
    from abmptools.udfcharge import restore_formal_charge

    S = 1
    udf, _ = _uniform_neutral_udf(tmp_path, RESTORE_A, S)
    res = restore_formal_charge(udf, S, tmp_path / "u_restored.udf", mode="uniform")

    assert res.mode == "uniform"
    assert res.lam is None
    assert abs(res.shift - S / len(RESTORE_A)) < 1e-12
    assert abs(res.output_total - 1.0) < 1e-9
    a_hat = read_molecule_charges(res.out_path).charges
    for got, exp in zip(a_hat, RESTORE_A):
        assert abs(got - exp) < 1e-9        # 均等分配は厳密 (符号問題なし)


def test_uniform_and_proportional_disagree_on_uniform_data(tmp_path):
    """uniform で中和したデータを proportional で復元すると元と一致しない
    (= mode を中和方法に合わせる必要があることの確認)。"""
    from abmptools.udfcharge import restore_formal_charge

    S = 1
    udf, _ = _uniform_neutral_udf(tmp_path, RESTORE_A, S)
    a_uni = read_molecule_charges(
        restore_formal_charge(udf, S, tmp_path / "a.udf", mode="uniform").out_path
    ).charges
    a_prop = read_molecule_charges(
        restore_formal_charge(udf, S, tmp_path / "b.udf", mode="proportional").out_path
    ).charges
    # uniform は元と一致、 proportional はズレる
    assert max(abs(g - e) for g, e in zip(a_uni, RESTORE_A)) < 1e-9
    assert max(abs(g - e) for g, e in zip(a_prop, RESTORE_A)) > 1e-3


def test_restore_unknown_mode_raises(tmp_path):
    from abmptools.udfcharge import restore_formal_charge

    udf, _ = _uniform_neutral_udf(tmp_path, RESTORE_A, 1)
    with pytest.raises(ValueError, match="mode"):
        restore_formal_charge(udf, 1, tmp_path / "x.udf", mode="bogus")


def test_cli_restore_uniform(tmp_path):
    from abmptools.udfcharge.__main__ import main

    udf, _ = _uniform_neutral_udf(tmp_path, RESTORE_A, 1)
    out = tmp_path / "cli_u.udf"
    rc = main(["restore", "--udf", str(udf), "--formal-charge", "1",
               "--mode", "uniform", "--out", str(out)])
    assert rc == 0
    a_hat = read_molecule_charges(out).charges
    for got, exp in zip(a_hat, RESTORE_A):
        assert abs(got - exp) < 1e-9
