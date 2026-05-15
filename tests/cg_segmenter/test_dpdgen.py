# -*- coding: utf-8 -*-
"""tests/cg_segmenter/test_dpdgen.py — DPDgen exporter の検証。"""
from __future__ import annotations

from pathlib import Path

import pytest


@pytest.fixture(scope="session")
def naphthalene_segmenter(tmp_path_factory):
    """fused 2 ring naphthalene。ring-ring 0.6 / 200 を検証する用。"""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from abmptools.fragmenter.cg_segmenter import CGSegmenter, CGSegmenterConfig

    pdb = tmp_path_factory.mktemp("dpd_naph") / "naph.pdb"
    mol = Chem.AddHs(Chem.MolFromSmiles("c1ccc2ccccc2c1"))
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol)
    Chem.MolToPDBFile(mol, str(pdb))

    config = CGSegmenterConfig(
        pdb_path=str(pdb),
        output_dir=str(tmp_path_factory.mktemp("dpd_naph_out")),
        target_mw=200,
    )
    return CGSegmenter.from_pdb(config)


@pytest.fixture(scope="session")
def cyclohexyl_octyl_segmenter(tmp_path_factory):
    """ring + chain。ring-chain は 0.86 / 50 を検証する用。"""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from abmptools.fragmenter.cg_segmenter import CGSegmenter, CGSegmenterConfig

    pdb = tmp_path_factory.mktemp("dpd_co") / "co.pdb"
    mol = Chem.AddHs(Chem.MolFromSmiles("CCCCCCCCC1CCCCC1"))
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol)
    Chem.MolToPDBFile(mol, str(pdb))

    config = CGSegmenterConfig(
        pdb_path=str(pdb),
        output_dir=str(tmp_path_factory.mktemp("dpd_co_out")),
        target_mw=200,
    )
    return CGSegmenter.from_pdb(config)


def test_export_dpdgen_writes_two_files(naphthalene_segmenter, tmp_path):
    """monomer file + calc_sett file が生成される。"""
    sg = naphthalene_segmenter
    mono, calc = sg.export_dpdgen(
        output_dir=str(tmp_path),
        monomer_name="naph",
        box_size=(10.0, 10.0, 10.0),
    )
    assert Path(mono).exists()
    assert Path(calc).exists()
    assert Path(mono).name == "naph_monomer"
    assert Path(calc).name == "naph_calc_sett"


def test_ring_ring_bond_is_0_60_stiffness_200(naphthalene_segmenter, tmp_path):
    """fused naphthalene の ring-ring bond は dist=0.6, stiff=200。"""
    sg = naphthalene_segmenter
    mono, _ = sg.export_dpdgen(output_dir=str(tmp_path), monomer_name="naph")
    content = Path(mono).read_text()
    # bond12h に [1, 0, 1, 0.6, 200, 0] が含まれる
    assert "0.6, 200" in content, f"expected ring-ring 0.6/200 in {content[:500]}"


def test_chain_chain_bond_is_0_86_stiffness_50(cyclohexyl_octyl_segmenter, tmp_path):
    """ring + chain では ring-chain bond は 0.86/50 (chain default)。"""
    sg = cyclohexyl_octyl_segmenter
    mono, _ = sg.export_dpdgen(output_dir=str(tmp_path), monomer_name="co")
    content = Path(mono).read_text()
    # ring と chain のどちらかが non-ring kind を持つなら 0.86 / 50
    assert "0.86, 50" in content, f"expected ring-chain 0.86/50 in {content[:500]}"


def test_calc_sett_includes_monomer_filename(naphthalene_segmenter, tmp_path):
    sg = naphthalene_segmenter
    _, calc = sg.export_dpdgen(
        output_dir=str(tmp_path), monomer_name="myname", box_size=(15, 15, 15),
    )
    content = Path(calc).read_text()
    assert "monomer_file = 'myname_monomer'" in content
    assert "'x': 15" in content


def test_dpdgen_exec_loads_in_python(naphthalene_segmenter, tmp_path):
    """生成された monomer / calc_sett は Python script として valid (exec できる)。"""
    sg = naphthalene_segmenter
    mono, calc = sg.export_dpdgen(output_dir=str(tmp_path), monomer_name="naph")
    # monomer: bond12, bond12h, bond13_150, bond14_150, angle13, angle13data を定義
    ns: dict = {}
    exec(Path(mono).read_text(), ns)
    for k in ("bond12", "bond12h", "bond13_150", "bond13_150h",
              "bond14_150", "bond14_150h", "angle13", "angle13data"):
        assert k in ns, f"{k} missing in monomer"
    assert isinstance(ns["bond12"], list)
    assert all(isinstance(b, list) and len(b) == 2 for b in ns["bond12"])
    assert all(len(b) == 6 for b in ns["bond12h"])
    assert all(isinstance(t, list) and len(t) == 3 for t in ns["angle13"])
    # angle13data: [a, b, c, eq_余角, stiffness] = 5 要素
    assert all(len(a) == 5 for a in ns["angle13data"])
    # calc_sett: total_num_list, phys_param 等
    ns2: dict = {}
    exec(Path(calc).read_text(), ns2)
    assert "total_num_list" in ns2
    assert "phys_param" in ns2
    assert "monomer_file" in ns2


def test_cli_dpdgen_subcommand(cyclohexyl_octyl_segmenter, tmp_path):
    """`python -m abmptools.fragmenter.cg_segmenter dpdgen ...` が動く。"""
    from abmptools.fragmenter.cg_segmenter.__main__ import main
    out = tmp_path / "cli_out"
    code = main([
        "dpdgen",
        "--pdb", cyclohexyl_octyl_segmenter.config.pdb_path,
        "--output-dir", str(out),
        "--monomer-name", "test",
        "--box", "8.0",
    ])
    assert code == 0
    assert (out / "test_monomer").exists()
    assert (out / "test_calc_sett").exists()


# --- Fused-ring full-connect + bond13_180 / bond13_120 ----------------------

@pytest.fixture(scope="session")
def cholesterol_segmenter_for_dpd(tmp_path_factory):
    """sample/cg_segmenter/cholesterol_rdkit.pdb から 4 ring fused 構造を読む。"""
    import os
    from abmptools.fragmenter.cg_segmenter import CGSegmenter, CGSegmenterConfig

    pdb = os.path.abspath("sample/cg_segmenter/cholesterol_rdkit.pdb")
    if not Path(pdb).exists():
        pytest.skip("sample/cg_segmenter/cholesterol_rdkit.pdb not found")
    config = CGSegmenterConfig(
        pdb_path=pdb,
        output_dir=str(tmp_path_factory.mktemp("dpd_chol_out")),
        target_mw=200,
    )
    return CGSegmenter.from_pdb(config)


@pytest.fixture(scope="session")
def decane_segmenter(tmp_path_factory):
    """C10 alkane → 複数 chain segment で bond13_180 を検証する用。"""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from abmptools.fragmenter.cg_segmenter import CGSegmenter, CGSegmenterConfig

    pdb = tmp_path_factory.mktemp("dpd_dec") / "decane.pdb"
    mol = Chem.AddHs(Chem.MolFromSmiles("CCCCCCCCCC"))
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol)
    Chem.MolToPDBFile(mol, str(pdb))
    config = CGSegmenterConfig(
        pdb_path=str(pdb),
        output_dir=str(tmp_path_factory.mktemp("dpd_dec_out")),
        target_mw=45,
    )
    return CGSegmenter.from_pdb(config)


@pytest.fixture(scope="session")
def hexene_segmenter(tmp_path_factory):
    """hex-3-ene (cis): CC/C=C\\CC → 中央の C=C で bond13_120 を検証する用。"""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from abmptools.fragmenter.cg_segmenter import CGSegmenter, CGSegmenterConfig

    pdb = tmp_path_factory.mktemp("dpd_hex") / "hexene.pdb"
    mol = Chem.AddHs(Chem.MolFromSmiles("CC/C=C\\CC"))
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol)
    Chem.MolToPDBFile(mol, str(pdb))
    config = CGSegmenterConfig(
        pdb_path=str(pdb),
        output_dir=str(tmp_path_factory.mktemp("dpd_hex_out")),
        target_mw=30,
    )
    return CGSegmenter.from_pdb(config)


def test_cholesterol_path_hierarchy(cholesterol_segmenter_for_dpd, tmp_path):
    """4 fused rings (cholesterol) で path-based hierarchy を検証:
    bond12 = 隣接 (A-B, B-C, C-D, A-tail) — 4 pair
    bond13_150 = 1-skip ring (A-C, B-D) — 2 pair、 distance 1.661 / stiff 200
    bond14_150 = 2-skip ring (A-D) — 1 pair、 distance 2.502 / stiff 200
    """
    sg = cholesterol_segmenter_for_dpd
    mono, _ = sg.export_dpdgen(output_dir=str(tmp_path), monomer_name="chol")
    ns: dict = {}
    exec(Path(mono).read_text(), ns)
    bond12 = ns["bond12"]
    bond13_150 = ns["bond13_150"]
    bond13_150h = ns["bond13_150h"]
    bond14_150 = ns["bond14_150"]
    bond14_150h = ns["bond14_150h"]

    ring_ids = sorted(s.segment_id for s in sg.segments if s.kind.startswith("ring"))
    assert len(ring_ids) == 4, f"expected 4 ring segments, got {ring_ids}"

    # bond12: 隣接 3 + ring-tail 1 = 4 pair
    assert len(bond12) == 4, f"expected 4 bond12 pairs, got {bond12}"
    # bond13_150: 1-skip 2 pair
    assert len(bond13_150) == 2, f"expected 2 bond13_150 pairs, got {bond13_150}"
    assert all(b[3] == 1.661 and b[4] == 200 for b in bond13_150h), bond13_150h
    # bond14_150: 2-skip 1 pair (A-D)
    assert len(bond14_150) == 1, f"expected 1 bond14_150 pair (A-D), got {bond14_150}"
    assert all(b[3] == 2.502 and b[4] == 200 for b in bond14_150h), bond14_150h


def test_angle13_chain_linear_eq_0(decane_segmenter, tmp_path):
    """3 segment 以上の chain → angle13 で eq=0 (180° 直線想定)、 stiff 5.0。"""
    sg = decane_segmenter
    mono, _ = sg.export_dpdgen(output_dir=str(tmp_path), monomer_name="dec")
    ns: dict = {}
    exec(Path(mono).read_text(), ns)
    angle13 = ns["angle13"]
    angle13data = ns["angle13data"]
    if len(sg.segments) >= 3:
        assert len(angle13) > 0, (
            f"expected angle13 entries (n_seg={len(sg.segments)}), got {angle13}"
        )
        # chain (= 全 non-ring) → eq 0
        for a in angle13data:
            assert a[3] == 0, f"chain angle eq should be 0, got {a}"
            assert a[4] == 5.0, f"angle stiff should be 5.0, got {a}"


def test_angle13_ring_bend_eq_30(cholesterol_segmenter_for_dpd, tmp_path):
    """cholesterol で ring-ring-ring の angle → eq=30 (= 150° ring bend)。"""
    sg = cholesterol_segmenter_for_dpd
    mono, _ = sg.export_dpdgen(output_dir=str(tmp_path), monomer_name="chol")
    ns: dict = {}
    exec(Path(mono).read_text(), ns)
    angle13data = ns["angle13data"]

    seg_kind = {s.segment_id: s.kind for s in sg.segments}
    # 全 ring の triple は eq=30、 chain を含む triple は eq=0
    n_ring_bend = 0
    for a in angle13data:
        all_ring = (
            seg_kind[a[0]].startswith("ring") and
            seg_kind[a[1]].startswith("ring") and
            seg_kind[a[2]].startswith("ring")
        )
        if all_ring:
            assert a[3] == 30, f"ring-ring-ring angle should have eq=30, got {a}"
            n_ring_bend += 1
        else:
            assert a[3] == 0, f"non-ring-ring-ring should have eq=0, got {a}"
        assert a[4] == 5.0, f"angle stiff should be 5.0, got {a}"
    assert n_ring_bend >= 2, f"expected >=2 ring-bend angles in cholesterol, got {n_ring_bend}"


def test_angle13_cis_double_bond_eq_60(hexene_segmenter, tmp_path):
    """hex-3-ene の中央 C=C → angle13 で該当 triple の eq=60 (= 120°)。"""
    sg = hexene_segmenter
    mono, _ = sg.export_dpdgen(output_dir=str(tmp_path), monomer_name="hex")
    ns: dict = {}
    exec(Path(mono).read_text(), ns)
    angle13data = ns["angle13data"]
    if len(sg.segments) >= 3:
        # cis double bond の周辺 angle が少なくとも 1 つ eq=60
        cis_count = sum(1 for a in angle13data if a[3] == 60)
        assert cis_count > 0, (
            f"expected angle13 with eq=60 for cis C=C "
            f"(n_seg={len(sg.segments)}), got: {angle13data}"
        )
        for a in angle13data:
            assert a[4] == 5.0, f"angle stiff should be 5.0, got {a}"


def test_bond13_180_120_default_commented(naphthalene_segmenter, tmp_path):
    """default では bond13_180 / bond13_120 は **コメントアウト** で出る
    (angle ポテンシャルで代替するため、 距離制約は user 手動 uncomment)。"""
    sg = naphthalene_segmenter
    mono, _ = sg.export_dpdgen(output_dir=str(tmp_path), monomer_name="naph")
    content = Path(mono).read_text()
    # 実行コードに bond13_180 / bond13_120 が含まれていない
    ns: dict = {}
    exec(content, ns)
    assert "bond13_180" not in ns, "bond13_180 should not be active (commented out)"
    assert "bond13_120" not in ns
    # コメントとしては書かれている (uncomment 用 template)
    assert "# bond13_180 = []" in content
    assert "# bond13_120 = []" in content
