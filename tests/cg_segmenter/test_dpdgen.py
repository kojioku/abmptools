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
    # monomer: bond12, bond12h を定義する純粋な Python
    ns: dict = {}
    exec(Path(mono).read_text(), ns)
    assert "bond12" in ns
    assert "bond12h" in ns
    assert isinstance(ns["bond12"], list)
    assert all(isinstance(b, list) and len(b) == 2 for b in ns["bond12"])
    assert all(len(b) == 6 for b in ns["bond12h"])
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


def test_fused_ring_full_connect_cholesterol(cholesterol_segmenter_for_dpd, tmp_path):
    """4 fused rings (cholesterol) → A-B-C-D 全 6 ring pair が bond12 に含まれる。"""
    import itertools
    sg = cholesterol_segmenter_for_dpd
    mono, _ = sg.export_dpdgen(output_dir=str(tmp_path), monomer_name="chol")
    ns: dict = {}
    exec(Path(mono).read_text(), ns)
    bond12 = ns["bond12"]
    ring_ids = sorted(s.segment_id for s in sg.segments if s.kind.startswith("ring"))
    assert len(ring_ids) == 4, f"expected 4 ring segments for cholesterol, got {ring_ids}"
    for a, b in itertools.combinations(ring_ids, 2):
        pair = [a, b]
        assert pair in bond12, (
            f"fused-ring pair {pair} missing from bond12 = {bond12}"
        )


def test_bond13_180_path_length_2(decane_segmenter, tmp_path):
    """3 segment 以上の chain → bond13_180 (path-length-2) entry が出る。"""
    sg = decane_segmenter
    mono, _ = sg.export_dpdgen(output_dir=str(tmp_path), monomer_name="dec")
    ns: dict = {}
    exec(Path(mono).read_text(), ns)
    bond13_180 = ns["bond13_180"]
    bond13_180h = ns["bond13_180h"]
    if len(sg.segments) >= 3:
        assert len(bond13_180) > 0, (
            f"expected bond13_180 entries (n_seg={len(sg.segments)}), got {bond13_180}"
        )
        # chain は stiff 80
        assert all(b[3] == 1.72 for b in bond13_180h), bond13_180h
        assert all(b[4] == 80 for b in bond13_180h), bond13_180h


def test_bond13_180_excludes_bond12_pairs(cholesterol_segmenter_for_dpd, tmp_path):
    """ring group 内部の 1-3 (例: A-C) は既に bond12 にあるので bond13_180 には出ない。"""
    sg = cholesterol_segmenter_for_dpd
    mono, _ = sg.export_dpdgen(output_dir=str(tmp_path), monomer_name="chol")
    ns: dict = {}
    exec(Path(mono).read_text(), ns)
    bond12 = [tuple(p) for p in ns["bond12"]]
    bond13_180 = [tuple(p) for p in ns["bond13_180"]]
    for p in bond13_180:
        assert p not in bond12, f"bond13_180 {p} also appears in bond12"


def test_bond13_120_cis_double_bond(hexene_segmenter, tmp_path):
    """中央に C=C を持つ hex-3-ene → bond13_120 entry が出る (3+ segments)。"""
    sg = hexene_segmenter
    mono, _ = sg.export_dpdgen(output_dir=str(tmp_path), monomer_name="hex")
    ns: dict = {}
    exec(Path(mono).read_text(), ns)
    bond13_120 = ns["bond13_120"]
    bond13_120h = ns["bond13_120h"]
    if len(sg.segments) >= 3:
        assert len(bond13_120) > 0, (
            f"expected bond13_120 entries for cis C=C (n_seg={len(sg.segments)}), "
            f"got {bond13_120}"
        )
        for b in bond13_120h:
            assert b[3] == 1.49
            assert b[4] == 70


def test_no_bond13_for_single_or_double_segment(naphthalene_segmenter, tmp_path):
    """naphthalene (2 ring segments) → bond13_* は空 (path length 2 ペアなし)。"""
    sg = naphthalene_segmenter
    mono, _ = sg.export_dpdgen(output_dir=str(tmp_path), monomer_name="naph")
    ns: dict = {}
    exec(Path(mono).read_text(), ns)
    if len(sg.segments) <= 2:
        assert ns["bond13_180"] == []
        assert ns["bond13_120"] == []
