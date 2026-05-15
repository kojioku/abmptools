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
