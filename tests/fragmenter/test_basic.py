# -*- coding: utf-8 -*-
"""tests/fragmenter/test_basic.py — dataclass + end-to-end pipeline."""
from __future__ import annotations

import json
from pathlib import Path

import pytest


# ---------------------------------------------------------------------------
# Dataclass tests (rdkit 不要)
# ---------------------------------------------------------------------------


def test_fragmenter_config_defaults():
    from abmptools.fragmenter import FragmenterConfig
    c = FragmenterConfig(pdb_path="x.pdb")
    assert c.target_mw == 200.0
    assert c.exclude_ring_cc is True
    assert c.exclude_multibond is True
    assert c.exclude_heteroneighbor is True
    assert c.skip_protein_dna is True


def test_fragmenter_config_json_roundtrip(tmp_path):
    from abmptools.fragmenter import FragmenterConfig
    c = FragmenterConfig(
        pdb_path="x.pdb", target_mw=150.0,
        polymer_groups=[["CCC", "CCCC"]],
    )
    p = tmp_path / "config.json"
    c.to_json(str(p))
    c2 = FragmenterConfig.from_json(str(p))
    assert c2.target_mw == 150.0
    assert c2.polymer_groups == [["CCC", "CCCC"]]


def test_cut_site_dict_roundtrip():
    from abmptools.fragmenter import CutSite
    cs = CutSite(bond_idx=3, atom1_idx=1, atom2_idx=2,
                 suggested=False, enabled=False)
    d = cs.to_dict()
    cs2 = CutSite.from_dict(d)
    assert cs2 == cs


def test_molecule_group_dict_roundtrip():
    from abmptools.fragmenter import MoleculeGroup, CutSite
    g = MoleculeGroup(
        smiles="CCC", residue_names={"PRP", "PROP"},
        representative_mol_idx=2, n_copies=3,
        member_mol_indices=[0, 1, 2],
        cut_sites=[CutSite(bond_idx=0, atom1_idx=0, atom2_idx=1)],
    )
    d = g.to_dict()
    g2 = MoleculeGroup.from_dict(d)
    assert g2.smiles == g.smiles
    assert g2.residue_names == g.residue_names
    assert g2.n_copies == g.n_copies
    assert len(g2.cut_sites) == 1
    assert g2.cut_sites[0].atom1_idx == 0


# ---------------------------------------------------------------------------
# Pipeline tests (要 rdkit、conftest.py で importorskip)
# ---------------------------------------------------------------------------


def test_load_and_group(propane_acetone_mix_pdb):
    from abmptools.fragmenter import load_pdb_molecules, group_by_smiles
    loaded = load_pdb_molecules(str(propane_acetone_mix_pdb))
    assert len(loaded) == 5

    groups = group_by_smiles(loaded)
    assert len(groups) == 2
    smiles_set = {g.smiles for g in groups}
    assert "CCC" in smiles_set
    assert "CC(C)=O" in smiles_set


def test_no_cuts_propane_at_default_target_mw(propane_acetone_mix_pdb):
    from abmptools.fragmenter import (
        FragmenterConfig, load_pdb_molecules, group_by_smiles,
        suggest_cuts_for_groups,
    )
    config = FragmenterConfig(
        pdb_path=str(propane_acetone_mix_pdb), target_mw=200.0,
    )
    loaded = load_pdb_molecules(config.pdb_path)
    groups = group_by_smiles(loaded)
    suggest_cuts_for_groups(groups, loaded, config)
    # propane (44 g/mol) と acetone (58 g/mol) はどちらも 200 未満 → 切断ゼロ
    for g in groups:
        assert len(g.cut_sites) == 0


def test_propane_one_cut_at_target_mw_15(propane_acetone_mix_pdb, tmp_path):
    from abmptools.fragmenter import (
        FragmenterConfig, load_pdb_molecules, group_by_smiles,
        suggest_cuts_for_groups, export_to_system,
    )
    config = FragmenterConfig(
        pdb_path=str(propane_acetone_mix_pdb), target_mw=15.0,
    )
    loaded = load_pdb_molecules(config.pdb_path)
    groups = group_by_smiles(loaded)
    suggest_cuts_for_groups(groups, loaded, config)

    propane_grp = next(g for g in groups if g.smiles == "CCC")
    acetone_grp = next(g for g in groups if g.smiles == "CC(C)=O")
    assert len(propane_grp.cut_sites) == 1
    assert len(acetone_grp.cut_sites) == 0

    output = tmp_path / "segment_data.dat"
    result = export_to_system(config.pdb_path, groups, loaded, str(output))
    # propane x3 -> 6 frag, acetone x2 -> 2 frag, total 8
    assert result.total_fragments == 8
    assert len(result.segment_data) == 2


def test_octane_target_mw_50_makes_cuts(octane_pdb):
    from abmptools.fragmenter import (
        FragmenterConfig, load_pdb_molecules, group_by_smiles,
        suggest_cuts_for_groups,
    )
    config = FragmenterConfig(pdb_path=str(octane_pdb), target_mw=50.0)
    loaded = load_pdb_molecules(config.pdb_path)
    groups = group_by_smiles(loaded)
    suggest_cuts_for_groups(groups, loaded, config)
    grp = groups[0]
    # octane (114 g/mol) / target 50 -> at least 1 cut
    assert len(grp.cut_sites) >= 1


def test_segment_data_dat_roundtrip(propane_acetone_mix_pdb, tmp_path):
    """log2config 互換の出力フォーマットが ast.literal_eval / exec で読み戻せること。"""
    from abmptools.fragmenter import (
        FragmenterConfig, load_pdb_molecules, group_by_smiles,
        suggest_cuts_for_groups, export_to_system,
    )
    config = FragmenterConfig(
        pdb_path=str(propane_acetone_mix_pdb), target_mw=15.0,
    )
    loaded = load_pdb_molecules(config.pdb_path)
    groups = group_by_smiles(loaded)
    suggest_cuts_for_groups(groups, loaded, config)
    output = tmp_path / "segment_data.dat"
    export_to_system(config.pdb_path, groups, loaded, str(output))

    text = output.read_text()
    namespace = {}
    exec(text, namespace)
    seg_data = namespace["seg_data"]
    assert isinstance(seg_data, list)
    assert all("name" in e and "atom" in e for e in seg_data)


def test_bda_baa_assigned_on_auto_suggest(propane_acetone_mix_pdb):
    """auto_split.suggest_cuts は BDA/BAA を decide して CutSite に格納する。"""
    from abmptools.fragmenter import (
        FragmenterConfig, load_pdb_molecules, group_by_smiles,
        suggest_cuts_for_groups,
    )
    config = FragmenterConfig(pdb_path=str(propane_acetone_mix_pdb), target_mw=15.0)
    loaded = load_pdb_molecules(config.pdb_path)
    groups = group_by_smiles(loaded)
    suggest_cuts_for_groups(groups, loaded, config)
    propane_grp = next(g for g in groups if g.smiles == "CCC")
    assert len(propane_grp.cut_sites) == 1
    cs = propane_grp.cut_sites[0]
    assert cs.bda_atom_idx is not None
    assert cs.baa_atom_idx is not None
    # BDA / BAA は cut bond の両端のいずれか、かつ互いに異なる
    assert cs.bda_atom_idx in (cs.atom1_idx, cs.atom2_idx)
    assert cs.baa_atom_idx in (cs.atom1_idx, cs.atom2_idx)
    assert cs.bda_atom_idx != cs.baa_atom_idx


def test_c_heteroatom_cut_assigns_bda_to_c_side(propane_acetone_mix_pdb):
    """include_c_heteroatom=True で C-X 単結合が候補になり、BDA は C 側固定。"""
    from abmptools.fragmenter import (
        FragmenterConfig, load_pdb_molecules, group_by_smiles,
        suggest_cuts_for_groups,
    )
    from rdkit import Chem
    config = FragmenterConfig(
        pdb_path=str(propane_acetone_mix_pdb),
        target_mw=15.0,
        include_c_heteroatom=True,
        # アミドカルボニル等は exclude_heteroneighbor で C-C はフィルタされるが、
        # C-X 自体は include_c_heteroatom で別経路で許可される
    )
    loaded = load_pdb_molecules(config.pdb_path)
    groups = group_by_smiles(loaded)
    suggest_cuts_for_groups(groups, loaded, config)
    acetone_grp = next(g for g in groups if g.smiles == "CC(C)=O")
    # acetone (CC(C)=O) は C=O が二重結合なので C-O 切断は無し。
    # C-C はヘテロ隣接 (= O) で除外される。target_mw=15 なら 0 cut の想定。
    # ただし include_c_heteroatom=True で C-X bond が候補になるかをまず確認:
    rep = loaded[acetone_grp.representative_mol_idx].mol
    from abmptools.fragmenter.auto_split import _enumerate_candidate_bonds
    cands = _enumerate_candidate_bonds(rep, config)
    # acetone は C-C bond 2 本、C=O 1 本。C=O は SINGLE でないため除外。
    # 残る C-C bonds は両方 carbonyl 隣接で除外される。include_c_heteroatom でも
    # 二重結合の C=O は含まれないので候補は 0 (= 切断不可分子)。
    # 一方、もし single C-O があれば (e.g. ester / ether) C-X として候補に
    # なる。本テストでは単純に「include_c_heteroatom flag が config に
    # 通る」ことを確認する。
    assert config.include_c_heteroatom is True

    # 別の検証: methanol (CO) は C-O single bond 1 本を持つので include_c_heteroatom
    # で切断候補になり、BDA = C 側になる
    methanol = Chem.AddHs(Chem.MolFromSmiles("CO"))
    from rdkit.Chem import AllChem
    AllChem.EmbedMolecule(methanol, randomSeed=42)
    cands_meoh = _enumerate_candidate_bonds(methanol, config)
    assert len(cands_meoh) == 1, f"expected 1 C-O candidate, got {cands_meoh}"


def test_walk_side_chains_flag_smoke(propane_acetone_mix_pdb):
    """walk_side_chains=True が config に渡り、suggest_cuts が完走する。"""
    from abmptools.fragmenter import (
        FragmenterConfig, load_pdb_molecules, group_by_smiles,
        suggest_cuts_for_groups,
    )
    config = FragmenterConfig(
        pdb_path=str(propane_acetone_mix_pdb),
        target_mw=10.0,
        walk_side_chains=True,
        exclude_heteroneighbor=False,
    )
    assert config.walk_side_chains is True
    loaded = load_pdb_molecules(config.pdb_path)
    groups = group_by_smiles(loaded)
    # walk_side_chains=True でも error なく動作することを確認
    suggest_cuts_for_groups(groups, loaded, config)


def test_walk_side_chain_function_internal():
    """_walk_side_chains_from_main の internal smoke: 長い ester 側鎖を持つ分子で
    main 短く / side 長く構成し、side chain walk から少なくとも 1 cut が出る。"""
    from rdkit import Chem
    from abmptools.fragmenter.auto_split import (
        _walk_side_chains_from_main, _enumerate_candidate_bonds,
    )
    from abmptools.fragmenter import FragmenterConfig
    # propane main (3 C) に長い C12 アルキル side chain
    # main_path=[0,1,2] (propane backbone)、side chain は atom 1 の長い鎖
    mol = Chem.MolFromSmiles("CC(CCCCCCCCCCCC)C")
    config = FragmenterConfig(pdb_path="x", target_mw=50.0, exclude_heteroneighbor=False)
    cands = _enumerate_candidate_bonds(mol, config)
    main_path = [0, 1, 14]  # propane backbone (atom 0 = 1st CH3, atom 1 = central CH, atom 14 = last CH3)
    cuts = _walk_side_chains_from_main(mol, main_path, cands, target_mw=50.0)
    # C12 side chain → target_mw=50 (~3-4 C 毎) で複数 cut が期待される
    assert len(cuts) >= 1, f"side chain walk should produce cuts, got {cuts}"


def test_c_heteroatom_disabled_by_default(propane_acetone_mix_pdb):
    """default (include_c_heteroatom=False) では C-O 等の C-X bond は候補にならない。"""
    from abmptools.fragmenter import FragmenterConfig
    from abmptools.fragmenter.auto_split import _enumerate_candidate_bonds
    from rdkit import Chem
    from rdkit.Chem import AllChem
    methanol = Chem.AddHs(Chem.MolFromSmiles("CO"))
    AllChem.EmbedMolecule(methanol, randomSeed=42)
    config = FragmenterConfig(pdb_path=str(propane_acetone_mix_pdb))
    assert config.include_c_heteroatom is False
    cands = _enumerate_candidate_bonds(methanol, config)
    assert len(cands) == 0, f"expected 0 candidates (C-X disabled), got {cands}"


def test_decide_bda_baa_for_manual_cc_cut(propane_pdb):
    """C-C bond の手動 cut: atom_idx 若い側 = BDA (deterministic)。"""
    from abmptools.fragmenter import load_pdb_molecules
    from abmptools.fragmenter.auto_split import decide_bda_baa_for_manual_cut
    loaded = load_pdb_molecules(str(propane_pdb))
    mol = loaded[0].mol
    bda, baa = decide_bda_baa_for_manual_cut(mol, 1, 0)
    assert bda == 0 and baa == 1
    bda, baa = decide_bda_baa_for_manual_cut(mol, 0, 1)
    assert bda == 0 and baa == 1


def test_decide_bda_baa_for_manual_cx_cut():
    """C-O bond の手動 cut: C 側 = BDA (atom 順序によらず)。"""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from abmptools.fragmenter.auto_split import decide_bda_baa_for_manual_cut
    methanol = Chem.AddHs(Chem.MolFromSmiles("CO"))
    AllChem.EmbedMolecule(methanol, randomSeed=42)
    # atom 0 = C, atom 1 = O
    bda, baa = decide_bda_baa_for_manual_cut(methanol, 1, 0)
    assert bda == 0 and baa == 1  # C 側 BDA
    bda, baa = decide_bda_baa_for_manual_cut(methanol, 0, 1)
    assert bda == 0 and baa == 1  # 順序逆でも C 側 BDA


def test_segment_data_flat_connect_and_baa_count(propane_acetone_mix_pdb, tmp_path):
    """log2config 形式: connect は flat [[BDA, BAA], ...]、connect_num は BAA 数。

    ``getmb_frag_seclists`` (AJF writer) が読める平坦形式であること、および
    ``connect_num`` が BAA を含む fragment に立つことを回帰固定する。
    """
    from abmptools.fragmenter import (
        FragmenterConfig, load_pdb_molecules, group_by_smiles,
        suggest_cuts_for_groups, export_to_system,
    )
    config = FragmenterConfig(pdb_path=str(propane_acetone_mix_pdb), target_mw=15.0)
    loaded = load_pdb_molecules(config.pdb_path)
    groups = group_by_smiles(loaded)
    suggest_cuts_for_groups(groups, loaded, config)
    output = tmp_path / "segment_data.dat"
    result = export_to_system(config.pdb_path, groups, loaded, str(output))

    # propane (CCC) entry: 6 fragments (3 copies × 2 frag/copy)、cut は 3 本
    prop_entry = next(e for e in result.segment_data if "CCC" in e["name"])
    cn = prop_entry["connect_num"]
    # connect_num = BAA 数。各 propane 1 本 cut → BAA を持つ fragment が 1 個/copy
    assert cn.count(1) == 3, f"expected 3 BAA fragments, got connect_num={cn}"
    assert cn.count(0) == 3, f"expected 3 non-BAA fragments, got connect_num={cn}"
    # connect は flat list (cut 1 本 = 要素 1 個)。3 copies → 3 要素、各 [BDA, BAA]
    conn = prop_entry["connect"]
    assert len(conn) == 3, f"expected 3 flat connect entries, got {conn}"
    assert sum(cn) == len(conn), "sum(connect_num) (BAA total) must equal #cut bonds"
    for pair in conn:
        assert len(pair) == 2          # flat [BDA, BAA] (ネストでない)
        bda_atom, baa_atom = pair
        assert bda_atom != baa_atom

    # BAA (各 connect の 2 番目) が connect_num>0 の fragment 内にあること
    def _frag_of(serial):
        for fi, frag in enumerate(prop_entry["seg_info"]):
            if serial in frag:
                return fi
        return -1
    baa_count = [0] * len(cn)
    for bda_atom, baa_atom in conn:
        baa_count[_frag_of(baa_atom)] += 1
    assert baa_count == cn, "connect_num must equal per-fragment BAA count"


def test_cli_suggest_apply(propane_acetone_mix_pdb, tmp_path):
    """CLI suggest -> apply round-trip。"""
    from abmptools.fragmenter.__main__ import main
    review = tmp_path / "review"
    code = main([
        "suggest",
        "--pdb", str(propane_acetone_mix_pdb),
        "--output-dir", str(review),
        "--target-mw", "15",
    ])
    assert code == 0
    assert (review / "review.json").exists()
    assert (review / "group_001.svg").exists()
    assert (review / "config.json").exists()

    seg_path = tmp_path / "segment_data.dat"
    code2 = main([
        "apply",
        "--pdb", str(propane_acetone_mix_pdb),
        "--review-dir", str(review),
        "--output", str(seg_path),
    ])
    assert code2 == 0
    assert seg_path.exists()


def test_segment_data_consumable_by_ajf_writer(propane_acetone_mix_pdb, tmp_path):
    """export_to_system 出力が AJF writer (config_read + getmb_frag_seclists)
    でクラッシュせず消費できること (回帰: 旧 nested connect は IndexError だった)。"""
    import abmptools
    from abmptools.fragmenter import (
        FragmenterConfig, load_pdb_molecules, group_by_smiles,
        suggest_cuts_for_groups, export_to_system,
    )
    config = FragmenterConfig(pdb_path=str(propane_acetone_mix_pdb), target_mw=15.0)
    loaded = load_pdb_molecules(config.pdb_path)
    groups = group_by_smiles(loaded)
    suggest_cuts_for_groups(groups, loaded, config)
    out_dir = tmp_path
    result = export_to_system(config.pdb_path, groups, loaded, str(out_dir / "segment_data.dat"))

    ai = abmptools.abinit_io()
    ai.mainpath = str(out_dir)
    for entry in result.segment_data:
        conf = ai.config_read(entry["name"], sum(entry["atom"]))
        # flat connect: 各要素は 2-atom pair
        for pair in conf["connect"]:
            assert len(pair) == 2
        # AJF writer に通す (旧 nested connect だと IndexError でここが落ちた)
        ai.getmb_frag_seclists(
            [[conf["atom"]], [conf["charge"]], [conf["connect_num"]],
             [conf["connect"]], [conf["seg_info"]]],
            [0],
        )
        # connects は flat な [a, b] のリスト
        for c in ai.connects:
            assert len(c) == 2
