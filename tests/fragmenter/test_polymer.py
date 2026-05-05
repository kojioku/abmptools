# -*- coding: utf-8 -*-
"""tests/fragmenter/test_polymer.py — γ 経路 (declare_same_pattern)。"""
from __future__ import annotations


def test_pe10_pe11_default_grouping(pe_n10_n11_mix_pdb):
    """同一視せずに group_by_smiles するだけだと PE N=10 と N=11 は別 group。"""
    from abmptools.fragmenter import load_pdb_molecules, group_by_smiles
    loaded = load_pdb_molecules(str(pe_n10_n11_mix_pdb))
    assert len(loaded) == 2
    groups = group_by_smiles(loaded)
    assert len(groups) == 2
    smis = sorted(g.smiles for g in groups)
    # PE N=10 と PE N=11 の canonical SMILES は鎖長が違うので別
    assert smis[0].count("C") + 1 != smis[1].count("C") + 1 or smis[0] != smis[1]


def test_pe10_pe11_pattern_transfer(pe_n10_n11_mix_pdb):
    """γ 経路で PE N=10 / N=11 を同一視 → 短い方も同じパターンで切断される。"""
    from abmptools.fragmenter import (
        FragmenterConfig, load_pdb_molecules, group_by_smiles,
        suggest_cuts_for_groups,
    )
    from abmptools.fragmenter.polymer import declare_same_pattern

    config = FragmenterConfig(pdb_path=str(pe_n10_n11_mix_pdb), target_mw=80.0)
    loaded = load_pdb_molecules(config.pdb_path)
    groups = group_by_smiles(loaded)
    assert len(groups) == 2

    suggest_cuts_for_groups(groups, loaded, config)

    # 鎖長が違う 2 つの group の SMILES を取得
    smiles_list = [g.smiles for g in groups]
    declare_same_pattern(groups, loaded, [smiles_list])

    # γ 経路適用後、両 group とも cut_sites が ≥ 1 になっているはず
    for g in groups:
        assert len(g.cut_sites) >= 1, (
            f"group {g.smiles} has no cut_sites after pattern transfer"
        )


def test_declare_same_pattern_no_match():
    """指定した SMILES が groups に存在しない場合は warning だけで例外を投げない。"""
    from abmptools.fragmenter import MoleculeGroup
    from abmptools.fragmenter.polymer import declare_same_pattern

    groups = [MoleculeGroup(smiles="CCC", representative_mol_idx=0, n_copies=1)]
    # 該当 SMILES なし → そのまま (例外なし)
    out = declare_same_pattern(groups, molecules=[], same_pattern_smiles=[["XYZ", "ABC"]])
    assert out is groups


def test_declare_same_pattern_empty():
    """same_pattern_smiles が空なら何もしない。"""
    from abmptools.fragmenter import MoleculeGroup
    from abmptools.fragmenter.polymer import declare_same_pattern

    groups = [MoleculeGroup(smiles="CCC", representative_mol_idx=0, n_copies=1)]
    out = declare_same_pattern(groups, molecules=[], same_pattern_smiles=[])
    assert out is groups
