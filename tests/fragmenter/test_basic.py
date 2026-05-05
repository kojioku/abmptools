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
