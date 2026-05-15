# -*- coding: utf-8 -*-
"""tests/cg_segmenter/test_basic.py — CG segmentation end-to-end and edge cases."""
from __future__ import annotations

import json
from pathlib import Path

import pytest


# ---------------------------------------------------------------------------
# Dataclass roundtrip
# ---------------------------------------------------------------------------


def test_cg_segmenter_config_json_roundtrip(tmp_path):
    from abmptools.fragmenter.cg_segmenter import CGSegmenterConfig
    c = CGSegmenterConfig(
        pdb_path="x.pdb", target_mw=150.0,
        hetero_cap_methyl_elements=["N", "O"],
    )
    p = tmp_path / "config.json"
    c.to_json(str(p))
    c2 = CGSegmenterConfig.from_json(str(p))
    assert c2.target_mw == 150.0
    assert c2.hetero_cap_methyl_elements == ["N", "O"]


def test_segment_dict_roundtrip():
    from abmptools.fragmenter.cg_segmenter import Segment, CapAtom
    seg = Segment(
        segment_id=2, atom_indices=[1, 2, 3], kind="ring",
        ring_indices=[0],
        cap_atoms=[CapAtom(parent_atom_idx=1, element="H",
                           position=(0.0, 1.0, 0.0))],
    )
    d = seg.to_dict()
    seg2 = Segment.from_dict(d)
    assert seg2.segment_id == 2
    assert seg2.kind == "ring"
    assert len(seg2.cap_atoms) == 1
    assert seg2.cap_atoms[0].element == "H"
    assert seg2.n_atoms == 4  # 3 ring + 1 H cap


# ---------------------------------------------------------------------------
# End-to-end via CGSegmenter
# ---------------------------------------------------------------------------


def test_propane_single_chain_no_caps(propane_pdb, tmp_path):
    from abmptools.fragmenter.cg_segmenter import CGSegmenter, CGSegmenterConfig
    config = CGSegmenterConfig(
        pdb_path=str(propane_pdb), target_mw=200,
        output_dir=str(tmp_path / "out"),
    )
    sg = CGSegmenter.from_pdb(config)
    result = sg.export()
    assert len(result.segments) == 1
    seg = result.segments[0]
    assert seg.kind == "chain"
    assert len(seg.atom_indices) == 3
    assert len(seg.cap_atoms) == 0


def test_octane_two_chain_segments_with_h_caps(octane_pdb, tmp_path):
    from abmptools.fragmenter.cg_segmenter import CGSegmenter, CGSegmenterConfig
    config = CGSegmenterConfig(
        pdb_path=str(octane_pdb), target_mw=50,
        output_dir=str(tmp_path / "out"),
    )
    sg = CGSegmenter.from_pdb(config)
    result = sg.export()
    assert len(result.segments) == 2
    for seg in result.segments:
        assert seg.kind == "chain"
        assert len(seg.cap_atoms) == 1
        assert seg.cap_atoms[0].element == "H"
        assert seg.cap_atoms[0].is_methyl_cap is False
    assert len(result.shared_atom_pairs) == 0


def test_methyl_acetate_o_side_gets_methyl_cap(methyl_acetate_pdb, tmp_path):
    from abmptools.fragmenter.cg_segmenter import CGSegmenter, CGSegmenterConfig
    config = CGSegmenterConfig(
        pdb_path=str(methyl_acetate_pdb), target_mw=30,
        output_dir=str(tmp_path / "out"),
    )
    sg = CGSegmenter.from_pdb(config)
    result = sg.export()
    assert len(result.segments) >= 2
    # ester O が atom_indices に含まれる segment は CH3 cap を持つ (boundary に O が来る)
    # 検証: 全 cap の中で少なくとも 1 つは CH3 cap で、その parent atom は O
    from rdkit import Chem
    mol = Chem.MolFromPDBFile(str(methyl_acetate_pdb), removeHs=False)
    methyl_caps = [
        c for s in result.segments for c in s.cap_atoms
        if c.is_methyl_cap
    ]
    assert len(methyl_caps) >= 1
    o_methyl = [
        c for c in methyl_caps
        if mol.GetAtomWithIdx(c.parent_atom_idx).GetSymbol() == "O"
    ]
    assert len(o_methyl) >= 1, "expected at least 1 CH3 cap on O side"


def test_benzene_single_ring_no_caps(benzene_pdb, tmp_path):
    from abmptools.fragmenter.cg_segmenter import CGSegmenter, CGSegmenterConfig
    config = CGSegmenterConfig(
        pdb_path=str(benzene_pdb), output_dir=str(tmp_path / "out"),
    )
    sg = CGSegmenter.from_pdb(config)
    result = sg.export()
    assert len(result.segments) == 1
    seg = result.segments[0]
    assert seg.kind == "ring"
    assert len(seg.atom_indices) == 6
    assert len(seg.cap_atoms) == 0  # 環内のみ、境界なし


def test_naphthalene_fused_rings_share_atoms(naphthalene_pdb, tmp_path):
    from abmptools.fragmenter.cg_segmenter import CGSegmenter, CGSegmenterConfig
    config = CGSegmenterConfig(
        pdb_path=str(naphthalene_pdb), output_dir=str(tmp_path / "out"),
        allow_atom_sharing=True,
    )
    sg = CGSegmenter.from_pdb(config)
    result = sg.export()
    assert len(result.segments) == 2
    assert all(s.kind == "ring" for s in result.segments)
    # 2 shared atoms (fused ring boundary)
    assert len(result.shared_atom_pairs) == 2
    shared_atoms = {p[0] for p in result.shared_atom_pairs}
    assert len(shared_atoms) == 2


def test_cyclohexanol_absorbs_oh(cyclohexanol_pdb, tmp_path):
    from abmptools.fragmenter.cg_segmenter import CGSegmenter, CGSegmenterConfig
    config = CGSegmenterConfig(
        pdb_path=str(cyclohexanol_pdb), output_dir=str(tmp_path / "out"),
        absorb_single_substituent=True,
    )
    sg = CGSegmenter.from_pdb(config)
    result = sg.export()
    assert len(result.segments) == 1
    seg = result.segments[0]
    assert seg.kind == "ring_with_substituent"
    assert len(seg.atom_indices) == 7  # 6 ring + 1 OH O


def test_cyclohexyl_octyl_ring_plus_chain(cyclohexyl_octyl_pdb, tmp_path):
    from abmptools.fragmenter.cg_segmenter import CGSegmenter, CGSegmenterConfig
    config = CGSegmenterConfig(
        pdb_path=str(cyclohexyl_octyl_pdb), output_dir=str(tmp_path / "out"),
        target_mw=200,
    )
    sg = CGSegmenter.from_pdb(config)
    result = sg.export()
    assert len(result.segments) == 2
    kinds = sorted(s.kind for s in result.segments)
    assert kinds == ["chain", "ring"]
    # ring-chain boundary に各 1 H cap
    for s in result.segments:
        assert len(s.cap_atoms) == 1
        assert s.cap_atoms[0].element == "H"


def test_output_files_pdb_xyz_summary(naphthalene_pdb, tmp_path):
    """PDB + XYZ + segments.json の出力ファイル確認。"""
    from abmptools.fragmenter.cg_segmenter import CGSegmenter, CGSegmenterConfig
    out_dir = tmp_path / "naph_out"
    config = CGSegmenterConfig(pdb_path=str(naphthalene_pdb), output_dir=str(out_dir))
    sg = CGSegmenter.from_pdb(config)
    sg.export()
    assert (out_dir / "seg_000.pdb").exists()
    assert (out_dir / "seg_000.xyz").exists()
    assert (out_dir / "seg_001.pdb").exists()
    assert (out_dir / "seg_001.xyz").exists()
    summary = json.loads((out_dir / "segments.json").read_text())
    assert summary["n_segments"] == 2
    assert len(summary["shared_atom_pairs"]) == 2


def test_cli_build(propane_pdb, tmp_path):
    from abmptools.fragmenter.cg_segmenter.__main__ import main
    out_dir = tmp_path / "cli_out"
    code = main([
        "build", "--pdb", str(propane_pdb),
        "--output-dir", str(out_dir),
        "--target-mw", "200",
    ])
    assert code == 0
    assert (out_dir / "segments.json").exists()
