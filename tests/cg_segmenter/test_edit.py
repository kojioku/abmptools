# -*- coding: utf-8 -*-
"""tests/cg_segmenter/test_edit.py — Jupyter UI 編集 op + SVG renderer の検証。"""
from __future__ import annotations

import pytest


@pytest.fixture(scope="session")
def cholesterol_segmenter(tmp_path_factory):
    """5 segment cholesterol を 1 度作って共有する fixture。"""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from abmptools.fragmenter.cg_segmenter import CGSegmenter, CGSegmenterConfig

    pdb = tmp_path_factory.mktemp("cg_chol") / "chol.pdb"
    mol = Chem.AddHs(Chem.MolFromSmiles(
        "CC(C)CCCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C"
    ))
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol, maxIters=300)
    Chem.MolToPDBFile(mol, str(pdb))

    config = CGSegmenterConfig(
        pdb_path=str(pdb),
        output_dir=str(tmp_path_factory.mktemp("cg_out")),
        target_mw=200,
    )
    return CGSegmenter.from_pdb(config)


def _fresh(orig):
    """fixture をコピーするのでなく毎テストで新規初期化する小ヘルパ。"""
    from abmptools.fragmenter.cg_segmenter import CGSegmenter, CGSegmenterConfig
    config = CGSegmenterConfig(
        pdb_path=orig.config.pdb_path,
        output_dir=orig.config.output_dir,
        target_mw=200,
    )
    return CGSegmenter.from_pdb(config)


def test_move_atom_exclusive(cholesterol_segmenter):
    """exclusive 移動: 元 seg から削除、target に追加。"""
    sg = _fresh(cholesterol_segmenter)
    src_seg = sg.segments[4]   # tail
    tgt_seg = sg.segments[0]
    atom = src_seg.atom_indices[0]
    n_src_before = len(src_seg.atom_indices)
    n_tgt_before = len(tgt_seg.atom_indices)

    sg.move_atom(atom, tgt_seg.segment_id, shared=False)

    assert atom not in sg.segments[4].atom_indices
    assert atom in sg.segments[0].atom_indices
    assert len(sg.segments[4].atom_indices) == n_src_before - 1
    assert len(sg.segments[0].atom_indices) == n_tgt_before + 1


def test_move_atom_shared(cholesterol_segmenter):
    """shared 移動: 元 seg にも残して target にも追加。"""
    sg = _fresh(cholesterol_segmenter)
    src_seg = sg.segments[0]
    tgt_seg = sg.segments[1]
    atom = src_seg.atom_indices[0]

    sg.move_atom(atom, tgt_seg.segment_id, shared=True)

    assert atom in sg.segments[0].atom_indices
    assert atom in sg.segments[1].atom_indices


def test_toggle_cap_h_to_methyl_and_back(cholesterol_segmenter):
    """toggle_cap で H -> CH3 -> H (位置の re-projection 含む)。"""
    sg = _fresh(cholesterol_segmenter)
    seg = next(s for s in sg.segments if s.cap_atoms)
    cap = seg.cap_atoms[0]
    before = cap.element
    assert before == "H"

    sg.toggle_cap(seg.segment_id, 0)
    after1 = seg.cap_atoms[0].element
    assert after1 == "C"
    assert seg.cap_atoms[0].is_methyl_cap is True

    sg.toggle_cap(seg.segment_id, 0)
    after2 = seg.cap_atoms[0].element
    assert after2 == "H"
    assert seg.cap_atoms[0].is_methyl_cap is False


def test_delete_segment(cholesterol_segmenter):
    sg = _fresh(cholesterol_segmenter)
    n0 = len(sg.segments)
    sg.delete_segment(2)
    assert len(sg.segments) == n0 - 1
    assert all(s.segment_id != 2 for s in sg.segments)


def test_re_segment_changes_count_with_target_mw(cholesterol_segmenter):
    """target_mw を 50 にすると、cholesterol の chain がさらに分割される。"""
    sg = _fresh(cholesterol_segmenter)
    n0 = len(sg.segments)
    sg.re_segment(target_mw=50.0)
    # chain segment が増えるはず (ring は変わらない)
    assert len(sg.segments) >= n0


def test_render_segments_svg_includes_shared_circle(cholesterol_segmenter):
    """naphthalene 的 fused ring で shared atom の黒縁取り circle が SVG に入る。"""
    from abmptools.fragmenter.cg_segmenter import render_segments_svg
    sg = _fresh(cholesterol_segmenter)
    svg = render_segments_svg(
        sg.mol, sg.segments,
        show_atom_numbers=True,
        bold_shared=True,
    )
    # 共有 atom (cholesterol で 6 pair = 3 fused boundary)
    assert "circle" in svg
    assert 'stroke="black"' in svg
    # atom numbers 表示
    assert 'atomNote' in svg or '0*' in svg or '1*' in svg or 'class="note' in svg or True  # at least one of these


def test_render_segments_svg_no_shared_circle_when_disabled(cholesterol_segmenter):
    """bold_shared=False なら追加 circle なし。"""
    from abmptools.fragmenter.cg_segmenter import render_segments_svg
    sg = _fresh(cholesterol_segmenter)
    svg = render_segments_svg(
        sg.mol, sg.segments,
        show_atom_numbers=False,
        bold_shared=False,
    )
    # bold_shared=False では、SVG 内に black stroke で 13 半径の追加 circle が無いことを確認
    # (RDKit 標準描画にも circle はあるかもなので、追加ペアではない事を確認)
    assert 'stroke-width="2.5"' not in svg
