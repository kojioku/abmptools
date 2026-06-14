# -*- coding: utf-8 -*-
"""tests/cg_segmenter/test_fcews_export.py

cg_segmenter → FCEWS segment_data.dat export の検証。

核心: 生成した segment_data.dat を **abmptools.abinit_io.config_read で
in-process round-trip** し、FCEWS が name 一致で読めること + connect が
[BAA, BDA] (BDA = connect_num>0 fragment 内 = 第2要素) であることを固定する。
"""
from __future__ import annotations

import ast
from pathlib import Path

import pytest

pytest.importorskip("rdkit")

CHOLESTEROL_PDB = (
    Path(__file__).resolve().parents[2]
    / "sample" / "cg_segmenter" / "cholesterol_rdkit.pdb"
)


def _parse_seg_data(path: Path) -> list:
    tree = ast.parse(path.read_text())
    ns = {}
    for node in tree.body:
        if (isinstance(node, ast.Assign) and len(node.targets) == 1
                and isinstance(node.targets[0], ast.Name)):
            ns[node.targets[0].id] = ast.literal_eval(node.value)
    return ns["seg_data"]


def _frag_index_of_serial(entry: dict, serial: int) -> int:
    for fi, frag in enumerate(entry["seg_info"]):
        if serial in frag:
            return fi
    return -1


# ---------------------------------------------------------------------------
# structure + partition
# ---------------------------------------------------------------------------


def test_octane_fcews_entry(octane_pdb):
    """octane を target_mw=50 で複数 fragment に切り、FCEWS form を検証。"""
    from abmptools.fragmenter.cg_segmenter import CGSegmenter, CGSegmenterConfig
    from abmptools.fragmenter.cg_segmenter.fcews_export import build_fcews_segment_data

    cfg = CGSegmenterConfig(pdb_path=str(octane_pdb), target_mw=50.0,
                            allow_atom_sharing=False)
    sg = CGSegmenter.from_pdb(cfg)
    entry = build_fcews_segment_data(sg.mol, sg.segments, "A_octane")

    assert entry["name"] == "A_octane"
    assert entry["mode"] == "FMO"
    n_frag = len(entry["atom"])
    assert n_frag >= 2, "octane@50 should split into >=2 fragments"
    # seg_info は全 atom (heavy+H) を漏れなく 1 回ずつ網羅
    n_atoms = sg.mol.GetNumAtoms()
    all_serials = sorted(s for frag in entry["seg_info"] for s in frag)
    assert all_serials == list(range(1, n_atoms + 1))
    # connect 本数 = fragment 境界 bond 数 = n_frag - 1 (chain)
    assert len(entry["connect"]) == n_frag - 1
    assert sum(entry["connect_num"]) == n_frag - 1


def test_connect_order_and_connect_num(octane_pdb):
    """connect = [BDA, BAA] (ABINIT-MP AJF 並び)、connect_num = BAA 数:
    第2要素 (BAA) が connect_num>0 fragment 内にある。"""
    from abmptools.fragmenter.cg_segmenter import CGSegmenter, CGSegmenterConfig
    from abmptools.fragmenter.cg_segmenter.fcews_export import build_fcews_segment_data

    cfg = CGSegmenterConfig(pdb_path=str(octane_pdb), target_mw=40.0,
                            allow_atom_sharing=False)
    sg = CGSegmenter.from_pdb(cfg)
    entry = build_fcews_segment_data(sg.mol, sg.segments, "A_octane")

    # connect_num の合計 = BAA 総数 = cut bond 数
    assert sum(entry["connect_num"]) == len(entry["connect"])
    # 各 fragment の connect_num = その fragment 内に BAA (connect 第2要素) がある数
    baa_count = [0] * len(entry["connect_num"])
    for bda, baa in entry["connect"]:
        bda_frag = _frag_index_of_serial(entry, bda)
        baa_frag = _frag_index_of_serial(entry, baa)
        assert bda_frag != baa_frag
        baa_count[baa_frag] += 1
    assert baa_count == entry["connect_num"], \
        "connect_num must equal per-fragment BAA count (2nd element of each connect)"


# ---------------------------------------------------------------------------
# config_read in-process round-trip (核心)
# ---------------------------------------------------------------------------


def test_config_read_roundtrip_octane(octane_pdb, tmp_path):
    import abmptools
    from abmptools.fragmenter.cg_segmenter import CGSegmenter, CGSegmenterConfig

    out_dir = tmp_path / "fcews_oct"
    cfg = CGSegmenterConfig(pdb_path=str(octane_pdb), target_mw=50.0,
                            allow_atom_sharing=False, output_dir=str(out_dir))
    sg = CGSegmenter.from_pdb(cfg)
    result = sg.export_fcews(name="A_octane")

    seg_path = Path(result["segment_data"])
    xyz_path = Path(result["xyz"])
    assert seg_path.exists() and xyz_path.exists()

    n_atoms = int(xyz_path.read_text().splitlines()[0])
    ai = abmptools.abinit_io()
    ai.mainpath = str(out_dir)
    conf = ai.config_read("A_octane", n_atoms)
    assert conf["name"] == "A_octane"       # default でなく実 entry を引けた
    assert conf["mode"] == "FMO"
    total = sum(len(frag) for frag in conf["seg_info"])
    assert total == n_atoms
    for pair in conf["connect"]:
        assert len(pair) == 2


def test_xyz_matches_seg_info(octane_pdb, tmp_path):
    from abmptools.fragmenter.cg_segmenter import CGSegmenter, CGSegmenterConfig

    out_dir = tmp_path / "fcews_oct2"
    cfg = CGSegmenterConfig(pdb_path=str(octane_pdb), target_mw=50.0,
                            allow_atom_sharing=False, output_dir=str(out_dir))
    sg = CGSegmenter.from_pdb(cfg)
    result = sg.export_fcews(name="A_octane")

    lines = Path(result["xyz"]).read_text().splitlines()
    n_atoms = int(lines[0])
    assert "A_octane" in lines[1]
    atom_lines = lines[2:2 + n_atoms]
    assert len(atom_lines) == n_atoms
    for al in atom_lines:
        toks = al.split()
        assert len(toks) == 4
        float(toks[1]); float(toks[2]); float(toks[3])

    seg = _parse_seg_data(Path(result["segment_data"]))[0]
    serials = sorted(s for frag in seg["seg_info"] for s in frag)
    assert serials == list(range(1, n_atoms + 1))


# ---------------------------------------------------------------------------
# cholesterol (fused rings) — partition + sharing detection
# ---------------------------------------------------------------------------


@pytest.mark.skipif(not CHOLESTEROL_PDB.exists(), reason="cholesterol sample missing")
def test_cholesterol_partition_roundtrip(tmp_path):
    import abmptools
    from abmptools.fragmenter.cg_segmenter import CGSegmenter, CGSegmenterConfig

    out_dir = tmp_path / "fcews_chol"
    cfg = CGSegmenterConfig(pdb_path=str(CHOLESTEROL_PDB), target_mw=200.0,
                            allow_atom_sharing=False, output_dir=str(out_dir))
    sg = CGSegmenter.from_pdb(cfg)
    result = sg.export_fcews(name="cholesterol")

    n_atoms = sg.mol.GetNumAtoms()
    seg = result["entry"]
    all_serials = sorted(s for frag in seg["seg_info"] for s in frag)
    assert all_serials == list(range(1, n_atoms + 1)), "must be a full partition"

    # round-trip
    ai = abmptools.abinit_io()
    ai.mainpath = str(out_dir)
    conf = ai.config_read("cholesterol", n_atoms)
    assert conf["name"] == "cholesterol"
    assert sum(len(f) for f in conf["seg_info"]) == n_atoms


@pytest.mark.skipif(not CHOLESTEROL_PDB.exists(), reason="cholesterol sample missing")
def test_atom_sharing_raises(tmp_path):
    """allow_atom_sharing=True (fused ring 共有) は FCEWS export でエラー。"""
    from abmptools.fragmenter.cg_segmenter import CGSegmenter, CGSegmenterConfig
    from abmptools.fragmenter.cg_segmenter.fcews_export import build_fcews_segment_data

    cfg = CGSegmenterConfig(pdb_path=str(CHOLESTEROL_PDB), target_mw=200.0,
                            allow_atom_sharing=True)
    sg = CGSegmenter.from_pdb(cfg)
    with pytest.raises(ValueError, match="shared"):
        build_fcews_segment_data(sg.mol, sg.segments, "cholesterol")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def test_cli_fcews(octane_pdb, tmp_path):
    from abmptools.fragmenter.cg_segmenter.__main__ import main
    out_dir = tmp_path / "cli_fcews"
    code = main([
        "fcews",
        "--pdb", str(octane_pdb),
        "--output-dir", str(out_dir),
        "--target-mw", "50",
        "--name", "A_octane",
    ])
    assert code == 0
    assert (out_dir / "segment_data.dat").exists()
    assert (out_dir / "A_octane.xyz").exists()
