# -*- coding: utf-8 -*-
"""
abmptools.fragmenter.cg_segmenter.fcews_export
----------------------------------------------
cg_segmenter の segment 分割結果を **FCEWS** がそのまま読める
``segment_data.dat`` + monomer ``.xyz`` に変換する。

cg_segmenter は元々 CG 粒子用に「物理的に分割して H/CH3 cap を付ける」ツールだが、
その **segment 分割 (ring 検出 + chain の MW walk)** はそのまま FMO フラグメント
定義に使える。本モジュールは cap を使わず、segment の **atom partition** から
FCEWS の segment_data.dat (``config_read`` が読む形式) を組み立てる。

FCEWS form の規約 (``abmptools.abinit_io.config_read`` / FCEWS ``gen_coord`` が
そのまま AJF へ書く):

  - 1 化学種 = 1 entry。``mode='FMO'``
  - atom 番号は **monomer 内 local 1-origin** (= mol atom idx + 1)。同じ番号体系で
    monomer ``.xyz`` を mol idx 順に書くので seg_info と xyz 行が一致する
  - ``seg_info`` : per-fragment の atom serial (heavy + H、1-origin)
  - ``connect``  : 各 cut bond を ``[BDA_atom, BAA_atom]`` で記録した **flat list**
  - ``connect_num`` : 各 fragment が持つ **BAA atom 数**

**connect の atom 順序は ``[BDA, BAA]``** (ABINIT-MP の AJF fragment 接続情報
セクションの並び)。**``connect_num`` は BAA 数** = 各 cut bond の BAA を含む
fragment にのみ加算される (LOGManager の ``fbaas`` と同じ。AJF の ``Frag.`` 列 =
BAA の所属 fragment)。したがって ``connect`` の **第2要素 (BAA) が
``connect_num>0`` の fragment 内**にある。FCEWS の手書き test-data (nafion:
``connect=[[4,7]]``, ``connect_num=[0,1]`` で BAA=7 が frag2 内) に一致する。

BDA / BAA の **役割** は abmptools の ``auto_split.decide_bda_baa_for_manual_cut``
(C-X → C 側 BDA / C-C → 若い atom idx を BDA) を流用する。

FMO フラグメントは atom を共有できないため、cg_segmenter は
``allow_atom_sharing=False`` (= partition) で動かす必要がある。共有が検出された
場合はエラーにする。
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from ..auto_split import decide_bda_baa_for_manual_cut
from .models import Segment

logger = logging.getLogger(__name__)


def build_fcews_segment_data(
    mol: Any,
    segments: List[Segment],
    name: str,
    mode: str = "FMO",
) -> Dict[str, Any]:
    """1 分子の cg_segmenter 分割結果を FCEWS form の segment_data entry にする。

    Parameters
    ----------
    mol
        RDKit Mol (H 含む、conformer 不要だがあれば xyz に流用)。
    segments
        CGSegmenter.segments。**atom を共有していないこと** (FMO fragment は
        partition 必須。``allow_atom_sharing=False`` で生成する)。
    name
        entry 名 (= 対応する monomer ``<name>.xyz`` の basename と一致させる)。
    mode
        ABINIT-MP segment mode。通常 'FMO'。

    Returns
    -------
    Dict
        config_read が読む FCEWS form の entry。

    Raises
    ------
    ValueError
        segment が atom を共有している / heavy atom が未カバー / segment が空。
    """
    heavy_atoms = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() != 1]

    # heavy atom -> fragment index (= segments list の位置) を作りつつ共有検出
    atom_to_frag: Dict[int, int] = {}
    for fi, seg in enumerate(segments):
        for a in seg.atom_indices:
            if mol.GetAtomWithIdx(a).GetAtomicNum() == 1:
                continue  # 念のため (segment は heavy のみのはず)
            if a in atom_to_frag:
                raise ValueError(
                    f"atom {a} is shared between fragment {atom_to_frag[a]} and "
                    f"{fi}. FCEWS FMO fragments cannot share atoms; run "
                    "CGSegmenter with allow_atom_sharing=False."
                )
            atom_to_frag[a] = fi

    missing = [a for a in heavy_atoms if a not in atom_to_frag]
    if missing:
        raise ValueError(
            f"heavy atoms not covered by any segment: {missing}. "
            "FCEWS export requires a full partition."
        )

    n_frag = len(segments)
    if n_frag == 0:
        raise ValueError("no segments to export")

    # H atom を親 heavy atom の fragment に割り当てる
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 1:
            continue
        h_idx = atom.GetIdx()
        heavy_parent = None
        for nb in atom.GetNeighbors():
            if nb.GetAtomicNum() != 1:
                heavy_parent = nb.GetIdx()
                break
        if heavy_parent is None:
            logger.warning("H atom %d has no heavy neighbor; assigned to frag 0.", h_idx)
            atom_to_frag[h_idx] = 0
        else:
            atom_to_frag[h_idx] = atom_to_frag[heavy_parent]

    # seg_info (per-fragment、1-origin local、heavy + H) / charge / atom
    seg_info: List[List[int]] = [[] for _ in range(n_frag)]
    for atom_idx, fi in atom_to_frag.items():
        seg_info[fi].append(atom_idx + 1)
    for fi in range(n_frag):
        seg_info[fi].sort()

    charge_pf: List[int] = [0] * n_frag
    for atom_idx, fi in atom_to_frag.items():
        charge_pf[fi] += mol.GetAtomWithIdx(atom_idx).GetFormalCharge()

    atom_pf = [len(s) for s in seg_info]

    # cut bond: heavy-heavy bond で fragment 境界をまたぐもの
    connect_num_pf = [0] * n_frag
    connect_flat: List[List[int]] = []
    seen: set = set()
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        if a1.GetAtomicNum() == 1 or a2.GetAtomicNum() == 1:
            continue
        i, j = a1.GetIdx(), a2.GetIdx()
        if atom_to_frag[i] == atom_to_frag[j]:
            continue  # 同一 fragment 内の bond は cut でない
        key = (min(i, j), max(i, j))
        if key in seen:
            continue
        seen.add(key)
        # BDA / BAA の役割を abmptools の decider で決定
        bda, baa = decide_bda_baa_for_manual_cut(mol, i, j)
        # connect は [BDA, BAA] (ABINIT-MP AJF の fragment 接続情報の並び)。
        connect_flat.append([bda + 1, baa + 1])
        # connect_num は **BAA 数** = BAA を含む fragment に加算 (FCEWS / nafion
        # 規約。LOGManager の fbaas と同じ。AJF の "Frag." 列 = BAA の所属 fragment)
        connect_num_pf[atom_to_frag[baa]] += 1

    # connect は deterministic 順 (BDA atom の 1-origin 番号で sort)
    connect_flat.sort(key=lambda p: p[0])

    entry = {
        "name": name,
        "mode": mode,
        "atom": atom_pf,
        "charge": charge_pf,
        "connect_num": connect_num_pf,
        "connect": connect_flat,
        "seg_info": seg_info,
    }
    logger.info(
        "FCEWS entry '%s': %d fragment(s), atom=%s, connect=%s",
        name, n_frag, atom_pf, connect_flat,
    )
    return entry


def write_fcews_segment_data(
    entries: List[Dict[str, Any]],
    output_path: str,
) -> None:
    """FCEWS form の segment_data.dat を書き出す (config_read 互換リテラル)。

    複数 entry (= 複数化学種) を 1 ファイルに並べられる。
    """
    out_path = Path(output_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    header = (
        "# segment_data.dat for FCEWS\n"
        "# generated by abmptools.fragmenter.cg_segmenter.fcews_export\n"
        "#\n"
        "# 'mode': 'FMO'         -> FMO, charge, connect_num, connect, seg_info\n"
        "# 'mode': 'solv'        -> no FMO, nummol_seg\n"
        "# 'mode': 'charge_only' -> no FMO, charge only\n"
        "#\n"
        "# atom        : per-fragment atom count (1 monomer)\n"
        "# charge      : per-fragment formal charge\n"
        "# connect_num : per-fragment BAA atom count\n"
        "# connect     : flat list of [BDA_atom, BAA_atom] (1-origin, local)\n"
        "#               BAA (2nd elem) is in the connect_num>0 fragment\n"
        "# seg_info    : per-fragment atom serials (1-origin, local)\n"
        "# NOTE: each entry 'name' must match its monomer <name>.xyz basename.\n"
        "#\n"
    )
    with out_path.open("w") as f:
        f.write(header)
        print("seg_data = [", file=f)
        for entry in entries:
            print("    {", file=f)
            print(f"    'name': '{entry['name']}',", file=f)
            print(f"    'mode': '{entry['mode']}',", file=f)
            print(f"    'atom': {entry['atom']},", file=f)
            print(f"    'charge': {entry['charge']},", file=f)
            print(f"    'connect_num': {entry['connect_num']},", file=f)
            print(f"    'connect': {entry['connect']},", file=f)
            print(f"    'seg_info': {entry['seg_info']},", file=f)
            print("    },", file=f)
        print("]", file=f)
    logger.info(
        "Wrote FCEWS segment_data.dat -> %s (%d entries)", out_path, len(entries)
    )


def write_molecule_xyz(
    mol: Any,
    name: str,
    output_path: str,
    embed_if_missing: bool = True,
) -> None:
    """分子全体の座標を FCEWS 互換 .xyz で書き出す (mol atom idx 順)。

    atom 行は **mol idx 順** (0..N-1) なので、segment_data の ``seg_info``
    (= idx+1) の atom 番号と xyz 行番号が一致する。

    形式::

        <n_atoms>
                       @@<name>     @@
        <elem>   <x>   <y>   <z>
        ...
    """
    work = mol
    conf = None
    if work.GetNumConformers() > 0:
        conf = work.GetConformer()
    elif embed_if_missing:
        work = _embed_3d(mol)
        conf = work.GetConformer()
    else:
        raise ValueError(f"mol for '{name}' has no conformer")

    out_path = Path(output_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    n = work.GetNumAtoms()
    lines = [str(n), f"               @@{name}     @@"]
    for atom in work.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        lines.append(
            f"{atom.GetSymbol():<2s} {pos.x:14.5f} {pos.y:14.5f} {pos.z:14.5f}"
        )
    out_path.write_text("\n".join(lines) + "\n")
    logger.info("Wrote molecule xyz -> %s (%d atoms)", out_path, n)


def _embed_3d(mol: Any) -> Any:
    """conformer の無い Mol を 3D 化する (ETKDGv3 + MMFF/UFF)。"""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    work = Chem.Mol(mol)
    if not any(a.GetAtomicNum() == 1 for a in work.GetAtoms()):
        work = Chem.AddHs(work)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    if AllChem.EmbedMolecule(work, params) != 0:
        params.useRandomCoords = True
        AllChem.EmbedMolecule(work, params)
    try:
        AllChem.MMFFOptimizeMolecule(work)
    except Exception:
        try:
            AllChem.UFFOptimizeMolecule(work)
        except Exception:
            logger.warning("MMFF/UFF optimization failed; using raw embed.")
    return work


def export_fcews(
    mol: Any,
    segments: List[Segment],
    output_dir: str,
    name: str,
    mode: str = "FMO",
) -> Dict[str, Any]:
    """1 分子の FCEWS 入力 (segment_data.dat + <name>.xyz) を出力する。

    Returns
    -------
    Dict
        {'segment_data': path, 'xyz': path, 'entry': dict}
    """
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    entry = build_fcews_segment_data(mol, segments, name, mode=mode)
    seg_path = out_dir / "segment_data.dat"
    write_fcews_segment_data([entry], str(seg_path))

    xyz_path = out_dir / f"{name}.xyz"
    write_molecule_xyz(mol, name, str(xyz_path))

    return {
        "segment_data": str(seg_path),
        "xyz": str(xyz_path),
        "entry": entry,
    }
