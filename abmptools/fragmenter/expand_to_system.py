# -*- coding: utf-8 -*-
"""
abmptools.fragmenter.expand_to_system
-------------------------------------
1 group の cut_sites を全コピーへ展開し、log2config 互換の segment_data.dat
を出力する。

log2config の出力形式 (abmptools.log2config.main):
    seg_data = [
        {
            'name': '<basename>',
            'atom':       [n_atoms_per_fragment, ...],
            'charge':     [charge_per_fragment, ...],
            'connect_num':[BAA_count_per_fragment, ...],
            'seg_info':   [[atom_idx_1origin, ...], ...],
            'connect':    [[(this_atom, other_atom), ...], ...],
            'nummol_seg': [1],
            'repeat':     [1],
            'pair_file':  [],
            'multi_xyz':  'none',
        },
    ]

新ツールでは「同じ SMILES の分子グループの全コピー」を 1 segment_data エントリ
にまとめる。各 fragment が seg_data 内 list の 1 要素になる。
"""
from __future__ import annotations

import logging
from dataclasses import asdict
from pathlib import Path
from typing import Any, Dict, List, Tuple

from .cut_apply import apply_cuts
from .models import FragmentResult, MoleculeGroup
from .pdb_loader import LoadedMolecule

logger = logging.getLogger(__name__)


def build_segment_data(
    pdb_basename: str,
    groups: List[MoleculeGroup],
    molecules: List[LoadedMolecule],
) -> Tuple[List[Dict[str, Any]], int]:
    """全 group の代表分子の cut_sites を、各メンバー分子に適用して系全体の
    segment_data 構造を組み立てる。

    Returns
    -------
    seg_data : List[Dict]
        log2config 互換 segment_data 構造。
    total_fragments : int
        全分子合計の fragment 数。
    """
    seg_data: List[Dict[str, Any]] = []
    total_fragments = 0

    for group in groups:
        atoms_pf: List[int] = []
        charges_pf: List[int] = []
        connect_num_pf: List[int] = []
        seg_info_pf: List[List[int]] = []
        connect_pf: List[List[Tuple[int, int]]] = []

        for mol_idx in group.member_mol_indices:
            lm = molecules[mol_idx]
            fragments = apply_cuts(lm.mol, group.cut_sites)
            for frag in fragments:
                # 各 fragment の atom 数 (heavy + H 全て、PDB 内に存在するもの)
                atoms_pf.append(len(frag.atom_indices))
                charges_pf.append(frag.charge)
                connect_num_pf.append(len(frag.baa_atom_pairs))
                # seg_info: PDB 内 atom index (1-origin)
                seg_info_pf.append([
                    lm.atom_indices_in_pdb[i] + 1 for i in frag.atom_indices
                ])
                # connect: (this_atom_pdb_1origin, other_atom_pdb_1origin) のリスト
                connect_list = [
                    (
                        lm.atom_indices_in_pdb[a_this] + 1,
                        lm.atom_indices_in_pdb[a_other] + 1,
                    )
                    for a_this, a_other in frag.baa_atom_pairs
                ]
                connect_pf.append(connect_list)
                total_fragments += 1

        # 1 group = 1 seg_data entry。fragmention済 1 分子分の繰り返しを atoms_pf
        # にすると n_copies × n_fragment_per_mol 個の要素が並ぶ。
        seg_entry = {
            "name": f"{pdb_basename}_grp{len(seg_data) + 1}_{_short_smiles(group.smiles)}",
            "atom": atoms_pf,
            "charge": charges_pf,
            "connect_num": connect_num_pf,
            "seg_info": seg_info_pf,
            "connect": connect_pf,
            "nummol_seg": [1],
            "repeat": [1],
            "pair_file": [],
            "multi_xyz": "none",
        }
        seg_data.append(seg_entry)

    return seg_data, total_fragments


def write_segment_data(seg_data: List[Dict[str, Any]], output_path: str) -> None:
    """segment_data を log2config と同等の Python リテラル形式で書き出す。

    log2config.py の出力フォーマットを踏襲し、pdb2fmo 等が ast.literal_eval で
    そのまま読み込める形にする。
    """
    out_path = Path(output_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as f:
        print("seg_data = [", file=f)
        for entry in seg_data:
            print("    {", file=f)
            print(f"    'name': '{entry['name']}',", file=f)
            print(f"    'atom': {entry['atom']},", file=f)
            print(f"    'charge': {entry['charge']},", file=f)
            print(f"    'connect_num': {entry['connect_num']},", file=f)
            print(f"    'seg_info': {entry['seg_info']},", file=f)
            print(f"    'connect': {entry['connect']},", file=f)
            print(f"    'nummol_seg': {entry['nummol_seg']},", file=f)
            print(f"    'repeat': {entry['repeat']},", file=f)
            print(f"    'pair_file': {entry['pair_file']},", file=f)
            print(f"    'multi_xyz': '{entry['multi_xyz']}'", file=f)
            print("    },", file=f)
        print("]", file=f)
    logger.info("Wrote segment_data.dat -> %s (%d entries)", out_path, len(seg_data))


def export_to_system(
    pdb_path: str,
    groups: List[MoleculeGroup],
    molecules: List[LoadedMolecule],
    output_path: str,
) -> FragmentResult:
    """suggest_cuts_for_groups 後に呼び出すエントリポイント。

    Parameters
    ----------
    pdb_path
        入力 PDB のパス (basename を name に使う)。
    groups
        cut_sites が埋まっている MoleculeGroup list。
    molecules
        load_pdb_molecules() の結果。
    output_path
        segment_data.dat の出力先。

    Returns
    -------
    FragmentResult
    """
    pdb_basename = Path(pdb_path).stem
    seg_data, total = build_segment_data(pdb_basename, groups, molecules)
    write_segment_data(seg_data, output_path)
    return FragmentResult(
        groups=groups,
        segment_data=seg_data,
        total_fragments=total,
    )


def _short_smiles(smiles: str, max_len: int = 16) -> str:
    """ファイル名安全に SMILES を短縮 (英数のみ、最大 max_len 文字)。"""
    safe = "".join(c for c in smiles if c.isalnum())
    return safe[:max_len] if safe else "mol"
