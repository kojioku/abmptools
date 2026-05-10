# -*- coding: utf-8 -*-
"""
abmptools.fragmenter.models
---------------------------
データクラス定義 (rdkit 非依存、JSON 往復可能)。
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Set


@dataclass
class FragmenterConfig:
    """Fragmenter の入力設定。

    Attributes
    ----------
    pdb_path
        入力 PDB ファイル。
    output_dir
        出力ディレクトリ (segment_data.dat / review bundle 等)。
    target_mw
        C-C 切断の目安分子量 (g/mol)。鎖を辿って累積分子量がこの値を超えるたび
        次の切断対象 C-C で分割する。
    exclude_ring_cc
        環内 C-C 結合を切断対象から除外する。
    exclude_multibond
        多重結合 (二重・三重) を切断対象から除外する。
    exclude_heteroneighbor
        両端 C 原子の隣接にヘテロ原子 (N/O/S/P/F/Cl/Br/I) があれば除外する。
        例: アミド C-N 隣の C-C、エステル C-O 隣の C-C は切らない。
    skip_protein_dna
        タンパク質・DNA は対象外として extract のみ行う (既存 log2config 経路へ
        受け渡す前提)。
    polymer_groups
        γ 経路: 同一視する SMILES グループのリスト (要素は SMILES 文字列の list)。
        例: [["CCCCCCCCCC", "CCCCCCCCCCC"]] とすると PE N=10 と N=11 を同一視。
    include_residue_name
        SMILES に加え residue_name も group key に含める。同一 SMILES でも
        residue_name が違えば別グループに分ける。
    include_c_heteroatom
        C-X (X = N/O/S/P/F/Cl/Br/I) の単結合も切断 candidate に含める。
        デフォルト False (C-C 単結合のみ切断対象)。True にすると
        ``exclude_heteroneighbor`` フィルタは C-C bond にのみ適用され、
        C-X bond は別途切断対象。BDA/BAA 役割は ABINIT-MP の慣習に従い
        **C 側を BDA、X 側を BAA** に固定 (walk 方向に依存しない)。
    """

    pdb_path: str
    output_dir: str = "./fragmenter_out"
    target_mw: float = 200.0
    exclude_ring_cc: bool = True
    exclude_multibond: bool = True
    exclude_heteroneighbor: bool = True
    skip_protein_dna: bool = True
    polymer_groups: List[List[str]] = field(default_factory=list)
    include_residue_name: bool = False
    include_c_heteroatom: bool = False

    def to_json(self, path: Optional[str] = None) -> str:
        text = json.dumps(asdict(self), indent=2, ensure_ascii=False)
        if path:
            Path(path).write_text(text)
        return text

    @classmethod
    def from_json(cls, path: str) -> "FragmenterConfig":
        data = json.loads(Path(path).read_text())
        return cls(**data)


@dataclass
class CutSite:
    """1 結合に対する切断状態。

    Attributes
    ----------
    bond_idx
        RDKit Mol 内での bond index。
    atom1_idx, atom2_idx
        結合両端の原子 index (RDKit Mol 内 0-origin)。
    suggested
        auto_split の自動提案で生成されたか (True) ユーザー手動追加か (False)。
    enabled
        切断するか (True) 切断しないか (False)。UI で toggle。
    bda_atom_idx
        BDA (Bond Detached Atom) = 切断後 electron pair を保持する側の atom
        index。`atom1_idx` か `atom2_idx` のどちらか。``None`` の場合は未決定
        (legacy データや手動 cut で方向が定まっていないケース)。
    baa_atom_idx
        BAA (Bond Attachment Atom) = H 擬似 capping される側の atom index。
        `bda_atom_idx` の反対側。``None`` の場合は未決定。

    BDA / BAA の決定ルール (auto_split で suggest される際):
        - C-C 単結合 (本サブパッケージの主対象): walk 起点側 (graph diameter
          path[0] を含む側) を BDA 側とする
        - C-X 単結合 (将来 ``--no-exclude-hetero`` 時): C 側を BDA
        - peptide bond (本サブパッケージ対象外): N 末側を BDA (log2config 経路)
    """

    bond_idx: int
    atom1_idx: int
    atom2_idx: int
    suggested: bool = True
    enabled: bool = True
    bda_atom_idx: Optional[int] = None
    baa_atom_idx: Optional[int] = None

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "CutSite":
        return cls(**d)


@dataclass
class MoleculeGroup:
    """canonical SMILES でグループ化された分子集合。

    Attributes
    ----------
    smiles
        canonical SMILES (group key)。
    residue_names
        系内に現れた残基名 (補助情報、UI 表示用)。
    representative_mol_idx
        系内での代表分子 ID (0-origin)。UI で表示し、cut_sites を編集する対象。
    n_copies
        同分子の数。
    member_mol_indices
        全コピーの分子 ID list (0-origin)。apply_to_system で展開時に使用。
    cut_sites
        代表分子に対する切断点リスト (auto_split.suggest で埋める)。
    """

    smiles: str
    residue_names: Set[str] = field(default_factory=set)
    representative_mol_idx: int = 0
    n_copies: int = 0
    member_mol_indices: List[int] = field(default_factory=list)
    cut_sites: List[CutSite] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "smiles": self.smiles,
            "residue_names": sorted(self.residue_names),
            "representative_mol_idx": self.representative_mol_idx,
            "n_copies": self.n_copies,
            "member_mol_indices": list(self.member_mol_indices),
            "cut_sites": [cs.to_dict() for cs in self.cut_sites],
        }

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "MoleculeGroup":
        return cls(
            smiles=d["smiles"],
            residue_names=set(d.get("residue_names", [])),
            representative_mol_idx=d.get("representative_mol_idx", 0),
            n_copies=d.get("n_copies", 0),
            member_mol_indices=list(d.get("member_mol_indices", [])),
            cut_sites=[CutSite.from_dict(cs) for cs in d.get("cut_sites", [])],
        )


@dataclass
class FragmentResult:
    """最終分割結果。

    Attributes
    ----------
    groups
        全 MoleculeGroup。
    segment_data
        log2config 互換構造の list。各要素は dict で keys = ['name', 'atom',
        'charge', 'connect_num', 'seg_info', 'connect', 'nummol_seg',
        'repeat', 'pair_file', 'multi_xyz']。
    total_fragments
        系全体のフラグメント総数。
    """

    groups: List[MoleculeGroup]
    segment_data: List[Dict[str, Any]] = field(default_factory=list)
    total_fragments: int = 0
