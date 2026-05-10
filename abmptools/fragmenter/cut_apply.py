# -*- coding: utf-8 -*-
"""
abmptools.fragmenter.cut_apply
------------------------------
CutSite list を RDKit Mol に適用し、connected components に分解して
fragment ごとの atom indices / formal charge / BAA を計算する。
"""
from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any, List, Tuple

from .models import CutSite

logger = logging.getLogger(__name__)


@dataclass
class FragmentInfo:
    """1 fragment 分の情報 (1 mol 内、cut 後の連結成分)。

    Attributes
    ----------
    fragment_id
        分子内 fragment ID (0-origin、cut 結果の出現順)。
    atom_indices
        この fragment に含まれる atom の mol 内 index (0-origin)。
    charge
        formal charge の sum。
    bda_pairs
        この fragment が **BDA holder** 側のペア list。各要素は ``(bda_atom_in_this_frag,
        baa_atom_in_other_frag)`` で、ABINIT-MP の log2config 形式の `connect`
        と一致する (BDA fragment のみが connect entry を持ち、`[BDA, BAA]` 順)。
    baa_pairs
        この fragment が **BAA 受け側** のペア list。各要素は ``(baa_atom_in_this_frag,
        bda_atom_in_other_frag)``。ABINIT-MP の正規出力には含まれないが、
        検証・UI 表示用に保持。
    baa_atom_pairs (property, deprecated alias)
        旧 API 互換。`bda_pairs + baa_pairs` を返す。CutSite に bda/baa が
        未設定の legacy ケースでは全ペアが `baa_pairs` 側に入っている。
    """

    fragment_id: int
    atom_indices: List[int] = field(default_factory=list)
    charge: int = 0
    bda_pairs: List[Tuple[int, int]] = field(default_factory=list)
    baa_pairs: List[Tuple[int, int]] = field(default_factory=list)

    @property
    def baa_atom_pairs(self) -> List[Tuple[int, int]]:
        """Deprecated alias: `bda_pairs + baa_pairs` を返す。

        旧 API (BDA/BAA を区別していなかった頃) との互換のため残してあるが、
        新しいコードでは `bda_pairs` / `baa_pairs` を直接使うこと。
        """
        return list(self.bda_pairs) + list(self.baa_pairs)


def apply_cuts(mol: Any, cut_sites: List[CutSite]) -> List[FragmentInfo]:
    """cut_sites の enabled=True のみ反映し、fragment list を返す。

    Parameters
    ----------
    mol
        RDKit Mol (元のまま、破壊しない)。
    cut_sites
        切断点リスト。

    Returns
    -------
    List[FragmentInfo]
        cut 後の各 fragment の情報。enabled=True の cut が無ければ全 atom が
        1 fragment として返る。
    """
    from rdkit import Chem

    enabled_cuts = [cs for cs in cut_sites if cs.enabled]

    rwmol = Chem.RWMol(mol)
    for cs in enabled_cuts:
        bond = rwmol.GetBondBetweenAtoms(cs.atom1_idx, cs.atom2_idx)
        if bond is not None:
            rwmol.RemoveBond(cs.atom1_idx, cs.atom2_idx)
    cut_mol = rwmol.GetMol()

    components = Chem.GetMolFrags(cut_mol, asMols=False, sanitizeFrags=False)

    # atom_idx -> fragment_id mapping
    atom_to_frag: dict = {}
    for frag_id, atom_indices in enumerate(components):
        for a in atom_indices:
            atom_to_frag[a] = frag_id

    fragments: List[FragmentInfo] = []
    for frag_id, atom_indices in enumerate(components):
        atoms_list = sorted(atom_indices)
        charge_sum = sum(
            cut_mol.GetAtomWithIdx(i).GetFormalCharge() for i in atoms_list
        )
        # BDA / BAA 役割を CutSite から取得して this fragment への振り分け:
        #   - cs.bda_atom_idx が this fragment 内 → bda_pairs に (bda, baa)
        #   - cs.baa_atom_idx が this fragment 内 → baa_pairs に (baa, bda)
        #   - bda/baa 未設定 (None) の legacy CutSite → baa_pairs に対称記録
        #     (旧 API 挙動の保持)
        bda_pairs: List[Tuple[int, int]] = []
        baa_pairs: List[Tuple[int, int]] = []
        for cs in enabled_cuts:
            bda = cs.bda_atom_idx
            baa = cs.baa_atom_idx
            if bda is not None and baa is not None:
                f_bda = atom_to_frag.get(bda, -1)
                f_baa = atom_to_frag.get(baa, -1)
                if f_bda == frag_id and f_baa != frag_id:
                    bda_pairs.append((bda, baa))
                elif f_baa == frag_id and f_bda != frag_id:
                    baa_pairs.append((baa, bda))
            else:
                # Legacy: BDA/BAA 未設定 → 旧挙動の対称記録
                a1, a2 = cs.atom1_idx, cs.atom2_idx
                f1 = atom_to_frag.get(a1, -1)
                f2 = atom_to_frag.get(a2, -1)
                if f1 == frag_id and f2 != frag_id:
                    baa_pairs.append((a1, a2))
                elif f2 == frag_id and f1 != frag_id:
                    baa_pairs.append((a2, a1))
        fragments.append(FragmentInfo(
            fragment_id=frag_id,
            atom_indices=atoms_list,
            charge=int(charge_sum),
            bda_pairs=sorted(bda_pairs),
            baa_pairs=sorted(baa_pairs),
        ))

    logger.debug(
        "apply_cuts: %d enabled cuts -> %d fragments (atoms per frag: %s)",
        len(enabled_cuts),
        len(fragments),
        [len(f.atom_indices) for f in fragments],
    )
    return fragments
