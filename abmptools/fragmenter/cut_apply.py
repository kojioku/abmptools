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
    baa_atom_pairs
        BAA (Bond Attachment Atom) ペア list。各要素は (this_frag_atom_idx,
        other_frag_atom_idx)。enabled=True の cut bond について、その bond の
        片端がこの fragment 内にあれば atom_idx (mol 内 0-origin) を記録する。
    """

    fragment_id: int
    atom_indices: List[int] = field(default_factory=list)
    charge: int = 0
    baa_atom_pairs: List[Tuple[int, int]] = field(default_factory=list)


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
        # BAA: enabled=True の cut bond の両端 atom について、片端がこの fragment
        # 内にあるなら (this_atom, other_atom) を記録
        baa_pairs: List[Tuple[int, int]] = []
        for cs in enabled_cuts:
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
            baa_atom_pairs=sorted(baa_pairs),
        ))

    logger.debug(
        "apply_cuts: %d enabled cuts -> %d fragments (atoms per frag: %s)",
        len(enabled_cuts),
        len(fragments),
        [len(f.atom_indices) for f in fragments],
    )
    return fragments
