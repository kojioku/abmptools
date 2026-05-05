# -*- coding: utf-8 -*-
"""
abmptools.fragmenter.polymer
----------------------------
γ 経路: 異なる SMILES のグループを「同じパターン」と明示的に同一視する。

ポリマーは鎖長が 1 残基違うだけで canonical SMILES が変わるため、
Auto-grouping だけでは N=10 と N=11 のポリエチレンが別々の MoleculeGroup
になってしまう。``declare_same_pattern`` でユーザーが明示的に「これらは
同じパターン」と宣言すると、cut_sites が最も多い (= 主鎖が最も長い)
group をマスターとして、他の group へ cut パターンを転送する。

転送ロジック:
    1. 各 group の代表分子の主鎖 (graph diameter = heavy-atom-only longest
       shortest path) を取得
    2. master の cut_site (atom1, atom2) が master_path のどの位置 (0-origin
       index) にあるかを調べる
    3. target_path で同じ index 位置にある atom ペアの bond で切断する
    4. target_path が短くて対応位置に達しない cut は捨てる
"""
from __future__ import annotations

import logging
from typing import List

from .auto_split import _find_main_chain_heavy
from .models import CutSite, MoleculeGroup
from .pdb_loader import LoadedMolecule

logger = logging.getLogger(__name__)


def declare_same_pattern(
    groups: List[MoleculeGroup],
    molecules: List[LoadedMolecule],
    same_pattern_smiles: List[List[str]],
) -> List[MoleculeGroup]:
    """指定された SMILES グループ間で cut_sites を同期する (in-place)。

    Parameters
    ----------
    groups
        suggest_cuts_for_groups 後の groups。
    molecules
        全分子リスト (representative_mol_idx の参照に使う)。
    same_pattern_smiles
        同一視するグループの SMILES リスト。
        例: [["CCCCCCCCCC", "CCCCCCCCCCC"]] = PE N=10 と N=11 を同一視。

    Returns
    -------
    List[MoleculeGroup]
        cut_sites が同期された groups (in-place 更新済)。
    """
    if not same_pattern_smiles:
        return groups

    smiles_to_group = {g.smiles: g for g in groups}
    for unify_set in same_pattern_smiles:
        targets = [smiles_to_group[s] for s in unify_set if s in smiles_to_group]
        if len(targets) < 2:
            logger.warning(
                "declare_same_pattern: less than 2 matching groups for %s "
                "(found groups: %s)",
                unify_set, [g.smiles for g in targets],
            )
            continue
        master = max(targets, key=lambda g: len(g.cut_sites))
        if not master.cut_sites:
            logger.warning(
                "declare_same_pattern: master group %s has no cut_sites "
                "(target_mw too high?). Skipping pattern transfer.",
                master.smiles,
            )
            continue
        for tgt in targets:
            if tgt is master:
                continue
            tgt.cut_sites = _transfer_cut_pattern(master, tgt, molecules)
            logger.info(
                "declare_same_pattern: transferred %d cut(s) from %s to %s",
                len(tgt.cut_sites), master.smiles, tgt.smiles,
            )
    return groups


def _transfer_cut_pattern(
    master: MoleculeGroup,
    target: MoleculeGroup,
    molecules: List[LoadedMolecule],
) -> List[CutSite]:
    """master の主鎖上の cut_sites を target の主鎖の対応位置に転送する。"""
    master_mol = molecules[master.representative_mol_idx].mol
    target_mol = molecules[target.representative_mol_idx].mol

    master_path = _find_main_chain_heavy(master_mol)
    target_path = _find_main_chain_heavy(target_mol)

    if len(master_path) < 2 or len(target_path) < 2:
        return []

    master_idx = {atom_idx: i for i, atom_idx in enumerate(master_path)}

    target_cuts: List[CutSite] = []
    for cs in master.cut_sites:
        if cs.atom1_idx not in master_idx or cs.atom2_idx not in master_idx:
            logger.debug(
                "Master cut (%d-%d) not on main chain; skipped.",
                cs.atom1_idx, cs.atom2_idx,
            )
            continue
        i1 = master_idx[cs.atom1_idx]
        i2 = master_idx[cs.atom2_idx]
        if abs(i1 - i2) != 1:
            logger.debug(
                "Master cut (%d-%d) not adjacent on main chain (path indices %d/%d); skipped.",
                cs.atom1_idx, cs.atom2_idx, i1, i2,
            )
            continue
        i_lo = min(i1, i2)
        if i_lo + 1 >= len(target_path):
            logger.debug(
                "Target chain too short (len=%d) for master cut at i=%d; skipped.",
                len(target_path), i_lo,
            )
            continue
        ta1 = target_path[i_lo]
        ta2 = target_path[i_lo + 1]
        bond = target_mol.GetBondBetweenAtoms(ta1, ta2)
        if bond is None:
            continue
        target_cuts.append(CutSite(
            bond_idx=bond.GetIdx(),
            atom1_idx=ta1,
            atom2_idx=ta2,
            suggested=True,
            enabled=cs.enabled,
        ))
    return target_cuts
