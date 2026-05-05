# -*- coding: utf-8 -*-
"""
abmptools.fragmenter.grouping
-----------------------------
canonical SMILES によるグループ化。

PDB の連結成分ごとに canonical SMILES を計算し、同じ SMILES の分子を
1 つの MoleculeGroup にまとめる。代表分子を 1 つだけ UI に表示し、
ユーザーの編集結果を全コピーへ展開できるようにする。
"""
from __future__ import annotations

import logging
from collections import OrderedDict
from typing import List, Tuple

from .models import MoleculeGroup
from .pdb_loader import LoadedMolecule

logger = logging.getLogger(__name__)


def group_by_smiles(
    molecules: List[LoadedMolecule],
    include_residue_name: bool = False,
) -> List[MoleculeGroup]:
    """canonical SMILES で分子をグループ化する。

    Parameters
    ----------
    molecules
        load_pdb_molecules() の結果。
    include_residue_name
        True なら (smiles, residue_name) のタプルを group key に使う。
        False なら smiles のみを key に使う (異なる残基名でも統合)。

    Returns
    -------
    List[MoleculeGroup]
        各グループの cut_sites は空リスト (auto_split.suggest で埋める)。
    """
    try:
        from rdkit import Chem
    except ImportError as e:
        raise ImportError(
            "rdkit is required for fragmenter. "
            "Install with: pip install 'abmptools[fragmenter]'"
        ) from e

    groups: "OrderedDict[Tuple, MoleculeGroup]" = OrderedDict()

    for mol_idx, lm in enumerate(molecules):
        # Heavy atom only の SMILES を group key に使う (可読性 + H 補完の差を吸収)。
        # cut_sites 計算で使う元 Mol (H 含む) は LoadedMolecule.mol に保持されたまま。
        # sanitize 不完全な Mol だと RemoveHs が valence エラーを投げる場合があるので
        # fallback で元 Mol からの SMILES に切り替える。
        try:
            heavy_mol = Chem.RemoveHs(lm.mol)
            smiles = Chem.MolToSmiles(heavy_mol, canonical=True)
        except Exception as e:
            logger.warning(
                "RemoveHs failed for mol %d (%s). Using original mol with explicit H.",
                mol_idx, e,
            )
            smiles = Chem.MolToSmiles(lm.mol, canonical=True)
        key = (smiles, lm.residue_name) if include_residue_name else (smiles,)

        if key not in groups:
            groups[key] = MoleculeGroup(
                smiles=smiles,
                representative_mol_idx=mol_idx,
                n_copies=0,
                member_mol_indices=[],
            )
        g = groups[key]
        g.member_mol_indices.append(mol_idx)
        g.n_copies += 1
        if lm.residue_name:
            g.residue_names.add(lm.residue_name)

    result = list(groups.values())
    logger.info(
        "Grouped %d molecules into %d groups (n_copies: %s)",
        len(molecules),
        len(result),
        [g.n_copies for g in result],
    )
    return result
