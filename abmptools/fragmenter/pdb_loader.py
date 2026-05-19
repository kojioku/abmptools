# -*- coding: utf-8 -*-
"""
abmptools.fragmenter.pdb_loader
-------------------------------
PDB → RDKit Mol 変換と連結成分分解。

Bond perception の戦略:
    1. PDB CONECT があれば優先 (RDKit が自動利用)
    2. CONECT 不足なら sanitize=False で再試行
    3. それでも失敗なら openbabel (`obabel`) で前処理 (subprocess)

rdkit が import できない環境ではこのモジュール自体が import できない
(abmptools.fragmenter.__init__ は try/except で握る)。
"""
from __future__ import annotations

import logging
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, List

logger = logging.getLogger(__name__)


@dataclass
class LoadedMolecule:
    """1 つの連結成分 (= PDB 内の独立分子)。

    Attributes
    ----------
    mol
        RDKit Mol オブジェクト (bond perception 済み)。
    residue_name
        PDB の残基名 (見つからない場合は空文字列)。
    atom_indices_in_pdb
        元 PDB ファイル内での 0-origin atom index list。
    """

    mol: Any                # rdkit.Chem.Mol (型は遅延読込)
    residue_name: str = ""
    atom_indices_in_pdb: List[int] = field(default_factory=list)


def load_pdb_molecules(
    pdb_path: str,
    fallback_obabel: bool = True,
) -> List[LoadedMolecule]:
    """PDB をロードし、連結成分単位の分子リストとして返す。

    Parameters
    ----------
    pdb_path
        入力 PDB ファイルパス。
    fallback_obabel
        True なら obabel フォールバックを許可する。False なら 1-2 のみ試す。

    Returns
    -------
    List[LoadedMolecule]
        各連結成分が 1 要素として並ぶ。

    Raises
    ------
    ImportError
        rdkit がインストールされていない。
    FileNotFoundError
        pdb_path が存在しない。
    RuntimeError
        全ての bond perception 戦略が失敗した。
    """
    try:
        from rdkit import Chem
    except ImportError as e:
        raise ImportError(
            "rdkit is required for fragmenter. "
            "Install with: pip install 'abmptools[fragmenter]'"
        ) from e

    path = Path(pdb_path)
    if not path.exists():
        raise FileNotFoundError(f"PDB not found: {pdb_path}")

    # Strategy 1: 標準 (proximityBonding=True, sanitize=True)
    # CONECT があれば優先、不足分は座標近接から推定。
    # 近接推定が valence を超える場合 RDKit が "Explicit valence ..." を stderr に
    # 吐くが、その後 Strategy 2 で正しく load されるため harmless。ユーザーが
    # 不安にならないよう、Strategy 1 試行時のみ rdkit log を抑制する。
    from rdkit import RDLogger
    _rdkit_logger = RDLogger.logger()
    _prev_level = _rdkit_logger.level if hasattr(_rdkit_logger, "level") else None
    RDLogger.DisableLog("rdApp.error")
    try:
        mol = Chem.MolFromPDBFile(str(path), removeHs=False, sanitize=True)
    finally:
        RDLogger.EnableLog("rdApp.error")

    # Strategy 2: sanitize 失敗 → proximity bonding を切って CONECT 厳密モードで再試行
    # 近接推定で誤結合 (例: 隣接分子間の C 同士をつないでしまう) が原因のことが多い。
    if mol is None:
        logger.info(
            "Sanitize failed. Retrying with proximityBonding=False "
            "(CONECT-only mode)."
        )
        mol = Chem.MolFromPDBFile(
            str(path), removeHs=False, sanitize=True, proximityBonding=False
        )

    # Strategy 3: 両方失敗 → sanitize=False で best effort
    if mol is None:
        logger.warning(
            "Both proximityBonding strategies failed with sanitize=True. "
            "Retrying with sanitize=False."
        )
        mol = Chem.MolFromPDBFile(str(path), removeHs=False, sanitize=False)

    # Strategy 4: obabel フォールバック
    if mol is None and fallback_obabel and shutil.which("obabel"):
        logger.info("Falling back to obabel preprocessing.")
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_pdb = Path(tmpdir) / "preprocessed.pdb"
            subprocess.run(
                ["obabel", str(path), "-O", str(tmp_pdb), "-h"],
                check=True,
                capture_output=True,
            )
            mol = Chem.MolFromPDBFile(str(tmp_pdb), removeHs=False, sanitize=True)

    if mol is None:
        raise RuntimeError(
            f"Failed to load PDB after all fallback strategies: {pdb_path}"
        )

    # 連結成分に分解
    components = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    indices_per_component = Chem.GetMolFrags(mol, asMols=False, sanitizeFrags=False)

    loaded: List[LoadedMolecule] = []
    for comp_mol, atom_indices in zip(components, indices_per_component):
        residue_name = _extract_residue_name(comp_mol)
        loaded.append(LoadedMolecule(
            mol=comp_mol,
            residue_name=residue_name,
            atom_indices_in_pdb=list(atom_indices),
        ))

    logger.info(
        f"Loaded {len(loaded)} molecules (connected components) from {pdb_path}"
    )
    return loaded


def _extract_residue_name(mol) -> str:
    """RDKit Mol から残基名を抽出 (最初の atom の PDBResidueInfo から)。

    残基名が見つからない場合は空文字列を返す。
    """
    if mol.GetNumAtoms() == 0:
        return ""
    atom = mol.GetAtomWithIdx(0)
    info = atom.GetPDBResidueInfo()
    if info is None:
        return ""
    return info.GetResidueName().strip()
