"""Small-molecule parameterizer (Phase 1 OpenFF route — Windows native).

`abmptools.formulation` の `acpype` / `sqm` / `antechamber` 依存を、 OpenFF
Sage 2.x (SMIRNOFF) + Interchange で置換するモジュール。 caprate (Na-C10) と
taurocholate を対象。

amorphous の :mod:`abmptools.amorphous.molecule_prep` の SMIRNOFF パターンを
formulation の small-molecule stage に薄く wrap した実装。 1 species ずつ
SMILES → 3D 配座 → OpenFF Molecule (charges pre-assigned) → PDB を生成する。

依存 (全 OS install 可):
- ``openff-toolkit`` (BSD-3)
- ``openff-interchange`` (BSD-3)
- ``rdkit`` (BSD-3)
- ``openff-nagl`` (optional、 ``charge_method="nagl"`` 時)
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional

from .models import BileSaltSpec, EnhancerSpec

logger = logging.getLogger(__name__)


__all__ = [
    "SmallMoleculeOpenFFResult",
    "parameterize_small_mol_openff",
]


@dataclass
class SmallMoleculeOpenFFResult:
    """Output for a single small-molecule species.

    Attributes
    ----------
    pdb_path
        単分子 3D PDB (packmol 入力)。
    molecule
        ``openff.toolkit.Molecule`` (conformer + charges 付き)。
        :func:`abmptools.formulation.topology_openff.merge_to_interchange`
        に渡す template。
    n_atoms_per_copy
        分子の atom 数。
    net_charge
        formal charge (cation = +1、 anion = -1 等)、 ions balance に使う。
    """

    pdb_path: Path
    molecule: Any
    n_atoms_per_copy: int
    net_charge: int


def parameterize_small_mol_openff(
    *,
    smiles: str,
    resname: str,
    net_charge: int,
    output_dir: Path,
    charge_method: str = "am1bcc",
) -> SmallMoleculeOpenFFResult:
    """SMILES → 3D 配座 → OpenFF Molecule (charged) + 1-frame PDB。

    Parameters
    ----------
    smiles
        構造 SMILES (charge 情報含む、 例: ``"CCCCCCCCCC(=O)[O-]"``)。
    resname
        GROMACS topology の resname (3-4 文字)、 例: ``"CPN"``、 ``"CPC"``。
    net_charge
        formal charge (acpype の ``-n`` と同義、 ions balance 用 metadata)。
    output_dir
        ``<resname>.pdb`` の出力 dir。
    charge_method
        - ``"am1bcc"`` (default): OpenFF の AM1-BCC、 small mol で安定、
          内部 ``sqm`` を要求するため Linux/macOS のみ
        - ``"nagl"``: NAGL の ML AM1-BCC (要 ``openff-nagl``)、 全 OS 対応
        - ``"gasteiger"``: 軽量 fallback (低精度)、 全 OS 対応
    """
    from ..amorphous.molecule_prep import prepare_molecule, write_single_mol_pdb

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info(
        "OpenFF route: parameterize %s (smiles=%s, charge_method=%s)",
        resname, smiles[:60] + ("..." if len(smiles) > 60 else ""), charge_method,
    )
    # amorphous.prepare_molecule は ``am1bcc`` を early return する
    # (Interchange downstream で sqm 計算する設計のため)。 formulation
    # OpenFF route では全 species を ``use_precomputed_charges`` で merge
    # するので、 明示的に pre-assign する必要がある。
    if charge_method == "am1bcc":
        # Build mol without pre-assign, then call assign_partial_charges('am1bcc')
        mol = prepare_molecule(smiles=smiles, name=resname, charge_method="")
        mol.assign_partial_charges(partial_charge_method="am1bcc")
    else:
        # gasteiger / nagl は prepare_molecule が pre-assign する
        mol = prepare_molecule(
            smiles=smiles, name=resname, charge_method=charge_method,
        )

    pdb_path = output_dir / f"{resname}.pdb"
    write_single_mol_pdb(mol, str(pdb_path))

    return SmallMoleculeOpenFFResult(
        pdb_path=pdb_path,
        molecule=mol,
        n_atoms_per_copy=int(mol.n_atoms),
        net_charge=net_charge,
    )
