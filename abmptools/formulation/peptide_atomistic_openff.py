"""Peptide builder (Phase 1 OpenFF route — Windows native).

`abmptools.formulation` の `tleap` 依存 peptide stage を、 OpenFF Toolkit
(`Molecule.from_polymer_pdb` + Sage 2.x SMIRNOFF) で置換するモジュール。

Phase 1 scope:
- **既存 PDB 入力前提** (``PeptideSpec.pdb_path`` 必須)
  - sequence からの 3D build (PDBFixer 経由) は Phase 2 で追加予定
- **whole-peptide SMIRNOFF** で typing (= 既存 D 体 GAFF mode と同じ思想)
  - 商用 OK (Apache-2.0 互換)
  - 真の ff14SB SMIRNOFF (``openff-amber-ff-ports``) を natural L-AA に
    使う経路は Phase 2 で追加

依存 (全 OS install 可):
- ``openff-toolkit`` (BSD-3)
- ``openff-interchange`` (BSD-3) — typing は :mod:`.topology_openff` 側
- ``rdkit`` (BSD-3)
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any, List, Optional, Sequence, Tuple

from .models import PeptideSpec

logger = logging.getLogger(__name__)


__all__ = [
    "PeptideOpenFFResult",
    "build_peptide_openff",
    "check_openff_amber_ff_ports_available",
]


@dataclass
class PeptideOpenFFResult:
    """Output of :func:`build_peptide_openff`.

    Attributes
    ----------
    pdb_path
        単一 peptide の 3D 構造 PDB (packmol 入力)。
    molecule
        ``openff.toolkit.Molecule`` (conformer + 任意で charges 付き)。
        :func:`abmptools.formulation.topology_openff.merge_to_interchange`
        に渡す template。
    n_atoms_per_copy
        単一 peptide の atom 数。
    """

    pdb_path: Path
    molecule: Any
    n_atoms_per_copy: int


def check_openff_amber_ff_ports_available() -> None:
    """Raise RuntimeError if ``openff-amber-ff-ports`` is not installed.

    Phase 2 (ff14SB SMIRNOFF for natural L-AA) で必須。 Phase 1 では
    whole-peptide SMIRNOFF (Sage) 経路を取るので **本関数を呼ばない経路でも
    Phase 1 build は通る**。
    """
    try:
        import openff.amber_ff_ports  # noqa: F401
    except ImportError as e:
        raise RuntimeError(
            "openff-amber-ff-ports is required for the ff14SB peptide route "
            "(Phase 2). Install: pip install openff-amber-ff-ports."
        ) from e


def build_peptide_openff(
    *,
    spec: PeptideSpec,
    output_dir: Path,
    charge_method: str = "gasteiger",
) -> PeptideOpenFFResult:
    """Whole-peptide SMIRNOFF route — Windows native compatible.

    Parameters
    ----------
    spec
        Peptide specification。 **OpenFF route では ``spec.pdb_path`` 必須**。
        ``spec.sequence`` 単独からの 3D build は Phase 2 で対応。
    output_dir
        ``build/input/`` 等の出力 dir。 ``<name>.pdb`` を書き出す。
    charge_method
        - ``"gasteiger"`` (default): 軽量 + 安定 (大 peptide で sqm 発散
          回避、 既存 D-octreotide GAFF mode と同じ慣習)、 全 OS 対応
        - ``"nagl"``: NAGL の ML AM1-BCC (要 ``openff-nagl`` install)、 全 OS 対応
        - ``"am1bcc"``: 内部 ``sqm`` を要求、 Linux/macOS のみ、 large peptide
          では SCF 発散リスクあり

    Returns
    -------
    PeptideOpenFFResult
        PDB path + Molecule + atom 数。

    Raises
    ------
    ValueError
        ``spec.pdb_path`` が未指定 (Phase 1 では sequence からの build 未対応)。
    RuntimeError
        OpenFF Toolkit や RDKit が未 install。
    """
    if not spec.pdb_path:
        raise ValueError(
            f"PeptideSpec[{spec.name}].pdb_path is required for the OpenFF "
            f"route (Phase 1). Sequence-only build is planned for Phase 2 "
            f"(see docs/platform_support.md)."
        )

    from ..amorphous.molecule_prep import prepare_molecule, write_single_mol_pdb

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info(
        "OpenFF route: build peptide '%s' from PDB=%s (charge_method=%s)",
        spec.name, spec.pdb_path, charge_method,
    )
    # amorphous.prepare_molecule の am1bcc は early return なので、 explicit
    # pre-assign が必要 (use_precomputed_charges=True と整合させるため)。
    if charge_method == "am1bcc":
        mol = prepare_molecule(
            pdb_path=str(spec.pdb_path), name=spec.name, charge_method="",
        )
        mol.assign_partial_charges(partial_charge_method="am1bcc")
    else:
        mol = prepare_molecule(
            pdb_path=str(spec.pdb_path),
            name=spec.name,
            charge_method=charge_method,
        )

    pdb_path = output_dir / f"{spec.name}.pdb"
    write_single_mol_pdb(mol, str(pdb_path))

    return PeptideOpenFFResult(
        pdb_path=pdb_path,
        molecule=mol,
        n_atoms_per_copy=int(mol.n_atoms),
    )
