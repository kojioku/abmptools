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


def _build_peptide_from_sequence_openff(
    *,
    sequence: str,
    name: str,
    output_dir: Path,
    cap_n: str = "",
    cap_c: str = "",
) -> Path:
    """Phase 2-B: PeptideBuilder で sequence → 3D extended chain PDB.

    Natural L-AA only (D-AA / 非標準残基は不可、 既存 PDB 入力経路で)。
    PeptideBuilder + Biopython は全 OS で pip install 可、 Windows native 動作。

    Cap (ACE/NME) 追加:
        現状実装は **cap なしの bare peptide PDB のみ生成**。 cap 必要なら
        ``spec.cap_n``/``cap_c`` を空にしておき、 既存の PDB 入力経路を使う
        か、 OpenFF Topology.add_residue で別途追加 (Phase 2-B+ で対応予定)。
    """
    try:
        import PeptideBuilder
        from PeptideBuilder import Geometry
        import Bio.PDB
    except ImportError as e:
        raise RuntimeError(
            "Phase 2-B sequence build requires PeptideBuilder + biopython.\n"
            "Install: pip install PeptideBuilder biopython"
        ) from e

    if cap_n or cap_c:
        logger.warning(
            "Phase 2-B: cap_n/cap_c (%s/%s) は PeptideBuilder 経路では "
            "未サポート、 bare peptide で build します。 cap 必要なら "
            "spec.pdb_path に cap 付き PDB を渡してください。",
            cap_n, cap_c,
        )

    # 各 1-letter code を Geometry に変換 (natural L-AA 20 種のみ)
    geo = Geometry.geometry(sequence[0])
    structure = PeptideBuilder.initialize_res(geo)
    for aa in sequence[1:]:
        structure = PeptideBuilder.add_residue(structure, Geometry.geometry(aa))
    # C 末 OXT
    PeptideBuilder.add_terminal_OXT(structure)

    pdb_path = output_dir / f"{name}_seq.pdb"
    io = Bio.PDB.PDBIO()
    io.set_structure(structure)
    io.save(str(pdb_path))
    return pdb_path


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
    from ..amorphous.molecule_prep import prepare_molecule, write_single_mol_pdb

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Phase 2-B: sequence → 3D PDB via PeptideBuilder if pdb_path 無し
    pdb_path = spec.pdb_path
    if not pdb_path:
        if not spec.sequence:
            raise ValueError(
                f"PeptideSpec[{spec.name}] requires either pdb_path or sequence."
            )
        pdb_path = _build_peptide_from_sequence_openff(
            sequence=spec.sequence,
            name=spec.name,
            output_dir=output_dir,
            cap_n=spec.cap_n,
            cap_c=spec.cap_c,
        )
        logger.info("Phase 2-B: built peptide '%s' from sequence '%s' → %s",
                    spec.name, spec.sequence, pdb_path)

    logger.info(
        "OpenFF route: build peptide '%s' from PDB=%s (charge_method=%s)",
        spec.name, spec.pdb_path, charge_method,
    )
    # amorphous.prepare_molecule の am1bcc は early return なので、 explicit
    # pre-assign が必要 (use_precomputed_charges=True と整合させるため)。
    if charge_method == "am1bcc":
        mol = prepare_molecule(
            pdb_path=str(pdb_path), name=spec.name, charge_method="",
        )
        mol.assign_partial_charges(partial_charge_method="am1bcc")
    else:
        mol = prepare_molecule(
            pdb_path=str(pdb_path),
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
