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
    is_protein: bool = True


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


def _pdbfix_protein(input_pdb: str, output_pdb: str, ph: float = 7.0) -> str:
    """Phase 2-C: PDBFixer で protein PDB を OpenFF が読める形に前処理。

    - crystal water (HETATM HOH) 除去
    - missing residue / atom 補完
    - **explicit hydrogen 付加** (OpenFF Toolkit は H 必須)

    全て pip install 可能な OpenMM/PDBFixer のみ使用 → Windows native 動作。
    """
    try:
        from pdbfixer import PDBFixer
        from openmm.app import PDBFile
    except ImportError as e:
        raise RuntimeError(
            "Phase 2-C protein route requires pdbfixer + openmm.\n"
            "Install: pip install pdbfixer openmm"
        ) from e
    fixer = PDBFixer(filename=str(input_pdb))
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(ph)
    with open(output_pdb, "w") as fh:
        PDBFile.writeFile(fixer.topology, fixer.positions, fh)
    return str(output_pdb)


def build_peptide_openff(
    *,
    spec: PeptideSpec,
    output_dir: Path,
    charge_method: str = "gasteiger",
) -> PeptideOpenFFResult:
    """Protein route (Phase 2-C) — Windows native、 multi-chain + disulfide 対応。

    全 OS で動く科学スタック (PDBFixer + OpenFF ``Topology.from_pdb`` +
    ff14SB SMIRNOFF) で標準アミノ酸ペプチド / タンパク質を組む。

    処理:
    1. peptide PDB を取得 (``spec.pdb_path`` か、 sequence から PeptideBuilder)
    2. :func:`_pdbfix_protein` で water 除去 + 欠損補完 + H 付加
    3. ``Topology.from_pdb`` で読み込み (multi-chain は 1 分子として認識、
       **disulfide bond を自動検出**)
    4. その単一 Molecule を template として返す。 charges は付与せず、
       後段 :func:`merge_to_interchange` で ff14SB **library charges** を使う

    Parameters
    ----------
    spec
        Peptide specification (``pdb_path`` か ``sequence``)。
    output_dir
        ``build/input/`` 等の出力 dir。 fixed PDB ``<name>.pdb`` を書き出す
        (packmol 入力 = H 付き)。
    charge_method
        後方互換のため残置。 protein route では使われない (ff14SB library
        charges を使用)。

    Returns
    -------
    PeptideOpenFFResult
        fixed PDB path + Molecule template (charges なし) + atom 数 +
        ``is_protein=True``。
    """
    from openff.toolkit import Topology

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # 1. peptide PDB を取得 (pdb_path か sequence build)
    raw_pdb = spec.pdb_path
    if not raw_pdb:
        if not spec.sequence:
            raise ValueError(
                f"PeptideSpec[{spec.name}] requires either pdb_path or sequence."
            )
        raw_pdb = _build_peptide_from_sequence_openff(
            sequence=spec.sequence,
            name=spec.name,
            output_dir=output_dir,
            cap_n=spec.cap_n,
            cap_c=spec.cap_c,
        )
        logger.info("Phase 2-B: built peptide '%s' from sequence '%s'",
                    spec.name, spec.sequence)

    # 2. PDBFixer 前処理 (water 除去 + H 付加) → packmol 入力 PDB
    fixed_pdb = output_dir / f"{spec.name}.pdb"
    _pdbfix_protein(str(raw_pdb), str(fixed_pdb))
    logger.info(
        "Phase 2-C: PDBFixer preprocessed protein '%s' (water strip + H add) → %s",
        spec.name, fixed_pdb,
    )

    # 3. Topology.from_pdb で multi-chain + disulfide を認識した Molecule template
    top = Topology.from_pdb(str(fixed_pdb))
    if top.n_molecules != 1:
        logger.warning(
            "Topology.from_pdb('%s') yielded %d molecules; using the first. "
            "(multi-chain proteins joined by disulfides should be 1 molecule.)",
            spec.name, top.n_molecules,
        )
    mol = top.molecule(0)
    n_ss = sum(
        1 for b in mol.bonds
        if b.atom1.symbol == "S" and b.atom2.symbol == "S"
    )
    logger.info(
        "Phase 2-C: '%s' loaded via Topology.from_pdb: %d atoms, %d disulfide(s)",
        spec.name, mol.n_atoms, n_ss,
    )

    return PeptideOpenFFResult(
        pdb_path=fixed_pdb,
        molecule=mol,
        n_atoms_per_copy=int(mol.n_atoms),
        is_protein=True,
    )
