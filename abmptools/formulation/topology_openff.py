"""Topology builder (Phase 1 OpenFF route — Windows native).

複数 species の OpenFF Molecule template と packmol で組まれた mixed PDB から、
OpenFF Interchange を作って GROMACS ``system.top`` / ``system.gro`` を export。

設計:
- amorphous.parameterizer.create_interchange + export_gromacs を呼ぶ薄い wrapper
- ``Topology.from_pdb(mixture_pdb, unique_molecules=[...])`` で全 species を
  1 system に combine、 ``Interchange.from_smirnoff`` で SMIRNOFF apply
- ``charge_from_molecules=[...]`` で NAGL / Gasteiger 事前 charges を再利用
  (Windows native で ``sqm`` 不要)
- TIP3P SMIRNOFF (``tip3p.offxml``) を後ろに stack、 water 部分を上書き

Phase 1 制限:
- water + ions は **mixture PDB 側に既に packmol で入れてある前提** (= tleap
  solvatebox に代わる「pre-solvated packmol」経路)
- OpenMM Modeller.addSolvent 経由の wet solvate は Phase 2 で追加

依存 (全 OS install 可):
- ``openff-toolkit``、 ``openff-interchange``、 ``openmm``
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

logger = logging.getLogger(__name__)


def _rename_water_moltype_to_sol(top_path: Path) -> None:
    """OpenFF Interchange の water moltype 名 (``MOLnn``) を gmx 慣習 ``SOL``
    に rename + ``[ molecules ]`` の対応 entry も更新 + 1 mol 既存行を削除
    (gmx solvate が新たに追加する SOL 行と重複しないため)。

    検出ロジック: ``[ moleculetype ]`` block 内の ``[ atoms ]`` section で
    resname が ``HOH`` (or ``WAT``/``TIP3``) の moltype を water と判定。
    """
    text = top_path.read_text()
    # [ moleculetype ] block を抽出
    # 各 block の先頭 1 行 (semicolon 以外) が "name nrexcl" の form
    blocks = re.split(r"(\[ moleculetype \][^\n]*\n)", text)
    # blocks = ['preamble', '[ moleculetype ]\n', 'block1', '[ moleculetype ]\n', 'block2', ...]
    water_moltype_name = None
    for i in range(2, len(blocks), 2):
        block = blocks[i]
        atoms_match = re.search(
            r"\[ atoms \](.*?)(\n\[|\Z)", block, re.DOTALL,
        )
        if atoms_match is None:
            continue
        atoms_section = atoms_match.group(1)
        # 判定 1: resname に HOH/WAT/TIP3 が明示されている
        # 判定 2: atom 数 = 3 (TIP3P water は 3 atom)
        is_water_by_resname = bool(
            re.search(r"\b(HOH|WAT|TIP3|TIP3P)\b", atoms_section)
        )
        n_atoms = sum(
            1 for line in atoms_section.splitlines()
            if line.strip()
            and not line.strip().startswith(";")
            and len(line.split()) >= 4
        )
        if not (is_water_by_resname or n_atoms == 3):
            continue
        # この block が water、 先頭の name を取得
        first_data = next(
            (line for line in block.splitlines()
             if line.strip() and not line.strip().startswith(";")),
            None,
        )
        if first_data is None:
            continue
        old_name = first_data.split()[0]
        logger.info(
            "Phase 2-A: water moltype detected: %s (atoms=%d) → renamed to SOL",
            old_name, n_atoms,
        )
        # name 行を rename
        new_first_data = re.sub(
            r"^\s*" + re.escape(old_name) + r"\s",
            f"SOL              ",
            first_data, count=1,
        )
        blocks[i] = block.replace(first_data, new_first_data, 1)
        water_moltype_name = old_name
        break

    text = "".join(blocks)
    if water_moltype_name is not None:
        # [ molecules ] section の対応 entry も SOL に rename
        # 既存 1 mol を削除せず残す (gmx solvate が後で既存 1 mol + 追加 N mol
        # = 全体 SOL count を [ molecules ] に再 append する)
        lines = text.splitlines()
        in_molecules = False
        new_lines = []
        for line in lines:
            if re.match(r"\[ molecules \]", line):
                in_molecules = True
                new_lines.append(line)
                continue
            if in_molecules and re.match(r"\[", line):
                in_molecules = False
            if in_molecules and re.match(
                r"^\s*" + re.escape(water_moltype_name) + r"\s",
                line,
            ):
                # rename to SOL
                new_lines.append(
                    re.sub(
                        r"^\s*" + re.escape(water_moltype_name),
                        "SOL              ",
                        line, count=1,
                    )
                )
                continue
            new_lines.append(line)
        text = "\n".join(new_lines) + "\n"
    top_path.write_text(text)


__all__ = [
    "MergedTopologyResult",
    "merge_to_interchange",
    "export_gromacs_files",
    "solvate_and_neutralize_gmx",
]


@dataclass
class MergedTopologyResult:
    """Output of :func:`merge_to_interchange` + :func:`export_gromacs_files`."""

    gro_path: Path
    top_path: Path
    interchange: Any
    n_atoms_total: int


def merge_to_interchange(
    *,
    species_molecules: Sequence[Any],
    species_counts: Sequence[int],
    mixture_pdb: Path,
    box_size_nm: float,
    forcefield_offxmls: Sequence[str] = (
        "openff_unconstrained-2.1.0.offxml",
        "tip3p.offxml",
    ),
    use_precomputed_charges: bool = True,
) -> Any:
    """Merge per-species OpenFF Molecule templates + mixture PDB into 1 Interchange.

    Parameters
    ----------
    species_molecules
        Per-species ``openff.toolkit.Molecule`` の list。 順序は packmol
        input と一致させる ([peptide, enhancer_neu, enhancer_chg, bile_salt,
        ...] の順)。
    species_counts
        対応する ``n_copies``。 informational only (sanity check 用途)。
    mixture_pdb
        packmol で組まれた全 species の mixed PDB。 water/ions が同 packmol
        run で入っている場合も含む。
    box_size_nm
        cubic box edge (nm)。
    forcefield_offxmls
        SMIRNOFF OFFXML の stack。 default = Sage 2.1 + TIP3P (= amorphous
        と同じ慣習)。
    use_precomputed_charges
        ``True`` (default) なら each Molecule の ``partial_charges`` を Interchange
        に伝えて ``sqm`` を呼ばない (Windows native 必須)。

    Returns
    -------
    interchange
        Parameterized ``openff.interchange.Interchange``。
    """
    from ..amorphous.parameterizer import create_interchange

    logger.info(
        "OpenFF route: merge %d species (counts=%s) → Interchange "
        "(FF stack=%s, precomputed_charges=%s)",
        len(species_molecules), list(species_counts), list(forcefield_offxmls),
        use_precomputed_charges,
    )
    interchange = create_interchange(
        molecules=list(species_molecules),
        counts=list(species_counts),
        box_size_nm=float(box_size_nm),
        mixture_pdb=str(mixture_pdb),
        forcefield_name=list(forcefield_offxmls),
        use_precomputed_charges=use_precomputed_charges,
    )
    return interchange


def export_gromacs_files(
    *,
    interchange: Any,
    gro_path: Path,
    top_path: Path,
) -> MergedTopologyResult:
    """``Interchange`` → GROMACS ``system.gro`` / ``system.top``。

    Parameters
    ----------
    interchange
        :func:`merge_to_interchange` の戻り値。
    gro_path, top_path
        出力 path。

    Returns
    -------
    MergedTopologyResult
    """
    from ..amorphous.parameterizer import export_gromacs

    gro_path = Path(gro_path)
    top_path = Path(top_path)
    out = export_gromacs(
        interchange=interchange,
        gro_path=str(gro_path),
        top_path=str(top_path),
    )
    # atom 数を gro から数える (Interchange API に直接 accessor が無い)
    n_atoms = 0
    try:
        with gro_path.open() as f:
            f.readline()  # title
            n_atoms = int(f.readline().strip())
    except (OSError, ValueError):
        n_atoms = 0

    logger.info("OpenFF route: GROMACS export done (%s, %s, %d atoms)",
                out["gro"], out["top"], n_atoms)
    return MergedTopologyResult(
        gro_path=Path(out["gro"]),
        top_path=Path(out["top"]),
        interchange=interchange,
        n_atoms_total=n_atoms,
    )


def solvate_and_neutralize_gmx(
    *,
    gro_path: Path,
    top_path: Path,
    box_size_nm: float,
    salt_concentration_M: float = 0.15,
    gmx: str = "gmx",
    workdir: Optional[Path] = None,
) -> MergedTopologyResult:
    """`gmx solvate` + `gmx genion` で TIP3P 充填 + Joung-Cheatham NaCl 追加。

    Phase 2-A 実装 — Windows native 動作。 OpenFF Interchange の dry build
    (peptide + small mol + water 1 mol を template として含む) を、 gmx
    subprocess で本番量の water + ions に拡張する。

    前提:
    - ``top_path`` の topology に **water moltype が既に登録されている**こと。
      :func:`merge_to_interchange` の ``forcefield_offxmls`` に ``tip3p.offxml``
      を含め、 packmol input に water 1 個を入れて typing しておく
      (builder._build_openff の Stage 2 で自動配置)。

    Parameters
    ----------
    gro_path
        dry build の ``system.gro``。 in-place で上書きされる。
    top_path
        dry build の ``system.top``。 in-place で `[ molecules ]` 更新。
    box_size_nm
        cubic box edge (nm)。
    salt_concentration_M
        Joung-Cheatham NaCl 濃度 (default 0.15 M、 physiological)。
    gmx
        ``gmx`` 実行 path (default は PATH 解決)。
    workdir
        intermediate file 置き場 (None → gro_path の parent)。

    Returns
    -------
    MergedTopologyResult
        水/ions 追加後の ``system.gro`` / ``system.top`` を再 stat した結果。

    Raises
    ------
    RuntimeError
        ``gmx`` が PATH に無い、 または subprocess が失敗。
    """
    from ..trajectory.postprocess import _resolve_gmx, GmxError
    import subprocess
    import shutil
    import re

    gmx_exec = _resolve_gmx(gmx)
    gro_path = Path(gro_path).resolve()
    top_path = Path(top_path).resolve()
    workdir = Path(workdir).resolve() if workdir is not None else gro_path.parent

    # Stage 0: system.top の water moltype を "SOL" に rename
    # (gmx solvate は water moltype 名を "SOL" hardcoded で `[ molecules ]`
    # に追加するため、 OpenFF Interchange の自動命名 ``MOLnn`` と mismatch する。)
    _rename_water_moltype_to_sol(top_path)

    # Stage 1: gmx solvate で水充填
    # SPC216 (216 water box) を template に使う (GROMACS 標準)。
    # gmx solvate は top_path に `[ molecules ]` への水 entry を append する。
    solvated_gro = workdir / "system_solvated.gro"
    logger.info("Phase 2-A: gmx solvate → %s (box=%.2f nm)", solvated_gro, box_size_nm)
    cmd = [
        gmx_exec, "solvate",
        "-cp", str(gro_path),
        "-cs", "spc216.gro",   # GROMACS 標準 water template
        "-o", str(solvated_gro),
        "-p", str(top_path),
    ]
    proc = subprocess.run(cmd, capture_output=True, check=False)
    if proc.returncode != 0:
        raise RuntimeError(
            f"gmx solvate failed (exit {proc.returncode}):\n"
            f"--- stderr ---\n{proc.stderr.decode(errors='replace')}\n"
            f"Note: top_path must already contain a water moltype "
            f"(forcefield_offxmls must include tip3p.offxml, and the dry "
            f"packmol must include 1 water as a typing template)."
        )
    # 上書き
    shutil.copy(str(solvated_gro), str(gro_path))

    # Stage 2: gmx genion で Na+/Cl- を追加 (charge neutral + 0.15 M)
    # まず簡易 mdp で grompp して .tpr を作る
    em_mdp = workdir / "_genion_em.mdp"
    em_mdp.write_text(
        "integrator = steep\nnsteps = 0\n"
        "cutoff-scheme = Verlet\nrcoulomb = 1.0\nrvdw = 1.0\nrlist = 1.0\n"
    )
    tpr = workdir / "_genion.tpr"
    logger.info("Phase 2-A: gmx grompp for genion")
    cmd = [
        gmx_exec, "grompp",
        "-f", str(em_mdp), "-c", str(gro_path), "-p", str(top_path),
        "-o", str(tpr), "-maxwarn", "10",
    ]
    proc = subprocess.run(cmd, capture_output=True, check=False)
    if proc.returncode != 0:
        raise RuntimeError(
            f"gmx grompp (genion) failed:\n"
            f"{proc.stderr.decode(errors='replace')}"
        )

    # genion 用の ndx を make_ndx で生成 (water group "SOL" を含む)
    genion_ndx = workdir / "_genion.ndx"
    cmd = [gmx_exec, "make_ndx", "-f", str(gro_path), "-o", str(genion_ndx)]
    proc = subprocess.run(cmd, input=b"q\n", capture_output=True, check=False)
    if proc.returncode != 0:
        raise RuntimeError(
            f"gmx make_ndx failed:\n{proc.stderr.decode(errors='replace')}"
        )

    logger.info("Phase 2-A: gmx genion → 0.15 M NaCl + neutralize")
    out_gro = workdir / "system_ionized.gro"
    cmd = [
        gmx_exec, "genion",
        "-s", str(tpr), "-o", str(out_gro), "-p", str(top_path),
        "-n", str(genion_ndx),
        "-pname", "NA", "-nname", "CL",
        "-neutral",
        "-conc", str(salt_concentration_M),
    ]
    # gmx genion は対話的に group を要求 → "SOL" (water group from make_ndx)
    proc = subprocess.run(
        cmd, input=b"SOL\n", capture_output=True, check=False,
    )
    if proc.returncode != 0:
        # Log stderr (some gmx versions emit OK on stderr) but only fail if
        # output gro doesn't exist
        stderr = proc.stderr.decode(errors="replace")
        logger.warning("gmx genion stderr: %s", stderr[-500:])
        if not out_gro.is_file():
            raise RuntimeError(f"gmx genion failed:\n{stderr}")
    shutil.copy(str(out_gro), str(gro_path))

    # atom 数を更新
    n_atoms = 0
    try:
        with gro_path.open() as f:
            f.readline()
            n_atoms = int(f.readline().strip())
    except (OSError, ValueError):
        pass

    logger.info("Phase 2-A complete: wet system %d atoms", n_atoms)
    return MergedTopologyResult(
        gro_path=gro_path,
        top_path=top_path,
        interchange=None,
        n_atoms_total=n_atoms,
    )
