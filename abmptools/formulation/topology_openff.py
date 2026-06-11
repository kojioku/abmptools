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
    "build_protein_route_topology",
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


# ff14SB SMIRNOFF port (openff-amber-ff-ports) — protein 用 library charges
FF14SB_OFFXML = "ff14sb_off_impropers_0.0.4.offxml"
SAGE_OFFXML = "openff_unconstrained-2.1.0.offxml"
TIP3P_OFFXML = "tip3p.offxml"


def merge_to_interchange(
    *,
    species_molecules: Sequence[Any],
    species_counts: Sequence[int],
    mixture_pdb: Path,
    box_size_nm: float,
    forcefield_offxmls: Sequence[str] = (
        SAGE_OFFXML,
        TIP3P_OFFXML,
    ),
    use_precomputed_charges: bool = True,
    protein_flags: Optional[Sequence[bool]] = None,
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
    protein_flags
        各 species が **protein か否か** の bool list (None なら全 False)。
        protein は ff14SB SMIRNOFF (openff-amber-ff-ports) の **library charges**
        を使い、 ``charge_from_molecules`` には含めない。 small molecule は
        precomputed charges (gasteiger / nagl) を使う。 protein が 1 つでも
        あれば FF stack の先頭に ``ff14sb`` を prepend する。

    Returns
    -------
    interchange
        Parameterized ``openff.interchange.Interchange``。
    """
    molecules = list(species_molecules)
    flags = list(protein_flags) if protein_flags is not None else [False] * len(molecules)
    has_protein = any(flags)

    if not has_protein:
        # 従来経路 (全 small molecule): amorphous の create_interchange に委譲
        from ..amorphous.parameterizer import create_interchange
        logger.info(
            "OpenFF route: merge %d species (counts=%s) → Interchange "
            "(FF stack=%s, precomputed_charges=%s)",
            len(molecules), list(species_counts), list(forcefield_offxmls),
            use_precomputed_charges,
        )
        return create_interchange(
            molecules=molecules,
            counts=list(species_counts),
            box_size_nm=float(box_size_nm),
            mixture_pdb=str(mixture_pdb),
            forcefield_name=list(forcefield_offxmls),
            use_precomputed_charges=use_precomputed_charges,
        )

    # protein を含む経路 (Phase 2-C): ff14SB + Sage + TIP3P stack、
    # mixed charges (protein=FF library, small mol=precomputed)
    from openff.toolkit import ForceField, Topology
    from openff.interchange import Interchange
    from openff.units import unit as off_unit
    import numpy as np

    offxmls = [FF14SB_OFFXML] + list(forcefield_offxmls)
    logger.info(
        "OpenFF route (protein): merge %d species (%d protein) → Interchange "
        "(FF stack=%s)",
        len(molecules), sum(flags), offxmls,
    )
    ff = ForceField(*offxmls)
    topology = Topology.from_pdb(str(mixture_pdb), unique_molecules=molecules)
    topology.box_vectors = np.eye(3) * float(box_size_nm) * off_unit.nanometer

    # small molecule (non-protein) のみ precomputed charges を渡す。
    # protein は ff14SB の LibraryCharges が自動適用される。
    charge_mols = [m for m, is_prot in zip(molecules, flags) if not is_prot]
    kwargs = {"force_field": ff, "topology": topology}
    if charge_mols:
        kwargs["charge_from_molecules"] = charge_mols
    return Interchange.from_smirnoff(**kwargs)


def build_protein_route_topology(
    *,
    species_molecules: Sequence[Any],
    species_counts: Sequence[int],
    protein_flags: Sequence[bool],
    mixture_pdb: Path,
    box_size_nm: float,
    gro_path: Path,
    top_path: Path,
    gmx: str = "gmx",
    forcefield_offxmls: Sequence[str] = (SAGE_OFFXML, TIP3P_OFFXML),
) -> MergedTopologyResult:
    """Protein 系の topology を **単一コピー parametrize + count 複製** で組む。

    巨大な protein 系で ``Interchange.from_smirnoff`` を full mixture に
    かけると nonbonded exception 生成が O(N²) で爆発する (insulin × 6 で
    2 時間+)。 対策として:

    1. **各 unique species 1 個ずつ**の topology を parametrize (~数秒)
    2. ``to_top`` で moleculetype 定義を得る
    3. ``[ molecules ]`` の count を実際の値に書き換え (GROMACS は同一
       moleculetype を count 参照する設計なので、 これで複製と等価)
    4. ``.gro`` は packmol mixture から ``gmx editconf`` で生成 (全コピーの
       座標)

    species_molecules / species_counts / protein_flags は packmol 順
    (peptide → enhancer → bile salt → water) に一致していること。
    """
    import subprocess
    from openff.toolkit import ForceField, Topology, Molecule
    from openff.interchange import Interchange
    from openff.units import unit as off_unit
    import numpy as np
    from ..trajectory.postprocess import _resolve_gmx

    molecules = list(species_molecules)
    counts = list(species_counts)
    flags = list(protein_flags)
    gro_path = Path(gro_path)
    top_path = Path(top_path)
    n_species = len(molecules)

    # ion moleculetype を生成するため Na+/Cl- を 1 個ずつ追加で parametrize
    # (mixture/.gro には含まれない → [ molecules ] には載せず、 定義だけ残す。
    # gmx genion が後で実 count を [ molecules ] に追記する。)
    na = Molecule.from_smiles("[Na+]"); na.generate_conformers(n_conformers=1); na.name = "NA"
    cl = Molecule.from_smiles("[Cl-]"); cl.generate_conformers(n_conformers=1); cl.name = "CL"
    parametrize_mols = molecules + [na, cl]

    # 1. 単一コピー topology
    single = Topology()
    for m in parametrize_mols:
        single.add_molecule(m)
    single.box_vectors = np.eye(3) * float(box_size_nm) * off_unit.nanometer

    offxmls = [FF14SB_OFFXML] + list(forcefield_offxmls)
    logger.info(
        "OpenFF protein route: single-copy parametrize (%d species + Na/Cl, "
        "%d protein, FF=%s)", n_species, sum(flags), offxmls,
    )
    ff = ForceField(*offxmls)
    # small molecule のみ precomputed charges。 protein は ff14SB library、
    # Na/Cl は ff (tip3p offxml) の ion charges を使う (charge_from_molecules 不要)
    charge_mols = [m for m, p in zip(molecules, flags) if not p]
    kwargs = {"force_field": ff, "topology": single}
    if charge_mols:
        kwargs["charge_from_molecules"] = charge_mols
    ic_single = Interchange.from_smirnoff(**kwargs)

    # 2. to_top → 3. [molecules] count 書き換え
    #    real species は実 count に、 末尾の Na/Cl 定義行は [molecules] から
    #    除外する (定義 [moleculetype] は top に残るので genion が参照可能)。
    ic_single.to_top(str(top_path))
    text = top_path.read_text()
    out_lines: List[str] = []
    in_mol = False
    idx = 0
    for ln in text.splitlines():
        if re.match(r"\[ molecules \]", ln):
            in_mol = True
            out_lines.append(ln)
            continue
        if in_mol and re.match(r"\[", ln):
            in_mol = False
        if in_mol and ln.strip() and not ln.strip().startswith(";"):
            if idx < n_species:
                name = ln.split()[0]
                out_lines.append(f"{name}  {counts[idx]}")
            # idx >= n_species は Na/Cl → [molecules] から除外 (genion が追加)
            idx += 1
            continue
        out_lines.append(ln)
    top_path.write_text("\n".join(out_lines) + "\n")
    logger.info(
        "OpenFF protein route: rewrote [molecules] for %d species "
        "(+ NA/CL moleculetypes defined for genion)", n_species,
    )

    # 4. mixture PDB → .gro (全コピー座標) via gmx editconf
    gmx_exec = _resolve_gmx(gmx)
    cmd = [
        gmx_exec, "editconf", "-f", str(mixture_pdb), "-o", str(gro_path),
        "-box", f"{box_size_nm:.3f}", f"{box_size_nm:.3f}", f"{box_size_nm:.3f}",
    ]
    proc = subprocess.run(cmd, capture_output=True, check=False)
    if proc.returncode != 0:
        raise RuntimeError(
            f"gmx editconf failed:\n{proc.stderr.decode(errors='replace')}"
        )

    n_atoms = 0
    try:
        with gro_path.open() as fh:
            fh.readline()
            n_atoms = int(fh.readline().strip())
    except (OSError, ValueError):
        pass
    logger.info("OpenFF protein route: GROMACS files done (%s, %s, %d atoms)",
                gro_path, top_path, n_atoms)
    return MergedTopologyResult(
        gro_path=gro_path,
        top_path=top_path,
        interchange=ic_single,
        n_atoms_total=n_atoms,
    )


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
