# -*- coding: utf-8 -*-
"""
abmptools.cg.peptide.system_packer
-----------------------------------
Box construction, multi-peptide insertion, solvation, and ion addition
for Martini 3 systems.

すべて GROMACS の ``gmx`` ツール (insert-molecules / solvate / grompp /
genion / make_ndx) を subprocess で呼び出す。自前の grid 配置や自前
solvent / ion 充填は持たない (Martini 3 公式推奨フローに沿う)。
"""
from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

from ._subprocess import ensure_dir, run_command, write_text

logger = logging.getLogger(__name__)


# ``gmx grompp`` requires a valid MDP even for ions.tpr; this is the
# minimal Martini-3-compatible template.
IONS_MDP = """\
; Minimal MDP for ions.tpr (used only by gmx grompp before genion).
integrator  = steep
nsteps      = 1
emtol       = 1000.0
emstep      = 0.01
nstlist     = 10
cutoff-scheme   = Verlet
coulombtype     = reaction-field
rcoulomb        = 1.1
epsilon_r       = 15
epsilon_rf      = 0
vdwtype         = cutoff
rvdw            = 1.1
vdw-modifier    = Potential-shift-Verlet
"""


@dataclass
class PackedSystem:
    """Output paths after pack + solvate + ions."""
    packed_gro: Path
    solvated_gro: Path
    ions_gro: Path
    ions_tpr: Path
    ions_mdp: Path


def insert_peptides(
    peptide_specs: Sequence[Tuple[str, Path, int]],
    output_dir: Path,
    box_size_nm: float = 0.0,
    box_lengths_nm: Optional[List[float]] = None,
    *,
    gmx_path: str = "gmx",
    seed: Optional[int] = None,
) -> Path:
    """Place peptide CG PDBs into a box via ``gmx insert-molecules``.

    Parameters
    ----------
    peptide_specs
        Sequence of ``(name, cg_pdb_path, count)`` tuples.
    output_dir
        Where to write packed.gro and intermediate files.
    box_size_nm
        Cubic box edge in nm (used if *box_lengths_nm* is None).
    box_lengths_nm
        Optional explicit ``[x, y, z]`` lengths in nm.
    gmx_path
        gmx CLI path (default ``"gmx"``).
    seed
        Optional seed for ``gmx insert-molecules -seed``.

    Returns
    -------
    Path to ``packed.gro``.
    """
    output_dir = ensure_dir(Path(output_dir)).resolve()

    if box_lengths_nm is not None:
        box = list(box_lengths_nm)
    else:
        if box_size_nm <= 0:
            raise ValueError(
                "Either box_size_nm > 0 or box_lengths_nm must be provided"
            )
        box = [box_size_nm] * 3
    box_args = ["-box", *(f"{L:.3f}" for L in box)]

    current_gro: Optional[Path] = None
    inserted_any = False
    for i, (name, cg_pdb, count) in enumerate(peptide_specs):
        if count <= 0:
            continue
        out_gro = output_dir / f"packed_{i:02d}_{name}.gro"
        cmd = [
            gmx_path, "insert-molecules",
            "-ci", str(Path(cg_pdb).resolve()),
            "-nmol", str(count),
            "-o", str(out_gro),
        ]
        if current_gro is None:
            cmd.extend(box_args)
        else:
            cmd.extend(["-f", str(current_gro)])
        if seed is not None:
            cmd.extend(["-seed", str(seed)])

        logger.info("gmx insert-molecules: %s x %d", name, count)
        run_command(cmd, cwd=output_dir)
        current_gro = out_gro
        inserted_any = True

    if not inserted_any or current_gro is None:
        raise ValueError("No peptides to insert (all counts were zero)")

    final = output_dir / "packed.gro"
    if final.exists():
        final.unlink()
    current_gro.rename(final)
    logger.info("Packed system written: %s", final)
    return final


def solvate(
    packed_gro: Path,
    topol_top: Path,
    output_dir: Path,
    martini_water_gro: Path,
    *,
    gmx_path: str = "gmx",
) -> Path:
    """Solvate with Martini W via ``gmx solvate``.

    Parameters
    ----------
    martini_water_gro
        Path to a Martini W solvent box (typically
        ``martini_v3.0.0_water.gro`` from cgmartini.nl). Not bundled with
        this package -- see :mod:`forcefield_check`.
    """
    output_dir = ensure_dir(Path(output_dir)).resolve()
    out = output_dir / "system_solv.gro"
    cmd = [
        gmx_path, "solvate",
        "-cp", str(Path(packed_gro).resolve()),
        "-cs", str(Path(martini_water_gro).resolve()),
        "-p", str(Path(topol_top).resolve()),
        "-o", str(out),
    ]
    logger.info("gmx solvate -> %s", out)
    run_command(cmd, cwd=output_dir)
    return out


def add_ions(
    solvated_gro: Path,
    topol_top: Path,
    output_dir: Path,
    *,
    gmx_path: str = "gmx",
    neutralize: bool = True,
    conc_molar: float = 0.15,
    pname: str = "NA",
    nname: str = "CL",
    solvent_group: str = "W",
) -> PackedSystem:
    """Generate ions.tpr and run ``gmx genion`` to add NaCl ions.

    The ``solvent_group`` identifier is fed to genion via stdin; for a
    Martini 3 build with the standard ``[ molecules ]`` ``W`` water,
    ``"W"`` is correct.
    """
    output_dir = ensure_dir(Path(output_dir)).resolve()
    ions_mdp = output_dir / "ions.mdp"
    write_text(ions_mdp, IONS_MDP)

    ions_tpr = output_dir / "ions.tpr"
    grompp_cmd = [
        gmx_path, "grompp",
        "-f", str(ions_mdp),
        "-c", str(Path(solvated_gro).resolve()),
        "-p", str(Path(topol_top).resolve()),
        "-o", str(ions_tpr),
        "-maxwarn", "10",
    ]
    logger.info("gmx grompp -> %s", ions_tpr)
    run_command(grompp_cmd, cwd=output_dir)

    out_gro = output_dir / "system_ions.gro"
    genion_cmd = [
        gmx_path, "genion",
        "-s", str(ions_tpr),
        "-o", str(out_gro),
        "-p", str(Path(topol_top).resolve()),
        "-pname", pname,
        "-nname", nname,
    ]
    if neutralize:
        genion_cmd.append("-neutral")
    if conc_molar > 0:
        genion_cmd.extend(["-conc", str(conc_molar)])

    logger.info(
        "gmx genion (neutralize=%s, conc=%s M) -> %s",
        neutralize, conc_molar, out_gro,
    )
    run_command(genion_cmd, cwd=output_dir, stdin_text=f"{solvent_group}\n")

    return PackedSystem(
        packed_gro=Path(solvated_gro).parent / "packed.gro",
        solvated_gro=Path(solvated_gro),
        ions_gro=out_gro,
        ions_tpr=ions_tpr,
        ions_mdp=ions_mdp,
    )


def make_index(
    gro_path: Path,
    output_dir: Path,
    *,
    gmx_path: str = "gmx",
) -> Path:
    """Build a default ``index.ndx`` via ``gmx make_ndx`` (auto-quit)."""
    output_dir = ensure_dir(Path(output_dir)).resolve()
    out = output_dir / "index.ndx"
    cmd = [
        gmx_path, "make_ndx",
        "-f", str(Path(gro_path).resolve()),
        "-o", str(out),
    ]
    logger.info("gmx make_ndx -> %s", out)
    run_command(cmd, cwd=output_dir, stdin_text="q\n")
    return out
