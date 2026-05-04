# -*- coding: utf-8 -*-
"""
abmptools.cg.membrane.insane_runner
------------------------------------
Wrapper around the `insane <https://github.com/Tsjerk/Insane>`_ CLI for
assembling Martini lipid bilayers (with optional embedded peptide).

Insane is distributed under GPL-2.0. This wrapper only invokes the
``insane`` command as an external subprocess -- insane source code is
**not bundled or modified**. abmptools itself remains MIT-licensed
(subprocess invocation is "mere aggregation" per GPL FAQ, identical to
the way :mod:`abmptools.cg.peptide.martinize_runner` handles vermouth).
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Iterable, List, Optional, Tuple

from ._subprocess import ensure_dir, run_command

logger = logging.getLogger(__name__)


def run_insane(
    output_dir: Path,
    *,
    lipid_resname: str,
    n_per_leaflet: int,
    box_d_nm: float,
    box_z_nm: float,
    peptide_pdb: Optional[Path] = None,
    peptide_z_offset_nm: float = 3.0,
    salt_concentration_M: float = 0.15,
    solvent: str = "W",
    pbc: str = "hexagonal",
    extra_args: Optional[Iterable[str]] = None,
    insane_path: str = "insane",
    output_gro_name: str = "bilayer.gro",
    output_top_name: str = "insane_topol.top",
) -> Tuple[Path, Path]:
    """Run insane to build a (peptide + ) Martini bilayer.

    Parameters
    ----------
    output_dir
        Directory to write ``output_gro_name`` and ``output_top_name`` into.
    lipid_resname
        Martini 3 lipid residue name (e.g. ``"POPC"``).
    n_per_leaflet
        Lipids per leaflet. Insane will place upper/lower symmetrically;
        some may be displaced when a peptide is inserted.
    box_d_nm
        ``-d`` argument: minimum image separation in xy. Insane sizes the
        xy bilayer patch so that the periodic image-to-image distance is
        ``box_d_nm`` after lipid placement.
    box_z_nm
        ``-dz`` argument: explicit box z dimension (nm).
    peptide_pdb
        If given, ``-f peptide_pdb`` is added so insane embeds the peptide
        in the bilayer. Coordinates from this PDB are placed at the bilayer
        centre and shifted by ``peptide_z_offset_nm`` along z.
    peptide_z_offset_nm
        ``-dm`` argument: peptide centre-of-mass z offset relative to
        bilayer centre.
    salt_concentration_M
        ``-salt`` argument (mol/L NaCl, charge-balanced when peptide carries
        a net charge).
    solvent
        ``-sol`` argument (default ``"W"`` = Martini 3 standard water bead).
    pbc
        ``-pbc`` argument (``"hexagonal"`` is insane's default).
    extra_args
        Additional argv tokens passed verbatim. Useful for power-user
        knobs (e.g. ``["-l", "DOPE", "-l", "POPC"]`` for mixtures, but
        v1 of the wrapper enforces single species at the dataclass layer).
    insane_path
        ``insane`` CLI command.
    output_gro_name / output_top_name
        Filenames inside ``output_dir``.

    Returns
    -------
    (system_gro, insane_topol)
        Absolute paths to insane's outputs. Topology is the *raw* insane
        output (``#include "martini.itp"``); post-processing happens in
        :mod:`topology_composer`.
    """
    output_dir = ensure_dir(Path(output_dir)).resolve()
    system_gro = output_dir / output_gro_name
    insane_top = output_dir / output_top_name

    # Build argv.
    cmd: List[str] = [
        insane_path,
        "-l", lipid_resname,
        "-u", lipid_resname,        # ensure both leaflets share the species
        "-d", f"{box_d_nm:.4f}",
        "-dz", f"{box_z_nm:.4f}",
        "-sol", solvent,
        "-salt", f"{salt_concentration_M:.4f}",
        "-pbc", pbc,
        "-o", str(system_gro),
        "-p", str(insane_top),
    ]

    if peptide_pdb is not None:
        peptide_pdb = Path(peptide_pdb).resolve()
        if not peptide_pdb.exists():
            raise FileNotFoundError(
                f"peptide_pdb does not exist: {peptide_pdb}"
            )
        cmd.extend(["-f", str(peptide_pdb)])
        cmd.extend(["-dm", f"{peptide_z_offset_nm:.4f}"])
        # ``-charge auto`` lets insane balance the salt around the peptide.
        cmd.extend(["-charge", "auto"])

    if extra_args:
        cmd.extend(list(extra_args))

    logger.info(
        "Running insane: lipid=%s n_per_leaflet=%d d=%.2f dz=%.2f peptide=%s",
        lipid_resname, n_per_leaflet, box_d_nm, box_z_nm,
        peptide_pdb.name if peptide_pdb else "(none)",
    )
    run_command(cmd, cwd=output_dir)

    if not system_gro.exists():
        raise RuntimeError(
            f"insane did not produce {system_gro}; check stderr above."
        )
    if not insane_top.exists():
        raise RuntimeError(
            f"insane did not produce {insane_top}; check stderr above."
        )

    logger.info(
        "insane complete: gro=%s top=%s", system_gro, insane_top,
    )
    return system_gro, insane_top
