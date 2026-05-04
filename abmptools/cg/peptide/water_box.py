# -*- coding: utf-8 -*-
"""
abmptools.cg.peptide.water_box
-------------------------------
Martini 3 W water box auto-generation helper.

cgmartini.nl does not distribute a Martini 3 ``water.gro`` directly
(martini_v300.zip ships only ITP files; the lipid tutorial archive is
the usual source). This module provides a self-contained generator that
seeds a single W bead and uses ``gmx insert-molecules`` to fill a small
cubic box at the standard W density. The result is suitable as a
``gmx solvate -cs`` template for the cg.peptide build pipeline.

Density rationale:
    Martini W bead = 4 H2O molecules (mass ~72 g/mol).
    Liquid water density 1 g/cm^3 -> ~8.36 W beads / nm^3.
"""
from __future__ import annotations

import logging
from pathlib import Path

from ._subprocess import ensure_dir, run_command

logger = logging.getLogger(__name__)


#: Standard Martini W number density (beads per nm^3) for liquid water.
MARTINI_W_DENSITY = 8.36


def make_martini_water_box(
    output_path: Path,
    *,
    gmx_path: str = "gmx",
    box_nm: float = 5.0,
    n_w_per_nm3: float = MARTINI_W_DENSITY,
) -> Path:
    """Generate a small Martini W cubic water box.

    Writes a 1-bead seed ``.gro``, then calls ``gmx insert-molecules`` to
    fill *box_nm*^3 cubic at the target W density.

    Parameters
    ----------
    output_path
        Destination ``.gro`` path. Parent dirs are created as needed.
    gmx_path
        gmx CLI command (default ``"gmx"``).
    box_nm
        Cubic box edge in nm. 5 nm gives ~1000 W beads which is plenty
        as a ``gmx solvate -cs`` template.
    n_w_per_nm3
        W bead number density. Default ~8.36 (Martini W standard).

    Returns
    -------
    Path
        The written water box ``.gro``.
    """
    output_path = Path(output_path).resolve()
    ensure_dir(output_path.parent)

    if box_nm <= 0:
        raise ValueError(f"box_nm must be > 0, got {box_nm}")
    if n_w_per_nm3 <= 0:
        raise ValueError(f"n_w_per_nm3 must be > 0, got {n_w_per_nm3}")

    n_w = max(1, int(round(box_nm ** 3 * n_w_per_nm3)))
    seed_gro = output_path.parent / f"_seed_W_{output_path.stem}.gro"

    # Single-W seed gro (W bead at origin in a tiny box).
    seed_gro.write_text(
        "Martini W single-bead seed\n"
        "    1\n"
        "    1W       W    1   0.000   0.000   0.000\n"
        f"   {box_nm:.5f}   {box_nm:.5f}   {box_nm:.5f}\n"
    )

    cmd = [
        gmx_path, "insert-molecules",
        "-ci", str(seed_gro),
        "-nmol", str(n_w),
        "-box", f"{box_nm:.3f}", f"{box_nm:.3f}", f"{box_nm:.3f}",
        "-o", str(output_path),
    ]
    logger.info(
        "Generating Martini W box: %d beads in %.1f nm cubic -> %s",
        n_w, box_nm, output_path,
    )
    try:
        run_command(cmd, cwd=output_path.parent)
    finally:
        seed_gro.unlink(missing_ok=True)

    if not output_path.exists():
        raise RuntimeError(
            f"gmx insert-molecules failed to produce {output_path}"
        )
    return output_path
