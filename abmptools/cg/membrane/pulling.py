# -*- coding: utf-8 -*-
"""
abmptools.cg.membrane.pulling
------------------------------
Reaction-coordinate generation: peptide pulled along z through the bilayer.

Generic helpers (``parse_pullx_xvg``, ``find_pbc_center_atom``,
``estimate_initial_pull_coord``, ``extract_window_frames``) are FF-agnostic
and re-exported from :mod:`abmptools.membrane.pulling`.

Only :func:`write_pulling_mdp_cg` is CG-specific (uses
:mod:`mdp_templates` for the Martini 3 NVT chassis with
``direction-periodic`` pull geometry, no Pcoupl).
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

# Re-export FF-agnostic helpers.
from abmptools.membrane.pulling import (  # noqa: F401  (re-export)
    _gmx_trjconv_dump,
    _override_mdp_field,
    _read_gro_positions,
    _read_ndx_groups,
    estimate_initial_pull_coord,
    extract_window_frames,
    find_pbc_center_atom,
    parse_pullx_xvg,
)

from . import mdp_templates
from ._subprocess import ensure_dir
from .models import MembraneCGBuildConfig

logger = logging.getLogger(__name__)


def write_pulling_mdp_cg(
    *, config: MembraneCGBuildConfig, pull_dir: Path,
    pull_init_nm: float,
    pbc_atom_g1: Optional[int] = None,
    pbc_atom_g2: Optional[int] = None,
    filename: str = "pull.mdp",
) -> Path:
    """Write ``pull.mdp`` for CG steered MD along z.

    Parameters
    ----------
    pull_init_nm
        Starting value of the pull coordinate (peptide_z - bilayer_z) at
        the start of pulling. Compute via :func:`estimate_initial_pull_coord`
        from the equilibrated NPT .gro and index file.
    pbc_atom_g1
        1-based atom index used for ``pull-group1-pbcatom``. REQUIRED when
        the Bilayer group's diameter exceeds half the shortest box vector.
        Compute with :func:`find_pbc_center_atom`.
    pbc_atom_g2
        1-based atom index for ``pull-group2-pbcatom`` (peptide group).
        Usually unnecessary because peptides are smaller than half the box.

    Returns
    -------
    Path
        Absolute path to the written ``pull.mdp``.
    """
    pd = ensure_dir(Path(pull_dir)).resolve()
    text = mdp_templates.render_pull_mdp(
        config=config,
        pull_init_nm=pull_init_nm,
        pbc_atom_g1=pbc_atom_g1,
        pbc_atom_g2=pbc_atom_g2,
    )
    out = pd / filename
    out.write_text(text)
    logger.info(
        "pull.mdp written: %s (init=%.3f nm, rate=%.5f nm/ps)",
        out, pull_init_nm, config.pulling.pull_rate_nm_per_ps,
    )
    return out
