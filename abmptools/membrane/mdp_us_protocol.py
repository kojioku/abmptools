# -*- coding: utf-8 -*-
"""
abmptools.membrane.mdp_us_protocol
----------------------------------
MDP file writers for the membrane-US pipeline.

Stages emitted
--------------
- ``em.mdp``      — energy minimisation (steep)
- ``nvt.mdp``     — NVT heating (V-rescale, no Pcoupl)
- ``npt.mdp``     — NPT relaxation, **semiisotropic** Pcoupl (Berendsen
                   during eq; switch to c-rescale for production)
- (the pull and window MDPs reuse :func:`render_pull_block` and the
  npt template; written by :mod:`pulling` and :mod:`umbrella`.)

Thermostat group convention
---------------------------
Two groups:
- ``Bilayer``     — all Lipid21 atoms
- ``Non_Bilayer`` — peptide + water + ions

This decouples the slower lipid relaxation timescale from the
faster solvent timescale and avoids small-group pathologies that
would arise with a separate Ions group (only ~10–20 atoms).

Pressure coupling convention
----------------------------
Semiisotropic — the bilayer plane (xy) couples independently from
the membrane normal (z). Both default to 1.0 bar; compressibility
4.5e-5 bar⁻¹ on both axes (≈ water value, sufficient for lipids).

License note
------------
All MDP options used here are part of GROMACS' free LGPL distribution.
No third-party MDP fragments are shipped or required.
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, Optional

from .models import MembraneConfig, USProtocol

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Public entry points
# ---------------------------------------------------------------------------

def write_equilibration_mdps(
    *, config: MembraneConfig, equil_dir: str,
) -> Dict[str, str]:
    """Write em / nvt / npt MDPs for pre-pulling equilibration.

    Returns
    -------
    dict
        Keys ``"em"``, ``"nvt"``, ``"npt"`` mapped to absolute file paths.
    """
    ed = Path(equil_dir).resolve()
    ed.mkdir(parents=True, exist_ok=True)

    paths = {
        "em":  str(ed / "em.mdp"),
        "nvt": str(ed / "nvt.mdp"),
        "npt": str(ed / "npt.mdp"),
    }
    Path(paths["em"]).write_text(render_em_mdp(config=config))
    Path(paths["nvt"]).write_text(render_nvt_mdp(config=config))
    Path(paths["npt"]).write_text(render_npt_mdp(config=config))
    logger.info("equilibration MDPs written: %s", list(paths.values()))
    return paths


# ---------------------------------------------------------------------------
# Renderers
# ---------------------------------------------------------------------------

def render_em_mdp(*, config: MembraneConfig) -> str:
    """Steepest-descent energy minimisation."""
    eq = config.equilibration
    return _strip("""
    ; em.mdp — steepest descent
    integrator      = steep
    emtol           = {emtol:.1f}
    emstep          = 0.01
    nsteps          = {emsteps:d}

    nstlist         = 10
    cutoff-scheme   = Verlet
    rlist           = 1.2
    coulombtype     = PME
    rcoulomb        = 1.2
    vdwtype         = cut-off
    vdw-modifier    = potential-shift
    rvdw            = 1.2

    constraints     = none
    pbc             = xyz
    """).format(
        emtol=eq.em_tol,
        emsteps=eq.em_steps,
    )


def render_nvt_mdp(*, config: MembraneConfig) -> str:
    """NVT heating (V-rescale, no Pcoupl)."""
    eq = config.equilibration
    seed = eq.gen_seed if eq.gen_seed is not None else -1
    return _strip("""
    ; nvt.mdp — V-rescale heating
    integrator      = md
    dt              = {dt:.3f}
    nsteps          = {nsteps:d}

    nstxout-compressed = {nstxtc:d}
    nstenergy          = {nstene:d}
    nstlog             = {nstene:d}

    nstlist         = 20
    cutoff-scheme   = Verlet
    rlist           = 1.2
    coulombtype     = PME
    rcoulomb        = 1.2
    vdwtype         = cut-off
    vdw-modifier    = potential-shift
    rvdw            = 1.2
    DispCorr        = EnerPres

    constraints       = h-bonds
    constraint-algorithm = lincs
    lincs-iter        = 1
    lincs-order       = 4

    tcoupl          = V-rescale
    tc-grps         = Bilayer Non_Bilayer
    tau-t           = {tau_t:.2f} {tau_t:.2f}
    ref-t           = {ref_t:.2f} {ref_t:.2f}

    pcoupl          = no

    gen-vel         = yes
    gen-temp        = {ref_t:.2f}
    gen-seed        = {seed:d}

    pbc             = xyz
    """).format(
        dt=eq.dt_ps,
        nsteps=eq.nvt_nsteps,
        nstxtc=eq.nstxout_compressed,
        nstene=eq.nstenergy,
        tau_t=eq.tau_t_ps,
        ref_t=eq.temperature_K,
        seed=seed,
    )


def render_npt_mdp(*, config: MembraneConfig,
                   pcoupl: str = "Berendsen") -> str:
    """NPT relaxation with semiisotropic Pcoupl.

    Default *pcoupl* is Berendsen for initial relaxation. Use
    ``"c-rescale"`` for production-like NPT (recommended for
    GROMACS 2021+) or ``"Parrinello-Rahman"`` for fluctuation-correct
    sampling once well-equilibrated.
    """
    eq = config.equilibration
    return _strip("""
    ; npt.mdp — semiisotropic NPT
    integrator      = md
    dt              = {dt:.3f}
    nsteps          = {nsteps:d}

    nstxout-compressed = {nstxtc:d}
    nstenergy          = {nstene:d}
    nstlog             = {nstene:d}

    nstlist         = 20
    cutoff-scheme   = Verlet
    rlist           = 1.2
    coulombtype     = PME
    rcoulomb        = 1.2
    vdwtype         = cut-off
    vdw-modifier    = potential-shift
    rvdw            = 1.2
    DispCorr        = EnerPres

    constraints       = h-bonds
    constraint-algorithm = lincs
    lincs-iter        = 1
    lincs-order       = 4

    tcoupl          = V-rescale
    tc-grps         = Bilayer Non_Bilayer
    tau-t           = {tau_t:.2f} {tau_t:.2f}
    ref-t           = {ref_t:.2f} {ref_t:.2f}

    pcoupl          = {pcoupl:s}
    pcoupltype      = semiisotropic
    tau-p           = {tau_p:.2f}
    ref-p           = {ref_p:.2f} {ref_p:.2f}
    compressibility = 4.5e-5 4.5e-5
    refcoord-scaling = com

    gen-vel         = no

    pbc             = xyz
    """).format(
        dt=eq.dt_ps,
        nsteps=eq.npt_nsteps,
        nstxtc=eq.nstxout_compressed,
        nstene=eq.nstenergy,
        tau_t=eq.tau_t_ps,
        ref_t=eq.temperature_K,
        tau_p=eq.tau_p_ps,
        ref_p=eq.pressure_bar,
        pcoupl=pcoupl,
    )


def render_pull_block(
    *, config: MembraneConfig,
    pull_init_nm: float,
    pull_rate_nm_per_ps: float,
    pull_k: Optional[float] = None,
    nst_pull_xf: int = 50,
) -> str:
    """Render the [pull] / [pull-coord*] block for membrane US.

    Parameters
    ----------
    pull_init_nm : float
        Reference COM-COM z distance (peptide minus bilayer).
        For a static umbrella window, this is the window centre.
        For pulling, this is the **starting** value (since
        ``pull-coord1-start = no``).
    pull_rate_nm_per_ps : float
        Constant velocity (nm/ps). Use 0.0 for static umbrella windows.
    pull_k : float or None
        Force constant (kJ/mol/nm²). If None, defaults to
        ``config.umbrella.force_constant_kj_mol_nm2``.
    nst_pull_xf : int
        Output stride for ``pullx.xvg`` / ``pullf.xvg``.
    """
    if pull_k is None:
        pull_k = config.umbrella.force_constant_kj_mol_nm2

    geom = config.umbrella.pull_geometry
    vec_line = (
        f"pull-coord1-vec            = {config.umbrella.pull_vec}\n"
        if geom.startswith("direction") else ""
    )
    return _strip("""
    ; pull-code (signed z-distance via {geom:s})
    pull                       = yes
    pull-ngroups               = 2
    pull-ncoords               = 1
    pull-group1-name           = {g1:s}
    pull-group2-name           = {g2:s}
    pull-coord1-type           = umbrella
    pull-coord1-geometry       = {geom:s}
    pull-coord1-groups         = 1 2
    pull-coord1-dim            = {dim:s}
    """).format(
        g1=config.umbrella.pull_group1,
        g2=config.umbrella.pull_group2,
        geom=geom,
        dim=config.umbrella.pull_dim,
    ) + vec_line + _strip("""
    pull-coord1-init           = {init:.4f}
    pull-coord1-rate           = {rate:.6f}
    pull-coord1-k              = {k:.2f}
    pull-coord1-start          = no
    pull-pbc-ref-prev-step-com = yes
    pull-print-com             = yes
    pull-nstxout               = {nst:d}
    pull-nstfout               = {nst:d}
    """).format(
        init=pull_init_nm,
        rate=pull_rate_nm_per_ps,
        k=pull_k,
        nst=nst_pull_xf,
    )


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _strip(template: str) -> str:
    """Strip the leading 4-space indent from a multi-line MDP template.

    Keeps the source code visually tidy while emitting GROMACS-friendly
    flush-left lines. Blank lines are preserved.
    """
    lines = template.splitlines()
    # drop leading blank line
    while lines and not lines[0].strip():
        lines.pop(0)
    out = []
    for line in lines:
        if line.startswith("    "):
            out.append(line[4:])
        else:
            out.append(line)
    return "\n".join(out) + "\n"
