# -*- coding: utf-8 -*-
"""
abmptools.cg.membrane.mdp_templates
------------------------------------
Martini 3 MDP renderers for peptide-membrane builds.

The CG flavour differs from the AA membrane (CHARMM36 / AMBER Lipid21)
in three places:

1. **Non-bonded**: Martini 3 uses ``coulombtype = reaction-field`` with
   ``rcoulomb = 1.1 nm`` and ``epsilon_r = 15``. AA uses PME with cutoff
   1.2 nm.
2. **Constraints**: Martini 3 has no h-bond constraints (CG beads have
   no explicit hydrogens). AA uses ``constraints = h-bonds`` (LINCS).
3. **Pressure coupling**: ``compressibility = 3e-4`` (Martini standard,
   ~6× water value because CG dampens density fluctuations) and
   ``tau_p = 12.0 ps`` (Martini standard). AA uses 4.5e-5 / 5.0 ps.

The 2-group thermostat (``Bilayer`` / ``Non_Bilayer``) and semiisotropic
P-coupling are CG/AA-common.

For the pull-code block we delegate to
:func:`abmptools.membrane.mdp_us_protocol.render_pull_block`, which only
reads ``config.umbrella.{pull_*, force_constant_kj_mol_nm2}`` (field
names match between :class:`UmbrellaCGProtocol` and ``USProtocol`` for
the relevant fields, so duck-typing works).
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, Optional

from abmptools.membrane.mdp_us_protocol import _strip
from abmptools.membrane.mdp_us_protocol import (
    render_pull_block as _aa_render_pull_block,
)

from ._subprocess import ensure_dir, write_text
from .models import MembraneCGBuildConfig

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Building blocks (Martini 3 specific)
# ---------------------------------------------------------------------------

# Non-bonded: M3 standard. Verlet + reaction-field (epsilon_r=15) + cutoff vdw
# at 1.1 nm with potential-shift modifier.
_CG_NB_BLOCK = _strip("""
    nstlist                  = 20
    cutoff-scheme            = Verlet
    pbc                      = xyz
    verlet-buffer-tolerance  = 0.005
    coulombtype              = reaction-field
    rcoulomb                 = 1.1
    epsilon_r                = 15
    epsilon_rf               = 0
    vdwtype                  = cutoff
    rvdw                     = 1.1
    vdw-modifier             = Potential-shift-Verlet
""")


def _tc_2group_block(ref_t_K: float, tau_t_ps: float) -> str:
    """V-rescale thermostat with two groups (Bilayer / Non_Bilayer).

    Decouples slower lipid relaxation from solvent-side timescales and
    avoids small-group pathologies that a separate ``Ions`` group would
    introduce (only ~10–20 ions in a typical box).
    """
    return _strip("""
    tcoupl                   = v-rescale
    tc-grps                  = Bilayer Non_Bilayer
    tau_t                    = {tau:.2f} {tau:.2f}
    ref_t                    = {ref:.1f} {ref:.1f}
    """).format(tau=tau_t_ps, ref=ref_t_K)


def _pc_semiiso_block(ref_p_bar: float, tau_p_ps: float) -> str:
    """Semiisotropic c-rescale P-coupling for Martini 3 bilayer."""
    return _strip("""
    pcoupl                   = c-rescale
    pcoupltype               = semiisotropic
    tau_p                    = {tau:.2f}
    ref_p                    = {ref:.2f} {ref:.2f}
    compressibility          = 3e-4 3e-4
    refcoord-scaling         = com
    """).format(tau=tau_p_ps, ref=ref_p_bar)


# ---------------------------------------------------------------------------
# Public renderers: equilibration (em / nvt / npt)
# ---------------------------------------------------------------------------

def render_em_mdp(*, config: MembraneCGBuildConfig) -> str:
    """Energy minimisation (steepest descent)."""
    eq = config.equilibration
    return _strip("""
    ; em.mdp -- Martini 3 steepest-descent energy minimisation
    integrator               = steep
    emtol                    = {emtol:.1f}
    emstep                   = 0.01
    nsteps                   = {emsteps:d}

    """).format(
        emtol=eq.em_tol,
        emsteps=eq.em_steps,
    ) + _CG_NB_BLOCK + "\n" + _strip("""
    constraints              = none
    """)


def render_nvt_mdp(*, config: MembraneCGBuildConfig) -> str:
    """NVT heating (V-rescale 2-group, no Pcoupl)."""
    eq = config.equilibration
    seed = eq.gen_seed if eq.gen_seed is not None else -1
    return _strip("""
    ; nvt.mdp -- Martini 3 NVT, 2-group V-rescale
    integrator               = md
    dt                       = {dt:.4f}
    nsteps                   = {nsteps:d}
    nstxout-compressed       = {nstxtc:d}
    nstenergy                = {nstene:d}
    nstlog                   = {nstene:d}

    """).format(
        dt=eq.dt_ps,
        nsteps=eq.nvt_nsteps,
        nstxtc=eq.nstxout_compressed,
        nstene=eq.nstenergy,
    ) + _CG_NB_BLOCK + "\n" + _tc_2group_block(
        ref_t_K=eq.temperature_K, tau_t_ps=eq.tau_t_ps,
    ) + _strip("""
    pcoupl                   = no

    gen-vel                  = yes
    gen-temp                 = {ref:.1f}
    gen-seed                 = {seed:d}

    constraints              = none
    """).format(ref=eq.temperature_K, seed=seed)


def render_npt_mdp(*, config: MembraneCGBuildConfig) -> str:
    """NPT relaxation with semiisotropic c-rescale Pcoupl."""
    eq = config.equilibration
    return _strip("""
    ; npt.mdp -- Martini 3 NPT, semiisotropic, 2-group V-rescale
    integrator               = md
    dt                       = {dt:.4f}
    nsteps                   = {nsteps:d}
    nstxout-compressed       = {nstxtc:d}
    nstenergy                = {nstene:d}
    nstlog                   = {nstene:d}

    """).format(
        dt=eq.dt_ps,
        nsteps=eq.npt_nsteps,
        nstxtc=eq.nstxout_compressed,
        nstene=eq.nstenergy,
    ) + _CG_NB_BLOCK + "\n" + _tc_2group_block(
        ref_t_K=eq.temperature_K, tau_t_ps=eq.tau_t_ps,
    ) + _pc_semiiso_block(
        ref_p_bar=eq.pressure_bar, tau_p_ps=eq.tau_p_ps,
    ) + _strip("""

    gen-vel                  = no

    constraints              = none
    """)


# ---------------------------------------------------------------------------
# Pull MDP (constant-rate pulling along z) and per-window MDP (rate=0)
# ---------------------------------------------------------------------------

def render_pull_mdp(
    *, config: MembraneCGBuildConfig,
    pull_init_nm: float,
    pbc_atom_g1: Optional[int] = None,
    pbc_atom_g2: Optional[int] = None,
) -> str:
    """Constant-rate pulling MDP (NVT base + pull block, no Pcoupl).

    The pulling stage must run with **no pressure coupling** because
    GROMACS rejects ``direction-periodic`` pull geometry under a dynamic
    box (semi/isotropic Pcoupl). The peptide's z-distance to bilayer COM
    transiently exceeds ``0.49 × shortest-box-vector`` during traversal,
    which only ``direction-periodic`` handles cleanly. The bilayer barely
    deforms over the pulling timescale (a few ns) without Pcoupl, so this
    NVT-pull / NPT-window split is a standard membrane US convention.

    The window stages re-introduce semiisotropic Pcoupl and geometry=
    ``direction`` (see :func:`render_window_mdp`).
    """
    pull = config.pulling
    return _replace_chassis_with_pull(
        chassis_text=render_nvt_mdp(config=config),
        pull_block=_aa_render_pull_block(
            config=config,
            pull_init_nm=pull_init_nm,
            pull_rate_nm_per_ps=pull.pull_rate_nm_per_ps,
            pull_k=pull.pull_force_constant,
            nst_pull_xf=max(50, pull.nstxout_compressed // 10),
            pbc_atom_g1=pbc_atom_g1,
            pbc_atom_g2=pbc_atom_g2,
            geometry_override="direction-periodic",
        ),
        nsteps=pull.nsteps,
        nstxout_compressed=pull.nstxout_compressed,
        regenerate_velocities=False,   # NVT chassis already has gen-vel=yes;
                                        # for pulling we want continuation
    )


def render_window_mdp(
    *, config: MembraneCGBuildConfig,
    window_z_nm: float,
    pbc_atom_g1: Optional[int] = None,
    pbc_atom_g2: Optional[int] = None,
) -> str:
    """Per-window umbrella MDP (NPT base + static pull at *window_z_nm*).

    Window MDPs use ``geometry=direction`` (signed z projection),
    compatible with semiisotropic dynamic-box. The window NPT base also
    preserves the post-pull box dimensions through Pcoupl relaxation.
    """
    umb = config.umbrella
    return _replace_chassis_with_pull(
        chassis_text=render_npt_mdp(config=config),
        pull_block=_aa_render_pull_block(
            config=config,
            pull_init_nm=window_z_nm,
            pull_rate_nm_per_ps=0.0,
            pull_k=umb.force_constant_kj_mol_nm2,
            nst_pull_xf=umb.window_nstxout_compressed,
            pbc_atom_g1=pbc_atom_g1,
            pbc_atom_g2=pbc_atom_g2,
            # Windows use ``direction`` (signed z, dynamic-box compatible)
        ),
        nsteps=umb.window_nsteps,
        nstxout_compressed=umb.window_nstxout_compressed,
        regenerate_velocities=False,
    )


def _replace_chassis_with_pull(
    *, chassis_text: str, pull_block: str,
    nsteps: int, nstxout_compressed: int,
    regenerate_velocities: bool,
) -> str:
    """Splice a pull block onto a (NVT or NPT) MDP chassis.

    Patches ``nsteps`` / ``nstxout-compressed`` to the new values, and
    flips ``gen-vel = yes`` to ``no`` when *regenerate_velocities* is
    False (typically for pull/window stages that continue from the
    previous .gro/.cpt instead of regenerating velocities).
    """
    out_lines = []
    for line in chassis_text.splitlines(keepends=True):
        stripped = line.strip()
        if stripped.startswith("nsteps") and "=" in stripped:
            out_lines.append(f"nsteps                   = {nsteps:d}\n")
        elif stripped.startswith("nstxout-compressed"):
            out_lines.append(
                f"nstxout-compressed       = {nstxout_compressed:d}\n"
            )
        elif (
            not regenerate_velocities
            and stripped.startswith("gen-vel")
            and "yes" in stripped
        ):
            out_lines.append("gen-vel                  = no\n")
        else:
            out_lines.append(line)
    return "".join(out_lines) + "\n" + pull_block


# ---------------------------------------------------------------------------
# Bulk writers
# ---------------------------------------------------------------------------

def write_equilibration_mdps(
    *, config: MembraneCGBuildConfig, equil_dir: Path,
) -> Dict[str, Path]:
    """Write em / nvt / npt MDPs into *equil_dir*. Returns ``{name: path}``."""
    ed = ensure_dir(Path(equil_dir)).resolve()
    paths = {
        "em": ed / "em.mdp",
        "nvt": ed / "nvt.mdp",
        "npt": ed / "npt.mdp",
    }
    write_text(paths["em"], render_em_mdp(config=config))
    write_text(paths["nvt"], render_nvt_mdp(config=config))
    write_text(paths["npt"], render_npt_mdp(config=config))
    logger.info(
        "equilibration MDPs written: %s",
        [str(p) for p in paths.values()],
    )
    return paths
