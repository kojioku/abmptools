# -*- coding: utf-8 -*-
"""
abmptools.genesis.grest.inp_writer
----------------------------------
Render GENESIS ``.inp`` control files for the gREST_SSCR pipeline.

Four renderers:

- :func:`render_minimize_inp`     -- atdyn minimisation (step 1)
- :func:`render_equilibrate_inp`  -- spdyn NPT equilibration (step 2)
- :func:`render_grest_inp`        -- spdyn gREST_SSCR production (step 3)
- :func:`render_remd_convert_inp` -- ``remd_convert`` parameter sort
                                     (post-MD, step 5)

Each returns a ``str`` (no I/O); the builder is responsible for writing
the result to disk. Templates derive directly from POC
``/home/okuwaki/llm-project/SI/grest-POC.md`` and GENESIS tutorial 12.3.

The :class:`SystemMetadata` dataclass carries runtime-only data
(box dimensions, residue counts) that the dataclass schema cannot know
at config-construction time -- it is filled by the builder after
:mod:`abmptools.genesis.grest.system_builder` runs tleap.
"""
from __future__ import annotations

from dataclasses import dataclass
from textwrap import dedent
from typing import List, Tuple

from .models import GrestBuildConfig
from .replica_temperatures import format_ladder_line


# ---------------------------------------------------------------------------
# Runtime metadata
# ---------------------------------------------------------------------------

@dataclass
class SystemMetadata:
    """Runtime data extracted from the prmtop/coor.

    Attributes
    ----------
    box_size_A
        ``(x, y, z)`` in Angstrom from prmtop ``%FLAG BOX_DIMENSIONS``.
    n_protein_residues
        Total non-solvent, non-ion residues. Used to render the default
        ``[SELECTION] group=rno:1-N and heavy`` for posres in the
        minimisation stage.
    rest_selection_string
        GENESIS ``rno:`` selection covering the REST solute, e.g.
        ``"rno:21,96,274-275"``.
    """
    box_size_A: Tuple[float, float, float]
    n_protein_residues: int
    rest_selection_string: str


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _on_off(b: bool) -> str:
    return "YES" if b else "NO"


def _bool_token(b: bool) -> str:
    """T/F as expected by some legacy GENESIS controls."""
    return "T" if b else "F"


# ---------------------------------------------------------------------------
# step1 -- minimize (atdyn)
# ---------------------------------------------------------------------------

def render_minimize_inp(
    cfg: GrestBuildConfig,
    meta: SystemMetadata,
    prmtop_name: str = "system.prmtop",
    coor_name: str = "system.coor",
    ref_coor_name: str = "system.coor",
) -> str:
    """Render ``step1_minimize.inp`` for atdyn.

    Posres group is ``rno:1-N and heavy`` over the protein (non-solvent)
    residues, with force constant
    :attr:`MinimizationStage.posres_force_kcal_per_mol`.
    """
    bx, by, bz = meta.box_size_A
    m = cfg.minimize
    return dedent(f"""\
        [INPUT]
        prmtopfile = {prmtop_name}                # AMBER topology
        ambcrdfile = {coor_name}                  # AMBER inpcrd
        ambreffile = {ref_coor_name}              # reference for posres

        [OUTPUT]
        dcdfile = step1.dcd
        rstfile = step1.rst

        [ENERGY]
        forcefield       = AMBER
        electrostatic    = PME
        switchdist       = {m.cutoffdist_A:.2f}     # AMBER: switch == cutoff
        cutoffdist       = {m.cutoffdist_A:.2f}
        pairlistdist     = {m.pairlistdist_A:.2f}
        vdw_force_switch = NO                       # AMBER FF
        contact_check    = YES

        [MINIMIZE]
        method        = {m.method}
        nsteps        = {m.nsteps}
        eneout_period = {m.eneout_period}
        crdout_period = {m.crdout_period}
        rstout_period = {m.rstout_period}

        [BOUNDARY]
        type       = PBC
        box_size_x = {bx:.6f}
        box_size_y = {by:.6f}
        box_size_z = {bz:.6f}

        [SELECTION]
        group1 = rno:1-{meta.n_protein_residues} and heavy

        [RESTRAINTS]
        nfunctions    = 1
        function1     = POSI
        direction1    = ALL
        constant1     = {m.posres_force_kcal_per_mol:.4f}
        select_index1 = 1
        """).rstrip() + "\n"


# ---------------------------------------------------------------------------
# step2 -- equilibrate (spdyn)
# ---------------------------------------------------------------------------

def render_equilibrate_inp(
    cfg: GrestBuildConfig,
    meta: SystemMetadata,
    prmtop_name: str = "system.prmtop",
    rst_in: str = "step1.rst",
) -> str:
    """Render ``step2_equilibrate.inp`` for spdyn (NPT).

    Reads coordinates from the minimisation restart and runs
    :attr:`EquilibrationStage.nsteps` of NPT to relax the box density
    before the gREST production phase. Position restraint on
    backbone atoms keeps the protein from drifting during heating.
    """
    e = cfg.equilibrate
    return dedent(f"""\
        [INPUT]
        prmtopfile = {prmtop_name}                # AMBER topology
        ambreffile = system.coor                  # original ref for posres
        rstfile    = {rst_in}                     # restart from step1

        [OUTPUT]
        dcdfile = step2.dcd
        rstfile = step2.rst

        [ENERGY]
        forcefield       = AMBER
        electrostatic    = PME
        switchdist       = {e.cutoffdist_A:.2f}
        cutoffdist       = {e.cutoffdist_A:.2f}
        pairlistdist     = {e.pairlistdist_A:.2f}
        vdw_force_switch = NO

        [DYNAMICS]
        integrator       = {e.integrator}
        timestep         = {e.timestep_ps:.4f}
        nsteps           = {e.nsteps}
        eneout_period    = {e.eneout_period}
        crdout_period    = {e.crdout_period}
        rstout_period    = {e.rstout_period}
        nbupdate_period  = {e.nbupdate_period}
        hydrogen_mr      = {_on_off(e.hydrogen_mr)}
        hmr_ratio        = {e.hmr_ratio:.2f}
        hmr_ratio_xh1    = {e.hmr_ratio_xh1:.2f}

        [CONSTRAINTS]
        rigid_bond = YES                          # SHAKE on H-bonds
        fast_water = YES                          # SETTLE for water
        water_model = WAT

        [ENSEMBLE]
        ensemble    = {e.ensemble}
        tpcontrol   = {e.tpcontrol}
        temperature = {e.temperature_K:.3f}
        pressure    = {e.pressure_bar:.3f}

        [BOUNDARY]
        type = PBC

        [SELECTION]
        group1 = rno:1-{meta.n_protein_residues} and (an:CA or an:N or an:C)

        [RESTRAINTS]
        nfunctions    = 1
        function1     = POSI
        direction1    = ALL
        constant1     = 1.0
        select_index1 = 1
        """).rstrip() + "\n"


# ---------------------------------------------------------------------------
# step3 -- gREST_SSCR production (spdyn)
# ---------------------------------------------------------------------------

def render_grest_inp(
    cfg: GrestBuildConfig,
    meta: SystemMetadata,
    ladder: List[float],
    prmtop_name: str = "system.prmtop",
    rst_in: str = "step2.rst",
) -> str:
    """Render ``step3_grest.inp`` for spdyn (gREST_SSCR production).

    Adds the ``[REMD]`` block with ``type1=REST`` over
    ``meta.rest_selection_string`` and a temperature ladder of length
    ``len(ladder)``. ``param_type1`` is rendered as a space-separated
    list of :attr:`RESTSelectionSpec.param_types` tokens (default
    ``"C L"`` for SSCR).
    """
    g = cfg.grest
    rest = cfg.rest_selection
    n_replicas = len(ladder)
    ladder_str = format_ladder_line(ladder)
    param_type_str = " ".join(rest.param_types)

    return dedent(f"""\
        [INPUT]
        prmtopfile = {prmtop_name}                # AMBER topology
        rstfile    = {rst_in}                     # restart from step2

        [OUTPUT]
        dcdfile = step3_rep{{}}.dcd               # one DCD per replica
        rstfile = step3_rep{{}}.rst
        logfile = step3_rep{{}}.log
        remfile = step3_rep{{}}.rem               # exchange history
        enefile = step3_rep{{}}.ene               # per-replica energies

        [ENERGY]
        forcefield       = AMBER
        electrostatic    = PME
        switchdist       = {g.cutoffdist_A:.2f}
        cutoffdist       = {g.cutoffdist_A:.2f}
        pairlistdist     = {g.pairlistdist_A:.2f}
        vdw_force_switch = NO

        [DYNAMICS]
        integrator         = {g.integrator}
        timestep           = {g.timestep_ps:.4f}
        nsteps             = {g.nsteps}
        eneout_period      = {g.eneout_period}
        crdout_period      = {g.crdout_period}
        rstout_period      = {g.rstout_period}
        nbupdate_period    = {g.nbupdate_period}
        elec_long_period   = {g.elec_long_period}
        thermostat_period  = {g.thermostat_period}
        barostat_period    = {g.barostat_period}

        [CONSTRAINTS]
        rigid_bond  = YES
        fast_water  = YES
        water_model = WAT

        [ENSEMBLE]
        ensemble    = {g.ensemble}
        tpcontrol   = {g.tpcontrol}
        temperature = {g.temperature_K:.3f}

        [BOUNDARY]
        type = PBC

        [REMD]
        dimension       = 1
        exchange_period = {g.exchange_period}
        type1           = REST
        nreplica1       = {n_replicas}
        parameters1     = {ladder_str}
        select_index1   = 1
        param_type1     = {param_type_str}
        analysis_grest  = {_on_off(g.analysis_grest)}

        [SELECTION]
        group1 = {meta.rest_selection_string}
        """).rstrip() + "\n"


# ---------------------------------------------------------------------------
# step5 -- remd_convert (parameter sort)
# ---------------------------------------------------------------------------

def render_remd_convert_inp(
    cfg: GrestBuildConfig,
    meta: SystemMetadata,
    ladder: List[float],
    prmtop_name: str = "system.prmtop",
    coor_name: str = "system.coor",
    convert_ids: List[int] = None,
) -> str:
    """Render ``step5_remd_convert.inp`` (post-MD parameter sort).

    Converts replica-indexed trajectories (``step3_rep1.dcd`` ...
    ``step3_repN.dcd``) into parameter-indexed trajectories
    (``param1.dcd`` ... ``paramN.dcd``).

    Parameters
    ----------
    convert_ids
        Empty list (default, POC convention) emits ``convert_ids = ``
        which makes ``remd_convert`` extract only the lowest-T replica.
        Provide e.g. ``[1]`` for explicit lowest-T sort, or list all
        replicas.
    """
    g = cfg.grest
    n_replicas = len(ladder)
    if convert_ids is None:
        convert_ids = []
    convert_ids_str = (
        " ".join(str(i) for i in convert_ids) if convert_ids else ""
    )

    return dedent(f"""\
        [INPUT]
        prmtopfile = {prmtop_name}
        ambcrdfile = {coor_name}
        ambreffile = {coor_name}
        dcdfile    = step3_rep{{}}.dcd
        remfile    = step3_rep{{}}.rem
        enefile    = step3_rep{{}}.ene

        [OUTPUT]
        pdbfile = param.pdb
        trjfile = param{{}}.dcd
        enefile = param{{}}.ene

        [SELECTION]
        group1 = all                              # complete system

        [FITTING]
        fitting_method = NO

        [OPTION]
        check_only      = NO
        convert_type    = PARAMETER
        convert_ids     = {convert_ids_str}
        num_replicas    = {n_replicas}
        nsteps          = {g.nsteps}
        exchange_period = {g.exchange_period}
        crdout_period   = {g.crdout_period}
        eneout_period   = {g.eneout_period}
        trjout_format   = DCD
        trjout_type     = COOR+BOX
        trjout_atom     = 1
        pbc_correct     = NO
        """).rstrip() + "\n"


__all__ = [
    "SystemMetadata",
    "render_minimize_inp",
    "render_equilibrate_inp",
    "render_grest_inp",
    "render_remd_convert_inp",
]
