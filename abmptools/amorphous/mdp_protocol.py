# -*- coding: utf-8 -*-
"""
abmptools.amorphous.mdp_protocol
---------------------------------
Generate 5-stage annealing MDP files for GROMACS.

Stages:
  1. Energy minimisation (steep)
  2. NVT at T_high
  3. NPT at T_high
  4. NPT simulated annealing T_high → T_low
  5. NPT equilibration at T_low
"""
from __future__ import annotations

import os
from pathlib import Path
from typing import List, Optional

from ..core.system_model import AnnealProtocol


def _common_mdp(protocol: AnnealProtocol) -> dict:
    """Shared MDP settings."""
    return {
        "dt": protocol.dt,
        "nstxout": 0,
        "nstvout": 0,
        "nstfout": 0,
        "nstxout-compressed": protocol.nstxout_compressed,
        "nstenergy": protocol.nstenergy,
        "nstlog": protocol.nstenergy,
        "pbc": "xyz",
        "coulombtype": "PME",
        "rcoulomb": 1.0,
        "vdwtype": "Cut-off",
        "rvdw": 1.0,
        "DispCorr": "EnerPres",
        "constraints": "h-bonds",
        "constraint-algorithm": "lincs",
        "lincs-order": 4,
        "lincs-iter": 1,
        "cutoff-scheme": "Verlet",
        "nstlist": 10,
    }


def _thermostat_block(protocol: AnnealProtocol, ref_t: float,
                      tc_grps: str = "System") -> dict:
    # GROMACS requires the cardinality of ref-t / tau-t to match the
    # number of tokens in tc-grps. With multiple components in the
    # system the builder emits e.g. ``tc-grps = A_methanol B_water``,
    # so the single-valued ref-t / tau-t below would trigger
    # ``Invalid T coupling input: 2 groups, 1 ref-t values``. Repeat
    # the scalars to match the group count — uniform thermostat across
    # components is what amorphous-build runs want anyway.
    n_groups = max(1, len(tc_grps.split()))
    return {
        "tcoupl": "V-rescale",
        "tc-grps": tc_grps,
        "tau-t": " ".join([str(protocol.tau_t)] * n_groups),
        "ref-t": " ".join([str(ref_t)] * n_groups),
    }


def _barostat_block(protocol: AnnealProtocol) -> dict:
    return {
        "pcoupl": "Parrinello-Rahman",
        "pcoupltype": "isotropic",
        "tau-p": protocol.tau_p,
        "ref-p": protocol.P_ref,
        "compressibility": 4.5e-5,
    }


def _gen_vel_block(ref_t: float, seed: Optional[int]) -> dict:
    d: dict = {
        "gen-vel": "yes",
        "gen-temp": ref_t,
    }
    if seed is not None:
        d["gen-seed"] = seed
    return d


def _format_mdp(params: dict, title: str = "") -> str:
    """Format a dict into MDP text."""
    lines = []
    if title:
        lines.append(f"title = {title}")
    for k, v in params.items():
        if isinstance(v, bool):
            v = "yes" if v else "no"
        lines.append(f"{k:<30s} = {v}")
    return "\n".join(lines) + "\n"


# ---- Public API ----

def generate_em_mdp(protocol: AnnealProtocol) -> str:
    """Stage 1: Steepest-descent energy minimisation."""
    params = {
        "integrator": "steep",
        "emtol": protocol.em_tol,
        "emstep": 0.01,
        "nsteps": protocol.em_steps,
        "nstxout-compressed": protocol.nstxout_compressed,
        "nstenergy": protocol.nstenergy,
        "nstlog": protocol.nstenergy,
        "pbc": "xyz",
        "coulombtype": "PME",
        "rcoulomb": 1.0,
        "vdwtype": "Cut-off",
        "rvdw": 1.0,
        "DispCorr": "EnerPres",
        "cutoff-scheme": "Verlet",
        "nstlist": 10,
    }
    return _format_mdp(params, title="Energy Minimisation")


def generate_nvt_high_mdp(protocol: AnnealProtocol,
                          tc_grps: str = "System") -> str:
    """Stage 2: NVT at T_high."""
    params = {"integrator": "md", "nsteps": protocol.nvt_high_nsteps}
    params.update(_common_mdp(protocol))
    params.update(_thermostat_block(protocol, protocol.T_high, tc_grps))
    params["pcoupl"] = "no"
    params.update(_gen_vel_block(protocol.T_high, protocol.gen_seed))
    return _format_mdp(params, title="NVT high-T equilibration")


def generate_npt_high_mdp(protocol: AnnealProtocol,
                          tc_grps: str = "System") -> str:
    """Stage 3: NPT at T_high."""
    params = {"integrator": "md", "nsteps": protocol.npt_high_nsteps}
    params.update(_common_mdp(protocol))
    params.update(_thermostat_block(protocol, protocol.T_high, tc_grps))
    params.update(_barostat_block(protocol))
    params["continuation"] = "yes"
    params["gen-vel"] = "no"
    return _format_mdp(params, title="NPT high-T equilibration")


def generate_anneal_mdp(protocol: AnnealProtocol,
                        tc_grps: str = "System") -> str:
    """Stage 4: NPT simulated annealing T_high → T_low."""
    anneal_time_ps = protocol.anneal_nsteps * protocol.dt
    params = {"integrator": "md", "nsteps": protocol.anneal_nsteps}
    params.update(_common_mdp(protocol))
    params.update(_thermostat_block(protocol, protocol.T_high, tc_grps))
    params.update(_barostat_block(protocol))
    params["continuation"] = "yes"
    params["gen-vel"] = "no"
    # annealing schedule
    params["annealing"] = "single"
    params["annealing-npoints"] = 2
    params["annealing-time"] = f"0 {anneal_time_ps:.1f}"
    params["annealing-temp"] = f"{protocol.T_high} {protocol.T_low}"
    return _format_mdp(params, title="Simulated annealing")


def generate_npt_final_mdp(protocol: AnnealProtocol,
                           tc_grps: str = "System") -> str:
    """Stage 5: NPT equilibration at T_low."""
    params = {"integrator": "md", "nsteps": protocol.npt_low_nsteps}
    params.update(_common_mdp(protocol))
    params.update(_thermostat_block(protocol, protocol.T_low, tc_grps))
    params.update(_barostat_block(protocol))
    params["continuation"] = "yes"
    params["gen-vel"] = "no"
    return _format_mdp(params, title="NPT final equilibration")


def write_all_mdp(protocol: AnnealProtocol,
                  output_dir: str,
                  tc_grps: str = "System") -> List[str]:
    """Write all 5 MDP files and return their paths.

    Parameters
    ----------
    protocol : AnnealProtocol
        Annealing protocol parameters.
    output_dir : str
        Directory to write MDP files into.
    tc_grps : str
        Temperature coupling groups string for GROMACS.

    Returns
    -------
    list of str
        Paths to the 5 generated MDP files.
    """
    os.makedirs(output_dir, exist_ok=True)
    generators = [
        ("01_em.mdp", generate_em_mdp),
        ("02_nvt_highT.mdp", generate_nvt_high_mdp),
        ("03_npt_highT.mdp", generate_npt_high_mdp),
        ("04_anneal.mdp", generate_anneal_mdp),
        ("05_npt_final.mdp", generate_npt_final_mdp),
    ]
    paths = []
    for fname, gen_func in generators:
        path = os.path.join(output_dir, fname)
        if gen_func is generate_em_mdp:
            text = gen_func(protocol)
        else:
            text = gen_func(protocol, tc_grps=tc_grps)
        Path(path).write_text(text)
        paths.append(path)
    return paths


def write_run_script(output_dir: str, gro: str = "system.gro",
                     top: str = "system.top",
                     ndx: Optional[str] = "system.ndx") -> str:
    """Write a GROMACS run_all.sh script.

    Parameters
    ----------
    output_dir : str
        Directory to write the script into.
    gro, top, ndx : str
        Relative paths to the build outputs (from md/ perspective).

    Returns
    -------
    str
        Path to the generated script.
    """
    build = "../build"
    ndx_flag = f" -n {build}/{ndx}" if ndx else ""
    # ``-ntmpi 1`` forces a single MPI rank → no domain decomposition.
    # Without this, the cooling phase (anneal / npt_final) can shrink
    # the box below ``rlist × DD-cells`` and crash with
    # ``box size in direction X is too small for a cut-off`` once the
    # liquid re-condenses. Single rank still gets full multi-thread
    # speed via OpenMP, so the slowdown is negligible for the small
    # systems amorphous build typically targets.
    lines = [
        "#!/bin/bash",
        "set -e",
        "",
        f'BUILD="{build}"',
        f'GRO="$BUILD/{gro}"',
        f'TOP="$BUILD/{top}"',
        'MDRUN_OPTS="${MDRUN_OPTS:--ntmpi 1}"',
        "",
        "# Stage 1: Energy minimisation",
        f'gmx grompp -f 01_em.mdp -c "$GRO" -p "$TOP"{ndx_flag} -o 01_em.tpr -maxwarn 2',
        'gmx mdrun $MDRUN_OPTS -deffnm 01_em',
        "",
        "# Stage 2: NVT high-T",
        f'gmx grompp -f 02_nvt_highT.mdp -c 01_em.gro -p "$TOP"{ndx_flag} -o 02_nvt_highT.tpr -maxwarn 2',
        'gmx mdrun $MDRUN_OPTS -deffnm 02_nvt_highT',
        "",
        "# Stage 3: NPT high-T",
        f'gmx grompp -f 03_npt_highT.mdp -c 02_nvt_highT.gro -p "$TOP"{ndx_flag} -o 03_npt_highT.tpr -maxwarn 2',
        'gmx mdrun $MDRUN_OPTS -deffnm 03_npt_highT',
        "",
        "# Stage 4: Simulated annealing",
        f'gmx grompp -f 04_anneal.mdp -c 03_npt_highT.gro -p "$TOP"{ndx_flag} -o 04_anneal.tpr -maxwarn 2',
        'gmx mdrun $MDRUN_OPTS -deffnm 04_anneal',
        "",
        "# Stage 5: NPT final equilibration",
        f'gmx grompp -f 05_npt_final.mdp -c 04_anneal.gro -p "$TOP"{ndx_flag} -o 05_npt_final.tpr -maxwarn 2',
        'gmx mdrun $MDRUN_OPTS -deffnm 05_npt_final',
        "",
        'echo "All stages completed successfully."',
    ]
    script_path = os.path.join(output_dir, "run_all.sh")
    Path(script_path).write_text("\n".join(lines) + "\n")
    os.chmod(script_path, 0o755)
    return script_path


_DEFAULT_OPENFF_STAGES = (
    "02_nvt_highT", "03_npt_highT", "04_anneal", "05_npt_final",
)
_DEFAULT_OPENFF_FINAL_STAGE = "05_npt_final"


def write_wrap_script(output_dir: str,
                      ndx: Optional[str] = "system.ndx",
                      stages: Optional[List[str]] = None,
                      final_stage: Optional[str] = None) -> str:
    """Write a PBC-wrap script (wrap_pbc.sh) for VMD-friendly trajectories.

    Applies ``gmx trjconv -pbc mol -ur compact`` to every produced .xtc and
    to the chosen final stage's .gro, writing ``<stage>_pbc.xtc`` /
    ``<final>_pbc.gro``.

    Parameters
    ----------
    output_dir : str
        Directory to write the script into (typically the md/ directory).
    ndx : str or None
        Relative path (from md/ perspective) to the system index file.
        If ``None``, no ``-n`` flag is used; group 0 (System) is selected.
    stages : list of str, optional
        Stage basenames (each yielding ``<stage>.tpr/.xtc``) to wrap.
        Defaults to the OpenFF amorphous protocol stages
        ``02_nvt_highT, 03_npt_highT, 04_anneal, 05_npt_final``. Pass the
        4-stage hybrid-route stages
        (``01_nvt_eq, 02_npt_high, 03_npt_low, 04_nvt_sampling``) when
        wrapping that protocol.
    final_stage : str, optional
        Stage whose ``.gro`` is also wrapped (for use as the VMD initial
        frame). Defaults to ``05_npt_final`` for the OpenFF protocol or
        the last entry of ``stages`` when ``stages`` is provided.

    Returns
    -------
    str
        Path to the generated script.
    """
    if stages is None:
        stage_list = list(_DEFAULT_OPENFF_STAGES)
    else:
        stage_list = list(stages)
    if final_stage is None:
        final_stage = (
            _DEFAULT_OPENFF_FINAL_STAGE if stages is None else stage_list[-1]
        )
    build = "../build"
    ndx_flag = f' -n "{build}/{ndx}"' if ndx else ""
    stage_array = " ".join(stage_list)
    lines = [
        "#!/bin/bash",
        "# Post-processing: wrap trajectories for visualization (e.g. VMD).",
        "#",
        "#   -pbc mol    : wrap by molecules (fix molecules broken at box edges)",
        "#   -ur compact : compact unit-cell representation",
        "#",
        "# Run this after run_all.sh finishes.",
        "set -e",
        "",
        f'STAGES=({stage_array})',
        "",
        'for stage in "${STAGES[@]}"; do',
        '    [ -f "${stage}.xtc" ] || continue',
        '    echo "Wrapping ${stage}.xtc ..."',
        f'    echo 0 | gmx trjconv -s "${{stage}}.tpr" -f "${{stage}}.xtc" -o "${{stage}}_pbc.xtc" -pbc mol -ur compact{ndx_flag}',
        "done",
        "",
        "# Wrap the final structure (useful as the VMD initial frame).",
        f'if [ -f "{final_stage}.gro" ]; then',
        f'    echo 0 | gmx trjconv -s {final_stage}.tpr -f {final_stage}.gro -o {final_stage}_pbc.gro -pbc mol -ur compact{ndx_flag}',
        "fi",
        "",
        'echo ""',
        'echo "PBC wrap complete. Example VMD usage:"',
        f'echo "  vmd {final_stage}_pbc.gro -xtc {final_stage}_pbc.xtc"',
    ]
    script_path = os.path.join(output_dir, "wrap_pbc.sh")
    Path(script_path).write_text("\n".join(lines) + "\n")
    os.chmod(script_path, 0o755)
    return script_path
