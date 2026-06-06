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
    # GROMACS annealing parameters are per-group, mirroring tc-grps
    # cardinality: ``annealing`` / ``annealing-npoints`` give one
    # token per group, ``annealing-time`` / ``annealing-temp`` give
    # (group × npoints) values total. With 1 group the legacy single
    # schedule still applies (the loop emits exactly one copy).
    n_groups = max(1, len(tc_grps.split()))
    params["annealing"] = " ".join(["single"] * n_groups)
    params["annealing-npoints"] = " ".join(["2"] * n_groups)
    params["annealing-time"] = " ".join([f"0 {anneal_time_ps:.1f}"] * n_groups)
    params["annealing-temp"] = " ".join(
        [f"{protocol.T_high} {protocol.T_low}"] * n_groups
    )
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
                  tc_grps: str = "System",
                  define_posres: Optional[str] = None) -> List[str]:
    """Write all 5 MDP files and return their paths.

    Parameters
    ----------
    protocol : AnnealProtocol
        Annealing protocol parameters.
    output_dir : str
        Directory to write MDP files into.
    tc_grps : str
        Temperature coupling groups string for GROMACS.
    define_posres : str, optional
        ``#define`` flag (e.g. ``"POSRES_TRIMER"``) added to 02_nvt_highT
        and later stages via mdp ``define = -DPOSRES_TRIMER``. Used to
        activate ``[ position_restraints ]`` blocks added inside specific
        moleculetypes in system.top (see builder.AmorphousBuilder.
        _add_trimer_posres_to_top). EM (01_em.mdp) is intentionally
        skipped because posres of pre-relaxed packmol overlaps yields
        ill-conditioned forces.

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
    define_block = ""
    if define_posres:
        define_block = (
            f"\n; Activate [ position_restraints ] in system.top\n"
            f"define                   = -D{define_posres}\n"
        )
    for fname, gen_func in generators:
        path = os.path.join(output_dir, fname)
        if gen_func is generate_em_mdp:
            text = gen_func(protocol)
        else:
            text = gen_func(protocol, tc_grps=tc_grps)
        if define_block and gen_func is not generate_em_mdp:
            text = text.rstrip() + "\n" + define_block
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
        f'gmx grompp -f 01_em.mdp -c "$GRO" -r "$GRO" -p "$TOP"{ndx_flag} -o 01_em.tpr -maxwarn 2',
        'gmx mdrun $MDRUN_OPTS -deffnm 01_em',
        "",
        "# Stage 2: NVT high-T",
        f'gmx grompp -f 02_nvt_highT.mdp -c 01_em.gro -r 01_em.gro -p "$TOP"{ndx_flag} -o 02_nvt_highT.tpr -maxwarn 2',
        'gmx mdrun $MDRUN_OPTS -deffnm 02_nvt_highT',
        "",
        "# Stage 3: NPT high-T",
        f'gmx grompp -f 03_npt_highT.mdp -c 02_nvt_highT.gro -r 02_nvt_highT.gro -p "$TOP"{ndx_flag} -o 03_npt_highT.tpr -maxwarn 2',
        'gmx mdrun $MDRUN_OPTS -deffnm 03_npt_highT',
        "",
        "# Stage 4: Simulated annealing",
        f'gmx grompp -f 04_anneal.mdp -c 03_npt_highT.gro -r 03_npt_highT.gro -p "$TOP"{ndx_flag} -o 04_anneal.tpr -maxwarn 2',
        'gmx mdrun $MDRUN_OPTS -deffnm 04_anneal',
        "",
        "# Stage 5: NPT final equilibration",
        f'gmx grompp -f 05_npt_final.mdp -c 04_anneal.gro -r 04_anneal.gro -p "$TOP"{ndx_flag} -o 05_npt_final.tpr -maxwarn 2',
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
    """Write a PBC-wrap Python script (``wrap_pbc.py``) for VMD-friendly trajectories.

    Generates a cross-platform Python script that calls
    :func:`abmptools.trajectory.wrap_pbc` (``gmx trjconv -pbc mol -ur compact``)
    on every produced .xtc and on the chosen final stage's .gro, writing
    ``<stage>_pbc.xtc`` / ``<final>_pbc.gro``.

    旧 ``wrap_pbc.sh`` (bash 専用) の Windows-compatible 置換。

    Parameters
    ----------
    output_dir : str
        Directory to write the script into (typically the md/ directory).
    ndx : str or None
        Relative path (from md/ perspective) to the system index file.
        If ``None``, no index file is used; group 0 (System) is selected.
    stages : list of str, optional
        Stage basenames (each yielding ``<stage>.tpr/.xtc``) to wrap.
        Defaults to the OpenFF amorphous protocol stages
        ``02_nvt_highT, 03_npt_highT, 04_anneal, 05_npt_final``.
    final_stage : str, optional
        Stage whose ``.gro`` is also wrapped (for use as the VMD initial
        frame). Defaults to ``05_npt_final`` for the OpenFF protocol or
        the last entry of ``stages`` when ``stages`` is provided.

    Returns
    -------
    str
        Path to the generated ``wrap_pbc.py``.
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
    ndx_repr = f'"{build}/{ndx}"' if ndx else "None"
    stages_repr = ", ".join(f'"{s}"' for s in stage_list)
    lines = [
        "#!/usr/bin/env python",
        '"""Post-processing: wrap trajectories for visualization (e.g. VMD).',
        "",
        "Generated by abmptools.amorphous.mdp_protocol.write_wrap_script.",
        "Cross-platform (Linux / macOS / Windows) — runs on any system with",
        "Python + gmx + abmptools installed.",
        "",
        "    -pbc mol    : wrap by molecules (fix molecules broken at box edges)",
        "    -ur compact : compact unit-cell representation",
        "",
        "Run after run_all.sh (mdrun pipeline) finishes:",
        "    python wrap_pbc.py",
        '"""',
        "from pathlib import Path",
        "",
        "from abmptools.trajectory import wrap_pbc",
        "",
        f"STAGES = [{stages_repr}]",
        f'FINAL_STAGE = "{final_stage}"',
        f"NDX = {ndx_repr}",
        "",
        "for stage in STAGES:",
        '    xtc = Path(f"{stage}.xtc")',
        '    tpr = Path(f"{stage}.tpr")',
        "    if not xtc.is_file():",
        "        continue",
        '    print(f"Wrapping {xtc.name} ...")',
        "    wrap_pbc(",
        "        trajectory=xtc, tpr=tpr,",
        '        output=Path(f"{stage}_pbc.xtc"),',
        '        group="0",  # System (group index 0)',
        "        ndx=NDX,",
        "    )",
        "",
        "# Wrap the final structure (useful as the VMD initial frame).",
        'final_gro = Path(f"{FINAL_STAGE}.gro")',
        "if final_gro.is_file():",
        "    wrap_pbc(",
        "        trajectory=final_gro,",
        '        tpr=Path(f"{FINAL_STAGE}.tpr"),',
        '        output=Path(f"{FINAL_STAGE}_pbc.gro"),',
        '        group="0", ndx=NDX,',
        "    )",
        "",
        'print()',
        'print("PBC wrap complete. Example VMD usage:")',
        'print(f"  vmd {FINAL_STAGE}_pbc.gro -xtc {FINAL_STAGE}_pbc.xtc")',
    ]
    script_path = os.path.join(output_dir, "wrap_pbc.py")
    Path(script_path).write_text("\n".join(lines) + "\n")
    os.chmod(script_path, 0o755)
    return script_path


def write_udf_export_script(output_dir: str,
                            ndx: Optional[str] = "system.ndx",
                            stage: Optional[str] = None,
                            n_energy_terms: int = 50) -> str:
    """Write a UDF / J-OCTA export Python script (``gen_for_udf.py``).

    Generates a cross-platform Python script that dumps two OCTA / J-OCTA
    compatible inputs from the production stage of the OpenFF amorphous
    protocol via :mod:`abmptools.trajectory`:

    * ``<stage>_energy.xvg`` — :func:`abmptools.trajectory.gmx_energy` with
      ``terms=range(1, n_energy_terms+1)``。 gmx は存在しない index を
      silently skip するため、 大きめ上限で全 term を一括取得。
    * ``<stage>_nojump.gro`` — :func:`abmptools.trajectory.nojump`
      (``gmx trjconv -pbc nojump``)。 OCTA / J-OCTA Viewer + gro2udf 下流
      用に分子を連続化。

    旧 ``gen_for_udf.sh`` (bash 専用) の Windows-compatible 置換。

    Parameters
    ----------
    output_dir : str
        Directory to write the script into (typically the md/ directory).
    ndx : str or None
        Relative path (from md/ perspective) to the system index file.
        If ``None``, no index file is used; group 0 (System) is selected.
    stage : str, optional
        Stage basename whose ``.edr`` / ``.trr`` / ``.tpr`` are exported.
        Defaults to ``05_npt_final`` (the OpenFF production stage).
    n_energy_terms : int, optional
        Upper bound for the energy term range (``range(1, N+1)``). Default 50.

    Returns
    -------
    str
        Path to the generated ``gen_for_udf.py``.
    """
    if stage is None:
        stage = _DEFAULT_OPENFF_FINAL_STAGE
    build = "../build"
    ndx_repr = f'"{build}/{ndx}"' if ndx else "None"
    lines = [
        "#!/usr/bin/env python",
        '"""Post-processing: export UDF / J-OCTA compatible inputs from MD outputs.',
        "",
        "Generated by abmptools.amorphous.mdp_protocol.write_udf_export_script.",
        "Cross-platform (Linux / macOS / Windows) — runs on any system with",
        "Python + gmx + abmptools installed.",
        "",
        "    gmx energy              : dump every energy term (1..N) to <stage>_energy.xvg",
        "    gmx trjconv -pbc nojump : keep molecules continuous across PBC for OCTA",
        "                              / J-OCTA Viewer and downstream UDF conversion",
        "                              (in contrast to wrap_pbc.py, which uses -pbc mol",
        "                              for VMD-compatible compact unit-cell rendering)",
        "",
        "Run after run_all.sh (mdrun pipeline) finishes:",
        "    python gen_for_udf.py",
        '"""',
        "from pathlib import Path",
        "",
        "from abmptools.trajectory import gmx_energy, nojump",
        "",
        f'STAGE = "{stage}"',
        f"N_ENERGY_TERMS = {n_energy_terms}",
        f"NDX = {ndx_repr}",
        "",
        "# 1. All energy terms -> <stage>_energy.xvg",
        'edr = Path(f"{STAGE}.edr")',
        "if edr.is_file():",
        '    print(f"Exporting energy terms to {STAGE}_energy.xvg ...")',
        "    gmx_energy(",
        "        edr=edr,",
        '        output=Path(f"{STAGE}_energy.xvg"),',
        "        terms=range(1, N_ENERGY_TERMS + 1),",
        "    )",
        "",
        "# 2. Trajectory with -pbc nojump -> <stage>_nojump.gro",
        "#    Prefer .trr (positions + velocities) when available, fall back to .xtc.",
        "input_path = None",
        'for ext in ("trr", "xtc"):',
        '    candidate = Path(f"{STAGE}.{ext}")',
        "    if candidate.is_file():",
        "        input_path = candidate",
        "        break",
        'tpr = Path(f"{STAGE}.tpr")',
        "if input_path is not None and tpr.is_file():",
        '    print(f"Exporting nojump trajectory to {STAGE}_nojump.gro (from {input_path.name}) ...")',
        "    nojump(",
        "        trajectory=input_path, tpr=tpr,",
        '        output=Path(f"{STAGE}_nojump.gro"),',
        '        group="0",  # System (group index 0)',
        "        ndx=NDX,",
        "    )",
        "",
        'print()',
        'print("UDF / J-OCTA export complete:")',
        'print(f"  {STAGE}_energy.xvg   (gmx energy)")',
        'print(f"  {STAGE}_nojump.gro   (gmx trjconv -pbc nojump)")',
    ]
    script_path = os.path.join(output_dir, "gen_for_udf.py")
    Path(script_path).write_text("\n".join(lines) + "\n")
    os.chmod(script_path, 0o755)
    return script_path
