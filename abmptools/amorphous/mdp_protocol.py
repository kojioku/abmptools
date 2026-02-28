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
    return {
        "tcoupl": "V-rescale",
        "tc-grps": tc_grps,
        "tau-t": protocol.tau_t,
        "ref-t": ref_t,
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
    lines = [
        "#!/bin/bash",
        "set -e",
        "",
        f'BUILD="{build}"',
        f'GRO="$BUILD/{gro}"',
        f'TOP="$BUILD/{top}"',
        "",
        "# Stage 1: Energy minimisation",
        f'gmx grompp -f 01_em.mdp -c "$GRO" -p "$TOP"{ndx_flag} -o 01_em.tpr -maxwarn 2',
        "gmx mdrun -deffnm 01_em",
        "",
        "# Stage 2: NVT high-T",
        f'gmx grompp -f 02_nvt_highT.mdp -c 01_em.gro -p "$TOP"{ndx_flag} -o 02_nvt_highT.tpr -maxwarn 2',
        "gmx mdrun -deffnm 02_nvt_highT",
        "",
        "# Stage 3: NPT high-T",
        f'gmx grompp -f 03_npt_highT.mdp -c 02_nvt_highT.gro -p "$TOP"{ndx_flag} -o 03_npt_highT.tpr -maxwarn 2',
        "gmx mdrun -deffnm 03_npt_highT",
        "",
        "# Stage 4: Simulated annealing",
        f'gmx grompp -f 04_anneal.mdp -c 03_npt_highT.gro -p "$TOP"{ndx_flag} -o 04_anneal.tpr -maxwarn 2',
        "gmx mdrun -deffnm 04_anneal",
        "",
        "# Stage 5: NPT final equilibration",
        f'gmx grompp -f 05_npt_final.mdp -c 04_anneal.gro -p "$TOP"{ndx_flag} -o 05_npt_final.tpr -maxwarn 2',
        "gmx mdrun -deffnm 05_npt_final",
        "",
        'echo "All stages completed successfully."',
    ]
    script_path = os.path.join(output_dir, "run_all.sh")
    Path(script_path).write_text("\n".join(lines) + "\n")
    os.chmod(script_path, 0o755)
    return script_path
