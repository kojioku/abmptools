# -*- coding: utf-8 -*-
"""
mdp_writer.py
-------------
Writes the GROMACS .mdp simulation parameter file from a SystemModel.
"""
from __future__ import annotations
from ....core.system_model import SystemModel, SimulationParams


def _strl(total: str, data: str, n: int) -> str:
    return total + data.ljust(n)[:n] + " "


class MdpWriter:
    """Generates the content of a .mdp file."""

    def write(self, model: SystemModel, filepath: str) -> None:
        content = self._build(model)
        with open(filepath, "w") as f:
            f.write(content)

    def _build(self, model: SystemModel) -> str:
        sp = model.sim_params
        s = ""

        s = _strl(s, "title", 20) + " = " + sp.title + "\n"
        s = _strl(s, "cpp",   20) + " = " + "/lib/cpp" + "\n"

        # vel_gen
        if sp.vel_gen:
            s = _strl(s, "gen_vel",  20) + " = yes\n"
            s = _strl(s, "gen_temp", 20) + " = " + str(sp.gen_temp) + "\n"

        # integrator
        s = _strl(s, "integrator", 20) + " = " + sp.integrator + "\n"
        if sp.ld_seed is not None:
            s = _strl(s, "ld_seed", 20) + " = " + str(sp.ld_seed) + "\n"

        # steps / output
        s = _strl(s, "nsteps",    20) + " = " + str(sp.nsteps) + "\n"
        s = _strl(s, "nstlist",   20) + " = 10\n"
        s = _strl(s, "nstxout",   20) + " = " + str(sp.outputinterval2) + "\n"
        s = _strl(s, "nstvout",   20) + " = " + str(sp.outputinterval2) + "\n"
        s = _strl(s, "nstxtcout", 20) + " = " + str(sp.outputinterval) + "\n"
        s = _strl(s, "nstlog",    20) + " = " + str(sp.outputinterval) + "\n"
        output_energy = max(1, int(sp.outputinterval / 1000))
        s = _strl(s, "nstenergy", 20) + " = " + str(output_energy) + "\n"
        s = _strl(s, "dt",        20) + " = " + "%.8f" % sp.dt + "\n"

        # constraints
        if sp.rattle_bond and not sp.rattle_angle:
            s = _strl(s, "constraints",         20) + " = all-bonds\n"
            s = _strl(s, "constraint-algorithm", 20) + " = Lincs\n"
        elif not sp.rattle_bond and sp.rattle_angle:
            s = _strl(s, "constraints",         20) + " = all-angles\n"
            s = _strl(s, "constraint-algorithm", 20) + " = Lincs\n"
        elif sp.rattle_bond and sp.rattle_angle:
            print(" Error: Please choose one constraint type (Bond or Angle)")
        else:
            s = _strl(s, "constraints", 20) + " = none\n"

        s = _strl(s, "ns_type", 20) + " = grid\n"

        # electrostatics
        if sp.calcQQ == 1:
            if sp.qq_algorithm == "Ewald":
                print(" Electrostatic algorithm was changed to PME.")
                s = _strl(s, "coulombtype",      20) + " = PME\n"
                s = _strl(s, "fourier-spacing",  20) + " = 0.1\n"
            else:
                print(" Electrostatic algorithm was changed to Cutoff_Coulomb.")
                s = _strl(s, "coulombtype", 20) + " = cut-off\n"
            s = _strl(s, "coulomb-modifier", 20) + " = None\n"
        else:
            s = _strl(s, "coulombtype",       20) + " = cut-off\n"
            s = _strl(s, "coulomb-modifier",  20) + " = None\n"

        # cutoffs
        s = _strl(s, "rlist",   20) + " = " + str(sp.lj_cutoff)      + "\n"
        s = _strl(s, "rcoulomb", 20) + " = " + str(sp.coulomb_cutoff) + "\n"
        s = _strl(s, "rvdw",    20) + " = " + str(sp.coulomb_cutoff) + "\n"

        if sp.tail_correction == 1:
            s = _strl(s, "dispcorr", 20) + " = EnerPres\n"

        s = _strl(s, "vdwtype",      20) + " = Cut-off\n"
        s = _strl(s, "vdw-modifier", 20) + " = Potential-shift\n"
        s = _strl(s, "cutoff-scheme", 20) + " = verlet\n"

        # temperature coupling
        s = _strl(s, "tcoupl", 20) + " = " + sp.t_coupl + "\n"
        if sp.t_coupl != "no":
            s = _strl(s, "nsttcouple",    20) + " = 1\n"
            s = _strl(s, "nh-chain-length", 20) + " = 1\n"
            s = _strl(s, "tc_grps",       20) + " = system\n"
            s = _strl(s, "tau_t",         20) + " = " + str(sp.tau_t) + " \n"
            s = _strl(s, "ref_t",         20) + " = " + str(sp.ref_t) + " \n"

        # pressure coupling
        s = _strl(s, "pcoupl", 20) + " = " + sp.p_coupl + "\n"
        if sp.p_coupl != "no":
            s = _strl(s, "nstpcouple",    20) + " = 1\n"
            s = _strl(s, "pcoupltype",    20) + " = " + sp.pcoupltype + "\n"

            s = _strl(s, "tau_p", 20) + " = " + str(sp.tau_p) + " \n"

            s = _strl(s, "ref_p", 20) + " = "
            if sp.ref_p_tensor is not None:
                for v in sp.ref_p_tensor:
                    s += str(v) + " "
            else:
                s += str(sp.ref_p) + " "
            s += "\n"

            s = _strl(s, "compressibility", 20) + " = "
            if sp.compressibility_tensor is not None:
                for v in sp.compressibility_tensor:
                    s += str(v) + " "
                s += "\n"
            else:
                s += str(sp.compressibility) + "\n"

        # pbc
        s = _strl(s, "pbc", 20) + " = " + sp.pbc + "\n"
        if sp.periodic_mol:
            s = _strl(s, "periodic-molecules", 20) + " = yes\n"

        # deformation
        if sp.deform_vel is not None:
            s = _strl(s, "deform", 20) + " = "
            for v in sp.deform_vel:
                s += "{} ".format(v)
            s += "\n"

        # constraint freeze
        if sp.freeze_grps is not None:
            s = _strl(s, "freezegrps", 20) + " = " + sp.freeze_grps + "\n"
            s = _strl(s, "freezedim",  20) + " = " + sp.freeze_dim  + "\n"

        return s
