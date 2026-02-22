# -*- coding: utf-8 -*-
"""
abmptools.mlopt
---------------
Structure optimization and minimization utilities.

Both classes are importable even when their optional runtime dependencies
are not installed; dependency errors are deferred to method call time.

Public API::

    from abmptools.mlopt import MacePdbOptimizer
    from abmptools.mlopt import OpenFFOpenMMMinimizer

    # MACE ML potential (requires ase, mace-torch, torch)
    opt = MacePdbOptimizer(model_name="small", device="auto")
    result = opt.optimize_pdb("in.pdb", "out.pdb")
    # {"energy": float[eV], "steps": int, "converged": bool, "out_pdb": str}

    # OpenFF + OpenMM force field (requires openmm, openff-toolkit)
    mini = OpenFFOpenMMMinimizer(forcefield="openff_unconstrained-2.1.0.offxml")
    result = mini.minimize_pdb("in.pdb", "out.pdb")
    # {"energy_before": float[kJ/mol], "energy_after": float[kJ/mol],
    #  "energy": float[kJ/mol], "converged": bool, "elapsed": float, "out_pdb": str}

    # PySCF DFT QM optimiser (requires pyscf, geometric or pyberny)
    opt = QMOptimizerPySCF(functional="B3LYP", basis="def2-SVP",
                            dispersion="d3bj")
    result = opt.optimize("in.xyz", "out.xyz")
    # {"energy": float[eV], "energy_hartree": float[Ha],
    #  "steps": int, "converged": bool, "out_xyz": str}
"""
from .mace_optimizer import MacePdbOptimizer
from .openff_openmm_minimizer import OpenFFOpenMMMinimizer
from .pyscf_optimizer import QMOptimizerPySCF

__all__ = ["MacePdbOptimizer", "OpenFFOpenMMMinimizer", "QMOptimizerPySCF"]
