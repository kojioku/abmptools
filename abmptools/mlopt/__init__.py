# -*- coding: utf-8 -*-
"""
abmptools.mlopt
---------------
Machine-learning-based structure optimization utilities.

This subpackage is importable even when optional dependencies
(ase, mace-torch, torch) are not installed.  The heavy imports
are deferred to method call time inside each class.

Public API::

    from abmptools.mlopt import MacePdbOptimizer

    opt = MacePdbOptimizer(model_name="small", device="auto")
    result = opt.optimize_pdb("in.pdb", "out.pdb")
"""
from .mace_optimizer import MacePdbOptimizer

__all__ = ["MacePdbOptimizer"]
