# -*- coding: utf-8 -*-
"""
abmptools.genesis.grest
-----------------------
GENESIS gREST_SSCR (generalized Replica-Exchange with Solute Tempering --
Solute Side-Chain Repartitioning) builder and post-MD analysis.

Subpackage layout (1.20.0):

- ``models``                — dataclass schema (``RESTSelectionSpec``,
                              ``ReplicaTemperatureSpec``,
                              ``MinimizationStage``,
                              ``EquilibrationStage``, ``GrestStage``,
                              ``GrestBuildConfig``).
- ``rest_selection``        — explicit / around-mode REST residue
                              resolution.
- ``replica_temperatures``  — geometric / manual ladder generators.
- ``inp_writer``            — render minimize / equilibrate / grest /
                              ``remd_convert`` ``.inp`` control files.
- ``system_builder``        — tleap-driven AMBER ``prmtop`` + ``coor``
                              build.
- ``forcefield_check``      — external tool resolution (``tleap``,
                              ``spdyn``, ``atdyn``, ``remd_convert``,
                              ``mpirun``).
- ``builder``               — :class:`GrestBuilder` orchestrating the
                              5-stage pipeline.
- ``grest_runner``          — ``run.sh`` template emission for
                              ``mpirun -np N spdyn``.
- ``analysis``              — replica transition, acceptance ratio,
                              1D distance PMF.
- ``cli``                   — argparse front-end for
                              ``python -m abmptools.genesis.grest``.

External dependencies (subprocess only, never bundled):

- GENESIS spdyn / atdyn / remd_convert >= 2.1 (LGPL-3.0-or-later)
- AmberTools tleap (recommended) and cpptraj (optional)
- mpirun (OpenMPI / MPICH) for parallel replica execution
- matplotlib (optional, for analysis plotting)
"""
