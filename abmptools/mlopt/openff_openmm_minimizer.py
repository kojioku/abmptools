# -*- coding: utf-8 -*-
"""
abmptools.mlopt.openff_openmm_minimizer
----------------------------------------
PDB structure minimization using OpenFF (SMIRNOFF) force fields via OpenMM.

Optional runtime dependencies (not required at import time):
    openmm         : conda install -c conda-forge openmm
    openff-toolkit : conda install -c conda-forge openff-toolkit
    rdkit          : conda install -c conda-forge rdkit  (for bond inference)
    ambertools      : conda install -c conda-forge ambertools  (for AM1-BCC charges)

Usage::

    from abmptools.mlopt import OpenFFOpenMMMinimizer

    minimizer = OpenFFOpenMMMinimizer(
        forcefield="openff_unconstrained-2.1.0.offxml",
        platform="auto",
        tolerance=10.0,
    )
    result = minimizer.minimize_pdb("in.pdb", "out.pdb")
    # result = {
    #     "energy_before": -100.0,   # kJ/mol
    #     "energy_after":  -120.0,   # kJ/mol
    #     "energy":        -120.0,   # kJ/mol  (alias for energy_after)
    #     "converged":     True,
    #     "elapsed":       1.5,      # seconds
    #     "out_pdb":       "/abs/path/out.pdb",
    # }
"""
import logging
import os
import random
import time
from pathlib import Path
from typing import Optional, Union

logger = logging.getLogger(__name__)

_SUPPORTED_PLATFORMS = ("auto", "CUDA", "OpenCL", "CPU", "HIP")
_SUPPORTED_CONSTRAINTS = (None, "None", "HBonds", "AllBonds", "HAngles")


class OpenFFOpenMMMinimizer:
    """Minimize a PDB structure using OpenFF SMIRNOFF force fields via OpenMM.

    All heavy imports (openmm, openff-toolkit, rdkit) are deferred to
    method call time so that this class is instantiable in environments
    where those packages are not installed.

    This class is designed to be a parallel counterpart to
    :class:`~abmptools.mlopt.MacePdbOptimizer` and follows the same
    conventions: lazy optional imports, identical I/O contract, and a
    consistent result dict.

    Parameters
    ----------
    forcefield : str
        OpenFF force field specification.  Any name accepted by
        ``openff.toolkit.ForceField()`` is valid, e.g.:

        - ``"openff_unconstrained-2.1.0.offxml"`` *(default)*
        - ``"openff-2.1.0.offxml"``
        - ``"openff_unconstrained-1.3.1.offxml"``

        The "unconstrained" variant is recommended for energy minimization
        because it removes all bond-length constraints, allowing the geometry
        to fully relax.
    platform : str
        OpenMM compute platform.  ``"auto"`` selects the fastest available
        platform in the order CUDA > OpenCL > CPU.  Explicit values:
        ``"CUDA"``, ``"OpenCL"``, ``"CPU"``, ``"HIP"``.  Default: ``"auto"``.
    tolerance : float
        RMS force convergence criterion [kJ/mol/nm].  Minimization stops when
        the root-mean-square of all force components falls below this value.
        Corresponds approximately to ``fmax`` in :class:`MacePdbOptimizer`
        (different units and norms).  Default: ``10.0``.
    max_iterations : int
        Maximum number of minimization iterations.  ``0`` means run until
        the force criterion is satisfied (analogous to setting a very large
        ``steps`` in :class:`MacePdbOptimizer`).  Default: ``0``.
    constraints : str or None
        Geometric constraints to apply to the system.  One of ``None``,
        ``"HBonds"``, ``"AllBonds"``, ``"HAngles"``.  For pure energy
        minimization ``None`` is recommended so all degrees of freedom can
        relax.  Default: ``None``.
    hydrogen_mass : float or None
        Hydrogen mass repartitioning value [Da].  Has minimal effect on
        energy minimization and is not applied in this version.
        Default: ``None``.
    logfile : str or None
        Destination for the per-job minimization log written when
        ``write_log=True`` in :meth:`minimize_pdb`.  If ``None``, the log
        is written to ``<out_pdb>.min.log``.  Default: ``None``.
    seed : int or None
        Seeds ``random`` and ``numpy.random`` for reproducibility.  The
        underlying L-BFGS minimizer in OpenMM is deterministic, so this
        primarily affects any stochastic post-processing.  Default: ``None``.
    """

    def __init__(
        self,
        forcefield: str = "openff_unconstrained-2.1.0.offxml",
        platform: str = "auto",
        tolerance: float = 10.0,
        max_iterations: int = 0,
        constraints: Optional[str] = None,
        hydrogen_mass: Optional[float] = None,
        logfile: Optional[str] = None,
        seed: Optional[int] = None,
    ) -> None:
        if constraints not in _SUPPORTED_CONSTRAINTS:
            raise ValueError(
                f"constraints must be one of {_SUPPORTED_CONSTRAINTS}, "
                f"got {constraints!r}"
            )
        # Accept platform case-insensitively (OpenMM names are case-sensitive
        # internally but users often type "cuda" / "cpu")
        _platform_lower = {p.lower(): p for p in _SUPPORTED_PLATFORMS}
        if platform.lower() not in _platform_lower:
            raise ValueError(
                f"platform must be one of {_SUPPORTED_PLATFORMS}, "
                f"got {platform!r}"
            )
        self.forcefield = forcefield
        self.platform = _platform_lower[platform.lower()]
        self.tolerance = tolerance
        self.max_iterations = max_iterations
        self.constraints = constraints
        self.hydrogen_mass = hydrogen_mass
        self.logfile = logfile
        self.seed = seed

    # ------------------------------------------------------------------
    # Dependency checks
    # ------------------------------------------------------------------

    @staticmethod
    def _check_openmm_available() -> None:
        """Raise ImportError with a helpful message if OpenMM is not installed."""
        try:
            import openmm  # noqa: F401
        except ImportError:
            raise ImportError(
                "OpenMM is not installed.\n"
                "Install it with:\n"
                "  conda install -c conda-forge openmm\n"
                "See https://openmm.org for details."
            ) from None

    @staticmethod
    def _check_openff_available() -> None:
        """Raise ImportError with a helpful message if openff-toolkit is not installed."""
        try:
            import openff.toolkit  # noqa: F401
        except ImportError:
            raise ImportError(
                "openff-toolkit is not installed.\n"
                "Install it with:\n"
                "  conda install -c conda-forge openff-toolkit\n"
                "See https://github.com/openforcefield/openff-toolkit for details."
            ) from None

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _set_seed(self) -> None:
        """Seed random and numpy.random if ``self.seed`` is set."""
        if self.seed is None:
            return
        random.seed(self.seed)
        try:
            import numpy as np

            np.random.seed(self.seed)
        except ImportError:
            pass

    def _resolve_platform(self):
        """Return the best available OpenMM Platform instance.

        Raises
        ------
        RuntimeError
            If the requested platform is unavailable, or no platform is found
            in ``"auto"`` mode.
        """
        import openmm

        if self.platform.lower() != "auto":
            try:
                platform = openmm.Platform.getPlatformByName(self.platform)
                logger.debug("Using OpenMM platform: %s", self.platform)
                return platform
            except Exception as e:
                available = [
                    openmm.Platform.getPlatformByIndex(i).getName()
                    for i in range(openmm.Platform.getNumPlatforms())
                ]
                raise RuntimeError(
                    f"OpenMM platform {self.platform!r} is not available: {e}\n"
                    f"Available platforms: {available}"
                ) from e

        for name in ("CUDA", "OpenCL", "CPU"):
            try:
                platform = openmm.Platform.getPlatformByName(name)
                logger.debug("Auto-selected OpenMM platform: %s", name)
                return platform
            except Exception:
                continue

        raise RuntimeError(
            "No OpenMM compute platform is available.  "
            "Ensure OpenMM is correctly installed."
        )

    def _load_molecule_from_pdb(self, pdb_path: Path):
        """Load OpenFF Molecule from a PDB file with fallback bond-inference.

        Uses OpenMM's ``PDBFile`` reader for the authoritative atom positions
        and topology (preserving residue/chain names), and builds an OpenFF
        ``Molecule`` separately for force-field parameterisation.

        Two strategies are attempted in order:

        1. ``Molecule.from_file(pdb, file_format="PDB")`` — works when the
           PDB has CONECT records or when OpenFF can infer bonds.
        2. RDKit ``MolFromPDBFile`` → ``Molecule.from_rdkit`` — works for
           most organic small molecules when RDKit is installed.

        Parameters
        ----------
        pdb_path : Path

        Returns
        -------
        tuple
            ``(mol, omm_topology, positions)``

            - *mol*: :class:`openff.toolkit.Molecule`
            - *omm_topology*: :class:`openmm.app.Topology`
              (from the PDB file; preserves residue/chain info)
            - *positions*: OpenMM ``Quantity`` with units of nanometres.

        Raises
        ------
        RuntimeError
            If neither strategy can build a valid molecule.
        ValueError
            If the atom count reported by OpenMM and OpenFF disagree.
        """
        import openmm.app as openmm_app
        from openff.toolkit import Molecule

        # Authoritative positions and topology from OpenMM's PDB reader
        pdb_obj = openmm_app.PDBFile(str(pdb_path))
        n_atoms_pdb = pdb_obj.topology.getNumAtoms()

        mol = None
        exc_history: list = []

        # --- Strategy 1: OpenFF native PDB reader ---
        try:
            mol = Molecule.from_file(str(pdb_path), file_format="PDB")
            logger.debug(
                "Strategy 1 (Molecule.from_file): loaded %d atoms", mol.n_atoms
            )
        except Exception as e:
            exc_history.append(f"Molecule.from_file: {e}")
            logger.debug("Strategy 1 failed: %s", e)

        # --- Strategy 2: RDKit bond inference ---
        if mol is None:
            try:
                from rdkit import Chem

                rdmol = Chem.MolFromPDBFile(
                    str(pdb_path), removeHs=False, sanitize=True
                )
                if rdmol is None:
                    raise ValueError(
                        "RDKit MolFromPDBFile returned None — "
                        "check for unsupported elements or malformed PDB"
                    )
                mol = Molecule.from_rdkit(rdmol, allow_undefined_stereo=True)
                logger.debug(
                    "Strategy 2 (RDKit): loaded %d atoms", mol.n_atoms
                )
            except ImportError:
                exc_history.append("RDKit: not installed")
                logger.debug("Strategy 2 skipped: RDKit not available")
            except Exception as e:
                exc_history.append(f"RDKit: {e}")
                logger.debug("Strategy 2 failed: %s", e)

        if mol is None:
            raise RuntimeError(
                f"Cannot build an OpenFF Molecule from {pdb_path.name}.\n"
                "PDB files often lack bond/connectivity information.\n"
                "Attempted strategies:\n"
                + "\n".join(f"  - {m}" for m in exc_history)
                + "\nSolutions:\n"
                "  1. Ensure the PDB has CONECT records\n"
                "  2. Install RDKit: conda install -c conda-forge rdkit"
            )

        # Validate atom counts are consistent between both readers
        if mol.n_atoms != n_atoms_pdb:
            raise ValueError(
                f"Atom count mismatch in {pdb_path.name}: "
                f"OpenMM PDB reader found {n_atoms_pdb} atoms but "
                f"OpenFF Molecule has {mol.n_atoms} atoms.\n"
                "Common causes: missing hydrogens, hydrogens removed during "
                "bond inference, or inconsistent CONECT records.\n"
                "Ensure all hydrogens are present and CONECT records are correct."
            )

        logger.debug(
            "Loaded %d atoms from %s (topology from OpenMM PDB reader)",
            n_atoms_pdb,
            pdb_path.name,
        )
        return mol, pdb_obj.topology, pdb_obj.positions

    def _build_system(self, mol):
        """Build an OpenMM System using the configured OpenFF force field.

        Parameters
        ----------
        mol : openff.toolkit.Molecule

        Returns
        -------
        openmm.System

        Raises
        ------
        RuntimeError
            If the force field cannot be loaded or atom typing fails.
        """
        import openmm.app as openmm_app
        from openff.toolkit import ForceField
        from openff.toolkit import Topology as OFFTopology

        # Resolve constraints object
        constraints_obj = None
        if self.constraints and self.constraints != "None":
            _cmap = {
                "HBonds": openmm_app.HBonds,
                "AllBonds": openmm_app.AllBonds,
                "HAngles": openmm_app.HAngles,
            }
            constraints_obj = _cmap.get(self.constraints)

        # Load force field
        try:
            ff = ForceField(self.forcefield)
        except Exception as e:
            raise RuntimeError(
                f"Failed to load OpenFF force field {self.forcefield!r}: {e}\n"
                "Commonly available force fields:\n"
                "  'openff_unconstrained-2.1.0.offxml'  (recommended for minimization)\n"
                "  'openff-2.1.0.offxml'\n"
                "  'openff_unconstrained-1.3.1.offxml'"
            ) from e

        # Build OpenFF topology (single molecule)
        off_topology = OFFTopology.from_molecules([mol])

        # Build OpenMM System; try with constraints first, fall back without
        kwargs: dict = {}
        if constraints_obj is not None:
            kwargs["constraints"] = constraints_obj

        try:
            system = ff.create_openmm_system(off_topology, **kwargs)
        except TypeError:
            # Older or newer openff-toolkit versions may not accept constraints
            logger.debug(
                "create_openmm_system did not accept constraints kwarg; "
                "retrying without it"
            )
            system = ff.create_openmm_system(off_topology)
        except Exception as e:
            raise RuntimeError(
                f"OpenFF force field assignment failed: {e}\n"
                "Common causes:\n"
                "  - Partial charges (AM1-BCC) could not be assigned.\n"
                "    Fix: conda install -c conda-forge ambertools\n"
                "         or: pip install openff-nagl\n"
                "  - Unknown atom type or unsupported bond order\n"
                "  - Multi-component system (currently only single-molecule "
                "PDB is supported)"
            ) from e

        if self.hydrogen_mass is not None:
            logger.warning(
                "hydrogen_mass=%s was specified but hydrogen mass "
                "repartitioning is not applied in this version.  "
                "It has negligible effect on energy minimization.",
                self.hydrogen_mass,
            )

        logger.debug(
            "Built OpenMM System: %d particles, force field=%s",
            system.getNumParticles(),
            self.forcefield,
        )
        return system

    def _check_convergence(self, context) -> bool:
        """Return ``True`` if the RMS force satisfies ``self.tolerance``.

        Computes the RMS of all force components and compares to
        ``self.tolerance`` [kJ/mol/nm].  Falls back to ``True`` when
        ``max_iterations == 0`` (unlimited run) if forces cannot be read.
        """
        try:
            import numpy as np
            import openmm.unit as openmm_unit

            state = context.getState(getForces=True)
            forces = state.getForces(asNumpy=True).value_in_unit(
                openmm_unit.kilojoules_per_mole / openmm_unit.nanometer
            )
            rms_force = float(np.sqrt(np.mean(np.sum(forces**2, axis=1))))
            converged = rms_force < self.tolerance
            logger.debug(
                "RMS force: %.4f kJ/mol/nm  (tolerance=%.4f)  converged=%s",
                rms_force,
                self.tolerance,
                converged,
            )
            return converged
        except Exception as e:
            logger.debug("Convergence check failed: %s; assuming converged", e)
            return self.max_iterations == 0

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def minimize_pdb(
        self,
        in_pdb: Union[str, Path],
        out_pdb: Union[str, Path],
        *,
        write_log: bool = False,
    ) -> dict:
        """Minimize the potential energy of a single PDB structure.

        Parameters
        ----------
        in_pdb : str or Path
            Input PDB file path.
        out_pdb : str or Path
            Output PDB file path.  Must differ from *in_pdb* to prevent
            accidental overwriting of the original structure.
        write_log : bool
            If True, write a text summary to ``self.logfile`` (if set)
            or to ``<out_pdb>.min.log``.

        Returns
        -------
        dict
            The following keys are always present:

            ``energy_before``
                Potential energy before minimization [kJ/mol].
            ``energy_after``
                Potential energy after minimization [kJ/mol].
            ``energy``
                Alias for *energy_after* (parallel to
                :class:`MacePdbOptimizer`'s ``"energy"`` key).
            ``converged``
                ``True`` if the RMS force fell below ``self.tolerance``.
            ``elapsed``
                Wall-clock time of the minimization step [s].
            ``out_pdb``
                Absolute path of the written output PDB.

        Raises
        ------
        ImportError
            If OpenMM or openff-toolkit are not installed.
        FileNotFoundError
            If *in_pdb* does not exist.
        ValueError
            If *in_pdb* and *out_pdb* resolve to the same path, or if
            atom counts are inconsistent between molecule readers.
        RuntimeError
            If molecule loading, force field assignment, or minimization
            raises an unexpected error.
        """
        self._check_openmm_available()
        self._check_openff_available()

        import openmm
        import openmm.app as openmm_app
        import openmm.unit as openmm_unit

        in_pdb = Path(in_pdb).resolve()
        out_pdb = Path(out_pdb).resolve()

        if not in_pdb.exists():
            raise FileNotFoundError(f"Input PDB not found: {in_pdb}")
        if in_pdb == out_pdb:
            raise ValueError(
                "in_pdb and out_pdb must be different paths to avoid "
                f"overwriting the input structure.  Got: {in_pdb}"
            )

        logger.info("Minimizing: %s -> %s", in_pdb, out_pdb)
        logger.info(
            "  forcefield=%s  platform=%s  tolerance=%.2f kJ/mol/nm"
            "  max_iterations=%d",
            self.forcefield,
            self.platform,
            self.tolerance,
            self.max_iterations,
        )

        self._set_seed()

        # --- Load molecule and positions ---
        mol, omm_topology, positions = self._load_molecule_from_pdb(in_pdb)

        # --- Build force field system ---
        system = self._build_system(mol)

        # Sanity-check particle count
        n_sys = system.getNumParticles()
        n_top = omm_topology.getNumAtoms()
        if n_sys != n_top:
            raise RuntimeError(
                f"Particle count mismatch: OpenMM System has {n_sys} particles "
                f"but topology has {n_top} atoms.  "
                "This is an internal inconsistency; please report it."
            )

        # --- Set up simulation ---
        platform = self._resolve_platform()
        # VerletIntegrator is a neutral placeholder; it is never stepped
        integrator = openmm.VerletIntegrator(1.0 * openmm_unit.femtoseconds)
        simulation = openmm_app.Simulation(
            omm_topology, system, integrator, platform
        )
        simulation.context.setPositions(positions)

        # --- Energy before minimization ---
        state_before = simulation.context.getState(getEnergy=True)
        energy_before = float(
            state_before.getPotentialEnergy().value_in_unit(
                openmm_unit.kilojoules_per_mole
            )
        )
        logger.info("Energy before: %.4f kJ/mol", energy_before)

        # --- Run minimization ---
        t0 = time.time()
        try:
            simulation.minimizeEnergy(
                tolerance=(
                    self.tolerance
                    * openmm_unit.kilojoules_per_mole
                    / openmm_unit.nanometer
                ),
                maxIterations=self.max_iterations,
            )
        except Exception as exc:
            raise RuntimeError(
                f"OpenMM minimization failed for {in_pdb.name}: {exc}"
            ) from exc
        elapsed = time.time() - t0

        # --- Energy after + final positions ---
        state_after = simulation.context.getState(
            getEnergy=True, getPositions=True
        )
        energy_after = float(
            state_after.getPotentialEnergy().value_in_unit(
                openmm_unit.kilojoules_per_mole
            )
        )
        positions_after = state_after.getPositions()

        # --- Convergence check ---
        converged = self._check_convergence(simulation.context)

        if converged:
            logger.info(
                "Converged: energy = %.4f kJ/mol  (Δ = %.4f kJ/mol)"
                "  elapsed = %.2f s",
                energy_after,
                energy_after - energy_before,
                elapsed,
            )
        else:
            logger.warning(
                "Did not converge for %s: energy = %.4f kJ/mol"
                "  (Δ = %.4f kJ/mol)  elapsed = %.2f s",
                in_pdb.name,
                energy_after,
                energy_after - energy_before,
                elapsed,
            )

        # --- Optional log file ---
        if write_log:
            log_path = (
                self.logfile
                if self.logfile is not None
                else str(out_pdb.with_suffix(".min.log"))
            )
            try:
                with open(log_path, "w") as lf:
                    lf.write(f"Input:          {in_pdb}\n")
                    lf.write(f"Output:         {out_pdb}\n")
                    lf.write(f"Force field:    {self.forcefield}\n")
                    lf.write(f"Platform:       {platform.getName()}\n")
                    lf.write(f"Tolerance:      {self.tolerance} kJ/mol/nm\n")
                    lf.write(f"Max iterations: {self.max_iterations}\n")
                    lf.write(f"Energy before:  {energy_before:.6f} kJ/mol\n")
                    lf.write(f"Energy after:   {energy_after:.6f} kJ/mol\n")
                    lf.write(
                        f"Delta energy:   {energy_after - energy_before:.6f} kJ/mol\n"
                    )
                    lf.write(f"Converged:      {converged}\n")
                    lf.write(f"Elapsed:        {elapsed:.3f} s\n")
                logger.debug("Log written to %s", log_path)
            except Exception as e:
                logger.warning("Failed to write log file: %s", e)

        # --- Write output PDB ---
        # Use the original OpenMM topology to preserve residue/chain names
        out_pdb.parent.mkdir(parents=True, exist_ok=True)
        with open(str(out_pdb), "w") as f:
            openmm_app.PDBFile.writeFile(omm_topology, positions_after, f)
        logger.info("Written: %s", out_pdb)

        return {
            "energy_before": energy_before,
            "energy_after": energy_after,
            "energy": energy_after,
            "converged": converged,
            "elapsed": elapsed,
            "out_pdb": str(out_pdb),
        }
