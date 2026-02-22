# -*- coding: utf-8 -*-
"""
abmptools.geomopt.mace_optimizer
-------------------------------
PDB structure optimization using a MACE machine-learning potential via ASE.

Optional runtime dependencies (not required at import time):
    ase       : pip install ase
    mace-torch: pip install mace-torch
    torch     : pip install torch

Usage::

    from abmptools.geomopt import MacePdbOptimizer

    opt = MacePdbOptimizer(model_name="small", device="auto", fmax=0.05)
    result = opt.optimize_pdb("in.pdb", "out.pdb")
    # result = {"energy": -123.4, "steps": 42, "converged": True, "out_pdb": "/abs/out.pdb"}
"""
import logging
import random
from pathlib import Path
from typing import Optional, Union

logger = logging.getLogger(__name__)

_SUPPORTED_OPTIMIZERS = ("BFGS", "FIRE", "LBFGS")


class MacePdbOptimizer:
    """Optimize a PDB structure using a MACE ML potential via ASE.

    All heavy imports (ase, mace, torch) are deferred to method call time so
    that instantiating this class in a mace-free environment does not raise.

    Parameters
    ----------
    model_path : str or None
        Path to a MACE model file (.model or .pt).
        Mutually exclusive with ``model_name``.  When both are given,
        ``model_path`` takes precedence.
    model_name : str or None
        Name of a pre-trained MACE-OFF model ("small", "medium", "large").
        Used when ``model_path`` is None.  Default: ``"small"``.
    device : str
        Compute device.  ``"auto"`` automatically selects CUDA > MPS > CPU.
        Explicit values: ``"cpu"``, ``"cuda"``, ``"mps"``.  Default: ``"auto"``.
    dtype : str
        Floating-point precision for MACE: ``"float32"`` or ``"float64"``.
        Default: ``"float64"``.
    optimizer : str
        ASE optimizer: ``"BFGS"``, ``"FIRE"``, or ``"LBFGS"``.
        Default: ``"BFGS"``.
    fmax : float
        Force convergence threshold [eV/Å].  Default: ``0.05``.
    steps : int
        Maximum number of optimization steps.  Default: ``500``.
    trajectory : str or None
        ASE trajectory file path written when ``write_traj=True`` in
        :meth:`optimize_pdb`.  If *None* the output PDB path with ``.traj``
        extension is used.
    logfile : str or None
        ASE optimizer log file.  ``"-"`` writes to stdout, *None* silences.
    seed : int or None
        If given, seeds ``random``, ``numpy.random``, and ``torch`` for
        reproducibility.
    """

    def __init__(
        self,
        model_path: Optional[str] = None,
        model_name: Optional[str] = "small",
        device: str = "auto",
        dtype: str = "float64",
        optimizer: str = "BFGS",
        fmax: float = 0.05,
        steps: int = 500,
        trajectory: Optional[str] = None,
        logfile: Optional[str] = None,
        seed: Optional[int] = None,
    ) -> None:
        if optimizer not in _SUPPORTED_OPTIMIZERS:
            raise ValueError(
                f"optimizer must be one of {_SUPPORTED_OPTIMIZERS}, "
                f"got {optimizer!r}"
            )
        self.model_path = model_path
        self.model_name = model_name
        self.device = device
        self.dtype = dtype
        self.optimizer = optimizer
        self.fmax = fmax
        self.steps = steps
        self.trajectory = trajectory
        self.logfile = logfile
        self.seed = seed

    # ------------------------------------------------------------------
    # Dependency checks (static; safe to call before heavy imports)
    # ------------------------------------------------------------------

    @staticmethod
    def _check_ase_available() -> None:
        """Raise ImportError with a helpful message if ASE is not installed."""
        try:
            import ase  # noqa: F401
        except ImportError:
            raise ImportError(
                "ASE (Atomic Simulation Environment) is not installed. "
                "Install it with: pip install ase"
            ) from None

    @staticmethod
    def _check_mace_available() -> None:
        """Raise ImportError with a helpful message if mace-torch is not installed."""
        try:
            import mace  # noqa: F401
        except ImportError:
            raise ImportError(
                "mace-torch is not installed. "
                "Install it with: pip install mace-torch\n"
                "See https://github.com/ACEsuit/mace for details."
            ) from None

    @staticmethod
    def _check_torch_available() -> None:
        """Raise ImportError with a helpful message if torch is not installed."""
        try:
            import torch  # noqa: F401
        except ImportError:
            raise ImportError(
                "PyTorch is not installed. "
                "Install it with: pip install torch\n"
                "See https://pytorch.org for platform-specific instructions."
            ) from None

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _resolve_device(self) -> str:
        """Resolve ``"auto"`` to a concrete device string."""
        if self.device != "auto":
            return self.device
        try:
            import torch

            if torch.cuda.is_available():
                return "cuda"
            if hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
                return "mps"
        except ImportError:
            pass
        return "cpu"

    def _set_seed(self) -> None:
        """Seed torch, numpy, and random if ``self.seed`` is set."""
        if self.seed is None:
            return
        random.seed(self.seed)
        try:
            import numpy as np

            np.random.seed(self.seed)
        except ImportError:
            pass
        try:
            import torch

            torch.manual_seed(self.seed)
            if torch.cuda.is_available():
                torch.cuda.manual_seed_all(self.seed)
        except ImportError:
            pass

    def _build_calculator(self, device: str):
        """Build and return a MACE ASE calculator instance."""
        if self.model_path is not None:
            from mace.calculators import MACECalculator

            logger.debug("Loading MACE model from file: %s", self.model_path)
            return MACECalculator(
                model_paths=self.model_path,
                device=device,
                default_dtype=self.dtype,
            )
        else:
            from mace.calculators import mace_off

            model_name = self.model_name or "small"
            logger.debug("Loading MACE-OFF pretrained model: %s", model_name)
            return mace_off(
                model=model_name,
                device=device,
                default_dtype=self.dtype,
            )

    def _get_optimizer_cls(self):
        """Return the ASE optimizer class matching ``self.optimizer``."""
        from ase.optimize import BFGS, FIRE, LBFGS

        mapping = {"BFGS": BFGS, "FIRE": FIRE, "LBFGS": LBFGS}
        return mapping[self.optimizer]

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def optimize_pdb(
        self,
        in_pdb: Union[str, Path],
        out_pdb: Union[str, Path],
        *,
        write_traj: bool = False,
    ) -> dict:
        """Run geometry optimization on a single PDB file.

        Parameters
        ----------
        in_pdb : str or Path
            Input PDB file path.
        out_pdb : str or Path
            Output PDB file path.  Must differ from *in_pdb* to prevent
            accidental overwriting of the original structure.
        write_traj : bool
            If True, write an ASE ``.traj`` trajectory alongside *out_pdb*.

        Returns
        -------
        dict
            ``energy`` (float, eV), ``steps`` (int), ``converged`` (bool),
            ``out_pdb`` (str, absolute path).

        Raises
        ------
        ImportError
            If ASE, mace-torch, or torch are not installed.
        FileNotFoundError
            If *in_pdb* does not exist.
        ValueError
            If *in_pdb* and *out_pdb* resolve to the same path.
        RuntimeError
            If the ASE optimization step raises unexpectedly.
        """
        self._check_ase_available()
        self._check_mace_available()
        self._check_torch_available()

        from ase.io import read
        from ase.io import write as ase_write

        in_pdb = Path(in_pdb).resolve()
        out_pdb = Path(out_pdb).resolve()

        if not in_pdb.exists():
            raise FileNotFoundError(f"Input PDB not found: {in_pdb}")
        if in_pdb == out_pdb:
            raise ValueError(
                "in_pdb and out_pdb must be different paths to avoid "
                f"overwriting the input structure.  Got: {in_pdb}"
            )

        logger.info("Optimizing: %s -> %s", in_pdb, out_pdb)
        logger.info(
            "  model=%s  device=%s  dtype=%s  optimizer=%s  fmax=%.4f  steps=%d",
            self.model_path or self.model_name,
            self.device,
            self.dtype,
            self.optimizer,
            self.fmax,
            self.steps,
        )

        self._set_seed()
        device = self._resolve_device()
        logger.debug("Resolved device: %s", device)

        # --- Read structure ---
        atoms = read(str(in_pdb), format="proteindatabank")
        logger.debug("Read %d atoms from %s", len(atoms), in_pdb.name)

        # --- Attach MACE calculator ---
        atoms.calc = self._build_calculator(device)

        # --- Prepare optimizer ---
        opt_cls = self._get_optimizer_cls()
        traj_path: Optional[str] = None
        if write_traj:
            traj_path = self.trajectory or str(out_pdb.with_suffix(".traj"))

        # ASE writes log to "/dev/null" (POSIX) or os.devnull when silenced
        import os

        logfile = self.logfile if self.logfile is not None else os.devnull

        # --- Run optimization ---
        try:
            dyn = opt_cls(atoms, logfile=logfile, trajectory=traj_path)
            converged = dyn.run(fmax=self.fmax, steps=self.steps)
        except Exception as exc:
            raise RuntimeError(
                f"ASE optimization failed for {in_pdb.name}: {exc}"
            ) from exc

        energy = float(atoms.get_potential_energy())
        n_steps = int(dyn.nsteps)

        if converged:
            logger.info(
                "Converged in %d steps, energy = %.6f eV", n_steps, energy
            )
        else:
            logger.warning(
                "Did not converge within %d steps for %s "
                "(last energy: %.6f eV)",
                self.steps,
                in_pdb.name,
                energy,
            )

        # --- Write output PDB ---
        out_pdb.parent.mkdir(parents=True, exist_ok=True)
        ase_write(str(out_pdb), atoms, format="proteindatabank")
        logger.info("Written: %s", out_pdb)

        return {
            "energy": energy,
            "steps": n_steps,
            "converged": converged,
            "out_pdb": str(out_pdb),
        }
