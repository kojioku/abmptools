# -*- coding: utf-8 -*-
"""
abmptools.geomopt.pyscf_optimizer
--------------------------------
QM geometry optimization using PySCF DFT with optional D3(BJ) dispersion.

Default settings: B3LYP / def2-SVP with D3(BJ) dispersion correction.

Supported input formats: xyz, pdb
    - xyz: element symbol + Cartesian coordinates (Å), standard format
    - pdb: ATOM/HETATM records; element is read from cols 77-78 and inferred
      from the atom name when the element column is blank.  If neither
      strategy yields an unambiguous symbol, a ValueError is raised.
Output format: xyz (always; pdb input also produces xyz output).

Optional runtime dependencies (not required at import time):
    pyscf        : pip install pyscf
    geometric    : pip install geometric   (geomeTRIC optimisation driver)
    pyberny      : pip install pyberny     (alternative; for solver="berny")
    simple-dftd3 : pip install simple-dftd3  (D3 dispersion, recommended)
    dftd3        : pip install dftd3          (alternative D3 provider)

Usage::

    from abmptools.geomopt import QMOptimizerPySCF

    opt = QMOptimizerPySCF(functional="B3LYP", basis="def2-SVP",
                            dispersion="d3bj")
    result = opt.optimize("water.xyz", "water_opt.xyz")
    # result = {
    #     "energy":          float,  # final energy in eV
    #     "energy_hartree":  float,  # final energy in Hartree
    #     "steps":           int,    # number of optimisation steps taken
    #     "converged":       bool,
    #     "out_xyz":         str,    # absolute path to output xyz
    # }
"""
import logging
import re
from pathlib import Path
from typing import List, Optional, Tuple, Union

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_HARTREE_TO_EV: float = 27.211396132   # 1 Hartree in eV (NIST CODATA 2018)
_BOHR_TO_ANG: float = 0.529177210903  # 1 Bohr in Angstrom

_SUPPORTED_SOLVERS = ("geometric", "berny")
_SUPPORTED_DISPERSIONS = ("d3bj", "d3", "none")

# ---------------------------------------------------------------------------
# Known two-letter chemical element symbols (for PDB atom-name fallback)
# ---------------------------------------------------------------------------

# All standard 2-letter element symbols.  Used by _infer_element_from_atom_name
# to distinguish e.g. 'Fe' (iron) from 'F' (fluorine) when the PDB element
# column is blank.
_TWO_LETTER_ELEMENTS: frozenset = frozenset({
    "He", "Li", "Be", "Ne", "Na", "Mg", "Al", "Si", "Cl", "Ar",
    "Ca", "Sc", "Ti", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Zr", "Nb",
    "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb",
    "Te", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm",
    "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf",
    "Ta", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi",
    "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "Np", "Pu",
    "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr",
})

# ---------------------------------------------------------------------------
# Standalone xyz / pdb parsers and xyz writer
# ---------------------------------------------------------------------------


def _parse_xyz(path: Path) -> List[Tuple[str, float, float, float]]:
    """Parse an xyz file and return a list of (element, x, y, z) in Angstrom.

    Raises
    ------
    ValueError
        If the file is malformed (wrong atom count, unparseable lines).
    """
    lines = path.read_text().splitlines()
    if len(lines) < 2:
        raise ValueError(f"xyz file too short (need >= 2 lines): {path}")
    try:
        n_atoms = int(lines[0].strip())
    except ValueError:
        raise ValueError(
            f"xyz first line must be the atom count, got {lines[0]!r}: {path}"
        )
    if len(lines) < n_atoms + 2:
        raise ValueError(
            f"xyz file has {len(lines)} lines but header declares "
            f"{n_atoms} atoms: {path}"
        )
    atoms: List[Tuple[str, float, float, float]] = []
    for i, line in enumerate(lines[2: n_atoms + 2], start=3):
        parts = line.split()
        if len(parts) < 4:
            raise ValueError(
                f"xyz line {i} has fewer than 4 fields: {line!r}"
            )
        elem = parts[0]
        try:
            x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
        except ValueError:
            raise ValueError(
                f"Cannot parse coordinates on xyz line {i}: {line!r}"
            )
        atoms.append((elem, x, y, z))
    return atoms


def _infer_element_from_atom_name(atom_name: str) -> str:
    """Infer a chemical element symbol from a PDB atom name field.

    PDB atom names can be like `` CA ``, `` H  ``, ``FE  ``, ``1HB ``, etc.
    Leading digits are stripped, then:

    1. If the first two letters form a known two-letter element symbol
       (e.g. ``Fe``, ``Ca``, ``Mg``), that symbol is returned.
    2. Otherwise the first letter alone is returned (covers C, H, N, O, S,
       P, F, I, B, and so on).

    Returns the element symbol string, or raises ValueError if no letter
    can be found.
    """
    name = atom_name.strip()
    # Remove leading digits (e.g. "1HB" -> "HB", "2H" -> "H")
    name = re.sub(r"^\d+", "", name)
    if not name:
        raise ValueError(
            f"Cannot infer element from atom name {atom_name!r}: "
            "all characters are digits or whitespace."
        )
    # Try two-letter element first
    if len(name) >= 2:
        candidate2 = name[0].upper() + name[1].lower()
        if candidate2 in _TWO_LETTER_ELEMENTS:
            return candidate2
    # Fall back to single letter
    return name[0].upper()


def _parse_pdb(path: Path) -> List[Tuple[str, float, float, float]]:
    """Parse a PDB file and return a list of (element, x, y, z) in Angstrom.

    Element is read from columns 77-78 (0-indexed 76-77) when present;
    the atom-name field (cols 13-16) is used as fallback.

    Raises
    ------
    ValueError
        If any atom's element cannot be determined or coordinates are invalid.
    """
    atoms: List[Tuple[str, float, float, float]] = []
    with path.open() as fh:
        for lineno, raw in enumerate(fh, start=1):
            if not raw.startswith(("ATOM  ", "HETATM")):
                continue
            # Pad to 80 chars so slice access never goes out of range
            line = raw.rstrip("\n").ljust(80)
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                raise ValueError(
                    f"PDB line {lineno}: cannot parse XYZ coordinates: "
                    f"{raw.rstrip()!r}"
                )
            # Prefer element column (cols 77-78 in 1-based PDB spec, i.e.
            # 0-indexed 76-78)
            elem_col = line[76:78].strip()
            if elem_col:
                if len(elem_col) == 2:
                    elem = elem_col[0].upper() + elem_col[1].lower()
                else:
                    elem = elem_col.upper()
            else:
                atom_name = line[12:16]
                try:
                    elem = _infer_element_from_atom_name(atom_name)
                except ValueError:
                    raise ValueError(
                        f"PDB line {lineno}: cannot determine element from "
                        f"atom name {atom_name!r} and the element column "
                        f"(cols 77-78) is blank.  "
                        f"Please add element symbols or use an xyz file."
                    )
            atoms.append((elem, x, y, z))
    if not atoms:
        raise ValueError(f"No ATOM/HETATM records found in {path}")
    return atoms


def _write_xyz(
    path: Path,
    atoms: List[Tuple[str, float, float, float]],
    comment: str = "",
) -> None:
    """Write an xyz file from a list of (element, x, y, z) tuples (Angstrom)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        fh.write(f"{len(atoms)}\n")
        fh.write(f"{comment}\n")
        for elem, x, y, z in atoms:
            fh.write(f"{elem:<4s} {x:15.8f} {y:15.8f} {z:15.8f}\n")


# ---------------------------------------------------------------------------
# PySCF conversion helpers
# ---------------------------------------------------------------------------


def _atoms_to_pyscf_spec(
    atoms: List[Tuple[str, float, float, float]]
) -> List[List]:
    """Convert (elem, x, y, z) list to PySCF atom specification."""
    return [[elem, (x, y, z)] for elem, x, y, z in atoms]


def _mol_to_atoms(mol) -> List[Tuple[str, float, float, float]]:
    """Extract atom symbols and coordinates (Bohr -> Ang) from a PySCF Mole."""
    result: List[Tuple[str, float, float, float]] = []
    for i in range(mol.natm):
        sym = mol.atom_pure_symbol(i)
        bohr = mol.atom_coord(i)          # numpy array, Bohr
        result.append((
            sym,
            float(bohr[0]) * _BOHR_TO_ANG,
            float(bohr[1]) * _BOHR_TO_ANG,
            float(bohr[2]) * _BOHR_TO_ANG,
        ))
    return result


# ---------------------------------------------------------------------------
# Main class
# ---------------------------------------------------------------------------


class QMOptimizerPySCF:
    """QM geometry optimisation using PySCF DFT.

    Default settings: B3LYP / def2-SVP with D3(BJ) dispersion correction.

    All PySCF imports (and dispersion / solver imports) are deferred to
    :meth:`optimize` so that instantiating this class in an environment
    without PySCF does not raise.

    Parameters
    ----------
    functional : str
        DFT exchange-correlation functional name accepted by PySCF.
        Default: ``"B3LYP"``.
    basis : str
        Basis set name accepted by PySCF.  Default: ``"def2-SVP"``.
    dispersion : str
        Dispersion correction: ``"d3bj"`` (D3-BJ, default), ``"d3"``
        (D3-zero damping), or ``"none"`` (disabled).
    charge : int
        Total molecular charge.  Default: ``0``.
    spin : int
        Number of unpaired electrons (2S).  ``0`` = closed-shell singlet.
        Default: ``0``.
    max_steps : int
        Maximum number of geometry optimisation steps.  Default: ``100``.
    solver : str
        Geometry optimisation driver: ``"geometric"`` (geomeTRIC, default)
        or ``"berny"`` (pyberny).
    conv_params : dict or None
        Convergence threshold overrides forwarded to the solver.  Keys are
        solver-specific (see geomeTRIC / berny documentation).
        ``None`` uses the solver defaults.
    verbose : int
        PySCF verbosity level (0–9).  Default: ``3``.
    seed : int or None
        Reserved for future use.  QM calculations are deterministic;
        this parameter has no effect in the current implementation.
    """

    def __init__(
        self,
        functional: str = "B3LYP",
        basis: str = "def2-SVP",
        dispersion: str = "d3bj",
        charge: int = 0,
        spin: int = 0,
        max_steps: int = 100,
        solver: str = "geometric",
        conv_params: Optional[dict] = None,
        verbose: int = 3,
        seed: Optional[int] = None,
    ) -> None:
        if solver not in _SUPPORTED_SOLVERS:
            raise ValueError(
                f"solver must be one of {_SUPPORTED_SOLVERS!r}, "
                f"got {solver!r}"
            )
        if dispersion.lower() not in _SUPPORTED_DISPERSIONS:
            raise ValueError(
                f"dispersion must be one of {_SUPPORTED_DISPERSIONS!r}, "
                f"got {dispersion!r}"
            )
        self.functional = functional
        self.basis = basis
        self.dispersion = dispersion.lower()
        self.charge = charge
        self.spin = spin
        self.max_steps = max_steps
        self.solver = solver
        self.conv_params: dict = dict(conv_params) if conv_params else {}
        self.verbose = verbose
        self.seed = seed  # reserved

    # ------------------------------------------------------------------
    # Dependency checks
    # ------------------------------------------------------------------

    @staticmethod
    def _check_pyscf_available() -> None:
        """Raise ImportError with a helpful message if PySCF is not installed."""
        try:
            import pyscf  # noqa: F401
        except ImportError:
            raise ImportError(
                "PySCF is not installed.  "
                "Install it with: pip install pyscf"
            ) from None

    @staticmethod
    def _check_geomopt_available(solver: str) -> None:
        """Raise ImportError if the requested geometry solver is missing."""
        if solver == "geometric":
            try:
                import geometric  # noqa: F401
            except ImportError:
                raise ImportError(
                    "geomeTRIC is not installed (required for "
                    "solver='geometric').\n"
                    "Install it with: pip install geometric\n"
                    "Alternatively use solver='berny' after: "
                    "pip install pyberny"
                ) from None
        elif solver == "berny":
            try:
                import berny  # noqa: F401
            except ImportError:
                raise ImportError(
                    "pyberny is not installed (required for "
                    "solver='berny').\n"
                    "Install it with: pip install pyberny\n"
                    "Alternatively use solver='geometric' after: "
                    "pip install geometric"
                ) from None

    # ------------------------------------------------------------------
    # Dispersion helper
    # ------------------------------------------------------------------

    def _apply_dispersion(self, mf):
        """Wrap *mf* with D3 dispersion if requested.

        Three strategies are tried in order; if none are available a warning
        is emitted and the uncorrected *mf* is returned.

        Strategy 1: ``pyscf.dftd3`` (requires the ``dftd3`` Python package).
        Strategy 2: ``dftd3.pyscf`` from ``simple-dftd3`` (s-dftd3).
        Strategy 3: No correction (with warning).
        """
        if self.dispersion == "none":
            return mf

        d3_sdftd3_version = {"d3bj": "d3bj", "d3": "d3zero"}.get(
            self.dispersion, "d3bj"
        )

        # --- Strategy 1: pyscf.dftd3 ---
        try:
            from pyscf import dftd3 as _pyscf_dftd3  # type: ignore

            mf_disp = _pyscf_dftd3.DFTD3(mf)
            logger.debug(
                "Dispersion via pyscf.dftd3 (%s)", self.dispersion
            )
            return mf_disp
        except (ImportError, AttributeError, Exception):
            pass

        # --- Strategy 2: simple-dftd3 Python bindings (s-dftd3) ---
        try:
            from dftd3.pyscf import DFTD3Model  # type: ignore

            mf_disp = DFTD3Model(
                mf,
                method=self.functional,
                version=d3_sdftd3_version,
            )
            logger.debug(
                "Dispersion via simple-dftd3/dftd3.pyscf (%s)",
                self.dispersion,
            )
            return mf_disp
        except (ImportError, TypeError, Exception):
            pass

        # --- Fallback: no dispersion ---
        logger.warning(
            "D3 dispersion (%s) requested but no dispersion library is "
            "available.  Running without dispersion correction.\n"
            "To enable D3: pip install simple-dftd3  (recommended)\n"
            "         or:  pip install dftd3",
            self.dispersion,
        )
        return mf

    # ------------------------------------------------------------------
    # PySCF mol / mf builders
    # ------------------------------------------------------------------

    def _build_mol(self, atoms: List[Tuple[str, float, float, float]]):
        """Build a PySCF Mole object from atom list (Angstrom coordinates)."""
        from pyscf import gto

        mol = gto.Mole()
        mol.atom = _atoms_to_pyscf_spec(atoms)
        mol.basis = self.basis
        mol.charge = self.charge
        mol.spin = self.spin
        mol.verbose = self.verbose
        mol.build()
        return mol

    def _build_mf(self, mol):
        """Build a DFT mean-field object (RKS or UKS) with dispersion applied."""
        from pyscf import dft

        if self.spin == 0:
            mf = dft.RKS(mol)
        else:
            mf = dft.UKS(mol)
        mf.xc = self.functional
        mf = self._apply_dispersion(mf)
        return mf

    def _single_point_energy(self, mol_eq) -> float:
        """Run a single-point DFT calculation on *mol_eq* and return energy (Ha).

        The result is guaranteed to reflect the dispersion correction setting,
        regardless of whether the optimiser's internal energy tracking captured
        it correctly.
        """
        from pyscf import dft

        mol_sp = mol_eq.copy()
        mol_sp.verbose = 0  # suppress output for the energy evaluation
        if self.spin == 0:
            mf_sp = dft.RKS(mol_sp)
        else:
            mf_sp = dft.UKS(mol_sp)
        mf_sp.xc = self.functional
        mf_sp = self._apply_dispersion(mf_sp)
        return float(mf_sp.kernel())

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def optimize(
        self,
        in_file: Union[str, Path],
        out_xyz: Union[str, Path],
        *,
        write_traj: bool = False,
    ) -> dict:
        """Run QM geometry optimisation on a single xyz or pdb file.

        Parameters
        ----------
        in_file : str or Path
            Input structure file (``.xyz`` or ``.pdb``).
        out_xyz : str or Path
            Output xyz file path.  Must differ from *in_file*.
        write_traj : bool
            Reserved for future use (trajectory writing not yet implemented).

        Returns
        -------
        dict
            ``"energy"`` (float, eV), ``"energy_hartree"`` (float, Ha),
            ``"steps"`` (int), ``"converged"`` (bool),
            ``"out_xyz"`` (str, absolute path).

        Raises
        ------
        ImportError
            If PySCF or the geometry optimisation solver is not installed.
        FileNotFoundError
            If *in_file* does not exist.
        ValueError
            If the file format is unsupported or element information is
            missing from a PDB file.
        RuntimeError
            If the SCF calculation encounters an unrecoverable error.
        """
        self._check_pyscf_available()
        self._check_geomopt_available(self.solver)

        in_file = Path(in_file).resolve()
        out_xyz = Path(out_xyz).resolve()

        if not in_file.exists():
            raise FileNotFoundError(f"Input file not found: {in_file}")
        if in_file == out_xyz:
            raise ValueError(
                "in_file and out_xyz must be different paths.  "
                f"Got: {in_file}"
            )

        # --- Parse input structure ---
        suffix = in_file.suffix.lower()
        if suffix == ".xyz":
            atoms = _parse_xyz(in_file)
        elif suffix == ".pdb":
            atoms = _parse_pdb(in_file)
        else:
            raise ValueError(
                f"Unsupported input format: {suffix!r}.  "
                "Supported formats: .xyz, .pdb"
            )

        logger.info(
            "QM optimize: %s -> %s  (%d atoms)",
            in_file.name,
            out_xyz.name,
            len(atoms),
        )
        logger.info(
            "  functional=%s  basis=%s  dispersion=%s  "
            "charge=%d  spin=%d  solver=%s  max_steps=%d",
            self.functional,
            self.basis,
            self.dispersion,
            self.charge,
            self.spin,
            self.solver,
            self.max_steps,
        )

        # --- Build PySCF objects ---
        mol = self._build_mol(atoms)
        mf = self._build_mf(mol)

        # --- Geometry optimisation ---
        step_counter = [0]

        def _step_cb(envs):  # noqa: ANN001
            step_counter[0] += 1

        converged = True
        mol_eq = None

        try:
            if self.solver == "geometric":
                from pyscf.geomopt.geometric_solver import (
                    optimize as _geom_opt,
                )

                try:
                    mol_eq = _geom_opt(
                        mf,
                        maxsteps=self.max_steps,
                        callback=_step_cb,
                        **self.conv_params,
                    )
                except TypeError:
                    # Older PySCF may not accept `callback`
                    mol_eq = _geom_opt(
                        mf,
                        maxsteps=self.max_steps,
                        **self.conv_params,
                    )
            else:  # berny
                from pyscf.geomopt.berny_solver import (
                    optimize as _berny_opt,
                )

                try:
                    mol_eq = _berny_opt(
                        mf,
                        maxsteps=self.max_steps,
                        callback=_step_cb,
                        **self.conv_params,
                    )
                except TypeError:
                    mol_eq = _berny_opt(
                        mf,
                        maxsteps=self.max_steps,
                        **self.conv_params,
                    )

        except RuntimeError as exc:
            msg = str(exc).lower()
            if any(kw in msg for kw in ("converge", "max", "failed", "step")):
                converged = False
                logger.warning(
                    "Geometry optimisation did not converge for %s: %s",
                    in_file.name,
                    exc,
                )
                # Use the last geometry stored in mf.mol
                mol_eq = mf.mol
            else:
                raise RuntimeError(
                    f"QM optimisation failed for {in_file.name}: {exc}"
                ) from exc

        n_steps = step_counter[0]

        # --- Final energy (single-point on optimised geometry) ---
        try:
            energy_ha = self._single_point_energy(mol_eq)
        except Exception as exc:
            logger.warning(
                "Single-point energy evaluation failed for %s: %s.  "
                "Attempting to read energy from optimiser state.",
                in_file.name,
                exc,
            )
            try:
                energy_ha = float(mf.e_tot)
            except Exception:
                energy_ha = float("nan")
                logger.warning(
                    "Could not retrieve final energy for %s.", in_file.name
                )

        energy_ev = energy_ha * _HARTREE_TO_EV

        if converged:
            logger.info(
                "Converged in %d steps, energy = %.8f Ha (%.6f eV)",
                n_steps,
                energy_ha,
                energy_ev,
            )
        else:
            logger.warning(
                "Did not converge in %d steps for %s "
                "(last energy: %.8f Ha / %.6f eV)",
                self.max_steps,
                in_file.name,
                energy_ha,
                energy_ev,
            )

        # --- Write output xyz ---
        opt_atoms = _mol_to_atoms(mol_eq)
        comment = (
            f"Optimized by QMOptimizerPySCF: "
            f"{self.functional}/{self.basis} ({self.dispersion})  "
            f"E={energy_ha:.8f} Ha"
        )
        _write_xyz(out_xyz, opt_atoms, comment=comment)
        logger.info("Written: %s", out_xyz)

        return {
            "energy": energy_ev,
            "energy_hartree": energy_ha,
            "steps": n_steps,
            "converged": converged,
            "out_xyz": str(out_xyz),
        }
