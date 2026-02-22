# -*- coding: utf-8 -*-
"""
mdp_parser.py
-------------
Minimal GROMACS .mdp file parser for gro2udf --from-top mode.

Only the subset of parameters needed for COGNAC UDF generation is extracted:
- ref_t          : reference temperature [K]
- tau_t          : Nose-Hoover relaxation time [ps]
- coulombtype    : electrostatics algorithm (PME, Ewald, ...)
- rcoulomb       : real-space Coulomb cutoff [nm]
- constraints    : constraint algorithm (none, h-bonds, ...)
"""
from __future__ import annotations

from typing import Optional


def parse_mdp(mdp_path: str) -> dict:
    """
    Parse a GROMACS .mdp file and return a dict of key → value strings.

    Keys are normalised: leading/trailing whitespace stripped,
    hyphens replaced by underscores (e.g. ``tau-t`` → ``tau_t``).
    Values have inline comments (``; …``) stripped.
    """
    params: dict = {}
    with open(mdp_path, "r") as fh:
        for raw_line in fh:
            line = raw_line.strip()
            if not line or line.startswith(";"):
                continue
            if "=" not in line:
                continue
            key, _, val = line.partition("=")
            key = key.strip().replace("-", "_")
            val = val.split(";")[0].strip()
            params[key] = val
    return params


class MdpParams:
    """
    Typed helper wrapping the dict returned by :func:`parse_mdp`.

    Provides convenient typed accessors with sensible defaults.
    """

    def __init__(self, raw: dict):
        self._raw = raw

    # ------------------------------------------------------------------
    # Accessors
    # ------------------------------------------------------------------

    @property
    def ref_t(self) -> float:
        """Reference temperature [K] (first value if space-separated)."""
        val = self._raw.get("ref_t", "300.0")
        return float(val.split()[0])

    @property
    def tau_t(self) -> float:
        """Nose-Hoover relaxation time [ps] (first value)."""
        val = self._raw.get("tau_t", "0.1")
        return float(val.split()[0])

    @property
    def coulombtype(self) -> str:
        return self._raw.get("coulombtype", "PME")

    @property
    def rcoulomb(self) -> float:
        """Real-space Coulomb cutoff [nm]."""
        return float(self._raw.get("rcoulomb", "0.9"))

    @property
    def constraints(self) -> str:
        return self._raw.get("constraints", "none")

    @property
    def integrator(self) -> str:
        return self._raw.get("integrator", "md")

    def __repr__(self) -> str:  # pragma: no cover
        return (
            f"MdpParams(ref_t={self.ref_t}, tau_t={self.tau_t}, "
            f"coulombtype={self.coulombtype!r}, rcoulomb={self.rcoulomb})"
        )


def load_mdp(mdp_path: Optional[str]) -> Optional[MdpParams]:
    """
    Load *mdp_path* and return :class:`MdpParams`, or ``None`` if path is None.
    """
    if mdp_path is None:
        return None
    return MdpParams(parse_mdp(mdp_path))
