# -*- coding: utf-8 -*-
"""
abmptools.genesis.mmgbsa.analysis
---------------------------------
Parse GENESIS atdyn ``[STEP4] Compute Single Point Energy`` output for
the three MM/GBSA systems, compute ΔG_bind, and emit a CSV + bar plot.

GENESIS atdyn STEP4 layout::

    [STEP4] Compute Single Point Energy for Molecules

                STEP          ENERGY            BOND  ...     SOLVATION
     --------------- --------------- --------------- ---------------
                   0      -8855.2124        147.8618  ...    -2029.5026

The 4th line below the ``[STEP4]`` header has 9 numeric tokens (POC
``4_analyse.py`` regex). We extract the **2nd token (ENERGY)** and the
**9th token (SOLVATION)**.

GENESIS ``ENERGY`` column convention (per GENESIS doc 05_Energy.rst:564):

    U = U_FF + ΔG_elec + ΔG_np

That is, the ``ENERGY`` column is **already the total** ``U_FF + SOLVATION``
(the solvation free energy is folded into the potential energy
function as an effective term). The ``SOLVATION`` column is the
``ΔG_elec + ΔG_np`` portion **as a sub-component of ENERGY**, not
something to add on top.

Binding free energy (standard MM/GBSA, equivalent forms):

    ΔG_bind = E_complex - E_ligand - E_receptor              (form A: total energy difference)
            = (egas_c + S_c) - (egas_l + S_l + egas_r + S_r) (form B: POC ``4_analyse.py``)

where ``egas = ENERGY - SOLVATION`` (gas-phase MM contribution).
Forms A and B are **algebraically identical** because
``egas + S = ENERGY``.

Units: GENESIS-native ``kcal/mol`` (matches AMBER force field convention).
"""
from __future__ import annotations

import csv
import logging
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

#: Number of lines below the ``[STEP4]`` header where the data row sits.
STEP4_DATA_OFFSET = 4

#: POC regex: signed floats and bare ``0``.
TOKEN_RE = re.compile(r"-?\d+\.\d+|(?<!\d)0(?!\d)")

#: Index of ENERGY column in the STEP4 data row (after STEP=0).
ENERGY_IDX = 1

#: Index of SOLVATION column. POC line shows 9 numeric tokens
#: ``[STEP, ENERGY, BOND, ANGLE, DIHEDRAL, IMPROPER, VDWAALS, ELECT, SOLVATION]``.
SOLVATION_IDX = 8

SYSTEM_NAMES: List[str] = ["complex", "ligand", "receptor"]


# ---------------------------------------------------------------------------
# Result containers
# ---------------------------------------------------------------------------

@dataclass
class SystemEnergies:
    """``(energy, solvation)`` pair for one system."""
    energy: float
    solvation: float


@dataclass
class TargetEnergies:
    """Aggregated energies for one target's three systems + ΔG_bind."""
    name: str
    complex_e: float
    complex_s: float
    ligand_e: float
    ligand_s: float
    receptor_e: float
    receptor_s: float
    dg_bind: float                 # kcal/mol


# ---------------------------------------------------------------------------
# STEP4 log parsing
# ---------------------------------------------------------------------------

def parse_step4_log(log_text: str) -> SystemEnergies:
    """Extract ``(energy, solvation)`` from a GENESIS atdyn log.

    Locates the line containing ``[STEP4] Compute Single Point Energy``,
    skips :data:`STEP4_DATA_OFFSET` lines forward, then extracts the
    2nd and 9th numeric tokens from the data row (POC convention).

    Raises
    ------
    ValueError
        If the STEP4 header is missing or the data row has too few tokens.
    """
    lines = log_text.splitlines()
    for i, line in enumerate(lines):
        if "[STEP4] Compute Single Point Energy" in line:
            data_idx = i + STEP4_DATA_OFFSET
            if data_idx >= len(lines):
                raise ValueError(
                    f"STEP4 found at line {i+1} but log has only "
                    f"{len(lines)} lines; data row at {data_idx+1} missing."
                )
            data_line = lines[data_idx].strip()
            tokens = TOKEN_RE.findall(data_line)
            if len(tokens) < SOLVATION_IDX + 1:
                raise ValueError(
                    f"STEP4 data row has {len(tokens)} tokens "
                    f"(need >= {SOLVATION_IDX + 1}): {data_line!r}"
                )
            return SystemEnergies(
                energy=float(tokens[ENERGY_IDX]),
                solvation=float(tokens[SOLVATION_IDX]),
            )
    raise ValueError(
        "No '[STEP4] Compute Single Point Energy' line in log."
    )


def parse_step4_log_file(log_path: Path) -> SystemEnergies:
    """Read a log file and parse it via :func:`parse_step4_log`."""
    return parse_step4_log(Path(log_path).read_text())


# ---------------------------------------------------------------------------
# ΔG_bind aggregation
# ---------------------------------------------------------------------------

def compute_dg_bind(
    complex_e: float, complex_s: float,
    ligand_e: float, ligand_s: float,
    receptor_e: float, receptor_s: float,
) -> float:
    """Standard MM/GBSA binding free energy.

    ``ΔG_bind = E_complex - E_ligand - E_receptor``  (kcal/mol)

    GENESIS' ``ENERGY`` column is the total potential ``U = U_FF +
    ΔG_solv`` (per doc 05_Energy.rst:564), so the difference of
    ``ENERGY`` columns alone is the full MM/GBSA ΔG_bind. The
    ``*_s`` (SOLVATION) arguments are accepted for API symmetry and
    component reporting but **must not** be added on top of ``E``
    (that would double-count solvation).

    Equivalent decomposition (POC ``4_analyse.py``):
    ``(egas + S)_c - (egas + S)_l - (egas + S)_r``
    where ``egas = ENERGY - SOLVATION``. Algebraically identical
    because ``egas + S = ENERGY``.
    """
    # Sanity check: silence unused-arg lint while documenting that
    # SOLVATION components are intentionally not added (already in E).
    _ = complex_s, ligand_s, receptor_s
    return complex_e - ligand_e - receptor_e


def compute_dg_components(
    complex_e: float, complex_s: float,
    ligand_e: float, ligand_s: float,
    receptor_e: float, receptor_s: float,
) -> dict:
    """Return ΔG_bind decomposed into MM (gas-phase) and solvation parts.

    Useful for reporting the relative contribution of GBSA solvation
    to binding without changing the overall sum.

    Returns
    -------
    dict
        ``{"dg_mm": float, "dg_solv": float, "dg_bind": float}``
        where ``dg_bind == dg_mm + dg_solv``.
    """
    egas_c = complex_e - complex_s
    egas_l = ligand_e - ligand_s
    egas_r = receptor_e - receptor_s
    dg_mm = egas_c - egas_l - egas_r
    dg_solv = complex_s - ligand_s - receptor_s
    return {
        "dg_mm": dg_mm,
        "dg_solv": dg_solv,
        "dg_bind": dg_mm + dg_solv,
    }


def aggregate_target(target_dir: Path, name: Optional[str] = None) -> TargetEnergies:
    """Read ``complex.log`` / ``ligand.log`` / ``receptor.log`` from
    *target_dir* and return one :class:`TargetEnergies`.

    Parameters
    ----------
    target_dir
        Per-target directory containing the three log files.
    name
        Display name (default = ``target_dir.name``).
    """
    target_dir = Path(target_dir)
    if name is None:
        name = target_dir.name

    energies: Dict[str, SystemEnergies] = {}
    for sys_name in SYSTEM_NAMES:
        log_path = target_dir / f"{sys_name}.log"
        if not log_path.is_file():
            raise FileNotFoundError(
                f"Missing {sys_name}.log in {target_dir}"
            )
        energies[sys_name] = parse_step4_log_file(log_path)

    dg = compute_dg_bind(
        complex_e=energies["complex"].energy,
        complex_s=energies["complex"].solvation,
        ligand_e=energies["ligand"].energy,
        ligand_s=energies["ligand"].solvation,
        receptor_e=energies["receptor"].energy,
        receptor_s=energies["receptor"].solvation,
    )
    return TargetEnergies(
        name=name,
        complex_e=energies["complex"].energy,
        complex_s=energies["complex"].solvation,
        ligand_e=energies["ligand"].energy,
        ligand_s=energies["ligand"].solvation,
        receptor_e=energies["receptor"].energy,
        receptor_s=energies["receptor"].solvation,
        dg_bind=dg,
    )


def aggregate_all(target_dirs: List[Path]) -> List[TargetEnergies]:
    """Run :func:`aggregate_target` for each entry."""
    return [aggregate_target(td) for td in target_dirs]


# ---------------------------------------------------------------------------
# CSV / plot output
# ---------------------------------------------------------------------------

CSV_HEADER: List[str] = [
    "Target",
    "Complex_Energy", "Complex_Solvation",
    "Ligand_Energy", "Ligand_Solvation",
    "Receptor_Energy", "Receptor_Solvation",
    "dG_bind_kcal_per_mol",
]


def write_results_csv(
    results: List[TargetEnergies],
    csv_path: Path,
) -> Path:
    """Write aggregated results to CSV (one row per target)."""
    csv_path = Path(csv_path)
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(CSV_HEADER)
        for r in results:
            w.writerow([
                r.name,
                f"{r.complex_e:.4f}", f"{r.complex_s:.4f}",
                f"{r.ligand_e:.4f}", f"{r.ligand_s:.4f}",
                f"{r.receptor_e:.4f}", f"{r.receptor_s:.4f}",
                f"{r.dg_bind:.4f}",
            ])
    logger.info("Wrote %s", csv_path)
    return csv_path


def plot_dg_bind(
    results: List[TargetEnergies],
    out_png: Path,
) -> Path:
    """Render a bar plot of ΔG_bind per target.

    matplotlib is loaded lazily (``[mmgbsa]`` extra). Uses the ``Agg``
    backend so the function works on headless / CI machines.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if not results:
        raise ValueError("plot_dg_bind: empty results list")

    out_png = Path(out_png)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    names = [r.name for r in results]
    dgs = [r.dg_bind for r in results]

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.bar(names, dgs, color="steelblue")
    ax.set_xlabel("Target")
    ax.set_ylabel("ΔG_bind (kcal/mol)")
    ax.set_title("MM/GBSA single-point ΔG_bind")
    ax.tick_params(axis="x", rotation=45)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
    logger.info("Wrote %s", out_png)
    return out_png


# ---------------------------------------------------------------------------
# Top-level entry
# ---------------------------------------------------------------------------

@dataclass
class AnalyzeResult:
    """Outputs of :func:`analyze`."""
    csv_path: Path
    png_path: Path
    targets: List[TargetEnergies]


def analyze(
    target_dirs: List[Path],
    out_csv: Path,
    out_png: Path,
) -> AnalyzeResult:
    """Aggregate logs across *target_dirs*, write CSV + bar plot."""
    results = aggregate_all([Path(p) for p in target_dirs])
    csv_path = write_results_csv(results, out_csv)
    png_path = plot_dg_bind(results, out_png)
    return AnalyzeResult(
        csv_path=csv_path, png_path=png_path, targets=results
    )


__all__ = [
    "CSV_HEADER",
    "ENERGY_IDX",
    "SOLVATION_IDX",
    "SYSTEM_NAMES",
    "STEP4_DATA_OFFSET",
    "AnalyzeResult",
    "SystemEnergies",
    "TargetEnergies",
    "aggregate_all",
    "aggregate_target",
    "analyze",
    "compute_dg_bind",
    "compute_dg_components",
    "parse_step4_log",
    "parse_step4_log_file",
    "plot_dg_bind",
    "write_results_csv",
]
