# -*- coding: utf-8 -*-
"""
abmptools.genesis.grest.analysis
--------------------------------
Post-MD analysis for the gREST_SSCR pipeline.

Four entry points:

- :func:`run_remd_convert`           -- run GENESIS ``remd_convert``
                                         to extract the lowest-T
                                         (parameter-sorted) trajectory.
- :func:`plot_replica_transition`    -- random-walk diagram of each
                                         replica through parameter
                                         space (cf. POC image3 / 19).
- :func:`plot_acceptance_ratio`      -- exchange-pair acceptance ratio
                                         vs MD step (rolling-window
                                         mean, default window=100).
- :func:`compute_distance_pmf`       -- 1D distance PMF
                                         ``-kT log P(r)`` from a
                                         lowest-T trajectory and a
                                         pair of atom masks.

Optional dependencies:

- ``matplotlib`` -- plot output. Imported inside each ``plot_*``
                    function so the module loads cleanly on a
                    matplotlib-free install.
- ``mdtraj``     -- distance time-series in pure Python. If
                    unavailable, :func:`compute_distance_pmf`
                    falls back to subprocess'd cpptraj.

GENESIS ``.rem`` file format (one line per exchange):
``# step replica_id parameter_id`` (POC convention). The parser
tolerates leading ``#`` comments and blank lines.
"""
from __future__ import annotations

import logging
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from ._subprocess import CommandError, run_command, write_text
from .models import GrestBuildConfig

logger = logging.getLogger(__name__)


KB = 8.314462618e-3  # kJ/(mol K)


# ---------------------------------------------------------------------------
# remd_convert orchestration
# ---------------------------------------------------------------------------

def run_remd_convert(
    cfg: GrestBuildConfig,
    inp_path: Path,
    workdir: Path,
) -> Path:
    """Invoke ``remd_convert`` on a previously rendered ``.inp``.

    Parameters
    ----------
    cfg
        Validated :class:`GrestBuildConfig` (provides
        ``remd_convert_path``).
    inp_path
        Path to ``step5_remd_convert.inp``.
    workdir
        Working directory (typically ``output_dir/build``). The
        ``param.pdb`` / ``param{}.dcd`` outputs land here.

    Returns
    -------
    Path
        Path to ``param1.dcd`` (lowest-T sorted trajectory).
    """
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    logger.info("remd_convert: %s", inp_path)
    run_command(
        [cfg.remd_convert_path, str(inp_path)],
        cwd=workdir,
        capture=True,
    )
    expected = workdir / "param1.dcd"
    if not expected.exists():
        raise CommandError(
            cmd=f"{cfg.remd_convert_path} {inp_path}",
            returncode=1,
            stderr=(
                f"remd_convert finished but {expected} is missing. "
                "Check the input file's [OUTPUT] section."
            ),
        )
    return expected


# ---------------------------------------------------------------------------
# .rem file parsing + replica transition plot
# ---------------------------------------------------------------------------

def parse_rem_file(path: Path) -> np.ndarray:
    """Parse a single GENESIS ``.rem`` file.

    Format (POC convention): each non-comment line has at minimum
    ``step replica_id parameter_id`` (whitespace-separated, integers).
    Returns an ``(n, 3)`` ndarray of int.
    """
    rows: List[Tuple[int, int, int]] = []
    text = Path(path).read_text()
    for line in text.splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        toks = s.split()
        if len(toks) < 3:
            continue
        try:
            step = int(toks[0])
            replica = int(toks[1])
            param = int(toks[2])
        except ValueError:
            continue
        rows.append((step, replica, param))
    if not rows:
        return np.empty((0, 3), dtype=int)
    return np.asarray(rows, dtype=int)


def parse_rem_files(rem_files: List[Path]) -> Dict[int, np.ndarray]:
    """Parse a list of per-replica ``.rem`` files into a dict.

    The dict key is the replica id read from each file; the value is
    a ``(n_step, 2)`` array of ``[step, parameter_id]``. If a file
    contains a single replica id consistently, the key is that id;
    otherwise the file is grouped by the most common replica id.
    """
    out: Dict[int, np.ndarray] = {}
    for f in rem_files:
        arr = parse_rem_file(f)
        if arr.size == 0:
            continue
        # Group by replica id; in the typical "one rep per file" layout
        # there's only one bucket, but we tolerate mixed.
        for rep_id in np.unique(arr[:, 1]):
            mask = arr[:, 1] == rep_id
            out.setdefault(int(rep_id), arr[mask][:, [0, 2]])
    return out


def plot_replica_transition(
    rem_files: List[Path],
    out_png: Path,
    out_csv: Optional[Path] = None,
) -> Dict[int, np.ndarray]:
    """Plot replica random-walks through parameter space.

    One line per replica, x = MD step, y = parameter id (1..N).
    """
    import matplotlib  # noqa: deferred import (optional dep)
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    series = parse_rem_files(rem_files)
    if not series:
        raise ValueError(f"No usable data in {len(rem_files)} .rem files.")

    fig, ax = plt.subplots(figsize=(8, 4))
    for rep_id in sorted(series):
        data = series[rep_id]
        ax.plot(data[:, 0], data[:, 1], lw=0.6, label=f"rep{rep_id}")
    ax.set_xlabel("MD step")
    ax.set_ylabel("parameter id")
    ax.set_title("Replica transition through parameter space")
    ax.legend(ncol=4, fontsize=8, loc="upper right")
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
    logger.info("Wrote %s", out_png)

    if out_csv is not None:
        rows = []
        for rep_id in sorted(series):
            data = series[rep_id]
            for step, param in data:
                rows.append((step, rep_id, param))
        np.savetxt(
            out_csv,
            np.asarray(rows, dtype=int),
            fmt="%d",
            delimiter=",",
            header="step,replica,parameter",
            comments="",
        )
        logger.info("Wrote %s", out_csv)

    return series


# ---------------------------------------------------------------------------
# Acceptance ratio plot (from spdyn log)
# ---------------------------------------------------------------------------

@dataclass
class ExchangeEvent:
    step: int
    pattern: int       # ExchangePattern: 1, 2, ...
    replica_a: int
    replica_b: int     # 0 -> no exchange attempted (boundary)
    accepted: bool
    accept_n: int
    accept_total: int
    T_before: float
    T_after: float


_REMD_HEADER_RE = re.compile(
    r"REMD>\s+Step:\s+(\d+)\s+Dimension:\s+\d+\s+ExchangePattern:\s+(\d+)"
)
_REMD_LINE_RE = re.compile(
    r"^\s*(\d+)\s+(\d+)\s*>\s*(\d+)\s+([NA])\s+(\d+)\s*/\s*(\d+)"
    r"\s+([\d.+-Ee]+)\s+([\d.+-Ee]+)"
)


def parse_acceptance_log(log_path: Path) -> List[ExchangeEvent]:
    """Parse ``REMD>`` events from a spdyn stdout log.

    POC format:
    ``REMD> Step: 80  Dimension: 1  ExchangePattern: 2``
    ``  1  1 > 0  N  0 / 0  300.000  300.000``
    """
    events: List[ExchangeEvent] = []
    current_step = None
    current_pattern = None
    for line in Path(log_path).read_text().splitlines():
        m_h = _REMD_HEADER_RE.search(line)
        if m_h:
            current_step = int(m_h.group(1))
            current_pattern = int(m_h.group(2))
            continue
        if current_step is None:
            continue
        m_l = _REMD_LINE_RE.match(line)
        if not m_l:
            continue
        rep_a = int(m_l.group(2))
        rep_b = int(m_l.group(3))
        accepted = m_l.group(4) == "A"
        accept_n = int(m_l.group(5))
        accept_total = int(m_l.group(6))
        events.append(ExchangeEvent(
            step=current_step,
            pattern=current_pattern,
            replica_a=rep_a,
            replica_b=rep_b,
            accepted=accepted,
            accept_n=accept_n,
            accept_total=accept_total,
            T_before=float(m_l.group(7)),
            T_after=float(m_l.group(8)),
        ))
    return events


def plot_acceptance_ratio(
    log_path: Path,
    out_png: Path,
    burn_in: int = 100,
    out_csv: Optional[Path] = None,
) -> Dict[Tuple[int, int], List[Tuple[int, float]]]:
    """Plot acceptance ratio per replica pair vs MD step.

    Excludes the first ``burn_in`` exchange events (POC notes
    inheritance from prior runs causes spurious early values).
    Returns a dict ``{(rep_a, rep_b): [(step, cumulative_ratio), ...]}``.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    events = parse_acceptance_log(log_path)
    events = events[burn_in:]
    if not events:
        raise ValueError(
            f"No exchange events after burn_in={burn_in} in {log_path}"
        )

    series: Dict[Tuple[int, int], List[Tuple[int, float]]] = {}
    for ev in events:
        if ev.replica_b == 0:
            continue  # boundary "no exchange" entries skipped
        key = (min(ev.replica_a, ev.replica_b), max(ev.replica_a, ev.replica_b))
        ratio = ev.accept_n / ev.accept_total if ev.accept_total else 0.0
        series.setdefault(key, []).append((ev.step, ratio))

    fig, ax = plt.subplots(figsize=(8, 4))
    for (a, b), pts in sorted(series.items()):
        steps, ratios = zip(*pts)
        ax.plot(steps, ratios, lw=1.0, label=f"{a}<->{b}")
    ax.set_xlabel("MD step")
    ax.set_ylabel("cumulative acceptance ratio")
    ax.set_title("Replica-pair exchange acceptance")
    ax.set_ylim(0.0, 1.0)
    ax.legend(ncol=2, fontsize=8, loc="upper right")
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
    logger.info("Wrote %s", out_png)

    if out_csv is not None:
        rows = []
        for (a, b), pts in sorted(series.items()):
            for step, ratio in pts:
                rows.append((step, a, b, ratio))
        np.savetxt(
            out_csv,
            np.asarray(rows, dtype=float),
            fmt=("%d", "%d", "%d", "%.6f"),
            delimiter=",",
            header="step,rep_a,rep_b,ratio",
            comments="",
        )
        logger.info("Wrote %s", out_csv)

    # Warn on poor pairs.
    for key, pts in series.items():
        last_ratio = pts[-1][1]
        if last_ratio < 0.1:
            logger.warning(
                "Replica pair %s acceptance %.3f < 0.1; "
                "consider re-tuning the temperature ladder.",
                key, last_ratio,
            )

    return series


# ---------------------------------------------------------------------------
# 1D distance PMF
# ---------------------------------------------------------------------------

def _genesis_to_cpptraj_mask(mask: str) -> str:
    """Translate a GENESIS-style mask (``rno:96 and an:NZ``) to cpptraj
    syntax (``:96@NZ``)."""
    m_res = re.search(r"rno:(\S+)", mask)
    m_atom = re.search(r"an:(\S+)", mask)
    if not m_res:
        raise ValueError(
            f"Mask {mask!r} must include 'rno:<resno>' for cpptraj translation."
        )
    res = m_res.group(1).strip(",;")
    atom = m_atom.group(1).strip(",;") if m_atom else None
    if atom:
        return f":{res}@{atom}"
    return f":{res}"


def compute_distance_timeseries_cpptraj(
    prmtop: Path,
    dcd: Path,
    mask1: str,
    mask2: str,
    workdir: Path,
    cpptraj_path: str = "cpptraj",
) -> np.ndarray:
    """Return a (n_frames,) ndarray of pairwise distances via cpptraj."""
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    cpptraj_in = workdir / "distance.cpptraj"
    out_dat = workdir / "distance.dat"
    script = (
        f"parm {prmtop}\n"
        f"trajin {dcd}\n"
        f"distance d1 {_genesis_to_cpptraj_mask(mask1)} "
        f"{_genesis_to_cpptraj_mask(mask2)} out {out_dat}\n"
        f"go\n"
        f"quit\n"
    )
    write_text(cpptraj_in, script)
    run_command([cpptraj_path, "-i", str(cpptraj_in)], cwd=workdir)
    # cpptraj emits a 2-column dat file: frame, distance.
    arr = np.loadtxt(out_dat, comments="#")
    if arr.ndim == 1:
        arr = arr.reshape(1, -1)
    return arr[:, 1]


def compute_distance_pmf(
    distances: np.ndarray,
    T_K: float,
    out_png: Path,
    out_xvg: Path,
    n_bins: int = 60,
    range_A: Optional[Tuple[float, float]] = None,
) -> np.ndarray:
    """1D PMF ``-kT log P(r)`` from a precomputed distance time-series.

    Returns ``(n_bins, 2)`` array of ``[r_centre_A, pmf_kJ_mol]``.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if distances.size == 0:
        raise ValueError("compute_distance_pmf got empty distance array.")
    counts, edges = np.histogram(distances, bins=n_bins, range=range_A)
    centres = 0.5 * (edges[:-1] + edges[1:])
    p = counts / counts.sum()
    # Avoid log(0): replace zeros with the smallest positive bin.
    p_safe = np.where(p > 0, p, p[p > 0].min() if (p > 0).any() else 1e-12)
    pmf = -KB * T_K * np.log(p_safe)
    pmf -= pmf.min()  # shift so min == 0
    out = np.column_stack([centres, pmf])
    np.savetxt(out_xvg, out, fmt="%.6f",
               header="r_A pmf_kJ_per_mol", comments="# ")
    logger.info("Wrote %s", out_xvg)

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(centres, pmf, "-", lw=1.5)
    ax.set_xlabel("r (Å)")
    ax.set_ylabel("PMF (kJ/mol)")
    ax.set_title(f"1D distance PMF (T={T_K:.1f} K)")
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
    logger.info("Wrote %s", out_png)
    return out


__all__ = [
    "ExchangeEvent",
    "compute_distance_pmf",
    "compute_distance_timeseries_cpptraj",
    "parse_acceptance_log",
    "parse_rem_file",
    "parse_rem_files",
    "plot_acceptance_ratio",
    "plot_replica_transition",
    "run_remd_convert",
]
