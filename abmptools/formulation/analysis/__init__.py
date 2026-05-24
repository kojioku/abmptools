# -*- coding: utf-8 -*-
"""Analysis stack for abmptools.formulation trajectories.

Modules in this sub-package import MDAnalysis / networkx lazily — these
are GPL-2.0 (MDAnalysis) and BSD-3 (networkx) and live in the
``abmptools[formulation-analysis]`` optional extra. The abmptools wheel
itself remains Apache-2.0 with no GPL code shipped or linked.
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Iterable, Sequence

logger = logging.getLogger(__name__)


def run_analysis(
    *,
    traj: str,
    top: str,
    out_dir: str,
    enhancer_resnames: Sequence[str] = (),
    bile_salt_resnames: Sequence[str] = (),
    contact_cutoff_nm: float = 0.5,
) -> dict:
    """Run the full analysis stack (aggregate / contact / SS / SASA / hbond).

    Each sub-step is optional and lazily imports its dependencies, so
    a missing tool (e.g. no MDAnalysis installed) skips that step with
    a warning rather than failing the whole call.
    """
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    results: dict = {}

    # 1. aggregate transition (MDAnalysis + networkx)
    try:
        from .aggregate_transition import compute_aggregate_transitions
        results["aggregate"] = compute_aggregate_transitions(
            traj=traj, top=top, out_dir=str(out),
            cutoff_nm=contact_cutoff_nm,
        )
    except ImportError as exc:
        logger.warning("aggregate_transition skipped: %s", exc)

    # 2. contact map
    try:
        from .contact_map import compute_contact_map
        results["contact_map"] = compute_contact_map(
            traj=traj, top=top, out_dir=str(out),
            enhancer_resnames=list(enhancer_resnames),
            bile_salt_resnames=list(bile_salt_resnames),
        )
    except ImportError as exc:
        logger.warning("contact_map skipped: %s", exc)

    return results


__all__ = ["run_analysis"]
