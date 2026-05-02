# -*- coding: utf-8 -*-
"""
abmptools.membrane.pmf
----------------------
PMF analysis: WHAM (gmx wham) and / or PyMBAR (optional).

Inputs (from MD output)
-----------------------
- ``windows/win_NNN/window.tpr`` — used by gmx wham for k / x0
- ``windows/win_NNN/pullx.xvg`` — pull-coord values per window
- ``windows/win_NNN/pullf.xvg`` — pull forces per window

Outputs (to ``analysis/``)
--------------------------
- ``pmf.xvg``                  — PMF(z) [kJ/mol]
- ``histo.xvg``                — per-window histograms
- ``bsResult.xvg`` (optional)  — bootstrap error estimate
- ``tpr.dat`` / ``pullx.dat``  — gmx-wham input lists (kept for audit)
"""
from __future__ import annotations

import logging
import os
import subprocess
from pathlib import Path
from typing import Any, Dict, List, Optional

from .models import MembraneConfig
from .bilayer import _resolve_amberhome

logger = logging.getLogger(__name__)


def run_wham(
    *, windows_dir: str, analysis_dir: str, config: MembraneConfig,
    bootstrap_n: int = 0, temperature_K: Optional[float] = None,
) -> Dict[str, Any]:
    """Run ``gmx wham`` on completed window output.

    Parameters
    ----------
    windows_dir : str
        Path to the ``windows/`` directory containing ``win_NNN/``.
    analysis_dir : str
        Output directory for PMF files.
    config : MembraneConfig
        Build configuration (used for temperature default).
    bootstrap_n : int
        Number of bootstrap samples for error estimation. 0 disables.
    temperature_K : float, optional
        Temperature override; defaults to
        ``config.equilibration.temperature_K``.

    Returns
    -------
    dict
        Keys ``"pmf"``, ``"histo"``, ``"tpr_list"``, ``"pullx_list"``,
        and (when ``bootstrap_n > 0``) ``"bs_result"``.
    """
    wd = Path(windows_dir).resolve()
    ad = Path(analysis_dir).resolve()
    ad.mkdir(parents=True, exist_ok=True)

    if temperature_K is None:
        temperature_K = config.equilibration.temperature_K

    tpr_files: List[str] = []
    pullx_files: List[str] = []
    for d in sorted(wd.glob("win_*")):
        tpr = d / "window.tpr"
        pullx = d / "pullx.xvg"
        if not tpr.is_file():
            logger.warning("missing tpr: %s — skipping", tpr)
            continue
        if not pullx.is_file():
            logger.warning("missing pullx: %s — skipping", pullx)
            continue
        tpr_files.append(str(tpr))
        pullx_files.append(str(pullx))

    if not tpr_files:
        raise RuntimeError(
            f"no window.tpr / pullx.xvg pairs found under {wd}; "
            f"did the per-window MD complete?"
        )

    tpr_dat = ad / "tpr.dat"
    pullx_dat = ad / "pullx.dat"
    tpr_dat.write_text("\n".join(tpr_files) + "\n")
    pullx_dat.write_text("\n".join(pullx_files) + "\n")

    pmf_path = str(ad / "pmf.xvg")
    histo_path = str(ad / "histo.xvg")
    cmd: List[str] = [
        config.gmx_path, "wham",
        "-it", str(tpr_dat),
        "-ix", str(pullx_dat),
        "-o",  pmf_path,
        "-hist", histo_path,
        "-temp", f"{temperature_K:.2f}",
        "-unit", "kJ",
    ]
    bs_path: Optional[str] = None
    if bootstrap_n > 0:
        bs_path = str(ad / "bsResult.xvg")
        cmd += [
            "-bsnum", str(bootstrap_n),
            "-bsres", bs_path,
        ]

    env = os.environ.copy()
    amberhome = _resolve_amberhome(config.gmx_path)
    if amberhome:
        env.setdefault("AMBERHOME", amberhome)

    logger.info("running gmx wham: %s", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True, env=env)
    log_path = ad / "gmx_wham.log"
    log_path.write_text(
        f"=== argv ===\n{' '.join(cmd)}\n\n"
        f"=== stdout ===\n{result.stdout}\n"
        f"=== stderr ===\n{result.stderr}\n"
        f"=== returncode ===\n{result.returncode}\n"
    )
    if result.returncode != 0 or not Path(pmf_path).is_file():
        raise RuntimeError(
            f"gmx wham failed (rc={result.returncode}). See {log_path}."
        )

    out: Dict[str, Any] = {
        "pmf":        pmf_path,
        "histo":      histo_path,
        "tpr_list":   str(tpr_dat),
        "pullx_list": str(pullx_dat),
    }
    if bs_path is not None:
        out["bs_result"] = bs_path
    return out


def run_pymbar(
    *, windows_dir: str, analysis_dir: str, config: MembraneConfig,
) -> Dict[str, Any]:
    """Optional MBAR analysis (cross-check vs WHAM).

    Not implemented yet (Phase B). Will use ``pymbar.MBAR`` from the
    PyMBAR package (MIT licensed) once added as an optional dep.
    """
    raise NotImplementedError(
        "run_pymbar: optional MBAR analysis is not yet implemented. "
        "Use run_wham for now; pymbar is planned as Phase B follow-up."
    )


# ---------------------------------------------------------------------------
# CLI entry point (used by run.sh)
# ---------------------------------------------------------------------------

def _main() -> None:
    """``python -m abmptools.membrane.pmf`` — run gmx wham."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Run gmx wham on completed umbrella-sampling windows.",
    )
    parser.add_argument("--windows-dir", required=True)
    parser.add_argument("--analysis-dir", required=True)
    parser.add_argument("--config", required=True,
                        help="Path to MembraneConfig JSON.")
    parser.add_argument("--bootstrap-n", type=int, default=0,
                        help="Number of bootstrap samples (0 = disable).")
    parser.add_argument("--temperature", type=float, default=None,
                        help="Override the temperature [K].")
    parser.add_argument("--gmx-path", default=None,
                        help="Override config.gmx_path (must match the gmx "
                             "that produced the per-window .tpr files).")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(name)s | %(message)s")
    cfg = MembraneConfig.from_json(args.config)
    if args.gmx_path is not None:
        cfg.gmx_path = args.gmx_path
    out = run_wham(
        windows_dir=args.windows_dir,
        analysis_dir=args.analysis_dir,
        config=cfg,
        bootstrap_n=args.bootstrap_n,
        temperature_K=args.temperature,
    )
    for k, v in out.items():
        print(f"  {k:12s} {v}")


if __name__ == "__main__":
    _main()
