# -*- coding: utf-8 -*-
"""Analysis stack for abmptools.formulation trajectories.

ユーザーが MD 後の 1 コマンドで論文 Hossain 2023 準拠の解析を回せる
ワークフロー。 個別解析の細粒度関数も export している。

論文 Fig 対応:
- **Fig 1c/2b** (max aggregate size + size distribution) — :mod:`.aggregate_transition`
- **Fig 1b** (% aggregated vs time) — :mod:`.aggregate_transition` の同 CSV
- **Fig 4** (per-residue contacts: Peptide↔Enhancer / BileSalt / Peptide) —
  :mod:`.contact_map.compute_per_residue_contacts`
- **Fig 5** (secondary structure DSSP) — :mod:`.secondary_structure`
- **plots** — :mod:`.plots` (matplotlib)

Modules in this sub-package import MDAnalysis / networkx lazily — these
are GPL-2.0 (MDAnalysis) and BSD-3 (networkx) and live in the
``abmptools[formulation-analysis]`` optional extra.
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional, Sequence

logger = logging.getLogger(__name__)


def run_analysis(
    *,
    traj: str,
    top: str,
    out_dir: str,
    n_peptides: int,
    enhancer_resnames: Sequence[str] = (),
    bile_salt_resnames: Sequence[str] = (),
    contact_cutoff_nm: float = 0.5,
    target_n_frames: int = 100,
    run_aggregate: bool = True,
    run_contacts: bool = True,
    run_dssp: bool = False,
    run_plots: bool = True,
    gmx: str = "gmx",
    ndx: Optional[str] = None,
    peptide_selector: str = "protein",
) -> dict:
    """End-to-end Hossain 2023 風解析を 1 コマンドで回す orchestrator。

    Parameters
    ----------
    traj
        Production trajectory (.xtc / .trr)。 **raw `prod.xtc` を推奨**
        (nojump xtc は cluster 解析に不適)。
    top
        Topology (.gro 推奨 — tpr v138 等の version mismatch を避ける)。
    out_dir
        全解析の出力 dir。 各 sub-module が CSV/PNG/JSON を書き込む。
    n_peptides
        peptide copy 数。 GROMACS の resid 1-N リセットに対応する atom-index
        分割のため必須。
    enhancer_resnames
        permeation enhancer の resname list (例: ``["CPN", "CPC"]`` for caprate)。
    bile_salt_resnames
        bile salt の resname list (例: ``["TCH"]`` for taurocholate)。
    contact_cutoff_nm
        論文準拠 0.5 nm。
    target_n_frames
        解析対象 frame 数 (stride は自動計算)。 100 → 1 ns 分解能 (100 ns run 時)。
    run_aggregate / run_contacts / run_dssp / run_plots
        個別 on/off。 DSSP は gmx dssp (GROMACS >= 2023) が必要。
    gmx, ndx
        DSSP で gmx subprocess を呼ぶ際の path / index。

    Returns
    -------
    dict
        各 sub-module の戻り値を ``aggregate``/``contacts``/``dssp``/``plots``
        キーで集約。 一部失敗しても他は実行される。
    """
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    results: dict = {}

    # stride 自動算出
    try:
        from MDAnalysis import Universe
        u = Universe(top, traj)
        n_frames_total = u.trajectory.n_frames
        stride = max(1, n_frames_total // target_n_frames)
        logger.info("run_analysis: %d frames → stride=%d (~%d frames analyzed)",
                    n_frames_total, stride, n_frames_total // stride)
    except Exception as exc:
        logger.warning("Could not detect n_frames; falling back to stride=1: %s", exc)
        stride = 1

    # 1. Aggregate transition (Fig 1b/1c + Fig 2)
    if run_aggregate:
        try:
            from .aggregate_transition import compute_aggregate_transitions
            results["aggregate"] = compute_aggregate_transitions(
                traj=traj, top=top, out_dir=str(out / "aggregate"),
                n_peptides=n_peptides,
                cutoff_nm=contact_cutoff_nm,
                stride=stride,
                use_heavy_atom=True,
                peptide_selector=peptide_selector,
            )
        except Exception as exc:
            logger.warning("aggregate_transition skipped: %s", exc)
            results["aggregate_error"] = str(exc)

    # 2. Per-residue contacts (Fig 4)
    if run_contacts:
        try:
            from .contact_map import compute_per_residue_contacts
            results["contacts"] = compute_per_residue_contacts(
                traj=traj, top=top, out_dir=str(out / "contacts"),
                n_peptides=n_peptides,
                enhancer_resnames=list(enhancer_resnames),
                bile_salt_resnames=list(bile_salt_resnames),
                cutoff_nm=contact_cutoff_nm,
                stride=stride,
                peptide_selector=peptide_selector,
            )
        except Exception as exc:
            logger.warning("contact_map skipped: %s", exc)
            results["contacts_error"] = str(exc)

    # 3. Secondary structure (Fig 5、 gmx dssp wrap)
    if run_dssp:
        try:
            from .secondary_structure import run_gmx_dssp
            results["dssp"] = run_gmx_dssp(
                traj=traj, tpr=top.replace(".gro", ".tpr") if top.endswith(".gro") else top,
                out_dir=str(out / "dssp"),
                gmx_path=gmx, ndx_path=ndx,
            )
        except Exception as exc:
            logger.warning("dssp skipped: %s", exc)
            results["dssp_error"] = str(exc)

    # 4. Plots (matplotlib)
    if run_plots:
        try:
            from .plots import plot_workflow_outputs
            results["plots"] = plot_workflow_outputs(
                out_dir=str(out), results=results,
            )
        except Exception as exc:
            logger.warning("plots skipped: %s", exc)
            results["plots_error"] = str(exc)

    return results


__all__ = ["run_analysis"]
