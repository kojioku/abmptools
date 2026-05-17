"""
analyzer.py
-----------
High-level orchestrator: BDF in → CSV/colored-BDF/plot out.

Coordinates bdf_reader + functional_groups + hbond_detector +
classifier + colorizer for one or more trajectory frames.
"""
from __future__ import annotations

import csv
import os
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence

import numpy as np

from .bdf_reader import BDFTrajectory
from .classifier import ClassificationResult, classify
from .colorizer import DEFAULT_COLORS, DrawAttribute, colorize_udf
from .functional_groups import (
    AmideGroup, CarboxylGroup, detect_amides, detect_carboxyls
)
from .hbond_detector import (
    AcceptorSite, DonorSite, HBond, HBondCriteria, detect_hbonds
)


@dataclass
class FrameResult:
    """Result for one record/frame."""
    record: int
    n_dual_mols: int
    n_single_mols: int
    n_free_mols: int
    n_hbonds_cc: int
    n_hbonds_ca: int
    classification: ClassificationResult
    hbonds_cc: List[HBond]
    hbonds_ca: List[HBond]


@dataclass
class AnalyzerConfig:
    """Top-level configuration."""
    bdf_path: str
    out_prefix: str
    criteria_mode: str = "luzar-chandler"   # or "strict" or "custom"
    custom_criteria: Optional[HBondCriteria] = None
    record_start: int = 0
    record_end: int = -1
    base_mol_name: str = "IMC"
    color_attrs: Optional[Dict[str, DrawAttribute]] = None
    do_colorize: bool = True
    do_plot: bool = True
    verbose: bool = True

    def get_criteria(self) -> HBondCriteria:
        if self.criteria_mode == "luzar-chandler":
            return HBondCriteria.luzar_chandler()
        if self.criteria_mode == "strict":
            return HBondCriteria.strict()
        if self.criteria_mode == "custom" and self.custom_criteria is not None:
            return self.custom_criteria
        raise ValueError(f"Unknown criteria_mode: {self.criteria_mode}")


class Analyzer:
    """Top-level analyzer."""

    def __init__(self, config: AnalyzerConfig):
        self.config = config
        self.traj: Optional[BDFTrajectory] = None
        self.carboxyls: List[CarboxylGroup] = []
        self.amides: List[AmideGroup] = []
        self.frame_results: List[FrameResult] = []

    def load(self) -> None:
        c = self.config
        self.traj = BDFTrajectory(c.bdf_path)
        self.traj.load_topology()
        self.carboxyls = detect_carboxyls(self.traj.molecules)
        self.amides = detect_amides(self.traj.molecules)
        if c.verbose:
            print(f"Loaded {c.bdf_path}")
            print(f"  {len(self.traj.molecules)} molecules, "
                  f"{self.traj.n_records} record(s)")
            print(f"  detected {len(self.carboxyls)} carboxyls, "
                  f"{len(self.amides)} amides")

    def run(self) -> List[FrameResult]:
        if self.traj is None:
            self.load()
        c = self.config
        end = c.record_end if c.record_end >= 0 else self.traj.n_records
        criteria = c.get_criteria()

        carb_donors = [
            DonorSite(g.mol_index, g.oh_atom, g.ho_atom)
            for g in self.carboxyls
        ]
        carb_acceptors = [
            AcceptorSite(g.mol_index, g.o_atom) for g in self.carboxyls
        ]
        amide_acceptors = [
            AcceptorSite(g.mol_index, g.o_atom) for g in self.amides
        ]

        n_mol = len(self.traj.molecules)
        results: List[FrameResult] = []
        for rec in range(c.record_start, end):
            frame = self.traj.get_frame(rec)
            hb_cc = detect_hbonds(
                carb_donors, carb_acceptors, frame.positions,
                frame.cell, criteria=criteria
            )
            hb_ca = detect_hbonds(
                carb_donors, amide_acceptors, frame.positions,
                frame.cell, criteria=criteria
            )
            cls = classify(n_mol, hb_cc, hb_ca, self.carboxyls, self.amides)
            fr = FrameResult(
                record=rec,
                n_dual_mols=cls.n_dual_mols,
                n_single_mols=cls.n_single_mols,
                n_free_mols=cls.n_free_mols,
                n_hbonds_cc=len(hb_cc),
                n_hbonds_ca=len(hb_ca),
                classification=cls,
                hbonds_cc=hb_cc,
                hbonds_ca=hb_ca,
            )
            results.append(fr)
            if c.verbose:
                print(f"  rec={rec}: dual={fr.n_dual_mols}, "
                      f"single={fr.n_single_mols}, free={fr.n_free_mols} "
                      f"(cc={fr.n_hbonds_cc}, ca={fr.n_hbonds_ca})")
        self.frame_results = results
        return results

    def write_outputs(self) -> Dict[str, str]:
        """Write all output artifacts. Returns {kind: path} dict."""
        c = self.config
        os.makedirs(os.path.dirname(os.path.abspath(c.out_prefix)) or ".",
                    exist_ok=True)
        out_paths = {}

        # summary CSV
        summary_path = f"{c.out_prefix}_summary.csv"
        with open(summary_path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["record", "n_dual_mols", "n_single_mols",
                        "n_free_mols", "n_hbonds_cc", "n_hbonds_ca"])
            for fr in self.frame_results:
                w.writerow([fr.record, fr.n_dual_mols, fr.n_single_mols,
                            fr.n_free_mols, fr.n_hbonds_cc, fr.n_hbonds_ca])
        out_paths["summary"] = summary_path

        # pairs CSV
        pairs_path = f"{c.out_prefix}_pairs.csv"
        with open(pairs_path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow([
                "record", "kind", "donor_mol", "donor_d", "donor_h",
                "acceptor_mol", "acceptor_a", "d_da", "d_ha", "angle"
            ])
            for fr in self.frame_results:
                for hb in fr.hbonds_cc:
                    w.writerow([
                        fr.record, "cc", hb.donor_mol, hb.donor_d, hb.donor_h,
                        hb.acceptor_mol, hb.acceptor_a,
                        f"{hb.d_da:.4f}", f"{hb.d_ha:.4f}", f"{hb.angle:.2f}"
                    ])
                for hb in fr.hbonds_ca:
                    w.writerow([
                        fr.record, "ca", hb.donor_mol, hb.donor_d, hb.donor_h,
                        hb.acceptor_mol, hb.acceptor_a,
                        f"{hb.d_da:.4f}", f"{hb.d_ha:.4f}", f"{hb.angle:.2f}"
                    ])
        out_paths["pairs"] = pairs_path

        # colorized BDF (uses last frame's classification by default)
        if c.do_colorize and self.frame_results:
            last_cls = self.frame_results[-1].classification
            colored_path = f"{c.out_prefix}_colored.bdf"
            colorize_udf(c.bdf_path, colored_path, last_cls,
                         base_mol_name=c.base_mol_name,
                         color_attrs=c.color_attrs)
            out_paths["colored"] = colored_path

        # plot
        if c.do_plot and len(self.frame_results) >= 1:
            plot_path = f"{c.out_prefix}_count.png"
            try:
                import matplotlib
                matplotlib.use("Agg")
                import matplotlib.pyplot as plt
                recs = [fr.record for fr in self.frame_results]
                fig, ax = plt.subplots(figsize=(7, 4))
                ax.plot(recs, [fr.n_dual_mols for fr in self.frame_results],
                        "o-", color="red", label="dual COOH-COOH mols")
                ax.plot(recs, [fr.n_single_mols for fr in self.frame_results],
                        "s-", color="blue", label="single COOH-amide mols")
                ax.plot(recs, [fr.n_free_mols for fr in self.frame_results],
                        "^-", color="gray", label="free mols")
                ax.set_xlabel("record")
                ax.set_ylabel("number of molecules")
                ax.set_title(f"H-bond classification ({c.criteria_mode})")
                ax.legend()
                ax.grid(True, alpha=0.3)
                fig.tight_layout()
                fig.savefig(plot_path, dpi=100)
                plt.close(fig)
                out_paths["plot"] = plot_path
            except ImportError:
                print("Warning: matplotlib not available; skipping plot")

        return out_paths
