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
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence

import numpy as np

from .bdf_reader import BDFTrajectory
from .classifier import ClassificationResult, classify
from .colorizer import DEFAULT_COLORS, DrawAttribute, colorize_udf
from .func_tags import FunctionalTagMapping, detect_force_field, get_mapping
from .functional_groups import (
    AmideGroup, AmineDonorGroup, CarboxylGroup,
    detect_amides, detect_amine_donors, detect_carboxyls,
)
from .hbond_detector import (
    AcceptorSite, DonorSite, HBond, HBondCriteria, detect_hbonds
)
from .lifetime import (
    PairLifetime, compute_autocorrelation, compute_lifetimes,
    integrate_autocorrelation, summarize_lifetimes,
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
    # Force field (None = auto-detect from atom types)
    force_field: Optional[str] = None
    # Donor / acceptor functional-group selection.
    # If both are None, defaults to {COOH→COOH, COOH→amide} as in v1.25.0.
    # When set, accepts subset of {"carboxyl", "amide_donor", "amine_donor", "hydroxyl"}
    # for donors and {"carboxyl_O", "amide_O", "hydroxyl_O", "ether_O"} for acceptors.
    donor_groups: Optional[List[str]] = None
    acceptor_groups: Optional[List[str]] = None
    # Lifetime analysis (multi-record only)
    compute_lifetime: bool = True
    gap_tolerance: int = 0      # for intermittent lifetime
    dt: float = 1.0             # time between consecutive records (user units, e.g. ps)
    autocorr_max_lag: Optional[int] = None   # None = N/2

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
        self.mapping: Optional[FunctionalTagMapping] = None
        self.carboxyls: List[CarboxylGroup] = []
        self.amides: List[AmideGroup] = []
        self.amine_donors: List[AmineDonorGroup] = []
        self.frame_results: List[FrameResult] = []
        self.lifetimes: List[PairLifetime] = []
        self.autocorr: Optional[np.ndarray] = None
        self.tau_hb: Optional[float] = None

    def load(self) -> None:
        c = self.config
        self.traj = BDFTrajectory(c.bdf_path)
        self.traj.load_topology()
        # Resolve force field
        if c.force_field is None:
            sample_types = {a.atom_type for m in self.traj.molecules for a in m.atoms}
            try:
                ff_name = detect_force_field(sample_types)
            except Exception:
                ff_name = "GAFF2"
        else:
            ff_name = c.force_field
        self.mapping = get_mapping(ff_name)
        self.carboxyls = detect_carboxyls(self.traj.molecules, self.mapping)
        self.amides = detect_amides(self.traj.molecules, self.mapping)
        self.amine_donors = detect_amine_donors(
            self.traj.molecules, mapping=self.mapping,
        )
        if c.verbose:
            print(f"Loaded {c.bdf_path}")
            print(f"  {len(self.traj.molecules)} molecules, "
                  f"{self.traj.n_records} record(s)")
            print(f"  force field: {self.mapping.force_field}")
            print(f"  detected {len(self.carboxyls)} carboxyls, "
                  f"{len(self.amides)} amides "
                  f"({sum(1 for a in self.amides if a.tert)} tert), "
                  f"{len(self.amine_donors)} N-H donors")

    def run(self) -> List[FrameResult]:
        if self.traj is None:
            self.load()
        c = self.config
        end = c.record_end if c.record_end >= 0 else self.traj.n_records
        criteria = c.get_criteria()

        # Build donor and acceptor sites from selected functional groups.
        donor_groups = c.donor_groups or ["carboxyl"]
        acceptor_groups = c.acceptor_groups or ["carboxyl_O", "amide_O"]

        donors_carb: List[DonorSite] = []
        donors_amine: List[DonorSite] = []
        donors_hydroxyl: List[DonorSite] = []
        if "carboxyl" in donor_groups:
            donors_carb = [
                DonorSite(g.mol_index, g.oh_atom, g.ho_atom)
                for g in self.carboxyls
            ]
        if "amine_donor" in donor_groups or "amide_donor" in donor_groups:
            # Filter amine_donors by amide vs amine if both keys present
            keep_amide = "amide_donor" in donor_groups
            keep_amine = "amine_donor" in donor_groups
            donors_amine = [
                DonorSite(g.mol_index, g.n_atom, g.h_atom)
                for g in self.amine_donors
                if (g.from_amide and keep_amide) or
                   (not g.from_amide and keep_amine)
            ]
        if "hydroxyl" in donor_groups:
            from .functional_groups import detect_hydroxyls
            for g in detect_hydroxyls(self.traj.molecules, mapping=self.mapping):
                donors_hydroxyl.append(
                    DonorSite(g.mol_index, g.oh_atom, g.ho_atom)
                )

        carb_donors_all = donors_carb + donors_amine + donors_hydroxyl

        accept_carb_o = (
            [AcceptorSite(g.mol_index, g.o_atom) for g in self.carboxyls]
            if "carboxyl_O" in acceptor_groups else []
        )
        accept_amide_o = (
            [AcceptorSite(g.mol_index, g.o_atom) for g in self.amides]
            if "amide_O" in acceptor_groups else []
        )
        # Legacy backward-compat names used downstream
        carb_donors = donors_carb if donors_carb else carb_donors_all
        carb_acceptors = accept_carb_o
        amide_acceptors = accept_amide_o

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

        # Lifetime analysis (multi-record only)
        if c.compute_lifetime and len(self.frame_results) >= 2:
            out_paths.update(self._write_lifetime_outputs())

        return out_paths

    def _write_lifetime_outputs(self) -> Dict[str, str]:
        """Compute and write lifetime stats + autocorrelation."""
        c = self.config
        out_paths: Dict[str, str] = {}

        # Combine cc + ca hbonds per record (any pair, regardless of acceptor type)
        per_rec_all: List[List[HBond]] = [
            list(fr.hbonds_cc) + list(fr.hbonds_ca) for fr in self.frame_results
        ]
        rec_indices = [fr.record for fr in self.frame_results]
        self.lifetimes = compute_lifetimes(
            per_rec_all, record_indices=rec_indices,
            gap_tolerance=c.gap_tolerance,
        )
        self.autocorr = compute_autocorrelation(
            per_rec_all, max_lag=c.autocorr_max_lag,
        )
        self.tau_hb = integrate_autocorrelation(self.autocorr, dt=c.dt)

        # CSV: per-pair lifetime
        lt_path = f"{c.out_prefix}_lifetime.csv"
        with open(lt_path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow([
                "donor_mol", "donor_d", "donor_h",
                "acceptor_mol", "acceptor_a",
                "total_present", "occupancy",
                "continuous_max", "continuous_mean",
                "intermittent_max", "intermittent_mean",
            ])
            for l in self.lifetimes:
                w.writerow([
                    *l.pair,
                    l.total_present,
                    f"{l.occupancy:.4f}",
                    l.continuous_max, f"{l.continuous_mean:.2f}",
                    l.intermittent_max, f"{l.intermittent_mean:.2f}",
                ])
        out_paths["lifetime"] = lt_path

        # CSV: autocorrelation
        ac_path = f"{c.out_prefix}_autocorr.csv"
        with open(ac_path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["lag", "time", "C(t)"])
            for lag, ct in enumerate(self.autocorr):
                w.writerow([lag, f"{lag * c.dt:.4f}", f"{ct:.6f}"])
        out_paths["autocorr"] = ac_path

        # PNG plot
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            t = np.arange(len(self.autocorr)) * c.dt
            fig, ax = plt.subplots(figsize=(7, 4))
            ax.plot(t, self.autocorr, "o-", color="purple", markersize=3)
            ax.set_xlabel(f"time")
            ax.set_ylabel("C(t)")
            ax.set_title(
                f"H-bond autocorrelation (τ_HB ≈ {self.tau_hb:.3f} "
                f"per dt={c.dt})"
            )
            ax.grid(True, alpha=0.3)
            ax.axhline(0, color="black", lw=0.5)
            fig.tight_layout()
            ac_png = f"{c.out_prefix}_autocorr.png"
            fig.savefig(ac_png, dpi=100)
            plt.close(fig)
            out_paths["autocorr_plot"] = ac_png
        except ImportError:
            pass

        if c.verbose:
            stats = summarize_lifetimes(self.lifetimes)
            print(f"\n--- Lifetime summary ---")
            print(f"  unique pairs: {stats['n_unique_pairs']}")
            print(f"  mean occupancy: {stats['mean_occupancy']:.3f}")
            print(f"  longest continuous run: {stats['max_continuous']} frames")
            print(f"  τ_HB (integral C(t) dt): {self.tau_hb:.4f} "
                  f"(dt={c.dt})")

        return out_paths
