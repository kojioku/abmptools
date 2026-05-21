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
import shutil
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from .bdf_reader import BDFTrajectory
from .classifier import (
    AmideRole, CarboxylRole, ClassificationResult,
    FunctionalGroupClassification, classify,
)
from .colorizer import (
    DEFAULT_ACTION_COLORS, DEFAULT_COLORS, DEFAULT_GENERIC_COLORS,
    DrawAttribute, colorize_udf, colorize_udf_action,
    colorize_udf_action_generic, write_hbond_attributes,
    write_hbond_attributes_generic, write_show_python_script,
    write_show_python_script_generic,
)
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
from .pair_type_stats import (
    GenericPairClassification, PairTypeStat, classify_generic,
    summarize_pair_stats,
)


@dataclass
class FrameResult:
    """Result for one record/frame.

    The ``n_*_mols`` fields are the mol-level representative role counts
    (used by ``colorize_udf``). The primary metrics are the per-functional
    -group counts inside ``classification`` (``n_carboxyls_dual``,
    ``ratio_amide_accept``, etc.).
    """
    record: int
    n_dual_mols: int
    n_single_mols: int
    n_free_mols: int
    n_hbonds_cc: int
    n_hbonds_ca: int
    classification: Optional[FunctionalGroupClassification]
    hbonds_cc: List[HBond]
    hbonds_ca: List[HBond]
    n_chain_mols: int = 0
    # generic mode (classify_mode="generic") fields
    generic_classification: Optional[GenericPairClassification] = None
    hbonds_by_pair_type: Dict[Tuple[str, str], List[HBond]] = field(
        default_factory=dict
    )


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
    do_copy_uncolored: bool = True   # writes ``<prefix>.bdf`` with Mol_Name kept
    # Classification mode:
    #   "imc"     — COOH-centric 4-species (dual/chain/single/free + amide accept/free)
    #               + per-COOH classification.csv. Defaults preserved for IMC use.
    #   "generic" — generic donor-type × acceptor-type pair statistics + per-atom
    #               role tags (Donor/Acceptor/Both/Candidate) for any system
    #               (PVA, peptide, alcohol, mixtures, ...). No dual/chain/single
    #               concept; output is donor->acceptor pair counts.
    classify_mode: str = "imc"
    # Element + bond-graph fallback (default ON): when ``Atom_Type_Name``
    # is None or per-atom unique (OpenFF SMIRNOFF case), tag O/H/N/C from
    # element + neighbours so detect_carboxyls / detect_amides /
    # detect_hydroxyls still work. Set to False to require FF-mapped tags
    # only (strict mode for debugging / partial mappings).
    use_element_fallback: bool = True
    # Append per-atom Attributes[] entries (Name=attribute, Value=role) on
    # functional-group atoms of the uncolored copy. Useful for J-OCTA Viewer
    # category filtering when sphere overlay rendering is unavailable.
    do_write_attributes: bool = True
    attributes_name: str = "hbond"
    # Coloring strategy:
    #   "molname" — rename Mol_Name + Draw_Attributes.Molecule[] (v1.25 legacy)
    #   "action"  — emit Python action .act that paints per functional group
    #               (Mol_Name preserved; J-OCTA pre-render compatible)
    #   "both"    — write both <prefix>_colored.bdf and <prefix>_action.bdf
    colorize_mode: str = "molname"
    action_colors: Optional[Dict[str, List[float]]] = None
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
        # Propagate the element-fallback flag to the module-scope switch
        # used inside detect_carboxyls / detect_amides / ...
        from . import functional_groups as _fg
        _fg.USE_ELEMENT_FALLBACK = bool(c.use_element_fallback)
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

    def _build_donor_sites_by_type(
        self, donor_groups: List[str],
    ) -> Tuple[Dict[str, List[DonorSite]], Dict[str, List[Tuple[int, int, int]]]]:
        """Per donor-type DonorSite lists + meta tuples (mi, d, h) for generic mode."""
        sites_by_type: Dict[str, List[DonorSite]] = {}
        meta_by_type: Dict[str, List[Tuple[int, int, int]]] = {}
        if "carboxyl" in donor_groups:
            sites = [
                DonorSite(g.mol_index, g.oh_atom, g.ho_atom)
                for g in self.carboxyls
            ]
            sites_by_type["carboxyl"] = sites
            meta_by_type["carboxyl"] = [
                (g.mol_index, g.oh_atom, g.ho_atom) for g in self.carboxyls
            ]
        if "amide_donor" in donor_groups:
            amide_sites = [
                DonorSite(g.mol_index, g.n_atom, g.h_atom)
                for g in self.amine_donors if g.from_amide
            ]
            if amide_sites:
                sites_by_type["amide_donor"] = amide_sites
                meta_by_type["amide_donor"] = [
                    (g.mol_index, g.n_atom, g.h_atom)
                    for g in self.amine_donors if g.from_amide
                ]
        if "amine_donor" in donor_groups:
            amine_sites = [
                DonorSite(g.mol_index, g.n_atom, g.h_atom)
                for g in self.amine_donors if not g.from_amide
            ]
            if amine_sites:
                sites_by_type["amine_donor"] = amine_sites
                meta_by_type["amine_donor"] = [
                    (g.mol_index, g.n_atom, g.h_atom)
                    for g in self.amine_donors if not g.from_amide
                ]
        if "hydroxyl" in donor_groups:
            from .functional_groups import detect_hydroxyls
            hydroxyls = detect_hydroxyls(
                self.traj.molecules, mapping=self.mapping,
            )
            sites = [
                DonorSite(g.mol_index, g.oh_atom, g.ho_atom) for g in hydroxyls
            ]
            if sites:
                sites_by_type["hydroxyl"] = sites
                meta_by_type["hydroxyl"] = [
                    (g.mol_index, g.oh_atom, g.ho_atom) for g in hydroxyls
                ]
        return sites_by_type, meta_by_type

    def _build_acceptor_sites_by_type(
        self, acceptor_groups: List[str],
    ) -> Tuple[Dict[str, List[AcceptorSite]], Dict[str, List[Tuple[int, int]]]]:
        """Per acceptor-type AcceptorSite lists + meta (mi, a) for generic mode."""
        sites_by_type: Dict[str, List[AcceptorSite]] = {}
        meta_by_type: Dict[str, List[Tuple[int, int]]] = {}
        if "carboxyl_O" in acceptor_groups:
            sites = [
                AcceptorSite(g.mol_index, g.o_atom) for g in self.carboxyls
            ]
            sites_by_type["carboxyl_O"] = sites
            meta_by_type["carboxyl_O"] = [
                (g.mol_index, g.o_atom) for g in self.carboxyls
            ]
        if "amide_O" in acceptor_groups:
            sites = [
                AcceptorSite(g.mol_index, g.o_atom) for g in self.amides
            ]
            sites_by_type["amide_O"] = sites
            meta_by_type["amide_O"] = [
                (g.mol_index, g.o_atom) for g in self.amides
            ]
        if "hydroxyl_O" in acceptor_groups:
            from .functional_groups import detect_hydroxyls
            hydroxyls = detect_hydroxyls(
                self.traj.molecules, mapping=self.mapping,
            )
            sites = [
                AcceptorSite(g.mol_index, g.oh_atom) for g in hydroxyls
            ]
            if sites:
                sites_by_type["hydroxyl_O"] = sites
                meta_by_type["hydroxyl_O"] = [
                    (g.mol_index, g.oh_atom) for g in hydroxyls
                ]
        if "ether_O" in acceptor_groups:
            # detect ether-type O via tag scan (hydroxyl_O w/o H counts as ether)
            from .func_tags import TAG_CARBONYL_O
            sites: List[AcceptorSite] = []
            for mi, topo in enumerate(self.traj.molecules):
                tags = [self.mapping.get_tag(a.atom_type) for a in topo.atoms]
                for ai, t in enumerate(tags):
                    if t in (TAG_CARBONYL_O,):
                        pass  # already covered by carboxyl_O / amide_O / etc.
            # Best-effort fallback: rely on user-provided detect_hydroxyls path
            # for now; users with ether systems can extend func_tags mapping.
            if sites:
                sites_by_type["ether_O"] = sites
                meta_by_type["ether_O"] = [(s.mol_index, s.acceptor_a) for s in sites]
        return sites_by_type, meta_by_type

    def run(self) -> List[FrameResult]:
        if self.traj is None:
            self.load()
        c = self.config
        end = c.record_end if c.record_end >= 0 else self.traj.n_records
        criteria = c.get_criteria()
        mode = (c.classify_mode or "imc").lower()
        if mode not in {"imc", "generic"}:
            raise ValueError(
                f"Unknown classify_mode={c.classify_mode!r}; "
                "must be 'imc' or 'generic'"
            )

        donor_groups = c.donor_groups or (
            ["carboxyl"] if mode == "imc" else []
        )
        acceptor_groups = c.acceptor_groups or (
            ["carboxyl_O", "amide_O"] if mode == "imc" else []
        )

        donor_sites_by_type, donor_meta_by_type = (
            self._build_donor_sites_by_type(donor_groups)
        )
        acceptor_sites_by_type, acceptor_meta_by_type = (
            self._build_acceptor_sites_by_type(acceptor_groups)
        )

        # IMC mode keeps backward-compat names for hb_cc/hb_ca path
        carb_donors = donor_sites_by_type.get("carboxyl", [])
        carb_acceptors = acceptor_sites_by_type.get("carboxyl_O", [])
        amide_acceptors = acceptor_sites_by_type.get("amide_O", [])

        n_mol = len(self.traj.molecules)
        results: List[FrameResult] = []
        for rec in range(c.record_start, end):
            frame = self.traj.get_frame(rec)

            if mode == "imc":
                hb_cc = detect_hbonds(
                    carb_donors, carb_acceptors, frame.positions,
                    frame.cell, criteria=criteria,
                )
                hb_ca = detect_hbonds(
                    carb_donors, amide_acceptors, frame.positions,
                    frame.cell, criteria=criteria,
                )
                cls = classify(
                    n_mol, hb_cc, hb_ca, self.carboxyls, self.amides,
                )
                fr = FrameResult(
                    record=rec,
                    n_dual_mols=cls.n_dual_mols,
                    n_chain_mols=cls.n_chain_mols,
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
                    print(
                        f"  rec={rec}: "
                        f"COOH dual/chain/single/free="
                        f"{cls.n_carboxyls_dual}/{cls.n_carboxyls_chain}/"
                        f"{cls.n_carboxyls_single}/{cls.n_carboxyls_free} "
                        f"({cls.ratio_carboxyl_dual*100:.0f}%/"
                        f"{cls.ratio_carboxyl_chain*100:.0f}%/"
                        f"{cls.ratio_carboxyl_single*100:.0f}%/"
                        f"{cls.ratio_carboxyl_free*100:.0f}%), "
                        f"amide accept/free="
                        f"{cls.n_amides_accept}/{cls.n_amides_free} "
                        f"({cls.ratio_amide_accept*100:.0f}%/"
                        f"{cls.ratio_amide_free*100:.0f}%), "
                        f"hb_cc={fr.n_hbonds_cc}, hb_ca={fr.n_hbonds_ca}"
                    )
            else:  # mode == "generic"
                hbonds_by_pair: Dict[Tuple[str, str], List[HBond]] = {}
                total_hb = 0
                for dt, donors in donor_sites_by_type.items():
                    for at, acceptors in acceptor_sites_by_type.items():
                        hbs = detect_hbonds(
                            donors, acceptors, frame.positions, frame.cell,
                            criteria=criteria,
                        )
                        hbonds_by_pair[(dt, at)] = hbs
                        total_hb += len(hbs)
                gcls = classify_generic(
                    donor_meta_by_type, acceptor_meta_by_type, hbonds_by_pair,
                )
                fr = FrameResult(
                    record=rec,
                    n_dual_mols=0, n_chain_mols=0,
                    n_single_mols=0, n_free_mols=0,
                    n_hbonds_cc=0, n_hbonds_ca=0,
                    classification=None,
                    hbonds_cc=[], hbonds_ca=[],
                    generic_classification=gcls,
                    hbonds_by_pair_type=hbonds_by_pair,
                )
                results.append(fr)
                if c.verbose:
                    print(f"  rec={rec}: total H-bonds = {total_hb}")
                    for (dt, at), st in sorted(gcls.pair_stats.items()):
                        denom_d = gcls.n_donor_candidates.get(dt, 0) or 1
                        denom_a = gcls.n_acceptor_candidates.get(at, 0) or 1
                        print(
                            f"    {dt:>13s} -> {at:<13s}  "
                            f"n={st.n_hbonds:4d}  "
                            f"donors={st.n_uniq_donors}/{denom_d} "
                            f"acceptors={st.n_uniq_acceptors}/{denom_a}"
                        )
        self.frame_results = results
        return results

    def write_outputs(self) -> Dict[str, str]:
        """Write all output artifacts. Returns {kind: path} dict."""
        c = self.config
        os.makedirs(os.path.dirname(os.path.abspath(c.out_prefix)) or ".",
                    exist_ok=True)
        out_paths: Dict[str, str] = {}
        mode = (c.classify_mode or "imc").lower()

        if mode == "imc":
            # IMC mode: per-COOH 4-species summary + classification CSV
            summary_path = f"{c.out_prefix}_summary.csv"
            with open(summary_path, "w", newline="") as f:
                w = csv.writer(f)
                w.writerow([
                    "record", "n_carboxyls",
                    "n_carb_dual", "n_carb_chain",
                    "n_carb_single", "n_carb_free",
                    "n_amides", "n_amide_accept", "n_amide_free",
                    "ratio_carb_dual", "ratio_carb_chain",
                    "ratio_carb_single", "ratio_carb_free",
                    "ratio_amide_accept", "ratio_amide_free",
                    "n_hbonds_cc", "n_hbonds_ca",
                    "n_dual_mols", "n_chain_mols",
                    "n_single_mols", "n_free_mols",
                ])
                for fr in self.frame_results:
                    cls = fr.classification
                    w.writerow([
                        fr.record, cls.n_carboxyls,
                        cls.n_carboxyls_dual, cls.n_carboxyls_chain,
                        cls.n_carboxyls_single, cls.n_carboxyls_free,
                        cls.n_amides, cls.n_amides_accept, cls.n_amides_free,
                        f"{cls.ratio_carboxyl_dual:.4f}",
                        f"{cls.ratio_carboxyl_chain:.4f}",
                        f"{cls.ratio_carboxyl_single:.4f}",
                        f"{cls.ratio_carboxyl_free:.4f}",
                        f"{cls.ratio_amide_accept:.4f}",
                        f"{cls.ratio_amide_free:.4f}",
                        fr.n_hbonds_cc, fr.n_hbonds_ca,
                        fr.n_dual_mols, fr.n_chain_mols,
                        fr.n_single_mols, fr.n_free_mols,
                    ])
            out_paths["summary"] = summary_path
        else:
            # generic mode: donor-type × acceptor-type pair statistics CSV
            pair_path = f"{c.out_prefix}_pair_stats.csv"
            with open(pair_path, "w", newline="") as f:
                w = csv.writer(f)
                w.writerow([
                    "record", "donor_type", "acceptor_type",
                    "n_hbonds", "n_uniq_donors", "n_uniq_acceptors",
                    "ratio_donor_busy", "ratio_acceptor_busy",
                ])
                for fr in self.frame_results:
                    gcls = fr.generic_classification
                    if gcls is None:
                        continue
                    for (dt, at), st in sorted(gcls.pair_stats.items()):
                        denom_d = gcls.n_donor_candidates.get(dt, 0)
                        denom_a = gcls.n_acceptor_candidates.get(at, 0)
                        rb_d = st.n_uniq_donors / denom_d if denom_d else 0.0
                        rb_a = st.n_uniq_acceptors / denom_a if denom_a else 0.0
                        w.writerow([
                            fr.record, dt, at,
                            st.n_hbonds, st.n_uniq_donors, st.n_uniq_acceptors,
                            f"{rb_d:.4f}", f"{rb_a:.4f}",
                        ])
            out_paths["pair_stats"] = pair_path

        # classification CSV (per-functional-group role table; imc mode only)
        if mode == "imc":
            cls_path = f"{c.out_prefix}_classification.csv"
            with open(cls_path, "w", newline="") as f:
                w = csv.writer(f)
                w.writerow([
                    "record", "group_type", "mol_index", "group_index",
                    "role", "partner_count", "partners",
                ])
                for fr in self.frame_results:
                    for r in fr.classification.carboxyl_roles:
                        partners = (r.dual_partners if r.role == "dual"
                                    else r.single_acceptors)
                        partner_str = ";".join(
                            f"{m}:{g}" for (m, g) in sorted(partners)
                        )
                        w.writerow([
                            fr.record, "carboxyl",
                            r.mol_index, r.carboxyl_index,
                            r.role, len(partners), partner_str,
                        ])
                    for r in fr.classification.amide_roles:
                        partner_str = ";".join(
                            f"{m}:{g}" for (m, g) in sorted(r.donor_carboxyls)
                        )
                        w.writerow([
                            fr.record, "amide",
                            r.mol_index, r.amide_index,
                            r.role, len(r.donor_carboxyls), partner_str,
                        ])
            out_paths["classification"] = cls_path

        # pairs CSV (works for both modes)
        pairs_path = f"{c.out_prefix}_pairs.csv"
        with open(pairs_path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow([
                "record", "kind", "donor_mol", "donor_d", "donor_h",
                "acceptor_mol", "acceptor_a", "d_da", "d_ha", "angle"
            ])
            for fr in self.frame_results:
                if mode == "imc":
                    for hb in fr.hbonds_cc:
                        w.writerow([
                            fr.record, "cc", hb.donor_mol, hb.donor_d, hb.donor_h,
                            hb.acceptor_mol, hb.acceptor_a,
                            f"{hb.d_da:.4f}", f"{hb.d_ha:.4f}",
                            f"{hb.angle:.2f}",
                        ])
                    for hb in fr.hbonds_ca:
                        w.writerow([
                            fr.record, "ca", hb.donor_mol, hb.donor_d, hb.donor_h,
                            hb.acceptor_mol, hb.acceptor_a,
                            f"{hb.d_da:.4f}", f"{hb.d_ha:.4f}",
                            f"{hb.angle:.2f}",
                        ])
                else:
                    for (dt, at), hbs in sorted(fr.hbonds_by_pair_type.items()):
                        kind = f"{dt}->{at}"
                        for hb in hbs:
                            w.writerow([
                                fr.record, kind,
                                hb.donor_mol, hb.donor_d, hb.donor_h,
                                hb.acceptor_mol, hb.acceptor_a,
                                f"{hb.d_da:.4f}", f"{hb.d_ha:.4f}",
                                f"{hb.angle:.2f}",
                            ])
        out_paths["pairs"] = pairs_path

        # Uncolored copy (Mol_Name preserved, J-OCTA pre-render compatible)
        uncolored_path = None
        if c.do_copy_uncolored:
            uncolored_path = f"{c.out_prefix}.bdf"
            if os.path.abspath(uncolored_path) != os.path.abspath(c.bdf_path):
                shutil.copy(c.bdf_path, uncolored_path)
                out_paths["uncolored"] = uncolored_path
            else:
                if c.verbose:
                    print(f"  (uncolored copy skipped: out path == input)")
                uncolored_path = None

        # Per-atom Attributes[] tagging on the uncolored copy (J-OCTA filter)
        if (c.do_write_attributes and uncolored_path and self.frame_results):
            if mode == "imc" and self.carboxyls:
                write_hbond_attributes(
                    uncolored_path, self.frame_results[-1].classification,
                    carboxyls=self.carboxyls, amides=self.amides,
                    attribute_name=c.attributes_name,
                )
            elif mode == "generic":
                gcls = self.frame_results[-1].generic_classification
                if gcls is not None:
                    write_hbond_attributes_generic(
                        uncolored_path, gcls.atom_role,
                        attribute_name=c.attributes_name,
                    )

        # colorized BDF (uses last frame's classification by default)
        if c.do_colorize and self.frame_results:
            cmode = (c.colorize_mode or "molname").lower()
            if cmode not in {"molname", "action", "both"}:
                raise ValueError(
                    f"Unknown colorize_mode {c.colorize_mode!r}; "
                    "must be one of: molname, action, both"
                )
            target_basename = (
                os.path.basename(f"{c.out_prefix}.bdf")
                if c.do_copy_uncolored
                else os.path.basename(c.bdf_path)
            )
            if mode == "imc":
                last_cls = self.frame_results[-1].classification
                if cmode in {"molname", "both"}:
                    colored_path = f"{c.out_prefix}_colored.bdf"
                    colorize_udf(c.bdf_path, colored_path, last_cls,
                                 base_mol_name=c.base_mol_name,
                                 color_attrs=c.color_attrs)
                    out_paths["colored"] = colored_path
                if cmode in {"action", "both"}:
                    action_bdf = f"{c.out_prefix}_action.bdf"
                    action_act = f"{c.out_prefix}_show.act"
                    colorize_udf_action(
                        c.bdf_path, action_bdf, action_act, last_cls,
                        carboxyls=self.carboxyls, amides=self.amides,
                        action_colors=c.action_colors,
                    )
                    out_paths["action_bdf"] = action_bdf
                    out_paths["action_act"] = action_act
                    script_py = f"{c.out_prefix}_show.py"
                    write_show_python_script(
                        script_py, last_cls,
                        carboxyls=self.carboxyls, amides=self.amides,
                        target_bdf_basename=target_basename,
                        action_colors=c.action_colors,
                    )
                    out_paths["action_script"] = script_py
            else:  # generic mode
                gcls = self.frame_results[-1].generic_classification
                if gcls is not None and cmode in {"action", "both"}:
                    action_bdf = f"{c.out_prefix}_action.bdf"
                    action_act = f"{c.out_prefix}_show.act"
                    colorize_udf_action_generic(
                        c.bdf_path, action_bdf, action_act, gcls.atom_role,
                        action_colors=c.action_colors,
                    )
                    out_paths["action_bdf"] = action_bdf
                    out_paths["action_act"] = action_act
                    script_py = f"{c.out_prefix}_show.py"
                    write_show_python_script_generic(
                        script_py, gcls.atom_role,
                        target_bdf_basename=target_basename,
                        action_colors=c.action_colors,
                    )
                    out_paths["action_script"] = script_py
                # generic mode does not write Mol_Name-renamed <prefix>_colored.bdf

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
                ax.plot(recs, [fr.n_chain_mols for fr in self.frame_results],
                        "D-", color="magenta", label="chain COOH-COOH mols")
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
