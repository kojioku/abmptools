"""
abmptools.hbond
---------------
Carboxyl / amide hydrogen-bond analysis for COGNAC UDF/BDF trajectories.

Goal: count, per frame, how many molecules participate in
- **dual** COOH-COOH (centrosymmetric carboxyl dimer, both directions)
- **single** COOH→amide C=O (one-direction H-bond to amide carbonyl)
- **free** (neither)

Designed for amorphous indomethacin (IMC) but reusable for any system with
GAFF2-typed carboxyl/amide groups (atomtypes ``c``, ``oh``, ``ho``, ``o``,
``n``).

Public API:
    AnalyzerConfig      Run config dataclass
    Analyzer            High-level orchestrator (load → run → write_outputs)
    HBondCriteria       Geometric thresholds (Luzar-Chandler, strict, custom)
    DonorSite, AcceptorSite, HBond
                        Low-level H-bond primitives
    detect_carboxyls, detect_amides, detect_hydroxyls
                        GAFF2-based functional group detector
    detect_hbonds       Geometric H-bond detector with PBC
    classify            dual/single/free classifier
    colorize_udf        Write 3-color-grouped UDF for gourmet visualization
    BDFTrajectory       Trajectory reader (UDFManager wrapper)

CLI:
    python -m abmptools.hbond <bdf_path> [--out-prefix PREFIX]
                              [--criteria luzar-chandler|strict|custom]
                              [--mol-name IMC] [--no-colorize] [--no-plot]

Dependencies:
    - UDFManager (OCTA)
    - numpy
    - matplotlib (optional, for count plot)
    - ipywidgets, rdkit (optional, for notebook_ui)
"""
from .analyzer import Analyzer, AnalyzerConfig, FrameResult
from .bdf_reader import (
    AtomInfo, BDFTrajectory, BondInfo, CellBox,
    MoleculeTopology, TrajectoryFrame
)
from .classifier import (
    AmideRole, CarboxylRole, ClassificationResult,
    FunctionalGroupClassification, MolRole, classify,
)
from .colorizer import (
    DEFAULT_ACTION_COLORS, DEFAULT_COLORS, DEFAULT_GENERIC_COLORS,
    DrawAttribute, VALID_COLORS,
    colorize_udf, colorize_udf_action, colorize_udf_action_generic,
    write_hbond_attributes, write_hbond_attributes_generic,
    write_show_python_script, write_show_python_script_generic,
)
from .pair_type_stats import (
    GenericPairClassification, PairTypeStat,
    classify_generic, summarize_pair_stats,
)
from .func_tags import (
    BUILTIN_MAPPINGS, CHARMM36, GAFF2, OPENFF_SAGE, OPLS_AA,
    FunctionalTagMapping, detect_force_field, fallback_tag_by_element,
    get_mapping, tag_atoms,
)
from .functional_groups import (
    AmideGroup, AmineDonorGroup, CarboxylGroup, EtherGroup, HydroxylGroup,
    detect_amides, detect_amine_donors, detect_carboxyls,
    detect_ethers, detect_hydroxyls, summarize_groups,
)
from .lifetime import (
    PairKey, PairLifetime,
    compute_autocorrelation, compute_lifetimes,
    integrate_autocorrelation, summarize_lifetimes,
)
from .distance_dist import (
    DEFAULT_COLORS_BY_LABEL, DistanceStats,
    aggregate_distance_angle, aggregate_distances_generic,
    aggregate_distances_imc, compute_distance_stats, default_bin_edges,
    plot_distance_angle_2d, plot_distance_histogram,
    plot_distance_histogram_classified, write_distance_histogram_csv,
    write_distance_stats_csv,
)
from .hbond_detector import (
    AcceptorSite, DonorSite, HBond, HBondCriteria, detect_hbonds,
    minimum_image_vector
)

__all__ = [
    # analyzer
    "Analyzer", "AnalyzerConfig", "FrameResult",
    # reader
    "BDFTrajectory", "TrajectoryFrame", "CellBox",
    "MoleculeTopology", "AtomInfo", "BondInfo",
    # functional tag mapping (v1.26+ FF abstraction)
    "FunctionalTagMapping", "GAFF2", "OPLS_AA", "CHARMM36", "OPENFF_SAGE",
    "BUILTIN_MAPPINGS", "detect_force_field", "fallback_tag_by_element",
    "get_mapping", "tag_atoms",
    # functional groups
    "CarboxylGroup", "AmideGroup", "HydroxylGroup", "AmineDonorGroup",
    "EtherGroup",
    "detect_carboxyls", "detect_amides", "detect_hydroxyls",
    "detect_amine_donors", "detect_ethers", "summarize_groups",
    # detector
    "HBond", "HBondCriteria", "DonorSite", "AcceptorSite",
    "detect_hbonds", "minimum_image_vector",
    # classifier (per-functional-group)
    "CarboxylRole", "AmideRole", "FunctionalGroupClassification",
    "ClassificationResult", "MolRole", "classify",
    # colorizer
    "DrawAttribute", "DEFAULT_COLORS", "VALID_COLORS", "colorize_udf",
    "DEFAULT_ACTION_COLORS", "colorize_udf_action", "write_show_python_script",
    "write_hbond_attributes",
    # generic mode
    "DEFAULT_GENERIC_COLORS", "colorize_udf_action_generic",
    "write_hbond_attributes_generic", "write_show_python_script_generic",
    "GenericPairClassification", "PairTypeStat",
    "classify_generic", "summarize_pair_stats",
    # lifetime (v1.26+ multi-record)
    "PairKey", "PairLifetime",
    "compute_lifetimes", "compute_autocorrelation",
    "integrate_autocorrelation", "summarize_lifetimes",
    # distance / angle distribution
    "DistanceStats", "DEFAULT_COLORS_BY_LABEL", "default_bin_edges",
    "compute_distance_stats",
    "aggregate_distances_imc", "aggregate_distances_generic",
    "aggregate_distance_angle",
    "write_distance_stats_csv", "write_distance_histogram_csv",
    "plot_distance_histogram", "plot_distance_histogram_classified",
    "plot_distance_angle_2d",
]

# Optional UI (lazy)
try:
    from .notebook_ui import open_panel  # noqa: F401
    __all__.append("open_panel")
except ImportError:
    pass
