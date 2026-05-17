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
from .classifier import ClassificationResult, MolRole, classify
from .colorizer import DEFAULT_COLORS, DrawAttribute, VALID_COLORS, colorize_udf
from .func_tags import (
    BUILTIN_MAPPINGS, CHARMM36, GAFF2, OPENFF_SAGE, OPLS_AA,
    FunctionalTagMapping, detect_force_field, get_mapping, tag_atoms,
)
from .functional_groups import (
    AmideGroup, AmineDonorGroup, CarboxylGroup, HydroxylGroup,
    detect_amides, detect_amine_donors, detect_carboxyls,
    detect_hydroxyls, summarize_groups,
)
from .lifetime import (
    PairKey, PairLifetime,
    compute_autocorrelation, compute_lifetimes,
    integrate_autocorrelation, summarize_lifetimes,
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
    "BUILTIN_MAPPINGS", "detect_force_field", "get_mapping", "tag_atoms",
    # functional groups
    "CarboxylGroup", "AmideGroup", "HydroxylGroup", "AmineDonorGroup",
    "detect_carboxyls", "detect_amides", "detect_hydroxyls",
    "detect_amine_donors", "summarize_groups",
    # detector
    "HBond", "HBondCriteria", "DonorSite", "AcceptorSite",
    "detect_hbonds", "minimum_image_vector",
    # classifier
    "ClassificationResult", "MolRole", "classify",
    # colorizer
    "DrawAttribute", "DEFAULT_COLORS", "VALID_COLORS", "colorize_udf",
    # lifetime (v1.26+ multi-record)
    "PairKey", "PairLifetime",
    "compute_lifetimes", "compute_autocorrelation",
    "integrate_autocorrelation", "summarize_lifetimes",
]

# Optional UI (lazy)
try:
    from .notebook_ui import open_panel  # noqa: F401
    __all__.append("open_panel")
except ImportError:
    pass
