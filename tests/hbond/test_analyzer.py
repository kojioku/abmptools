"""Analyzer.run()-level tests (frame iteration, record subsampling, auto mode)."""
from types import SimpleNamespace

from abmptools.hbond.analyzer import Analyzer, AnalyzerConfig
from abmptools.hbond.bdf_reader import TrajectoryFrame, CellBox


class _StubTraj:
    """Minimal trajectory: N empty frames, no molecules/groups."""

    def __init__(self, n_records):
        self.n_records = n_records
        self.molecules = []

    def get_frame(self, rec):
        return TrajectoryFrame(record=rec, cell=CellBox(10.0, 10.0, 10.0),
                               positions=[])


def _run_with_stride(n_records, start, end, stride):
    cfg = AnalyzerConfig(
        bdf_path="x", out_prefix="x", classify_mode="generic",
        donor_groups=[], acceptor_groups=[], verbose=False,
        record_start=start, record_end=end, record_stride=stride,
    )
    an = Analyzer.__new__(Analyzer)          # bypass __init__/load (stub traj)
    an.config = cfg
    an.traj = _StubTraj(n_records)
    an.carboxyls = []
    an.amides = []
    return [fr.record for fr in an.run()]


def test_record_stride_default_every_frame():
    assert _run_with_stride(5, 0, -1, 1) == [0, 1, 2, 3, 4]


def test_record_stride_subsamples():
    # 30-record trajectory thinned by 10 -> records 0, 10, 20
    assert _run_with_stride(30, 0, 30, 10) == [0, 10, 20]


def test_record_stride_with_start_offset():
    assert _run_with_stride(20, 3, 20, 5) == [3, 8, 13, 18]


def test_record_stride_zero_or_negative_falls_back_to_one():
    # Guard: stride <= 0 must behave as stride=1 (no infinite/empty loop)
    assert _run_with_stride(4, 0, -1, 0) == [0, 1, 2, 3]


# --- classify_mode='auto' (functional-group auto detection) ------------------

def _analyzer_with_groups(carboxyls, amides, amine_donors):
    """Build an Analyzer with mocked detected groups (no real trajectory)."""
    an = Analyzer.__new__(Analyzer)
    an.config = AnalyzerConfig(
        bdf_path="x", out_prefix="x", classify_mode="auto", verbose=False,
    )
    an.traj = _StubTraj(1)        # molecules=[] → detect_hydroxyls returns []
    an.mapping = None
    an.carboxyls = carboxyls
    an.amides = amides
    an.amine_donors = amine_donors
    return an


def test_auto_groups_carboxyl_only():
    an = _analyzer_with_groups([SimpleNamespace()], [], [])
    donors, acceptors = an._auto_groups()
    assert donors == ["carboxyl"]
    assert acceptors == ["carboxyl_O"]


def test_auto_groups_amide_only_apz_like():
    # APZ-like: lactam → amide_O acceptor + secondary-amide N-H donor, no COOH
    an = _analyzer_with_groups([], [SimpleNamespace()],
                               [SimpleNamespace(from_amide=True)])
    donors, acceptors = an._auto_groups()
    assert donors == ["amide_donor"]
    assert acceptors == ["amide_O"]


def test_auto_groups_mixed():
    an = _analyzer_with_groups(
        [SimpleNamespace()], [SimpleNamespace()],
        [SimpleNamespace(from_amide=True), SimpleNamespace(from_amide=False)],
    )
    donors, acceptors = an._auto_groups()
    assert "carboxyl" in donors and "amide_donor" in donors and "amine_donor" in donors
    assert "carboxyl_O" in acceptors and "amide_O" in acceptors


def test_auto_groups_none_detected():
    an = _analyzer_with_groups([], [], [])
    donors, acceptors = an._auto_groups()
    assert donors == [] and acceptors == []


def test_run_auto_resolves_to_generic():
    # classify_mode='auto' must be accepted and resolve to the generic path.
    an = _analyzer_with_groups([], [], [])
    an.config.record_start = 0
    an.config.record_end = -1
    an.config.record_stride = 1
    an.traj = _StubTraj(3)
    results = an.run()
    assert [fr.record for fr in results] == [0, 1, 2]
    assert an.config.classify_mode == "generic"   # auto resolved downstream
