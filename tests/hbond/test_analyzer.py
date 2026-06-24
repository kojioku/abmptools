"""Analyzer.run()-level tests (frame iteration, record subsampling)."""
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
