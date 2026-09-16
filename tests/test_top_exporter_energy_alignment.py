# -*- coding: utf-8 -*-
"""energy の行と trajectory の frame の対応付け。

**既定の設定では energy のほうが密になる。** amorphous の mdp は
``nstenergy = 1000`` / ``nstxout-compressed = 5000`` なので、500 ps の run で
**energy 501 行 / frame 101 枚**、つまり 5:1 になる (実測)。

ここを行番号で対応させると、frame 1 (t=5 ps) に t=1 ps の energy が入る。
**値としては妥当な範囲に収まるので、出力を見ても気付けない。** 時刻で
最近傍を取っていることを、この検査で固定する。
"""
from __future__ import annotations

from types import SimpleNamespace

from abmptools.gro2udf.top_exporter import _aggregate_statistics_per_frame

KEY = ("Pressure", "", "[bar]")


def _run(frame_times, energy_times, pressures):
    frames = [SimpleNamespace(time=t) for t in frame_times]
    return _aggregate_statistics_per_frame(
        frames, list(energy_times),
        {"Pressure": list(pressures), "Pres. DC": [0.0] * len(pressures)},
    )


class TestEnergyDenserThanTrajectory:
    # energy は 1 ps 刻み、 frame は 5 ps 刻み (実機と同じ 5:1)
    E_TIMES = [float(i) for i in range(11)]
    E_VALS = [100.0 + i for i in range(11)]        # t ps -> 100+t bar
    F_TIMES = [0.0, 5.0, 10.0]

    def test_instantaneous_matches_by_time_not_by_row(self):
        per_frame = _run(self.F_TIMES, self.E_TIMES, self.E_VALS)
        got = [per_frame[i][KEY][0] for i in range(3)]
        assert got == [100.0, 105.0, 110.0], got
        # 行番号で拾っていれば [100, 101, 102] になる
        assert got != [100.0, 101.0, 102.0]

    def test_batch_average_spans_the_rows_since_the_previous_frame(self):
        """frame 間に挟まった行を捨てず、 平均として残す。"""
        per_frame = _run(self.F_TIMES, self.E_TIMES, self.E_VALS)
        batch1 = per_frame[1][KEY][1]
        assert abs(batch1 - sum(self.E_VALS[1:6]) / 5) < 1e-9, batch1
        batch2 = per_frame[2][KEY][1]
        assert abs(batch2 - sum(self.E_VALS[6:11]) / 5) < 1e-9, batch2

    def test_total_average_runs_from_the_first_row(self):
        per_frame = _run(self.F_TIMES, self.E_TIMES, self.E_VALS)
        total2 = per_frame[2][KEY][2]
        assert abs(total2 - sum(self.E_VALS[0:11]) / 11) < 1e-9, total2


class TestOneToOneStillWorks:
    def test_aligned_series_is_unchanged(self):
        per_frame = _run([0.0, 1.0, 2.0], [0.0, 1.0, 2.0], [10.0, 20.0, 30.0])
        assert [per_frame[i][KEY][0] for i in range(3)] == [10.0, 20.0, 30.0]
        # 1:1 なら batch は 1 行分なので instantaneous と同じ
        assert per_frame[1][KEY][1] == 20.0


class TestTrajectoryDenserThanEnergy:
    """逆向き (frame のほうが密) でも、 一番近い行を使う。"""

    def test_nearest_row_is_reused(self):
        per_frame = _run([0.0, 1.0, 2.0, 3.0], [0.0, 2.0], [50.0, 70.0])
        got = [per_frame[i][KEY][0] for i in range(4)]
        # t=1.0 は 0.0 と 2.0 の中間。 どちらでもよいが、 t=3.0 は 2.0 側
        assert got[0] == 50.0 and got[2] == 70.0 and got[3] == 70.0, got
