# -*- coding: utf-8 -*-
"""熱浴・圧力浴の変換。 どれも「変換は通るが値が間違っている」類の欠陥。

COGNAC は熱浴・圧力浴を**質量**で持ち (``Q`` / ``Cell_Mass``)、 GROMACS は
**時間**で持つ (``tau_t`` / ``tau_p``)。 この橋渡しを間違えても .mdp は
文法的に正しく、 grompp も通り、 短い NVT なら結果もそれらしく出る。

ここで押さえるのは 4 点:

1. ``tau_t`` は **系のサイズに依存しない**。 ``Q = g k_B T tau^2`` の ``g``
   で割り忘れると ``sqrt(3N)`` 倍に膨らむ (J-OCTA の ``Export_GROMACS.py``
   がこれで、 3050 原子で 5.48 ps。 正しくは 0.63 ps)
2. GROMACS の Nose-Hoover の ``tau_t`` は緩和時間ではなく**振動の周期**な
   ので ``2*pi`` が要る
3. ``tau_p`` は **UDF が時間で持っているときだけ変換する** (Berendsen)。
   ``Cell_Mass`` は質量なので変換しない
4. ``--barostat`` は **UDF が NVT でも効く**。 C-rescale は COGNAC に対応
   概念が無く、 指定でしか選べないため
"""
from __future__ import annotations

import logging
import math

import pytest

from abmptools.udf2gro.udf_adapter import (
    BAROSTAT_NAMES,
    DEFAULT_TAU_P_PS,
    UdfAdapter,
    canonical_barostat,
)

_KB = 0.0083144626      # amu nm^2 / (ps^2 K)


class _StubUDF:
    """パスと値の辞書。 ``get`` の第 2 引数 (単位) は無視する。"""

    def __init__(self, values):
        self._values = values

    def get(self, path, *_a, **_kw):
        return self._values.get(path)


class _Cell:
    def __init__(self, a=3.0, b=3.0, c=3.0):
        self.a, self.b, self.c = a, b, c


def _adapter(values, tau_t=None, tau_p=None, barostat=None):
    """__init__ を通さずに、 抽出メソッドだけ試せる最小の adapter。"""
    a = UdfAdapter.__new__(UdfAdapter)
    a._udf = _StubUDF(values)
    a._tau_t_override = tau_t
    a._tau_p_override = tau_p
    a._barostat_override = barostat
    return a


def _nvt_values(n_atoms, tau_ps, T=300.0, comm_linear=False):
    """狙った ``tau`` がちょうど出る ``Q`` を組み立てる。"""
    g = 3 * n_atoms - (3 if comm_linear else 0)
    Q = g * _KB * T * tau_ps ** 2            # unit_Mass = unit_L = 1 とする
    v = {
        "Simulation_Conditions.Dynamics_Conditions.Temperature.Temperature": T,
        "Unit_Parameter.Length": 1.0,
        "Unit_Parameter.Mass": 1.0,
        "Simulation_Conditions.Solver.Dynamics.NVT_Nose_Hoover.Q": Q,
    }
    if comm_linear:
        base = "Simulation_Conditions.Dynamics_Conditions.Moment."
        v[base + "Calc_Moment"] = 1
        v[base + "Stop_Translation"] = 1
    return v


# ---------------------------------------------------------------------------
# tau_t
# ---------------------------------------------------------------------------

def test_tau_t_is_two_pi_times_the_relaxation_time():
    """GROMACS の Nose-Hoover の tau_t は周期。 2*pi を落とすと 6.28 倍ずれる。"""
    a = _adapter(_nvt_values(1000, 0.1))
    _, tau_t, _ = a._extract_temperature("NVT_Nose_Hoover", 1000)

    assert tau_t == pytest.approx(2.0 * math.pi * 0.1, rel=1e-9)


def test_tau_t_does_not_change_with_system_size():
    """同じ tau なら 100 原子でも 10000 原子でも同じ値。

    ``g`` で割り忘れた式は ``sqrt(3N)`` 倍になるので、 ここで 10 倍ずれる。
    """
    small = _adapter(_nvt_values(100, 0.2))
    large = _adapter(_nvt_values(10000, 0.2))

    _, t_small, _ = small._extract_temperature("NVT_Nose_Hoover", 100)
    _, t_large, _ = large._extract_temperature("NVT_Nose_Hoover", 10000)

    assert t_small == pytest.approx(t_large, rel=1e-9)


def test_the_jocta_formula_would_differ_and_grow_with_n():
    """``2*pi*sqrt(Q_d/T)`` (g*k_B 抜け) との差を明示的に残す。"""
    n, tau, T = 3050, 0.1, 300.0
    values = _nvt_values(n, tau, T=T)
    Q = values["Simulation_Conditions.Solver.Dynamics.NVT_Nose_Hoover.Q"]

    _, ours, _ = _adapter(values)._extract_temperature("NVT_Nose_Hoover", n)
    jocta = 2.0 * math.pi * math.sqrt(Q / T)

    assert ours == pytest.approx(2.0 * math.pi * tau, rel=1e-9)
    # 比は sqrt(3N k_B) = 8.7。 実機 (3050 原子) の 5.48 ps と一致する
    assert jocta == pytest.approx(ours * math.sqrt(3 * n * _KB), rel=1e-9)
    assert jocta == pytest.approx(5.48, abs=0.01)
    assert ours == pytest.approx(0.628, abs=0.001)

    # そして N とともに増える。 これが「系を大きくすると熱浴が鈍る」の正体
    small = _nvt_values(100, tau, T=T)
    q_small = small["Simulation_Conditions.Solver.Dynamics.NVT_Nose_Hoover.Q"]
    assert 2.0 * math.pi * math.sqrt(q_small / T) < jocta


def test_stopping_com_motion_uses_3n_minus_3():
    """comm-mode = Linear なら自由度が 3 減る。"""
    n = 500
    free = _adapter(_nvt_values(n, 0.1, comm_linear=False))
    fixed_q = free._udf._values[
        "Simulation_Conditions.Solver.Dynamics.NVT_Nose_Hoover.Q"]

    # 同じ Q を comm-mode = Linear の UDF に置くと、 g が小さいぶん tau_t は上がる
    linear_values = _nvt_values(n, 0.1, comm_linear=True)
    linear_values[
        "Simulation_Conditions.Solver.Dynamics.NVT_Nose_Hoover.Q"] = fixed_q
    linear = _adapter(linear_values)

    _, t_free, _ = free._extract_temperature("NVT_Nose_Hoover", n)
    _, t_linear, _ = linear._extract_temperature("NVT_Nose_Hoover", n)

    assert t_linear == pytest.approx(
        t_free * math.sqrt(3 * n / float(3 * n - 3)), rel=1e-9)


def test_tau_t_override_wins():
    a = _adapter(_nvt_values(1000, 0.1), tau_t=1.5)
    _, tau_t, _ = a._extract_temperature("NVT_Nose_Hoover", 1000)
    assert tau_t == 1.5


def test_a_missing_q_does_not_crash_the_conversion():
    values = _nvt_values(1000, 0.1)
    values["Simulation_Conditions.Solver.Dynamics.NVT_Nose_Hoover.Q"] = 0.0
    _, tau_t, _ = _adapter(values)._extract_temperature("NVT_Nose_Hoover", 1000)
    assert tau_t > 0.0


# ---------------------------------------------------------------------------
# tau_p
# ---------------------------------------------------------------------------

def _npt_values(cell_mass=1.0e4):
    v = {
        "Unit_Parameter.Length": 1.0,
        "Unit_Parameter.Mass": 1.0,
        "Simulation_Conditions.Dynamics_Conditions.Pressure_Stress.Pressure": 1.0,
        "Simulation_Conditions.Solver.Dynamics."
        "NPT_Andersen_Nose_Hoover.Cell_Mass": cell_mass,
    }
    for e in ("xx", "yy", "zz", "yz", "zx", "xy"):
        v["Simulation_Conditions.Dynamics_Conditions."
          "Pressure_Stress.Stress." + e] = 0.0
    return v


def _pressure(adapter, algorithm="NPT_Andersen_Nose_Hoover"):
    return adapter._extract_pressure(algorithm, "", 0, False, None, _Cell())


def test_cell_mass_is_not_turned_into_tau_p():
    """Cell_Mass は [mass]。 時間に読み替えると系ごとに違う tau_p が出る。"""
    light = _pressure(_adapter(_npt_values(cell_mass=1.0e3)))
    heavy = _pressure(_adapter(_npt_values(cell_mass=1.0e6)))

    assert light[2] == DEFAULT_TAU_P_PS
    assert heavy[2] == DEFAULT_TAU_P_PS


def test_tau_p_override_wins():
    out = _pressure(_adapter(_npt_values(), tau_p=5.0))
    assert out[2] == 5.0


def test_the_cell_mass_reference_value_is_only_logged(caplog):
    """採用はしないが、 元の値を知りたいことはあるので INFO に残す。"""
    with caplog.at_level(logging.INFO, logger="abmptools.udf2gro.udf_adapter"):
        _pressure(_adapter(_npt_values()))
    assert "Cell_Mass" in caplog.text


def test_berendsen_tau_p_is_converted_because_the_udf_holds_a_time():
    """COGNAC の NPT_Berendsen は tau_P を時間で持つ。 ここは単位換算するだけ。"""
    v = _npt_values()
    v["Simulation_Conditions.Dynamics_Conditions.Pressure_Stress.Pressure"] = 1.0
    v["Simulation_Conditions.Solver.Dynamics.NPT_Berendsen.tau_P"] = 2.0e5
    a = _adapter(v)
    # unit_P = Pressure[bar] / Pressure[P] = 1.0 になるよう、 同じ値を返す stub
    out = a._extract_pressure("NPT_Berendsen", "", 0, False, None, _Cell())

    assert out[0] == "berendsen"
    assert out[2] == pytest.approx(2.0e5 * 0.000045, rel=1e-9)
    assert out[2] != DEFAULT_TAU_P_PS


def test_tau_p_override_beats_the_berendsen_conversion():
    v = _npt_values()
    v["Simulation_Conditions.Solver.Dynamics.NPT_Berendsen.tau_P"] = 2.0e5
    out = _adapter(v, tau_p=3.0)._extract_pressure(
        "NPT_Berendsen", "", 0, False, None, _Cell())
    assert out[2] == 3.0


# ---------------------------------------------------------------------------
# --barostat
# ---------------------------------------------------------------------------

def test_barostat_applies_to_an_nvt_udf():
    """C-rescale は COGNAC に無いので、 NVT の UDF から始めるしかない。

    以前は ``pcoupl == "no"`` の早期 return が上書きより前にあり、 指定が
    黙って捨てられて NVT の .mdp が出ていた。
    """
    out = _adapter(_nvt_values(100, 0.1),
                   barostat="C-rescale")._extract_pressure(
        "NVT_Nose_Hoover", "", 0, False, None, _Cell())

    assert out[0] == "C-rescale"
    assert out[1] == "isotropic"
    assert out[2] == DEFAULT_TAU_P_PS


def test_changing_the_ensemble_is_announced(caplog):
    with caplog.at_level(logging.WARNING,
                         logger="abmptools.udf2gro.udf_adapter"):
        _adapter(_nvt_values(100, 0.1),
                 barostat="C-rescale")._extract_pressure(
            "NVT_Nose_Hoover", "", 0, False, None, _Cell())
    assert "NPT" in caplog.text


def test_barostat_replaces_the_one_the_udf_asks_for():
    out = _pressure(_adapter(_npt_values(), barostat="C-rescale"))
    assert out[0] == "C-rescale"


def test_barostat_no_turns_pressure_coupling_off():
    out = _pressure(_adapter(_npt_values(), barostat="no"))
    assert out[0] == "no"


def test_an_overridden_barostat_is_isotropic():
    """Parrinello-Rahman の anisotropic は J-OCTA 由来の既定。 指定時は外す。"""
    out = _adapter(_npt_values(), barostat="Parrinello-Rahman")._extract_pressure(
        "NPT_Parrinello_Rahman_Nose_Hoover", "", 0, False, None, _Cell())
    assert out[1] == "isotropic"


def test_names_are_normalised_to_the_gromacs_spelling():
    assert canonical_barostat("c-rescale") == "C-rescale"
    assert canonical_barostat("  PARRINELLO_RAHMAN ") == "Parrinello-Rahman"
    assert canonical_barostat("berendsen") == "Berendsen"


def test_a_misspelled_barostat_is_refused_not_passed_through():
    """綴り違いをそのまま .mdp に書くと grompp で初めて落ちる。"""
    with pytest.raises(ValueError):
        canonical_barostat("crescal")


def test_every_name_maps_to_something_gromacs_accepts():
    accepted = {"C-rescale", "Parrinello-Rahman", "Berendsen", "MTTK", "no"}
    assert set(BAROSTAT_NAMES.values()) <= accepted
