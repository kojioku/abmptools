# -*- coding: utf-8 -*-
"""熱浴・圧力浴の変換。 どれも「変換は通るが値が間違っている」類の欠陥。

COGNAC は熱浴・圧力浴を**質量**で持ち (``Q`` / ``Cell_Mass``)、 GROMACS は
**時間**で持つ (``tau_t`` / ``tau_p``)。 この橋渡しを間違えても .mdp は
文法的に正しく、 grompp も通り、 短い NVT なら結果もそれらしく出る。

ここで押さえるのは 4 点:

1. ``tau_t`` は **系のサイズに依存しない**。 ``Q = g k_B T tau^2`` の ``g``
   で割り忘れると ``sqrt(3N)`` 倍に膨らむ (下流の ``下流の GROMACS 変換器``
   がこれで、 3050 原子で 5.48 ps。 正しくは 0.63 ps)
2. GROMACS の Nose-Hoover の ``tau_t`` は緩和時間ではなく**振動の周期**な
   ので ``2*pi`` が要る
3. ``tau_p`` も **UDF から変換する**。 Berendsen は時間なので単位換算、
   Andersen / PR は ``Cell_Mass`` から運動方程式経由
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
    tau_t = a._build_thermostat("NVT_Nose_Hoover", 1000).tau_t

    assert tau_t == pytest.approx(2.0 * math.pi * 0.1, rel=1e-9)


def test_tau_t_does_not_change_with_system_size():
    """同じ tau なら 100 原子でも 10000 原子でも同じ値。

    ``g`` で割り忘れた式は ``sqrt(3N)`` 倍になるので、 ここで 10 倍ずれる。
    """
    small = _adapter(_nvt_values(100, 0.2))
    large = _adapter(_nvt_values(10000, 0.2))

    t_small = small._build_thermostat("NVT_Nose_Hoover", 100).tau_t
    t_large = large._build_thermostat("NVT_Nose_Hoover", 10000).tau_t

    assert t_small == pytest.approx(t_large, rel=1e-9)


def test_dropping_g_kb_would_differ_and_grow_with_n():
    """``2*pi*sqrt(Q_d/T)`` (g*k_B 抜け) との差を明示的に残す。"""
    n, tau, T = 3050, 0.1, 300.0
    values = _nvt_values(n, tau, T=T)
    Q = values["Simulation_Conditions.Solver.Dynamics.NVT_Nose_Hoover.Q"]

    ours = _adapter(values)._build_thermostat("NVT_Nose_Hoover", n).tau_t
    nodof = 2.0 * math.pi * math.sqrt(Q / T)

    assert ours == pytest.approx(2.0 * math.pi * tau, rel=1e-9)
    # 比は sqrt(3N k_B) = 8.7。 実機 (3050 原子) の 5.48 ps と一致する
    assert nodof == pytest.approx(ours * math.sqrt(3 * n * _KB), rel=1e-9)
    assert nodof == pytest.approx(5.48, abs=0.01)
    assert ours == pytest.approx(0.628, abs=0.001)

    # そして N とともに増える。 これが「系を大きくすると熱浴が鈍る」の正体
    small = _nvt_values(100, tau, T=T)
    q_small = small["Simulation_Conditions.Solver.Dynamics.NVT_Nose_Hoover.Q"]
    assert 2.0 * math.pi * math.sqrt(q_small / T) < nodof


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

    t_free = free._build_thermostat("NVT_Nose_Hoover", n).tau_t
    t_linear = linear._build_thermostat("NVT_Nose_Hoover", n).tau_t

    assert t_linear == pytest.approx(
        t_free * math.sqrt(3 * n / float(3 * n - 3)), rel=1e-9)


def test_tau_t_override_wins():
    a = _adapter(_nvt_values(1000, 0.1), tau_t=1.5)
    tau_t = a._build_thermostat("NVT_Nose_Hoover", 1000).tau_t
    assert tau_t == 1.5


def test_a_missing_q_does_not_crash_the_conversion():
    values = _nvt_values(1000, 0.1)
    values["Simulation_Conditions.Solver.Dynamics.NVT_Nose_Hoover.Q"] = 0.0
    tau_t = _adapter(values)._build_thermostat("NVT_Nose_Hoover", 1000).tau_t
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


def _pr_values(cell_mass):
    v = _npt_values(cell_mass=cell_mass)
    v["Simulation_Conditions.Solver.Dynamics."
      "NPT_Parrinello_Rahman_Nose_Hoover.Cell_Mass"] = cell_mass
    return v


def _pressure(adapter, algorithm="NPT_Andersen_Nose_Hoover"):
    return adapter._build_barostat(algorithm, "", 0, False, None, _Cell())


def _expected_pr_tau_p(C, max_L):
    beta = 0.000045 / 0.06022140762          # nm^3 mol / kJ
    return 2.0 * math.pi * math.sqrt(C * beta / (3.0 * max_L))


def test_cell_mass_is_converted_through_the_equation_of_motion():
    """PR: ``tau_p = 2*pi*sqrt(C*beta/(3*L))``、 ``L`` は最長のセル辺。

    GROMACS は ``W^-1 = 4*pi^2*beta/(3*tau_p^2*L)`` (``L`` = 最長の箱要素) で
    バロスタット質量を決めるので、 ``tau_p`` は箱の振動の周期そのもの。
    COGNAC PR は ``W h_ddot = dP h^-1 V`` (``PRsystem.cpp:7,58``)。 周期どうし
    を等置すればこの式になる。
    """
    C = 1.0e4
    out = _adapter(_pr_values(C))._build_barostat(
        "NPT_Parrinello_Rahman_Nose_Hoover", "", 0, False, None, _Cell())

    assert out.tau_p == pytest.approx(_expected_pr_tau_p(C, 3.0), rel=1e-9)


def test_andersen_is_slower_than_pr_by_sqrt3():
    """Andersen は体積を座標にしていて ``W = C*V^(-4/3)`` (``Anphsystem.cpp:17``)。

    ``dL`` で書き直すと ``omega^2 = L/(C*beta)``、 PR の ``3L/(C*beta)`` に
    対して 1/3 なので、 周期は ``sqrt(3)`` 倍。 同じ周期を PR で出したければ
    ``Cell_Mass`` を **3 倍**する。
    """
    C = 1.0e4
    andersen = _pressure(_adapter(_npt_values(cell_mass=C))).tau_p
    pr = _adapter(_pr_values(C))._build_barostat(
        "NPT_Parrinello_Rahman_Nose_Hoover", "", 0, False, None, _Cell()).tau_p

    assert andersen == pytest.approx(pr * math.sqrt(3.0), rel=1e-9)


def test_the_longest_edge_is_used_not_the_cube_root():
    """GROMACS の ``W^-1`` は最長の箱要素で定義されている。

    立方体では区別がつかないので、 扁平なセルで固定する
    (2.3 x 2.3 x 5.2 では ``V^(1/3)`` = 3.02 に対し ``max_L`` = 5.20)。
    """
    C = 1.0e4
    flat = _Cell(2.29911, 2.29911, 5.20238)
    out = _adapter(_pr_values(C))._build_barostat(
        "NPT_Parrinello_Rahman_Nose_Hoover", "", 0, False, None, flat)

    assert out.tau_p == pytest.approx(_expected_pr_tau_p(C, 5.20238), rel=1e-9)
    assert out.tau_p != pytest.approx(_expected_pr_tau_p(C, 3.0180), rel=1e-3)


def test_a_heavier_cell_gives_a_slower_barostat():
    """質量を捨てて既定を返していた頃は、 ここが同じ値になっていた。"""
    light = _pressure(_adapter(_npt_values(cell_mass=1.0e3))).tau_p
    heavy = _pressure(_adapter(_npt_values(cell_mass=1.0e5))).tau_p

    assert heavy == pytest.approx(light * 10.0, rel=1e-9)


def test_sigma_squared_is_not_applied_to_cell_mass():
    """``Cell_Mass`` は [mass]。 ``Q`` 用の ``unit_L^2`` を掛けると 100 倍ずれる。

    all-atom の UDF は ``Unit_Parameter.Length = 0.1`` (Å) なので、 旧実装は
    ここで 0.01 倍になり、 下限 2.0 ps に丸められて値が見えなくなっていた。
    """
    v = _npt_values(cell_mass=1.0e4)
    v["Unit_Parameter.Length"] = 0.1        # Å 系
    out = _pressure(_adapter(v))

    # 長さの単位は体積 (cell) 側で既に nm になっているので、 tau_p は変わらない
    assert out.tau_p == pytest.approx(
        _pressure(_adapter(_npt_values(cell_mass=1.0e4))).tau_p, rel=1e-9)


def test_the_real_system_reproduces_the_measured_values():
    """実機 (3050 原子、 Cell_Mass = 14076.4 amu、 2.3 x 2.3 x 5.2 nm)。"""
    cell = _Cell(2.29911, 2.29911, 5.20238)

    andersen = _adapter(_npt_values(cell_mass=14076.4))._build_barostat(
        "NPT_Andersen_Nose_Hoover", "", 0, False, None, cell).tau_p
    pr = _adapter(_pr_values(14076.4))._build_barostat(
        "NPT_Parrinello_Rahman_Nose_Hoover", "", 0, False, None, cell).tau_p

    assert andersen == pytest.approx(8.93, abs=0.01)
    assert pr == pytest.approx(5.16, abs=0.01)


def test_a_missing_cell_mass_falls_back_and_says_so(caplog):
    v = _npt_values(cell_mass=0.0)
    with caplog.at_level(logging.WARNING,
                         logger="abmptools.udf2gro.udf_adapter"):
        out = _pressure(_adapter(v))

    assert out.tau_p == DEFAULT_TAU_P_PS
    assert "Cell_Mass" in caplog.text


def test_an_implausible_tau_p_is_flagged_but_kept(caplog):
    """Cell_Mass が既定 (系の全質量) のままだと大きな系で数十 ps に伸びる。"""
    with caplog.at_level(logging.WARNING,
                         logger="abmptools.udf2gro.udf_adapter"):
        out = _pressure(_adapter(_npt_values(cell_mass=1.0e6)))

    assert out.tau_p > 20.0
    assert "--tau-p" in caplog.text


def test_tau_p_override_wins():
    out = _pressure(_adapter(_npt_values(), tau_p=5.0))
    assert out.tau_p == 5.0


def test_berendsen_tau_p_is_converted_because_the_udf_holds_a_time():
    """COGNAC の NPT_Berendsen は tau_P を時間で持つ。 ここは単位換算するだけ。"""
    v = _npt_values()
    v["Simulation_Conditions.Dynamics_Conditions.Pressure_Stress.Pressure"] = 1.0
    v["Simulation_Conditions.Solver.Dynamics.NPT_Berendsen.tau_P"] = 2.0e5
    a = _adapter(v)
    # unit_P = Pressure[bar] / Pressure[P] = 1.0 になるよう、 同じ値を返す stub
    out = a._build_barostat("NPT_Berendsen", "", 0, False, None, _Cell())

    assert out.p_coupl == "berendsen"
    assert out.tau_p == pytest.approx(2.0e5 * 0.000045, rel=1e-9)
    assert out.tau_p != DEFAULT_TAU_P_PS


def test_tau_p_override_beats_the_berendsen_conversion():
    v = _npt_values()
    v["Simulation_Conditions.Solver.Dynamics.NPT_Berendsen.tau_P"] = 2.0e5
    out = _adapter(v, tau_p=3.0)._build_barostat(
        "NPT_Berendsen", "", 0, False, None, _Cell())
    assert out.tau_p == 3.0


# ---------------------------------------------------------------------------
# --barostat
# ---------------------------------------------------------------------------

def test_barostat_applies_to_an_nvt_udf():
    """C-rescale は COGNAC に無いので、 NVT の UDF から始めるしかない。

    以前は ``pcoupl == "no"`` の早期 return が上書きより前にあり、 指定が
    黙って捨てられて NVT の .mdp が出ていた。
    """
    out = _adapter(_nvt_values(100, 0.1),
                   barostat="C-rescale")._build_barostat(
        "NVT_Nose_Hoover", "", 0, False, None, _Cell())

    assert out.p_coupl == "C-rescale"
    assert out.pcoupltype == "isotropic"
    assert out.tau_p == DEFAULT_TAU_P_PS


def test_changing_the_ensemble_is_announced(caplog):
    with caplog.at_level(logging.WARNING,
                         logger="abmptools.udf2gro.udf_adapter"):
        _adapter(_nvt_values(100, 0.1),
                 barostat="C-rescale")._build_barostat(
            "NVT_Nose_Hoover", "", 0, False, None, _Cell())
    assert "NPT" in caplog.text


def test_barostat_replaces_the_one_the_udf_asks_for():
    out = _pressure(_adapter(_npt_values(), barostat="C-rescale"))
    assert out.p_coupl == "C-rescale"


def test_barostat_no_turns_pressure_coupling_off():
    out = _pressure(_adapter(_npt_values(), barostat="no"))
    assert out.p_coupl == "no"


def test_an_overridden_barostat_is_isotropic():
    """Parrinello-Rahman の anisotropic は 下流の変換器 由来の既定。 指定時は外す。"""
    out = _adapter(_npt_values(), barostat="Parrinello-Rahman")._build_barostat(
        "NPT_Parrinello_Rahman_Nose_Hoover", "", 0, False, None, _Cell())
    assert out.pcoupltype == "isotropic"


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


# ---------------------------------------------------------------------------
# tau_t と tau_p の組合せ
#
# UDF から出した値は「その UDF が記述している計算」であって、 GROMACS で
# 普通に選ぶ値とは限らない。 grompp が言うことは grompp が言うが、 実務の
# 範囲から外れていることは誰も言わないので、 ここで言う。
# ---------------------------------------------------------------------------

from abmptools.udf2gro.udf_adapter import (          # noqa: E402
    PRACTICAL_TAU_P_PS,
    _check_tau_pair,
)


def test_a_slow_barostat_is_flagged(caplog):
    """Cell_Mass = 全質量 だと tau_p は箱の一辺に比例して伸びる。

    5 nm 立方で PR 11 ps / Andersen 19 ps。 値としては UDF に忠実だが、
    密度の緩和には数 x tau_p かかるので、 短い NPT では効かない。
    """
    with caplog.at_level(logging.WARNING,
                         logger="abmptools.udf2gro.udf_adapter"):
        _check_tau_pair("nose-hoover", 0.628, "Parrinello-Rahman", 8.93)

    assert "--tau-p" in caplog.text
    assert "2-5 ps" in caplog.text


def test_a_barostat_in_the_usual_range_is_quiet(caplog):
    with caplog.at_level(logging.WARNING,
                         logger="abmptools.udf2gro.udf_adapter"):
        _check_tau_pair("nose-hoover", 0.628, "Parrinello-Rahman", 3.0)

    assert caplog.text == ""


def test_the_resonance_condition_is_checked_before_grompp(caplog):
    """``tau_p < 2*tau_t`` は grompp も言うが、 こちらは両方の値を持っている。"""
    with caplog.at_level(logging.WARNING,
                         logger="abmptools.udf2gro.udf_adapter"):
        _check_tau_pair("nose-hoover", 2.0, "Parrinello-Rahman", 3.0)

    assert "twice" in caplog.text


def test_nothing_is_said_when_there_is_no_barostat(caplog):
    with caplog.at_level(logging.WARNING,
                         logger="abmptools.udf2gro.udf_adapter"):
        _check_tau_pair("nose-hoover", 5.0, "no", 2.0)

    assert caplog.text == ""


def test_the_practical_range_matches_what_the_docs_say():
    assert PRACTICAL_TAU_P_PS == (2.0, 5.0)
