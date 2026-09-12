# -*- coding: utf-8 -*-
"""Nose-Hoover の Q と Cell_Mass は J-OCTA の規約に合わせる。

`Export_GROMACS.py` はアルゴリズムごとに**別のフィールド**から `Q` を読み、
`tau_t = 2*pi*sqrt(Q*Mass*Length^2/T)` を作る。`NVT_Nose_Hoover.Q` しか
書かないと、J-OCTA 側で NPT に切り替えた瞬間に 0 が読まれ **tau_t = 0** に
なり、Nose-Hoover が 0 除算して箱が nan に飛ぶ (2026-09-12 実機)。
NVT で露見しないのはこのため。

`Cell_Mass` (バロスタット質量) も同様に未記入だと `tau_p` の式が 0 になり、
下限の 2.0 に丸められる。J-OCTA は `CognacSystemUtil.setCellMass` で
**系の全質量**を入れる。

## 自由度を 3N にする根拠

J-OCTA が書いた UDF 4 系 (原子数 23 / 80 / 110 / 3050) の `Q` は、いずれも
`Q = 3N * kB * T * (0.1 ps)^2` にぴったり乗る (逆算した tau は 4 系とも
0.099999972 ps)。**3N-3 では乗らない** —— 80 原子で 0.1006 とずれる。

物理的にも、J-OCTA が書き出す mdp は `comm-mode = None` なので重心運動が
除かれず、3N が正しい (GROMACS も "degrees of freedom ... is 9150.00" =
3*3050 と報告した)。`comm-mode = Linear` で流すなら 3N-3 が正しいので、
そちらは `--nh-dof 3N-3` で選ぶ。
"""
from __future__ import annotations

import math

import pytest

KB = 0.83144626      # amu A^2 / (ps^2 K) — top_exporter の KB_AMU_A2_PS2_K


def _q(n_atoms, dof, tau=0.1, T=300.0):
    g = 3 * n_atoms - 3 if dof == "3N-3" else 3 * n_atoms
    return g * KB * T * tau ** 2


def _tau_t(q, mass=1.0, length=0.1, T=300.0):
    """J-OCTA が Q から tau_t を戻す式。"""
    return 2.0 * math.pi * math.sqrt(q * mass * length * length / T)


# J-OCTA が実際に書いた 4 系 (原子数, Q)。regression の実測値。
JOCTA_FILES = [
    (23, 172.10927951148415),        # SI/udfcheck/testallatom.bdf
    (80, 598.64097),                 # udf_Q/sys2.bdf
    (110, 823.13134),                # udf_Q/sys1.bdf と sys3.bdf
    (3050, 22823.1870656533),        # 力場取り直し後
]


@pytest.mark.parametrize("n_atoms,q_jocta", JOCTA_FILES)
def test_3n_reproduces_what_jocta_writes(n_atoms, q_jocta):
    assert _q(n_atoms, "3N") == pytest.approx(q_jocta, rel=2e-5)


@pytest.mark.parametrize("n_atoms,q_jocta", JOCTA_FILES)
def test_3n_minus_3_does_not(n_atoms, q_jocta):
    """80 原子だと 3.8% ずれる。ここで規約を判別できた。"""
    if n_atoms > 1000:
        pytest.skip("大きい系では差が 0.03% で判別に使えない")
    assert _q(n_atoms, "3N-3") != pytest.approx(q_jocta, rel=2e-5)


def test_the_two_modes_differ_by_three_degrees_of_freedom():
    assert _q(80, "3N") / _q(80, "3N-3") == pytest.approx(240.0 / 237.0)


def test_the_difference_is_small_for_a_large_system():
    """3050 原子では 0.03%。物理的にはほぼ無意味な差。"""
    ratio = _q(3050, "3N") / _q(3050, "3N-3")
    assert 1.0 < ratio < 1.0005


def test_tau_t_is_not_zero_when_q_is_written():
    """これが 0 になると Nose-Hoover が 0 除算して箱が nan に飛ぶ。"""
    assert _tau_t(_q(3050, "3N")) == pytest.approx(5.4803, abs=1e-3)


def test_tau_t_is_zero_when_q_is_missing():
    """NPT 側の Q を書かなかったときに起きていたこと。"""
    assert _tau_t(0.0) == 0.0


# ---------------------------------------------------------------------------
# comm-mode の連動
# ---------------------------------------------------------------------------

class _MomentUDF:
    """Moment まわりの put だけ覚える stand-in。"""

    def __init__(self):
        self.puts = {}

    def jump(self, _r):
        pass

    def put(self, value, path, *_a, **_kw):
        self.puts[path] = value

    def get(self, *_a, **_kw):
        return None

    def size(self, *_a, **_kw):
        return 0


def _moment(nh_dof):
    """_set_default_condition が Moment に何を書くかだけ取り出す。"""
    from abmptools.gro2udf.top_exporter import TopExporter
    from abmptools.gro2udf.top_model import TopModel
    u = _MomentUDF()
    model = TopModel(
        comb_rule=2, fudge_lj=0.5, fudge_qq=0.8333,
        atom_type_specs=[], bond_type_specs=[], angle_type_specs=[],
        torsion_type_specs=[], mass_dict={}, mol_specs=[],
        mol_type_names=[], mol_instance_list=[], n_atoms_total=100,
    )
    try:
        TopExporter._set_default_condition(u, model, nh_dof=nh_dof)
    except Exception:                                    # noqa: BLE001
        pass          # Moment 以外の put で落ちても、そこまでの記録は使える
    pre = "Simulation_Conditions.Dynamics_Conditions.Moment."
    return {k[len(pre):]: v for k, v in u.puts.items() if k.startswith(pre)}


def test_3n_leaves_centre_of_mass_motion_alone():
    """comm-mode = None。 J-OCTA の既定と同じで、3N と整合する。"""
    m = _moment("3N")
    assert m.get("Calc_Moment") == 0 and m.get("Stop_Translation") == 0


def test_3n_minus_3_turns_on_removal():
    """comm-mode = Linear。 これを書かないと GROMACS は 3N で積分してしまい、
    --nh-dof 3N-3 が説明どおりに動かない。"""
    m = _moment("3N-3")
    assert m.get("Calc_Moment") == 1 and m.get("Stop_Translation") == 1


def test_rotation_is_never_stopped():
    """GROMACS は回転のみの除去を受け付けない (J-OCTA が例外を投げる)。"""
    for dof in ("3N", "3N-3"):
        assert _moment(dof).get("Stop_Rotation") == 0
