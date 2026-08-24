# -*- coding: utf-8 -*-
"""`2_optmask-frame.sh --check` が「途中まで」を見落とさないこと。

途中で落ちた実行を再開するとき、危ないのは**完成していないものを完成と
みなす**ことである。判定の根拠がファイルの存在だけだと、

- `mkdir` は作業の前に走るので、**ディレクトリがあっても中身は無い**
- 殺されて切り詰められた `.gro` / `.pdb` も「存在する」し「空でない」

の 2 つを取り違える。スクリプトは完了マーカーだけを済んだ証拠として使い、
マーカーを書く前に中身を検める。ここではその分類を、gmx も cpptraj も
使わずに確かめる (`--check` は何も実行せず状態を出すだけなので、
prmtop も mdp も要らない)。
"""
from __future__ import annotations

import os
import shutil
import subprocess

import pytest

SCRIPT = os.path.join(
    os.path.dirname(__file__), os.pardir,
    "tips", "md", "gmx", "md-fmo", "2_optmask-frame.sh",
)

pytestmark = pytest.mark.skipif(
    shutil.which("bash") is None, reason="bash が無い")


def _complete_gro(path, natoms=2):
    """2 行目の原子数と行数 (N+3) が合い、箱ベクトルで終わる gro。"""
    lines = ["toy", f"{natoms:5d}"]
    for i in range(natoms):
        lines.append(
            f"    1SOL     OW{i + 1:5d}   1.000   1.000   1.000")
    lines.append("   3.00000   3.00000   3.00000")
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def _truncated_gro(path, natoms=100):
    """原子数は 100 と宣言しているのに 2 行しか書けていない gro。

    途中で殺された gro はこうなる。空ではないので `-s` は通ってしまう。
    """
    with open(path, "w") as fh:
        fh.write("toy\n%5d\n" % natoms)
        fh.write("    1SOL     OW    1   1.000   1.000   1.000\n")


def _complete_pdb(path):
    with open(path, "w") as fh:
        fh.write(
            "ATOM      1  OW  SOL     1      10.000  10.000  10.000  1.00  0.00\n"
            "END\n")


def _truncated_pdb(path):
    """座標欄の途中で切れた原子行で終わる pdb。"""
    with open(path, "w") as fh:
        fh.write(
            "ATOM      1  OW  SOL     1      10.000  10.000  10.000  1.00  0.00\n"
            "ATOM      2  HW  SOL     1      10.5\n")


def _frame_dir(root, i):
    d = root / f"traj_{i}_fmo"
    d.mkdir()
    return d


def _run_check(root):
    proc = subprocess.run(
        ["bash", SCRIPT, "-n", "index.ndx", "-p", "system.top",
         "-f", "traj.xtc", "--check"],
        cwd=str(root), capture_output=True, text=True,
    )
    assert proc.returncode == 0, proc.stderr
    return proc.stdout


def _status_of(out, frame):
    """--check の出力から 1 フレームの (minimize, arrange) を返す。"""
    for line in out.splitlines():
        if not line.startswith(f"frame {frame} |"):
            continue
        _, m, a = (part.strip() for part in line.split("|"))
        return m.split(":", 1)[1].strip(), a.split(":", 1)[1].strip()
    raise AssertionError(f"frame {frame} の行が出力に無い:\n{out}")


@pytest.fixture
def scene(tmp_path):
    """フレーム 0-4 に、再開時に起こりうる 5 通りの状態を作る。"""
    # 0: 両方とも完了マーカーあり
    d = _frame_dir(tmp_path, 0)
    _complete_gro(d / "traj_0_fmo.gro")
    _complete_pdb(d / "traj_0_fmo_mask.pdb")
    (d / ".done.minimize").write_text("x")
    (d / ".done.arrange").write_text("x")

    # 1: minimize は済み、液滴はあるがマーカーが無い (arrange の直後で落ちた)
    d = _frame_dir(tmp_path, 1)
    _complete_gro(d / "traj_1_fmo.gro")
    (d / ".done.minimize").write_text("x")
    _complete_pdb(d / "traj_1_fmo_mask.pdb")

    # 2: gro が書きかけで切れている (mdrun の最中に殺された)
    d = _frame_dir(tmp_path, 2)
    _truncated_gro(d / "traj_2_fmo.gro")

    # 3: ディレクトリだけ (mkdir の直後に落ちた)
    _frame_dir(tmp_path, 3)

    # 4: 何も無い
    return tmp_path


def test_a_directory_alone_is_not_progress(scene):
    """mkdir は作業の前に走る。ディレクトリの存在を済んだ証拠にしない。"""
    m, a = _status_of(_run_check(scene), 3)
    assert m.startswith("未") and "ディレクトリのみ" in m
    assert a.startswith("未")


def test_nothing_at_all_is_reported_as_未着手(scene):
    m, a = _status_of(_run_check(scene), 4)
    assert m.startswith("未") and a.startswith("未")


def test_a_truncated_gro_is_not_mistaken_for_done(scene):
    """空ではないので `-s` は通る。行数まで見て初めて途中だと分かる。"""
    m, _ = _status_of(_run_check(scene), 2)
    assert "途中" in m and "壊れている" in m


def test_a_finished_stage_is_skipped(scene):
    m, a = _status_of(_run_check(scene), 0)
    assert m == "済" and a == "済"


def test_an_artifact_without_a_marker_is_half_done(scene):
    """成果物は完全でもマーカーが無いなら、まだ済んだことにはしない。"""
    m, a = _status_of(_run_check(scene), 1)
    assert m == "済"
    assert "途中" in a and "マーカー無し" in a


def test_a_truncated_droplet_is_not_mistaken_for_done(tmp_path):
    d = _frame_dir(tmp_path, 0)
    _complete_gro(d / "traj_0_fmo.gro")
    (d / ".done.minimize").write_text("x")
    _truncated_pdb(d / "traj_0_fmo_mask.pdb")
    _, a = _status_of(_run_check(tmp_path), 0)
    assert "途中" in a and "壊れている" in a


def test_check_runs_before_the_prmtop_is_built(tmp_path):
    """--check は何も実行しない。prmtop も mdp も無い段階で状態を見たい。"""
    assert not (tmp_path / "system.prmtop").exists()
    out = _run_check(tmp_path)
    assert "frame 0 |" in out and "minimize:" in out and "arrange:" in out


def test_the_script_no_longer_forces_overwrites():
    """`mv -f` は既存を黙って消す。park() でどける方針に置き換えてある。"""
    with open(SCRIPT, encoding="utf-8") as fh:
        body = [l for l in fh if not l.lstrip().startswith("#")]
    assert not [l for l in body if "mv -f" in l]
    assert not [l for l in body if l.lstrip().startswith("rm ")]
