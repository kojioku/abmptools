# -*- coding: utf-8 -*-
"""静的セルを書かないと、座標と箱が別物の UDF が黙って出る。

``gro2udf`` はセルを**レコードにしか書かなかった**ので、静的
``Structure.Unit_Cell`` と ``Initial_Structure.Initial_Unit_Cell`` には
テンプレートの値が残っていた。同梱テンプレートのそれは **20 Å 立方**と
**100 Å 立方**で、系とは何の関係も無い既定値である。

レコード側は正しいので、変換も下流のスキーマ検証も何も言わない。
ところが静的側を読む実装に渡すと、**座標は新しい箱・セルは古い箱**という
組み合わせになる。下流の GROMACS コンバータがそれで、変換自体は通り、
**実行時にセルが inf になった** (2026-09-11 に実機で報告)。

NVT では箱が変わらないので露見しない。**NPT を通した軌跡でだけ壊れる。**

UDFManager は CI の必須依存ではないので stub で差し替える。
"""
from __future__ import annotations

import logging

import pytest

from abmptools.gro2udf.top_exporter import (
    TopExporter,
    _warn_if_template_box_differs,
)
from abmptools.gro2udf.udf_writer import ff_comment, set_force_field_comment
from abmptools.gro2udf.top_model import GROFrameData


class _RecordingUDF:
    """put されたパスと値を覚えるだけの stand-in。"""

    def __init__(self, static_cell=None):
        self.puts = {}
        self._static_cell = static_cell
        self._comment = None

    def jump(self, _rec):
        pass

    def put(self, value, path, *_a, **_kw):
        self.puts[path] = value

    def get(self, path, *_a, **_kw):
        if path == "Structure.Unit_Cell.Cell_Size":
            return self._static_cell
        if path == "Unit_Parameter.Comment":
            return self._comment
        return None


def _frame(a, b, c):
    """箱だけを持つ最小フレーム (cell は nm)。"""
    return GROFrameData(step=0, time=0.0, coord_list=[], cell=(a, b, c))


# ---------------------------------------------------------------------------
# 静的セルを書く側
# ---------------------------------------------------------------------------

def test_static_cell_comes_from_the_first_frame():
    u = _RecordingUDF()
    TopExporter._write_static_cell(u, _frame(2.29911, 2.29911, 5.20238))

    assert u.puts["Structure.Unit_Cell.Cell_Size.a"] == 2.29911
    assert u.puts["Structure.Unit_Cell.Cell_Size.c"] == 5.20238


def test_initial_unit_cell_is_written_too():
    """ここを忘れると 100 Å 立方が全レコードに残り続ける。"""
    u = _RecordingUDF()
    TopExporter._write_static_cell(u, _frame(2.29911, 2.29911, 5.20238))

    assert u.puts["Initial_Structure.Initial_Unit_Cell.Cell_Size.a"] == 2.29911
    assert u.puts["Initial_Structure.Initial_Unit_Cell.Cell_Size.c"] == 5.20238


def test_angles_are_set_on_both():
    u = _RecordingUDF()
    TopExporter._write_static_cell(u, _frame(1.0, 2.0, 3.0))

    for field in ("alpha", "beta", "gamma"):
        assert u.puts["Structure.Unit_Cell.Cell_Size.%s" % field] == 90.0
        assert u.puts[
            "Initial_Structure.Initial_Unit_Cell.Cell_Size.%s" % field] == 90.0


def test_the_first_frame_wins_not_the_last():
    """Initial_Unit_Cell が名前どおりであるために、最初のフレームを使う。"""
    u = _RecordingUDF()
    TopExporter._write_static_cell(u, _frame(2.0, 2.0, 2.0))

    assert u.puts["Initial_Structure.Initial_Unit_Cell.Cell_Size.a"] == 2.0


# ---------------------------------------------------------------------------
# 食い違いを知らせる側
# ---------------------------------------------------------------------------

def test_warns_when_the_template_box_is_not_the_data_box(caplog):
    """同梱テンプレートの 20 Å 立方で変換したときに出るべき警告。"""
    u = _RecordingUDF(static_cell=[20.0, 20.0, 20.0, 90.0, 90.0, 90.0])
    with caplog.at_level(logging.WARNING):
        _warn_if_template_box_differs(u, "default_template.udf",
                                      _frame(2.29911, 2.29911, 5.20238))

    assert "20.0000" in caplog.text and "52.0238" in caplog.text


def test_the_warning_names_the_template(caplog):
    """MD 前の UDF をテンプレートに渡した事故も、ここで気付ける。"""
    u = _RecordingUDF(static_cell=[20.0, 20.0, 20.0, 90.0, 90.0, 90.0])
    with caplog.at_level(logging.WARNING):
        _warn_if_template_box_differs(u, "pre_md.udf", _frame(1.0, 1.0, 1.0))

    assert "pre_md.udf" in caplog.text and "pre-MD" in caplog.text


def test_silent_when_the_boxes_agree(caplog):
    u = _RecordingUDF(static_cell=[22.9911, 22.9911, 52.0238, 90.0, 90.0, 90.0])
    with caplog.at_level(logging.WARNING):
        _warn_if_template_box_differs(u, "t.udf",
                                      _frame(2.29911, 2.29911, 5.20238))

    assert caplog.text == ""


def test_a_template_without_a_cell_is_not_an_error(caplog):
    """節が無いテンプレートでも変換は続ける。"""
    u = _RecordingUDF(static_cell=None)
    with caplog.at_level(logging.WARNING):
        _warn_if_template_box_differs(u, "t.udf", _frame(1.0, 1.0, 1.0))

    assert caplog.text == ""


# ---------------------------------------------------------------------------
# 力場 ID (Unit_Parameter.Comment)
# ---------------------------------------------------------------------------

def test_gaff_is_ff_2():
    """力場を取り直すと FF=2 が入る。 実測 3 件で対応を確認済み。"""
    assert ff_comment("gaff") == "FF=2"
    assert ff_comment("GAFF2") == "FF=3"
    assert ff_comment("dreiding") == "FF=4"


def test_a_number_or_a_ready_made_string_passes_through():
    assert ff_comment(2) == "FF=2"
    assert ff_comment("FF=7") == "FF=7"


def test_an_unknown_force_field_is_refused():
    with pytest.raises(ValueError):
        ff_comment("charmm36")


def test_empty_means_leave_it_alone():
    u = _RecordingUDF()
    set_force_field_comment(u, "")
    assert "Unit_Parameter.Comment" not in u.puts


def test_the_id_is_written_when_the_field_is_empty():
    """空のままだと 下流が「力場が分からない」扱いにして書き出しが通らない。"""
    u = _RecordingUDF()
    u._comment = ""
    set_force_field_comment(u, "gaff")
    assert u.puts["Unit_Parameter.Comment"] == "FF=2"


def test_a_template_that_already_names_one_is_kept():
    """テンプレート由来の FF=n を既定で潰さない。"""
    u = _RecordingUDF()
    u._comment = "FF=4"
    set_force_field_comment(u, "gaff")
    assert "Unit_Parameter.Comment" not in u.puts


def test_overwrite_is_available_when_asked_for():
    u = _RecordingUDF()
    u._comment = "FF=4"
    set_force_field_comment(u, "gaff", overwrite=True)
    assert u.puts["Unit_Parameter.Comment"] == "FF=2"
