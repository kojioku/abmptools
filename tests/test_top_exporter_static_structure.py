# -*- coding: utf-8 -*-
"""gro2udf は静的 Structure を書き換えないので、空でないなら警告する。

UDF は「静的セクション」と「動的レコード」に分かれ、``gro2udf`` は座標とセルを
**レコードにしか書かない**。静的 Structure はテンプレートの内容がそのまま残る。
同梱テンプレートは空なので実害は無いが、実在の UDF をテンプレートに使うと
(J-OCTA や COGNAC の系を往復させるときは自然にそうなる) **変換前の座標と箱が
静的側に残ったまま**、レコードには変換後が入る。どちらも形式としては正しいので
下流は何も言わない。

UDFManager は CI の必須依存ではないので stub で差し替える。
"""
from __future__ import annotations

import logging
import sys
import types
from unittest.mock import patch

import pytest

from abmptools.gro2udf.top_exporter import TopExporter, _static_structure_mol_count
from abmptools.gro2udf.top_model import GROFrameData, TopModel


class _StubUDFManager:
    """静的 Structure の中身だけを差し替えられる最小の stand-in。"""

    #: size("Structure.Position.mol[]") が返す値。None なら例外を投げる。
    static_mols = 0

    def __init__(self, path):
        self._path = path

    def totalRecord(self):
        return 0

    def eraseRecord(self, _a, _b):
        pass

    def jump(self, _rec):
        pass

    def size(self, path):
        if path != "Structure.Position.mol[]":
            return 0
        if self.static_mols is None:
            raise RuntimeError("no such node")
        return self.static_mols

    def newRecord(self, *_a, **_kw):
        pass

    def put(self, *_a, **_kw):
        pass

    def get(self, *_a, **_kw):
        return 0

    def write(self, *_a, **_kw):
        pass


def _stub_module(static_mols):
    cls = type("_Stub", (_StubUDFManager,), {"static_mols": static_mols})
    mod = types.ModuleType("UDFManager")
    mod.UDFManager = cls  # type: ignore[attr-defined]
    return mod


# ---------------------------------------------------------------------------
# 数える側
# ---------------------------------------------------------------------------

def test_counts_the_molecules_the_static_block_declares():
    assert _static_structure_mol_count(_stub_module(3).UDFManager("x")) == 3


def test_empty_static_block_counts_as_zero():
    assert _static_structure_mol_count(_stub_module(0).UDFManager("x")) == 0


def test_a_template_without_the_node_counts_as_zero():
    """節が無いテンプレートで例外が出ても、変換を止める理由にはしない。"""
    assert _static_structure_mol_count(_stub_module(None).UDFManager("x")) == 0


def test_a_non_numeric_answer_counts_as_zero():
    class _Weird(_StubUDFManager):
        def size(self, _path):
            return "many"

    assert _static_structure_mol_count(_Weird("x")) == 0


# ---------------------------------------------------------------------------
# 警告を出す側
# ---------------------------------------------------------------------------

def _minimal_model():
    """書き出しが通る最小のモデル (原子も分子も持たない)。"""
    return TopModel(
        comb_rule=2,
        fudge_lj=0.5,
        fudge_qq=0.8333,
        atom_type_specs=[],
        bond_type_specs=[],
        angle_type_specs=[],
        torsion_type_specs=[],
        mass_dict={},
        mol_type_names=[],
        mol_specs=[],
        mol_instance_list=[],
        frames=[GROFrameData(step=0, time=0.0, coord_list=[], cell=[1.0, 1.0, 1.0])],
    )


def _export(tmp_path, static_mols, caplog):
    template = tmp_path / "template.udf"
    template.write_text('\\include{"cognac112.udf"}\n', encoding="utf-8")
    out = tmp_path / "out.udf"
    with patch.dict(sys.modules, {"UDFManager": _stub_module(static_mols)}), \
            patch.object(TopExporter, "_append_structure", lambda *a, **k: None), \
            patch.object(TopExporter, "_write_set_of_molecules", lambda *a, **k: None), \
            patch.object(TopExporter, "_set_default_condition", lambda *a, **k: None), \
            patch.object(TopExporter, "_write_molecular_attributes", lambda *a, **k: None), \
            patch.object(TopExporter, "_write_interactions", lambda *a, **k: None), \
            caplog.at_level(logging.WARNING, logger="abmptools.gro2udf.top_exporter"):
        TopExporter().export_model(_minimal_model(), str(template), str(out))
    return caplog.text


def test_warns_when_the_template_carries_a_structure(tmp_path, caplog):
    text = _export(tmp_path, 12, caplog)
    assert "static Structure" in text
    assert "12 molecule(s)" in text
    # 何をすればよいかまで書く
    assert "default_template.udf" in text


def test_says_nothing_for_an_empty_template(tmp_path, caplog):
    assert "static Structure" not in _export(tmp_path, 0, caplog)


def test_the_bundled_templates_are_empty():
    """同梱テンプレートで警告が出るなら、警告そのものが無意味になる。"""
    import os

    from abmptools.gro2udf import top_exporter as te

    for name in ("default_template.udf", "default_template_cognac101.udf"):
        path = os.path.join(os.path.dirname(te.__file__), name)
        assert os.path.exists(path), path
        text = open(path, encoding="utf-8", errors="replace").read()
        static = text.split("\\begin{record}", 1)[0]
        start = static.find("\nStructure:{")
        if start == -1:
            # 静的 Structure を持たないテンプレートもある (cognac101 版)。
            # 「空である」の最も強い形なので、それでよい。
            continue
        depth = 0
        for j in range(start + 1, len(static)):
            if static[j] == "{":
                depth += 1
            elif static[j] == "}":
                depth -= 1
                if depth == 0:
                    block = static[start:j + 1]
                    break
        else:
            pytest.fail(f"{name} の Structure ブロックが閉じていない")
        import re
        assert sum(a.count("{") for a in re.findall(r"\[(.*?)\]", block, re.S)) == 0
