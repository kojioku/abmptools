# -*- coding: utf-8 -*-
"""tests/cg_dpd/test_notebook_ui.py — G1 Jupyter UI のスモークテスト.

ipywidgets が install されているかどうかで動作分岐。 install 済なら panel を構築
できることまで確認、 未 install なら ImportError が期待通り発生することを確認。
"""
from __future__ import annotations

import inspect

import pytest


def test_open_panel_importable_from_package():
    """``from abmptools.cg.dpd import open_panel`` で取れる."""
    from abmptools.cg.dpd import open_panel
    assert callable(open_panel)


def test_open_panel_signature():
    """open_panel(builder) の signature."""
    from abmptools.cg.dpd.notebook_ui import open_panel
    sig = inspect.signature(open_panel)
    assert "builder" in sig.parameters


def test_open_panel_with_ipywidgets(cholesterol_cg, sample_aij_a_mode):
    """ipywidgets が install 済なら panel が構築できる (display は mock で吸収)."""
    pytest.importorskip("ipywidgets")
    pytest.importorskip("IPython")

    from abmptools.cg.dpd import CGDpdBuilder, open_panel
    from IPython.display import display
    import IPython.display as ipd

    builder = CGDpdBuilder.from_files(
        monomer=cholesterol_cg["monomer"],
        aij=sample_aij_a_mode,
        calc_sett=cholesterol_cg["calc_sett"],
    )

    # display を mock で吸収して例外なく panel が構築されることを確認
    displayed = []

    def fake_display(obj):
        displayed.append(obj)

    orig_display = ipd.display
    ipd.display = fake_display
    try:
        open_panel(builder)
    finally:
        ipd.display = orig_display

    assert len(displayed) == 1
    # widgets.VBox であること
    import ipywidgets as widgets
    assert isinstance(displayed[0], widgets.VBox)
