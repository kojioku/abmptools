# -*- coding: utf-8 -*-
"""
abmptools.cg.dpd.notebook_ui
----------------------------
Jupyter notebook 上で :class:`CGDpdBuilder` を interactive に build / verify する UI。

使い方:

    from abmptools.cg.dpd import CGDpdBuilder, open_panel

    builder = CGDpdBuilder.from_files(
        monomer="chol_monomer", aij="aij.dat", calc_sett="chol_calc_sett",
    )
    open_panel(builder)

UI 構成:
- **Summary**: monomer 一覧 + aij segment 一覧 + `validate()` 結果 (warning または OK)
- **R1 (Cognac UDF)**: output path text + Build ボタン
- **R2 (J-OCTA DPM、 B 案)**: template path + Virtual.mom + output dir + dpm filename + Build ボタン
- **Verify**: validate() を再実行して summary を更新

依存:
- ipywidgets >= 8 (extras: ``[jupyter]``)
- IPython (Jupyter 標準)
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .orchestrator import CGDpdBuilder

logger = logging.getLogger(__name__)


def open_panel(builder: "CGDpdBuilder") -> None:
    """Jupyter notebook 上で ``CGDpdBuilder`` の interactive UI を起動する。

    Parameters
    ----------
    builder : CGDpdBuilder
        :meth:`CGDpdBuilder.from_files` or :meth:`CGDpdBuilder.from_multi_files` で
        構築済の builder インスタンス。

    Raises
    ------
    ImportError
        ipywidgets / IPython が install されていない場合。
        ``pip install 'abmptools[jupyter]'`` を実行する。
    """
    try:
        import ipywidgets as widgets
        from IPython.display import display
    except ImportError as e:
        raise ImportError(
            "ipywidgets / IPython is required for open_panel(). "
            "Install with: pip install 'abmptools[jupyter]'"
        ) from e

    # --- Summary HTML 表示 ---
    summary = widgets.HTML(value="")

    def _refresh_summary():
        warnings = builder.validate()
        html = "<h3>cg.dpd Builder Summary</h3>"
        html += f"<p><b>Monomers</b>: {len(builder.spec.monomers)}</p><ul>"
        for m in builder.spec.monomers:
            names_str = ", ".join(m.particle_names) if m.particle_names else "(empty)"
            html += (
                f"<li><b>{m.name}</b>: {m.n_particles} particles "
                f"[{names_str}], bond12={len(m.bond12)} / "
                f"bond13_150={len(m.bond13_150)} / bond14_150={len(m.bond14_150)} / "
                f"angle13={len(m.angle13)}</li>"
            )
        html += "</ul>"
        html += (
            f"<p><b>aij</b>: mode=<code>{builder.spec.aij.mode}</code>, "
            f"aii={builder.spec.aij.aii}, "
            f"{len(builder.spec.aij.pairs)} pairs, "
            f"{len(builder.spec.aij.segments)} segments "
            f"({', '.join(builder.spec.aij.segments)})</p>"
        )
        if warnings:
            html += (
                f'<p style="color:#dc3545;"><b>Warnings ({len(warnings)})</b>:</p>'
                "<ul>"
            )
            for w in warnings:
                html += f'<li style="color:#dc3545;">{w}</li>'
            html += "</ul>"
        else:
            html += '<p style="color:#28a745;"><b>✓ 整合性 OK</b></p>'
        summary.value = html

    # --- R1 (Cognac UDF) ---
    udf_path_input = widgets.Text(
        value="output_uin.udf",
        description="UDF path:",
        layout=widgets.Layout(width="400px"),
    )
    udf_include_input = widgets.Text(
        value="cognac112.udf",
        description="\\include:",
        layout=widgets.Layout(width="300px"),
    )
    udf_btn = widgets.Button(
        description="Build R1 UDF",
        button_style="primary",
        layout=widgets.Layout(width="180px"),
    )
    udf_status = widgets.HTML(value="")

    def on_udf_click(_b):
        try:
            out = builder.build_udf(
                udf_path_input.value, include_file=udf_include_input.value,
            )
            udf_status.value = (
                f'<span style="color:#28a745;">UDF generated: '
                f'<code>{out}</code> ({Path(out).stat().st_size:,} bytes)</span>'
            )
        except Exception as e:
            udf_status.value = (
                f'<span style="color:#dc3545;">Error: {type(e).__name__}: {e}</span>'
            )

    # --- R2 (J-OCTA DPM、 B 案) ---
    dpm_template = widgets.Text(
        value="",
        placeholder="J-OCTA で作成した空 dpm template path",
        description="Template:",
        layout=widgets.Layout(width="500px"),
    )
    vmom_template = widgets.Text(
        value="",
        placeholder="J-OCTA Monomer Modeler で作成した Virtual.mom path (省略可)",
        description="Virtual.mom:",
        layout=widgets.Layout(width="500px"),
    )
    dpm_outdir = widgets.Text(
        value="./dpm_out",
        description="Output dir:",
        layout=widgets.Layout(width="350px"),
    )
    dpm_filename = widgets.Text(
        value="system.dpm",
        description="dpm name:",
        layout=widgets.Layout(width="280px"),
    )
    dpm_btn = widgets.Button(
        description="Build R2 DPM",
        button_style="info",
        layout=widgets.Layout(width="180px"),
    )
    dpm_status = widgets.HTML(value="")

    def on_dpm_click(_b):
        try:
            if not dpm_template.value:
                raise ValueError(
                    "Template path 必須 (B 案運用: J-OCTA で空 dpm を 1 回作成して指定)"
                )
            out = builder.build_dpm(
                template=dpm_template.value,
                output_dir=dpm_outdir.value,
                virtual_mom_template=(vmom_template.value or None),
                dpm_filename=dpm_filename.value,
            )
            dpm_status.value = (
                f'<span style="color:#28a745;">DPM generated: '
                f'<code>{out}</code></span>'
            )
        except Exception as e:
            dpm_status.value = (
                f'<span style="color:#dc3545;">Error: {type(e).__name__}: {e}</span>'
            )

    # --- Verify button ---
    verify_btn = widgets.Button(
        description="Re-verify (整合性チェック)",
        button_style="warning",
        layout=widgets.Layout(width="240px"),
    )

    def on_verify_click(_b):
        _refresh_summary()

    udf_btn.on_click(on_udf_click)
    dpm_btn.on_click(on_dpm_click)
    verify_btn.on_click(on_verify_click)

    # 初期描画
    _refresh_summary()

    panel = widgets.VBox([
        widgets.HTML("<h2>abmptools.cg.dpd — interactive panel</h2>"),
        summary,
        widgets.HTML("<hr><b>R1: Cognac DPD 入力 UDF を生成</b>"),
        widgets.HBox([udf_path_input, udf_include_input]),
        udf_btn,
        udf_status,
        widgets.HTML(
            "<hr><b>R2: J-OCTA DPM 生成 (B 案 = user template patch)</b><br>"
            "<i>Template + Virtual.mom は user が J-OCTA で 1 回だけ作成して指定</i>"
        ),
        dpm_template,
        vmom_template,
        widgets.HBox([dpm_outdir, dpm_filename]),
        dpm_btn,
        dpm_status,
        widgets.HTML("<hr>"),
        verify_btn,
    ])
    display(panel)
