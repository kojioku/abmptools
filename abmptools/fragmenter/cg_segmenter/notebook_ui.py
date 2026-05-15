# -*- coding: utf-8 -*-
"""
abmptools.fragmenter.cg_segmenter.notebook_ui
---------------------------------------------
Jupyter notebook 上で CGSegmenter を interactive に表示・編集する UI。

使い方:
    from abmptools.fragmenter.cg_segmenter import (
        CGSegmenter, CGSegmenterConfig, open_panel,
    )

    config = CGSegmenterConfig(pdb_path="cholesterol.pdb", target_mw=200)
    sg = CGSegmenter.from_pdb(config)
    open_panel(sg)

UI 構成:
    - Show atom numbers / Bold shared atoms toggle
    - SVG output: 各 segment を別色で highlight、shared atom を黒縁取り
    - Segments list: 各 segment ごとに [Delete] + [cap_n toggle ...] ボタン
    - Move atom: atom_idx + target_seg + shared checkbox + Move ボタン
    - Re-segment: target_mw / config flags + Re-segment ボタン (全 segment 上書き)
    - Export: PDB + XYZ + segments.json 一括出力

依存:
    - ipywidgets >= 8 (extras: [jupyter])
    - rdkit (extras: [fragmenter])
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .orchestrator import CGSegmenter

from .exporter import render_segments_svg

logger = logging.getLogger(__name__)


def open_panel(sg: "CGSegmenter") -> None:
    """Jupyter notebook 上で CGSegmenter の interactive UI を起動する。"""
    try:
        import ipywidgets as widgets
        from IPython.display import display
    except ImportError as e:
        raise ImportError(
            "ipywidgets / IPython is required for open_panel(). "
            "Install with: pip install 'abmptools[jupyter]'"
        ) from e

    if sg.mol is None:
        raise ValueError(
            "CGSegmenter.mol is None. Call CGSegmenter.from_pdb(config) first."
        )

    from rdkit import Chem

    # --- Widgets ---
    show_atom_nums = widgets.Checkbox(
        value=False, description="Show atom numbers",
        indent=False, layout=widgets.Layout(width="240px"),
    )
    show_bold_shared = widgets.Checkbox(
        value=True, description="Bold shared atoms",
        indent=False, layout=widgets.Layout(width="240px"),
    )

    svg_output = widgets.HTML(value="")
    seg_list_box = widgets.VBox([])

    # heavy mol の atom idx list (atom 移動 dropdown 用)
    heavy_mol_local = Chem.RemoveHs(sg.mol)
    atom_options = [
        (f"atom {a.GetIdx()} ({a.GetSymbol()})", a.GetIdx())
        for a in heavy_mol_local.GetAtoms()
    ]
    atom_move_atom = widgets.Dropdown(
        options=atom_options, description="Move atom:",
        layout=widgets.Layout(width="300px"),
    )
    atom_move_to = widgets.Dropdown(
        options=[(f"seg {s.segment_id}", s.segment_id) for s in sg.segments],
        description="To seg:",
        layout=widgets.Layout(width="200px"),
    )
    atom_move_shared = widgets.Checkbox(
        value=False, description="Keep in source (shared)", indent=False,
        layout=widgets.Layout(width="240px"),
    )
    atom_move_btn = widgets.Button(
        description="Move", button_style="info",
        layout=widgets.Layout(width="100px"),
    )
    atom_move_status = widgets.HTML(value="")

    target_mw_input = widgets.FloatText(
        value=sg.config.target_mw, description="Target MW:",
        layout=widgets.Layout(width="220px"),
    )
    sep_rings_cb = widgets.Checkbox(
        value=sg.config.separate_rings, description="Separate rings", indent=False,
    )
    atom_share_cb = widgets.Checkbox(
        value=sg.config.allow_atom_sharing, description="Atom sharing", indent=False,
    )
    absorb_sub_cb = widgets.Checkbox(
        value=sg.config.absorb_single_substituent, description="Absorb substituent", indent=False,
    )
    re_seg_btn = widgets.Button(
        description="Re-segment (overwrite)", button_style="warning",
        layout=widgets.Layout(width="220px"),
    )
    re_seg_status = widgets.HTML(value="")

    export_btn = widgets.Button(
        description="Export PDB + XYZ + JSON", button_style="primary",
        layout=widgets.Layout(width="280px"),
    )
    export_status = widgets.HTML(value="")

    dpdgen_name = widgets.Text(
        value="cg", description="DPDgen name:",
        layout=widgets.Layout(width="220px"),
    )
    dpdgen_box = widgets.FloatText(
        value=10.0, description="Box (DPD):",
        layout=widgets.Layout(width="180px"),
    )
    dpdgen_btn = widgets.Button(
        description="Export DPDgen (monomer + calc_sett)", button_style="success",
        layout=widgets.Layout(width="320px"),
    )
    dpdgen_status = widgets.HTML(value="")

    # --- Refresh functions ---
    def _refresh_svg():
        svg = render_segments_svg(
            sg.mol, sg.segments,
            show_atom_numbers=show_atom_nums.value,
            bold_shared=show_bold_shared.value,
        )
        svg_output.value = (
            '<div style="border:1px solid #ccc;padding:10px;background:#fff;">'
            + svg + "</div>"
        )

    def _refresh_seg_options():
        atom_move_to.options = [
            (f"seg {s.segment_id}", s.segment_id) for s in sg.segments
        ]

    def _refresh_seg_list():
        rows = []
        for seg in sg.segments:
            n_methyl = sum(1 for c in seg.cap_atoms if c.is_methyl_cap)
            n_h = sum(1 for c in seg.cap_atoms if not c.is_methyl_cap)
            header = widgets.HTML(
                value=(
                    f"<b>seg {seg.segment_id}</b> "
                    f"({seg.kind}): {len(seg.atom_indices)} heavy, "
                    f"{n_methyl} CH3 + {n_h} H caps"
                )
            )
            info = widgets.HTML(
                value=f"&nbsp;&nbsp;atom_indices: {seg.atom_indices[:20]}"
                       f"{'...' if len(seg.atom_indices) > 20 else ''}"
            )

            del_btn = widgets.Button(
                description=f"Delete seg {seg.segment_id}",
                button_style="danger",
                layout=widgets.Layout(width="160px"),
            )
            del_btn.on_click(_make_delete_handler(seg.segment_id))

            cap_buttons = []
            for ci, cap in enumerate(seg.cap_atoms):
                kind = "CH3" if cap.is_methyl_cap else "H"
                btn = widgets.Button(
                    description=f"cap[{ci}] p={cap.parent_atom_idx} ({kind})",
                    layout=widgets.Layout(width="220px"),
                )
                btn.on_click(_make_toggle_handler(seg.segment_id, ci))
                cap_buttons.append(btn)

            action_row = widgets.HBox([del_btn] + cap_buttons)
            rows.append(widgets.VBox([header, info, action_row]))
        seg_list_box.children = rows

    def _make_delete_handler(sid):
        def _on_click(_btn):
            try:
                sg.delete_segment(sid)
                _refresh_seg_list()
                _refresh_seg_options()
                _refresh_svg()
            except Exception as e:
                logger.exception("delete_segment failed")
        return _on_click

    def _make_toggle_handler(sid, cap_idx):
        def _on_click(_btn):
            try:
                sg.toggle_cap(sid, cap_idx)
                _refresh_seg_list()
                _refresh_svg()
            except Exception as e:
                logger.exception("toggle_cap failed")
        return _on_click

    # --- Event handlers ---
    def on_move_click(_btn):
        try:
            sg.move_atom(
                atom_move_atom.value, atom_move_to.value,
                shared=atom_move_shared.value,
            )
            atom_move_status.value = (
                f'<span style="color:#28a745;">moved atom {atom_move_atom.value} '
                f'-> seg {atom_move_to.value} '
                f'(shared={atom_move_shared.value})</span>'
            )
            _refresh_seg_list()
            _refresh_svg()
        except Exception as e:
            atom_move_status.value = f'<span style="color:#dc3545;">Error: {e}</span>'

    def on_re_seg_click(_btn):
        try:
            sg.re_segment(
                target_mw=target_mw_input.value,
                separate_rings=sep_rings_cb.value,
                allow_atom_sharing=atom_share_cb.value,
                absorb_single_substituent=absorb_sub_cb.value,
            )
            re_seg_status.value = (
                f'<span style="color:#28a745;">Re-segmented: '
                f'{len(sg.segments)} segments (target_mw={target_mw_input.value})</span>'
            )
            _refresh_seg_list()
            _refresh_seg_options()
            _refresh_svg()
        except Exception as e:
            re_seg_status.value = f'<span style="color:#dc3545;">Error: {e}</span>'

    def on_export_click(_btn):
        try:
            out_dir = sg.config.output_dir
            result = sg.export(out_dir)
            export_status.value = (
                f'<span style="color:#28a745;">'
                f'Wrote {result.total_atoms_with_cap} atoms across '
                f'{len(result.segments)} segments to {out_dir}</span>'
            )
        except Exception as e:
            export_status.value = f'<span style="color:#dc3545;">Error: {e}</span>'

    def on_dpdgen_click(_btn):
        try:
            box = float(dpdgen_box.value)
            name = dpdgen_name.value or "cg"
            mono, calc = sg.export_dpdgen(
                output_dir=sg.config.output_dir,
                monomer_name=name,
                box_size=(box, box, box),
            )
            dpdgen_status.value = (
                f'<span style="color:#28a745;">'
                f'DPDgen monomer: <code>{mono}</code><br>'
                f'DPDgen calc_sett: <code>{calc}</code><br>'
                f'(edit ratio_list / aij_file / phys_param, then '
                f'<code>python makeudf_dpd.py -p {calc}</code>)</span>'
            )
        except Exception as e:
            dpdgen_status.value = f'<span style="color:#dc3545;">Error: {e}</span>'

    show_atom_nums.observe(lambda _: _refresh_svg(), names="value")
    show_bold_shared.observe(lambda _: _refresh_svg(), names="value")
    atom_move_btn.on_click(on_move_click)
    re_seg_btn.on_click(on_re_seg_click)
    export_btn.on_click(on_export_click)
    dpdgen_btn.on_click(on_dpdgen_click)

    # 初期描画
    _refresh_seg_list()
    _refresh_svg()

    panel = widgets.VBox([
        widgets.HTML("<h3>abmptools.fragmenter.cg_segmenter -- interactive panel</h3>"),
        widgets.HBox([show_atom_nums, show_bold_shared]),
        svg_output,
        widgets.HTML("<b>Segments (Delete = remove from list; cap[n] click = toggle H<->CH3):</b>"),
        seg_list_box,
        widgets.HTML("<hr><b>Move atom across segments:</b>"),
        widgets.HBox([atom_move_atom, atom_move_to, atom_move_shared, atom_move_btn]),
        atom_move_status,
        widgets.HTML(
            "<hr><b>Re-segment</b> "
            "(<i>overwrites all current segments, including manual edits</i>):"
        ),
        widgets.HBox([target_mw_input, sep_rings_cb, atom_share_cb, absorb_sub_cb]),
        re_seg_btn,
        re_seg_status,
        widgets.HTML("<hr>"),
        export_btn,
        export_status,
        widgets.HTML(
            "<hr><b>Export DPDgen input</b> "
            "(monomer + calc_sett for OCTA COGNAC):"
        ),
        widgets.HBox([dpdgen_name, dpdgen_box, dpdgen_btn]),
        dpdgen_status,
    ])
    display(panel)
