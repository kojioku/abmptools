# -*- coding: utf-8 -*-
"""
abmptools.fragmenter.notebook_ui
--------------------------------
A 経路: Jupyter + RDKit + ipywidgets による interactive UI。

使い方:
    from abmptools.fragmenter import FragmenterConfig, AutoFragmenter, open_panel

    config = FragmenterConfig(pdb_path="input.pdb", target_mw=200,
                              output_dir="./fragmenter_out")
    af = AutoFragmenter.from_pdb(config)
    af.suggest_cuts()
    open_panel(af)   # Jupyter 上で結合の on/off を toggle

操作:
    - Dropdown で分子グループ選択
    - SVG (heavy atom only) で代表分子を表示、提案 cut sites を赤線ハイライト
    - 各 cut site の Checkbox で enabled を toggle
    - "Export segment_data.dat" ボタンで全コピーへ展開して出力

依存:
    - ipywidgets >= 8 (extras: [jupyter])
    - rdkit-pypi      (extras: [fragmenter])
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import List, Optional

from .auto_split import suggest_cuts_for_groups
from .expand_to_system import export_to_system
from .grouping import group_by_smiles
from .headless_io import (
    _render_svg,
    export_review_bundle as _export_bundle,
    import_edited_review,
    resync_member_indices,
)
from .models import FragmenterConfig, FragmentResult, MoleculeGroup
from .pdb_loader import LoadedMolecule, load_pdb_molecules

logger = logging.getLogger(__name__)


class AutoFragmenter:
    """fragmenter の state 管理 orchestrator。

    CLI からも使えるが、主に Jupyter UI (open_panel) との連携用。
    """

    def __init__(self, config: FragmenterConfig):
        self.config = config
        self.molecules: List[LoadedMolecule] = []
        self.groups: List[MoleculeGroup] = []

    @classmethod
    def from_pdb(cls, config: FragmenterConfig) -> "AutoFragmenter":
        """PDB をロード → グループ化までを一括実行。"""
        af = cls(config)
        af.molecules = load_pdb_molecules(config.pdb_path)
        af.groups = group_by_smiles(
            af.molecules, include_residue_name=config.include_residue_name
        )
        return af

    def suggest_cuts(self) -> None:
        """各 group の代表分子に切断候補を提案 (group.cut_sites を上書き)。"""
        suggest_cuts_for_groups(self.groups, self.molecules, self.config)

    def export_segment_data(self, output_path: Optional[str] = None) -> FragmentResult:
        """全コピーへ展開して segment_data.dat を出力。"""
        if output_path is None:
            output_path = str(Path(self.config.output_dir) / "segment_data.dat")
        return export_to_system(
            self.config.pdb_path, self.groups, self.molecules, output_path
        )

    def export_review_bundle(self, output_dir: Optional[str] = None) -> None:
        """ヘッドレス経路用の SVG + JSON bundle を出力。"""
        if output_dir is None:
            output_dir = self.config.output_dir
        _export_bundle(output_dir, self.groups, self.molecules, self.config)

    def import_edited(self, review_dir: str) -> None:
        """編集後 review_dir を読み込み、self.groups を上書き + member_mol_indices 再計算。"""
        edited = import_edited_review(review_dir)
        resync_member_indices(edited, self.molecules)
        self.groups = edited


def open_panel(fragmenter: AutoFragmenter) -> None:
    """Jupyter notebook 上で interactive UI を起動する。

    UI 構成:
        - Group dropdown: 分子グループ選択
        - Info label: SMILES / n_copies / residue_names / cut sites 数
        - SVG output: 代表分子の構造 (cut sites 赤線ハイライト)
        - Cut site checkboxes: enabled を toggle
        - Export button: segment_data.dat 出力
    """
    try:
        import ipywidgets as widgets
        from IPython.display import display
    except ImportError as e:
        raise ImportError(
            "ipywidgets / IPython is required for open_panel(). "
            "Install with: pip install 'abmptools[jupyter]'"
        ) from e

    if not fragmenter.groups:
        raise ValueError(
            "AutoFragmenter has no groups. "
            "Call AutoFragmenter.from_pdb(config) first."
        )

    state = {"group_idx": 0}

    group_dropdown = widgets.Dropdown(
        options=[
            (
                f"Group {i + 1}: {g.smiles[:30]}{'...' if len(g.smiles) > 30 else ''} "
                f"(n_copies={g.n_copies}, cuts={len(g.cut_sites)})",
                i,
            )
            for i, g in enumerate(fragmenter.groups)
        ],
        value=0,
        description="Group:",
        layout=widgets.Layout(width="700px"),
    )

    svg_output = widgets.HTML(value="")
    cuts_box = widgets.VBox([])
    info_label = widgets.Label(value="")
    export_button = widgets.Button(
        description="Export segment_data.dat",
        button_style="primary",
        layout=widgets.Layout(width="280px"),
    )
    export_status = widgets.HTML(value="")

    def _refresh_svg(group: MoleculeGroup):
        rep = fragmenter.molecules[group.representative_mol_idx]
        svg = _render_svg(rep.mol, group.cut_sites)
        svg_output.value = (
            '<div style="border:1px solid #ccc;padding:10px;background:#fff;">'
            + svg
            + "</div>"
        )

    def _make_cb_handler(group: MoleculeGroup, cs_idx: int):
        def _on_change(change):
            if change.get("name") == "value":
                group.cut_sites[cs_idx].enabled = bool(change["new"])
                _refresh_svg(group)
        return _on_change

    def render(group_idx: int):
        group = fragmenter.groups[group_idx]
        info_label.value = (
            f"SMILES: {group.smiles} | n_copies: {group.n_copies} | "
            f"residues: {sorted(group.residue_names)} | "
            f"cut sites: {len(group.cut_sites)}"
        )
        _refresh_svg(group)

        checkboxes = []
        for k, cs in enumerate(group.cut_sites):
            cb = widgets.Checkbox(
                value=cs.enabled,
                description=f"bond {cs.bond_idx} (atoms {cs.atom1_idx}-{cs.atom2_idx})",
                indent=False,
            )
            cb.observe(_make_cb_handler(group, k))
            checkboxes.append(cb)
        cuts_box.children = checkboxes

    def on_group_change(change):
        if change.get("name") == "value":
            state["group_idx"] = change["new"]
            render(state["group_idx"])

    def on_export_click(_btn):
        try:
            output_path = str(Path(fragmenter.config.output_dir) / "segment_data.dat")
            result = fragmenter.export_segment_data(output_path)
            export_status.value = (
                f'<span style="color:#28a745;">'
                f'OK -- Wrote {output_path} ({result.total_fragments} fragments, '
                f'{len(result.segment_data)} groups)</span>'
            )
        except Exception as e:
            export_status.value = f'<span style="color:#dc3545;">Error: {e}</span>'

    group_dropdown.observe(on_group_change)
    export_button.on_click(on_export_click)

    render(0)

    panel = widgets.VBox(
        [
            widgets.HTML("<h3>abmptools.fragmenter -- interactive panel</h3>"),
            group_dropdown,
            info_label,
            svg_output,
            widgets.HTML("<b>Cut sites (toggle to enable/disable):</b>"),
            cuts_box,
            export_button,
            export_status,
        ]
    )
    display(panel)
