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
        - "Show bond numbers" checkbox: SVG 上に bond_idx を表示 (Add Cut で参照)
        - SVG output: 代表分子の構造 (cut bond=赤、BDA=青の点線円、BAA=オレンジ塗り)
        - Cut sites: 各行に [enabled checkbox] + [Remove ボタン]
        - Add cut: bond_idx (heavy_mol) を入力 + "Add cut" ボタン
        - Re-suggest: target_mw を変更 + "Re-suggest" ボタン (全 cut を再提案で上書き)
        - Export button: `segment_data.dat` 出力
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

    from rdkit import Chem

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

    show_bond_nums = widgets.Checkbox(
        value=True,
        description="Show bond numbers on SVG (heavy_mol indexing)",
        indent=False,
        layout=widgets.Layout(width="500px"),
    )

    svg_output = widgets.HTML(value="")
    cuts_box = widgets.VBox([])
    info_label = widgets.Label(value="")

    add_cut_input = widgets.IntText(
        value=0, description="Bond idx:",
        layout=widgets.Layout(width="220px"),
    )
    add_cut_bda_choice = widgets.Dropdown(
        options=[
            ("auto (rule-based: C-X→C side, C-C→smaller idx)", "auto"),
            ("atom1 = BDA", "atom1"),
            ("atom2 = BDA", "atom2"),
        ],
        value="auto",
        description="BDA:",
        layout=widgets.Layout(width="500px"),
    )
    add_cut_button = widgets.Button(
        description="Add cut",
        button_style="info",
        layout=widgets.Layout(width="120px"),
    )
    add_cut_status = widgets.HTML(value="")

    target_mw_input = widgets.FloatText(
        value=fragmenter.config.target_mw,
        description="Target MW:",
        layout=widgets.Layout(width="220px"),
    )
    cb_exclude_hetero = widgets.Checkbox(
        value=fragmenter.config.exclude_heteroneighbor,
        description="exclude_heteroneighbor (C-C bonds adjacent to N/O/S/... 除外)",
        indent=False,
        layout=widgets.Layout(width="600px"),
    )
    cb_include_c_hetero = widgets.Checkbox(
        value=fragmenter.config.include_c_heteroatom,
        description="include_c_heteroatom (C-X 単結合も切断対象に: X 側が BAA)",
        indent=False,
        layout=widgets.Layout(width="600px"),
    )
    re_suggest_button = widgets.Button(
        description="Re-suggest (overwrite)",
        button_style="warning",
        layout=widgets.Layout(width="200px"),
    )
    re_suggest_status = widgets.HTML(value="")

    export_button = widgets.Button(
        description="Export segment_data.dat",
        button_style="primary",
        layout=widgets.Layout(width="280px"),
    )
    export_status = widgets.HTML(value="")

    def _refresh_svg(group: MoleculeGroup):
        rep = fragmenter.molecules[group.representative_mol_idx]
        svg = _render_svg(
            rep.mol, group.cut_sites,
            show_bond_numbers=show_bond_nums.value,
        )
        svg_output.value = (
            '<div style="border:1px solid #ccc;padding:10px;background:#fff;">'
            + svg + "</div>"
        )

    def _refresh_dropdown():
        """Group dropdown のラベル (cuts=N) を最新に更新。"""
        group_dropdown.options = [
            (
                f"Group {i + 1}: {g.smiles[:30]}{'...' if len(g.smiles) > 30 else ''} "
                f"(n_copies={g.n_copies}, cuts={len(g.cut_sites)})",
                i,
            )
            for i, g in enumerate(fragmenter.groups)
        ]

    def _make_cb_handler(group: MoleculeGroup, cs_idx: int):
        def _on_change(change):
            if change.get("name") == "value":
                group.cut_sites[cs_idx].enabled = bool(change["new"])
                _refresh_svg(group)
        return _on_change

    def _make_remove_handler(group: MoleculeGroup, cs_idx: int):
        def _on_click(_btn):
            del group.cut_sites[cs_idx]
            _refresh_dropdown()
            render(state["group_idx"])
        return _on_click

    def render(group_idx: int):
        group = fragmenter.groups[group_idx]
        info_label.value = (
            f"SMILES: {group.smiles[:50]}{'...' if len(group.smiles) > 50 else ''} | "
            f"n_copies: {group.n_copies} | "
            f"residues: {sorted(group.residue_names)} | "
            f"cut sites: {len(group.cut_sites)}"
        )
        _refresh_svg(group)

        rows = []
        for k, cs in enumerate(group.cut_sites):
            label = (
                f"bond {cs.bond_idx} (atoms {cs.atom1_idx}-{cs.atom2_idx})  "
                f"BDA={cs.bda_atom_idx}, BAA={cs.baa_atom_idx}, suggested={cs.suggested}"
            )
            cb = widgets.Checkbox(
                value=cs.enabled,
                description=label,
                indent=False,
                layout=widgets.Layout(width="600px"),
            )
            cb.observe(_make_cb_handler(group, k))
            rm_btn = widgets.Button(
                description="Remove",
                button_style="danger",
                layout=widgets.Layout(width="100px"),
            )
            rm_btn.on_click(_make_remove_handler(group, k))
            rows.append(widgets.HBox([cb, rm_btn]))
        cuts_box.children = rows

    def on_group_change(change):
        if change.get("name") == "value":
            state["group_idx"] = change["new"]
            render(state["group_idx"])

    def on_show_bond_nums_change(change):
        if change.get("name") == "value":
            _refresh_svg(fragmenter.groups[state["group_idx"]])

    def on_add_cut_click(_btn):
        from .auto_split import decide_bda_baa_for_manual_cut
        from .models import CutSite
        try:
            target_bond_idx = int(add_cut_input.value)
            group = fragmenter.groups[state["group_idx"]]
            rep = fragmenter.molecules[group.representative_mol_idx]
            try:
                heavy_mol = Chem.RemoveHs(rep.mol)
            except Exception:
                heavy_mol = rep.mol

            if target_bond_idx < 0 or target_bond_idx >= heavy_mol.GetNumBonds():
                add_cut_status.value = (
                    f'<span style="color:#dc3545;">Invalid bond_idx '
                    f'{target_bond_idx} (range 0..{heavy_mol.GetNumBonds() - 1})</span>'
                )
                return

            heavy_bond = heavy_mol.GetBondWithIdx(target_bond_idx)
            a1_h = heavy_bond.GetBeginAtomIdx()
            a2_h = heavy_bond.GetEndAtomIdx()
            heavy_atom_origs = [
                a.GetIdx() for a in rep.mol.GetAtoms() if a.GetAtomicNum() != 1
            ]
            a1_orig = heavy_atom_origs[a1_h]
            a2_orig = heavy_atom_origs[a2_h]
            orig_bond = rep.mol.GetBondBetweenAtoms(a1_orig, a2_orig)
            if orig_bond is None:
                add_cut_status.value = (
                    '<span style="color:#dc3545;">Could not map bond to original mol</span>'
                )
                return
            orig_bond_idx = orig_bond.GetIdx()

            for cs in group.cut_sites:
                if cs.bond_idx == orig_bond_idx:
                    add_cut_status.value = (
                        f'<span style="color:#dc3545;">Bond {target_bond_idx} (orig '
                        f'{orig_bond_idx}) already has a cut. Use Remove first.</span>'
                    )
                    return

            # BDA / BAA 決定: dropdown が "auto" なら rule-based、それ以外は手動指定
            bda_choice = add_cut_bda_choice.value
            if bda_choice == "atom1":
                bda, baa = a1_orig, a2_orig
            elif bda_choice == "atom2":
                bda, baa = a2_orig, a1_orig
            else:
                bda, baa = decide_bda_baa_for_manual_cut(rep.mol, a1_orig, a2_orig)
            new_cut = CutSite(
                bond_idx=orig_bond_idx,
                atom1_idx=a1_orig,
                atom2_idx=a2_orig,
                suggested=False,
                enabled=True,
                bda_atom_idx=bda,
                baa_atom_idx=baa,
            )
            group.cut_sites.append(new_cut)
            add_cut_status.value = (
                f'<span style="color:#28a745;">Added cut: heavy bond '
                f'{target_bond_idx}, orig bond {orig_bond_idx}, atoms '
                f'{a1_orig}-{a2_orig}, BDA={bda}, BAA={baa} '
                f'(mode={bda_choice})</span>'
            )
            _refresh_dropdown()
            render(state["group_idx"])
        except Exception as e:
            add_cut_status.value = f'<span style="color:#dc3545;">Error: {e}</span>'

    def on_re_suggest_click(_btn):
        try:
            new_mw = float(target_mw_input.value)
            fragmenter.config.target_mw = new_mw
            fragmenter.config.exclude_heteroneighbor = bool(cb_exclude_hetero.value)
            fragmenter.config.include_c_heteroatom = bool(cb_include_c_hetero.value)
            fragmenter.suggest_cuts()
            total = sum(len(g.cut_sites) for g in fragmenter.groups)
            re_suggest_status.value = (
                f'<span style="color:#28a745;">Re-suggested with target_mw='
                f'{new_mw}, exclude_hetero={cb_exclude_hetero.value}, '
                f'include_c_hetero={cb_include_c_hetero.value}: '
                f'{total} cut(s) total across {len(fragmenter.groups)} '
                f'group(s)</span>'
            )
            _refresh_dropdown()
            render(state["group_idx"])
        except Exception as e:
            re_suggest_status.value = f'<span style="color:#dc3545;">Error: {e}</span>'

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
    show_bond_nums.observe(on_show_bond_nums_change)
    add_cut_button.on_click(on_add_cut_click)
    re_suggest_button.on_click(on_re_suggest_click)
    export_button.on_click(on_export_click)

    render(0)

    panel = widgets.VBox(
        [
            widgets.HTML("<h3>abmptools.fragmenter -- interactive panel</h3>"),
            group_dropdown,
            info_label,
            show_bond_nums,
            svg_output,
            widgets.HTML("<b>Cut sites (toggle = enable/disable, Remove = delete):</b>"),
            cuts_box,
            widgets.HTML(
                '<hr><b>Add cut by heavy-mol bond index</b> '
                '(see numbers on the SVG; BDA を "auto" 以外で固定可):'
            ),
            widgets.HBox([add_cut_input, add_cut_button]),
            add_cut_bda_choice,
            add_cut_status,
            widgets.HTML(
                '<hr><b>Re-suggest auto cuts</b> '
                '(<i>overwrites all current cuts in all groups</i>):'
            ),
            widgets.HBox([target_mw_input, re_suggest_button]),
            cb_exclude_hetero,
            cb_include_c_hetero,
            re_suggest_status,
            widgets.HTML("<hr>"),
            export_button,
            export_status,
        ]
    )
    display(panel)
