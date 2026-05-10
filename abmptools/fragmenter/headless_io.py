# -*- coding: utf-8 -*-
"""
abmptools.fragmenter.headless_io
--------------------------------
ヘッドレス経路 (C 系統) のサポート。

レビューバンドル (SVG + JSON + config.json) を出力し、ユーザーが手動で
*.json の cut_sites.enabled を編集後、別 CLI コマンドで segment_data.dat
を生成するワークフローを提供する。

Output layout (export_review_bundle):
    output_dir/
        config.json              # FragmenterConfig (再 apply 用)
        review.json              # 全 group のサマリ
        group_001.svg            # 代表分子の SVG (cut_sites を赤線ハイライト)
        group_001.json           # 編集可能 JSON (cut_sites.enabled を toggle)
        group_002.svg
        group_002.json
        ...
"""
from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any, List

from .models import CutSite, FragmenterConfig, MoleculeGroup
from .pdb_loader import LoadedMolecule

logger = logging.getLogger(__name__)


def export_review_bundle(
    output_dir: str,
    groups: List[MoleculeGroup],
    molecules: List[LoadedMolecule],
    config: FragmenterConfig,
) -> None:
    """各 group の SVG + JSON + 全体 review.json + config.json を出力する。"""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    config.to_json(str(out / "config.json"))

    summary: List[dict] = []
    for i, group in enumerate(groups, start=1):
        rep = molecules[group.representative_mol_idx]

        svg_text = _render_svg(rep.mol, group.cut_sites)
        svg_path = out / f"group_{i:03d}.svg"
        svg_path.write_text(svg_text)

        json_path = out / f"group_{i:03d}.json"
        with json_path.open("w") as f:
            json.dump(group.to_dict(), f, indent=2, ensure_ascii=False)

        summary.append({
            "group_idx": i,
            "smiles": group.smiles,
            "residue_names": sorted(group.residue_names),
            "n_copies": group.n_copies,
            "n_cut_sites": len(group.cut_sites),
            "svg": svg_path.name,
            "json": json_path.name,
        })

    with (out / "review.json").open("w") as f:
        json.dump({"groups": summary}, f, indent=2, ensure_ascii=False)

    logger.info("Wrote review bundle to %s (%d groups)", out, len(groups))


def import_edited_review(review_dir: str) -> List[MoleculeGroup]:
    """編集後の group_*.json を読み込み、MoleculeGroup list として返す。

    member_mol_indices / n_copies / representative_mol_idx は JSON のまま読まれるが、
    PDB を再ロードして apply するときには resync_member_indices() で更新する。
    """
    review_path = Path(review_dir)
    if not review_path.exists():
        raise FileNotFoundError(f"review_dir not found: {review_dir}")

    group_files = sorted(review_path.glob("group_*.json"))
    if not group_files:
        raise FileNotFoundError(f"No group_*.json found in {review_dir}")

    groups: List[MoleculeGroup] = []
    for jf in group_files:
        with jf.open() as f:
            d = json.load(f)
        groups.append(MoleculeGroup.from_dict(d))
    return groups


def resync_member_indices(
    edited_groups: List[MoleculeGroup],
    molecules: List[LoadedMolecule],
) -> List[MoleculeGroup]:
    """edited_groups の SMILES と molecules の SMILES を match させて
    member_mol_indices / n_copies / representative_mol_idx を再計算する。

    PDB を再ロードした後の apply フェーズで呼び出すこと。
    """
    from rdkit import Chem

    mol_smiles: List[str] = []
    for lm in molecules:
        try:
            heavy = Chem.RemoveHs(lm.mol)
            mol_smiles.append(Chem.MolToSmiles(heavy, canonical=True))
        except Exception:
            mol_smiles.append(Chem.MolToSmiles(lm.mol, canonical=True))

    for eg in edited_groups:
        matches = [i for i, s in enumerate(mol_smiles) if s == eg.smiles]
        eg.member_mol_indices = matches
        eg.n_copies = len(matches)
        if matches:
            eg.representative_mol_idx = matches[0]
    return edited_groups


# ---------------------------------------------------------------------------
# SVG rendering
# ---------------------------------------------------------------------------


def _render_svg(
    mol: Any,
    cut_sites: List[CutSite],
    width: int = 600,
    height: int = 400,
    show_bond_numbers: bool = False,
) -> str:
    """RDKit で heavy-atom-only SVG を生成、cut_sites の bond と BDA/BAA atom を
    識別表示。

    PDB の 3D 配座をそのまま投影すると化学者向けの構造式として読みにくいので、
    `Compute2DCoords` で 2D 座標を再生成してから描画する。これで PE / PP の
    ような長鎖でもまっすぐな zigzag に整えて表示される。

    Display scheme:
        - cut bond:  red line (R=1.0, G=0.2, B=0.2)
        - **BDA atom**: 点線縁取りの青い円 (塗りつぶしなし) + 'BDA' atom note
          — electron pair holder
        - **BAA atom**: オレンジ塗り highlight + 'BAA' atom note
          — H で擬似 capping される側

    BDA を点線中抜き、BAA を塗りで対比することで「電子は BDA 側に donate、
    BAA 側は H で capping される」という化学的方向性を視覚化する。

    `cs.bda_atom_idx` / `cs.baa_atom_idx` が ``None`` の legacy CutSite では
    atom 装飾なし、cut bond の赤線のみ。
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit.Chem.Draw import rdMolDraw2D

    # H 非表示で描画 (cut_sites も heavy atom 同士の bond なので問題なし)
    try:
        heavy_mol = Chem.RemoveHs(mol)
    except Exception:
        heavy_mol = mol  # sanitize 不完全な mol はそのまま使う

    # 2D 座標を再計算 (PDB 由来の 3D 配座をそのまま 2D 投影すると線が交差して
    # 視認性が悪くなるため、化学者向けの 2D 構造式として描画し直す)。
    try:
        AllChem.Compute2DCoords(heavy_mol)
    except Exception:
        pass  # 失敗しても元 conformer で続行

    # 元 mol の heavy atom idx を heavy_mol の atom idx に mapping
    heavy_atom_origs = [
        a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() != 1
    ]
    orig_to_heavy = {orig: i for i, orig in enumerate(heavy_atom_origs)}

    # bond 番号表示 (Jupyter UI の Add Cut で使う)
    if show_bond_numbers:
        for bond in heavy_mol.GetBonds():
            bond.SetProp("bondNote", str(bond.GetIdx()))

    hl_bonds: List[int] = []
    bond_colors: dict = {}
    hl_atoms: List[int] = []        # BAA only (BDA は SVG 後処理で点線円を描く)
    atom_colors: dict = {}
    bda_h_atoms: List[int] = []     # BDA atom の heavy_mol 内 idx (post-process 用)
    BDA_OUTLINE = "rgb(76,127,255)"   # blue (0.3, 0.5, 1.0) を 0-255 換算
    BAA_COLOR = (1.0, 0.6, 0.0)        # orange
    CUT_BOND_COLOR = (1.0, 0.2, 0.2)   # red

    for cs in cut_sites:
        if not cs.enabled:
            continue
        a1_h = orig_to_heavy.get(cs.atom1_idx)
        a2_h = orig_to_heavy.get(cs.atom2_idx)
        if a1_h is None or a2_h is None:
            continue
        bond = heavy_mol.GetBondBetweenAtoms(a1_h, a2_h)
        if bond is not None:
            hl_bonds.append(bond.GetIdx())
            bond_colors[bond.GetIdx()] = CUT_BOND_COLOR

        # BDA: ハイライトせず atomNote 付与 + post-process で点線円
        if cs.bda_atom_idx is not None:
            bda_h = orig_to_heavy.get(cs.bda_atom_idx)
            if bda_h is not None:
                heavy_mol.GetAtomWithIdx(bda_h).SetProp("atomNote", "BDA")
                bda_h_atoms.append(bda_h)
        # BAA: オレンジ塗り highlight + atomNote
        if cs.baa_atom_idx is not None:
            baa_h = orig_to_heavy.get(cs.baa_atom_idx)
            if baa_h is not None:
                hl_atoms.append(baa_h)
                atom_colors[baa_h] = BAA_COLOR
                heavy_mol.GetAtomWithIdx(baa_h).SetProp("atomNote", "BAA")

    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    rdMolDraw2D.PrepareAndDrawMolecule(
        drawer,
        heavy_mol,
        highlightAtoms=hl_atoms,
        highlightAtomColors=atom_colors,
        highlightBonds=hl_bonds,
        highlightBondColors=bond_colors,
    )

    # BDA atom 位置を drawer から取得して、SVG 後処理で点線円を挿入
    bda_circle_elems: List[str] = []
    for bda_h in bda_h_atoms:
        try:
            pt = drawer.GetDrawCoords(bda_h)
            bda_circle_elems.append(
                f'<circle cx="{pt.x:.2f}" cy="{pt.y:.2f}" r="14" '
                f'fill="none" stroke="{BDA_OUTLINE}" '
                f'stroke-width="2" stroke-dasharray="3,3"/>'
            )
        except Exception:
            pass

    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()

    if bda_circle_elems:
        svg = svg.replace("</svg>", "\n".join(bda_circle_elems) + "\n</svg>")
    return svg
