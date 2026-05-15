# -*- coding: utf-8 -*-
"""
abmptools.fragmenter.cg_segmenter.exporter
------------------------------------------
per-segment PDB + XYZ + summary JSON 出力。

各 segment について以下を生成:
- `seg_{NNN}.pdb` (PDB format、cap 込み)
- `seg_{NNN}.xyz` (XYZ format、cap 込み)
- `segments.json` (全 segment のサマリ + shared_atom_pairs)
"""
from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any, List, Tuple

from .cap_attach import methyl_hydrogen_positions
from .models import Segment, SegmentResult

logger = logging.getLogger(__name__)


def export_segments(
    mol_with_h: Any,
    segments: List[Segment],
    output_dir: str,
) -> SegmentResult:
    """各 segment を PDB + XYZ で個別出力 + summary JSON。"""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    # 共有 atom 検出
    atom_to_segs: dict = {}
    for seg in segments:
        for a in seg.atom_indices:
            atom_to_segs.setdefault(a, []).append(seg.segment_id)
    shared_atom_pairs: List[Tuple[int, int, int]] = []
    for a, sids in atom_to_segs.items():
        if len(sids) > 1:
            for i in range(len(sids)):
                for j in range(i + 1, len(sids)):
                    shared_atom_pairs.append((a, sids[i], sids[j]))

    total_atoms = 0
    for seg in segments:
        atoms = _collect_atoms(mol_with_h, seg)
        _write_pdb(atoms, seg, out / f"seg_{seg.segment_id:03d}.pdb")
        _write_xyz(atoms, seg, out / f"seg_{seg.segment_id:03d}.xyz")
        total_atoms += len(atoms)

    summary = {
        "n_segments": len(segments),
        "total_atoms_with_cap": total_atoms,
        "shared_atom_pairs": [list(p) for p in shared_atom_pairs],
        "segments": [s.to_dict() for s in segments],
    }
    (out / "segments.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False)
    )

    logger.info(
        "export_segments: %d segments to %s (total %d atoms incl. caps)",
        len(segments), out, total_atoms,
    )

    return SegmentResult(
        segments=segments,
        total_atoms_with_cap=total_atoms,
        shared_atom_pairs=shared_atom_pairs,
    )


def _collect_atoms(mol: Any, seg: Segment) -> List[Tuple[str, float, float, float, str]]:
    """segment 内の (heavy atoms + 直接 attach の H) + cap atoms を集める。

    Returns
    -------
    List[(symbol, x, y, z, label)]
    """
    conf = mol.GetConformer()
    atoms: List[Tuple[str, float, float, float, str]] = []

    # heavy atom + 隣接 H
    for a in seg.atom_indices:
        atom = mol.GetAtomWithIdx(a)
        pos = conf.GetAtomPosition(a)
        atoms.append((atom.GetSymbol(), pos.x, pos.y, pos.z, f"A{a}"))
        for nb in atom.GetNeighbors():
            if nb.GetAtomicNum() == 1:
                h_pos = conf.GetAtomPosition(nb.GetIdx())
                atoms.append(("H", h_pos.x, h_pos.y, h_pos.z, f"H{nb.GetIdx()}"))

    # cap atoms (CH3 cap は central C + 3 H に展開)
    for k, cap in enumerate(seg.cap_atoms):
        cx, cy, cz = cap.position
        if cap.is_methyl_cap:
            atoms.append(("C", cx, cy, cz, f"CC{k}"))
            parent_pos = conf.GetAtomPosition(cap.parent_atom_idx)
            for hi, (hx, hy, hz) in enumerate(methyl_hydrogen_positions(
                cap.position,
                (parent_pos.x, parent_pos.y, parent_pos.z),
                bond_len_ch=1.09,
            )):
                atoms.append(("H", hx, hy, hz, f"CH{k}{hi}"))
        else:
            atoms.append((cap.element, cx, cy, cz, f"Cap{k}"))

    return atoms


def _write_pdb(
    atoms: List[Tuple[str, float, float, float, str]],
    seg: Segment,
    path: Path,
) -> None:
    """PDB ATOM record で 1 segment を書き出す。"""
    res_name = "SEG"
    lines: List[str] = [
        f"REMARK    abmptools.fragmenter.cg_segmenter segment {seg.segment_id} ({seg.kind})",
    ]
    for i, (sym, x, y, z, label) in enumerate(atoms, start=1):
        atom_name = label[:4].ljust(4)
        # PDB ATOM/HETATM record (heuristic alignment)
        lines.append(
            f"ATOM  {i:>5} {atom_name} {res_name} A   1    "
            f"{x:>8.3f}{y:>8.3f}{z:>8.3f}  1.00  0.00          {sym:>2}"
        )
    lines.append("END")
    path.write_text("\n".join(lines) + "\n")


def _write_xyz(
    atoms: List[Tuple[str, float, float, float, str]],
    seg: Segment,
    path: Path,
) -> None:
    """XYZ format で 1 segment を書き出す。"""
    lines: List[str] = [
        str(len(atoms)),
        f"segment {seg.segment_id} ({seg.kind})",
    ]
    for sym, x, y, z, _ in atoms:
        lines.append(f"{sym} {x:.4f} {y:.4f} {z:.4f}")
    path.write_text("\n".join(lines) + "\n")


# ---------------------------------------------------------------------------
# Multi-color SVG renderer for Jupyter UI
# ---------------------------------------------------------------------------

SEG_COLOR_PALETTE = [
    (1.00, 0.55, 0.55),   # red
    (0.55, 1.00, 0.55),   # green
    (0.55, 0.65, 1.00),   # blue
    (1.00, 0.85, 0.40),   # orange
    (0.85, 0.55, 1.00),   # purple
    (0.40, 0.85, 1.00),   # cyan
    (1.00, 0.65, 0.85),   # pink
    (0.65, 0.85, 0.40),   # lime
    (0.85, 0.70, 0.55),   # tan
    (0.70, 0.70, 0.70),   # grey
]


def render_segments_svg(
    mol_with_h: Any,
    segments: List[Segment],
    width: int = 700,
    height: int = 500,
    show_atom_numbers: bool = False,
    bold_shared: bool = True,
) -> str:
    """各 segment を別色で塗った heavy-atom-only SVG を生成する。

    - segment 色: SEG_COLOR_PALETTE で順に割当 (segment_id mod len)
    - shared atom は太い黒縁取りで強調 (`bold_shared=True`)
    - `show_atom_numbers=True` で各 heavy atom に番号、shared なら "*" suffix
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit.Chem.Draw import rdMolDraw2D

    try:
        heavy_mol = Chem.RemoveHs(mol_with_h)
    except Exception:
        heavy_mol = mol_with_h

    try:
        AllChem.Compute2DCoords(heavy_mol)
    except Exception:
        pass

    heavy_atom_origs = [
        a.GetIdx() for a in mol_with_h.GetAtoms() if a.GetAtomicNum() != 1
    ]
    orig_to_heavy = {orig: i for i, orig in enumerate(heavy_atom_origs)}

    atom_to_segs: dict = {}
    for s in segments:
        for a in s.atom_indices:
            atom_to_segs.setdefault(a, []).append(s.segment_id)

    hl_atoms: List[int] = []
    atom_colors: dict = {}
    for a_orig, sids in atom_to_segs.items():
        a_h = orig_to_heavy.get(a_orig)
        if a_h is None:
            continue
        hl_atoms.append(a_h)
        atom_colors[a_h] = SEG_COLOR_PALETTE[sids[0] % len(SEG_COLOR_PALETTE)]

    # atom note: 数字 + shared "*"
    for a in heavy_mol.GetAtoms():
        if a.HasProp("atomNote"):
            a.ClearProp("atomNote")
    if show_atom_numbers:
        heavy_to_orig = {h: o for o, h in orig_to_heavy.items()}
        for a in heavy_mol.GetAtoms():
            h_idx = a.GetIdx()
            a_orig = heavy_to_orig.get(h_idx)
            note = str(h_idx)
            if a_orig is not None and len(atom_to_segs.get(a_orig, [])) > 1:
                note += "*"
            a.SetProp("atomNote", note)

    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    rdMolDraw2D.PrepareAndDrawMolecule(
        drawer, heavy_mol,
        highlightAtoms=hl_atoms,
        highlightAtomColors=atom_colors,
    )

    # shared atom の黒縁取り (SVG 後処理)
    shared_circles: List[str] = []
    if bold_shared:
        for a_orig, sids in atom_to_segs.items():
            if len(sids) < 2:
                continue
            a_h = orig_to_heavy.get(a_orig)
            if a_h is None:
                continue
            try:
                pt = drawer.GetDrawCoords(a_h)
                shared_circles.append(
                    f'<circle cx="{pt.x:.2f}" cy="{pt.y:.2f}" r="13" '
                    f'fill="none" stroke="black" stroke-width="2.5"/>'
                )
            except Exception:
                pass

    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    if shared_circles:
        svg = svg.replace("</svg>", "\n".join(shared_circles) + "\n</svg>")
    return svg
