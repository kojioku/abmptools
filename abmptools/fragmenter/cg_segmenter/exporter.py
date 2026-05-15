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
