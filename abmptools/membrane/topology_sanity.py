"""Topology sanity-check helpers for membrane builder output.

Parses GROMACS ``system.top`` (and chain itp files for the CHARMM
backend, where pdb2gmx produces per-chain ``system_<Other|Protein>_chain_X.itp``
includes) and reports counts and charges per molecule type.

This is *content* validation — checks that the output topology is
chemically sensible (electroneutral, expected lipid / peptide / water /
ion counts), independent of whether ``grompp`` accepts it. Used by the
sanity-check integration test (``tests/integration/test_topology_sanity.py``)
and as a CLI helper (``python -m abmptools.membrane.topology_sanity <build_dir>``).
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Tuple


@dataclass
class MolBlock:
    name: str
    count: int
    atoms_per_mol: int
    charge_per_mol: float
    residues: List[str] = field(default_factory=list)


@dataclass
class TopologySummary:
    blocks: List[MolBlock] = field(default_factory=list)
    total_atoms: int = 0
    total_charge: float = 0.0

    def by_residue(self) -> Dict[str, Tuple[int, float]]:
        """Return {resname: (n_molecules, total_charge)}."""
        out: Dict[str, Tuple[int, float]] = {}
        for b in self.blocks:
            for r in b.residues:
                n, c = out.get(r, (0, 0.0))
                out[r] = (n + b.count, c + b.charge_per_mol * b.count)
        return out


def _iter_atoms_block(text: str):
    """Yield (atomname, resname, charge) tuples from every [atoms] block in text."""
    in_atoms = False
    for line in text.splitlines():
        s = line.strip()
        if s.startswith("[ ") and s.endswith(" ]"):
            in_atoms = (s == "[ atoms ]")
            continue
        if not in_atoms:
            continue
        if s.startswith("#"):
            in_atoms = False  # preprocessor ends section
            continue
        if not s or s.startswith(";"):
            continue
        parts = s.split()
        if len(parts) < 7:
            continue
        try:
            yield parts[4], parts[3], float(parts[6])
        except (ValueError, IndexError):
            continue


def _parse_inline_top(top_text: str) -> Dict[str, List[Tuple[str, str, float]]]:
    """Find each [moleculetype]+[atoms] inline block (AMBER-style top)."""
    blocks: Dict[str, List[Tuple[str, str, float]]] = {}
    moltype_pat = re.compile(
        r"^\[ moleculetype \]\n(?:;.*\n)*\s*(\S+)\s+\d+", re.MULTILINE
    )
    starts = [(m.start(), m.group(1)) for m in moltype_pat.finditer(top_text)]
    for i, (start, name) in enumerate(starts):
        end = starts[i + 1][0] if i + 1 < len(starts) else len(top_text)
        section = top_text[start:end]
        am = re.search(r"^\[ atoms \](.*?)(?=^\[|\Z)", section, re.MULTILINE | re.DOTALL)
        if not am:
            continue
        atoms = []
        for line in am.group(1).splitlines():
            s = line.strip()
            if not s or s.startswith(";") or s.startswith("#"):
                continue
            parts = s.split()
            if len(parts) >= 7:
                try:
                    atoms.append((parts[4], parts[3], float(parts[6])))
                except (ValueError, IndexError):
                    continue
        blocks[name] = atoms
    return blocks


def _parse_molecules_section(top_text: str) -> List[Tuple[str, int]]:
    m = re.search(r"^\[ molecules \](.*?)\Z", top_text, re.MULTILINE | re.DOTALL)
    if not m:
        return []
    out = []
    for line in m.group(1).splitlines():
        s = line.strip()
        if not s or s.startswith(";"):
            continue
        parts = s.split()
        if len(parts) >= 2:
            try:
                out.append((parts[0], int(parts[1])))
            except ValueError:
                continue
    return out


def summarize_topology(build_dir: Path) -> TopologySummary:
    """Parse ``system.top`` (and chain itps if present) and return summary.

    Handles both layouts:
    - AMBER: monolithic ``system.top`` with inline [moleculetype]+[atoms] blocks
    - CHARMM (pdb2gmx): ``system.top`` references ``system_*_chain_X.itp``
      via ``#include``; per-chain itps hold the [atoms] blocks.
    """
    build_dir = Path(build_dir)
    top_path = build_dir / "system.top"
    top_text = top_path.read_text()
    mols = _parse_molecules_section(top_text)
    inline_blocks = _parse_inline_top(top_text)

    summary = TopologySummary()
    for name, count in mols:
        atoms: List[Tuple[str, str, float]] = []
        if name in inline_blocks:
            atoms = inline_blocks[name]
        else:
            chain_itp = build_dir / f"system_{name}.itp"
            if chain_itp.exists():
                atoms = list(_iter_atoms_block(chain_itp.read_text()))
        if not atoms:
            continue
        residues = sorted({a[1] for a in atoms})
        charge_per_mol = sum(c for _, _, c in atoms)
        summary.blocks.append(
            MolBlock(
                name=name,
                count=count,
                atoms_per_mol=len(atoms),
                charge_per_mol=charge_per_mol,
                residues=residues,
            )
        )
        summary.total_atoms += len(atoms) * count
        summary.total_charge += charge_per_mol * count
    return summary


def format_summary(summary: TopologySummary) -> str:
    lines = []
    for b in summary.blocks:
        chain_atoms = b.atoms_per_mol * b.count
        chain_charge = b.charge_per_mol * b.count
        residues = "/".join(b.residues)
        lines.append(
            f"  {b.name:25s} ×{b.count:5d}: "
            f"{chain_atoms:6d} atoms, charge={chain_charge:+7.3f}, residues={residues}"
        )
    lines.append(f"  TOTAL: {summary.total_atoms} atoms, charge={summary.total_charge:+.3f}")
    return "\n".join(lines)


def main(argv=None) -> int:
    import argparse
    p = argparse.ArgumentParser(
        description="Summarize a GROMACS system.top from membrane builder output."
    )
    p.add_argument("build_dir", type=Path, help="build/ directory containing system.top")
    args = p.parse_args(argv)
    s = summarize_topology(args.build_dir)
    print(format_summary(s))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
