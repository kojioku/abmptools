# -*- coding: utf-8 -*-
"""
abmptools.amorphous.ndx_writer
--------------------------------
Write GROMACS .ndx (index) files for multi-component amorphous systems.
"""
from __future__ import annotations

import os
from pathlib import Path
from typing import Dict, List


def _format_index_group(name: str, indices: List[int], columns: int = 15) -> str:
    """Format a single NDX group."""
    lines = [f"[ {name} ]"]
    for i in range(0, len(indices), columns):
        chunk = indices[i:i + columns]
        lines.append(" ".join(f"{idx:>6d}" for idx in chunk))
    return "\n".join(lines)


def write_ndx(
    component_names: List[str],
    atom_counts_per_mol: List[int],
    mol_counts: List[int],
    output_path: str,
) -> str:
    """Write a GROMACS .ndx file with System + per-component groups.

    Parameters
    ----------
    component_names : list of str
        Name for each component (e.g. ["API", "Polymer"]).
    atom_counts_per_mol : list of int
        Number of atoms per single molecule of each component.
    mol_counts : list of int
        Number of molecules for each component.
    output_path : str
        Path to write the .ndx file.

    Returns
    -------
    str
        Absolute path to the written file.
    """
    all_indices: List[int] = []
    groups: Dict[str, List[int]] = {}

    atom_id = 1
    for comp_name, n_atoms, n_mol in zip(component_names, atom_counts_per_mol, mol_counts):
        comp_indices: List[int] = []
        for _ in range(n_mol):
            for _ in range(n_atoms):
                comp_indices.append(atom_id)
                atom_id += 1
        groups[comp_name] = comp_indices
        all_indices.extend(comp_indices)

    sections = [_format_index_group("System", all_indices)]
    for comp_name in component_names:
        sections.append(_format_index_group(comp_name, groups[comp_name]))

    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    text = "\n".join(sections) + "\n"
    Path(output_path).write_text(text)
    return str(Path(output_path).resolve())
