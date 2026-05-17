"""
colorizer.py
------------
Write H-bond classification results back to a COGNAC UDF/BDF so that
gourmet renders the participating molecules with three distinct colors.

Strategy: rename ``Set_of_Molecules.molecule[i].Mol_Name`` per classification,
then add ``Draw_Attributes.Molecule[]`` entries colored for each new Mol_Name.

Schema (verified via queryDefine on Draw_Attributes.Molecule[]):
- ``color`` is a ``select`` type with 9 named colors only:
    Red, Green, Blue, Magenta, Cyan, Yellow, White, Black, Gray.
  (The RGBA float syntax `[r,g,b,a]` in GOURMET docs applies only to the
   Python viewer's drawing functions like ``sphere()``, NOT to the
   ``Draw_Attributes`` schema.)
- ``transparency`` is a float [0.0, 1.0]; **1.0 = NO transparency (opaque)**.
- ``radius`` is a float (sphere radius scale).

Output file: ``<bdf_stem>_colored.bdf`` (or user-supplied path).
"""
from __future__ import annotations

import os
import shutil
from dataclasses import dataclass
from typing import Dict, List, Sequence, Tuple

from .classifier import ClassificationResult, MolRole


VALID_COLORS = {
    "Red", "Green", "Blue", "Magenta", "Cyan", "Yellow", "White", "Black", "Gray",
}


@dataclass
class DrawAttribute:
    """One Mol_Name -> color/transparency/radius mapping.

    ``color`` must be one of the 9 GOURMET named colors (see VALID_COLORS).
    ``transparency`` 1.0 = opaque, 0.0 = fully transparent.
    """
    mol_name: str
    color: str
    transparency: float = 1.0
    radius: float = 1.0

    def __post_init__(self):
        if self.color not in VALID_COLORS:
            raise ValueError(
                f"Invalid color {self.color!r}. Must be one of {sorted(VALID_COLORS)}"
            )


DEFAULT_COLORS: Dict[str, DrawAttribute] = {
    "dual":   DrawAttribute("_DUAL",   "Red",   transparency=1.0, radius=1.0),
    "single": DrawAttribute("_SINGLE", "Blue",  transparency=1.0, radius=1.0),
    "free":   DrawAttribute("_FREE",   "Gray",  transparency=0.3, radius=1.0),
}


def _make_role_to_mol_name(
    base_mol_name: str,
    attrs: Dict[str, DrawAttribute],
) -> Dict[str, str]:
    """Build role -> renamed Mol_Name. E.g., 'dual' -> 'IMC_dual'."""
    out = {}
    for role, attr in attrs.items():
        if attr.mol_name.startswith("_"):
            out[role] = f"{base_mol_name}{attr.mol_name}"
        else:
            out[role] = attr.mol_name
    return out


def colorize_udf(
    input_path: str,
    output_path: str,
    result: ClassificationResult,
    base_mol_name: str = "IMC",
    color_attrs: Dict[str, DrawAttribute] = None,
) -> str:
    """Rewrite UDF with renamed Mol_Names and Draw_Attributes.

    Parameters
    ----------
    input_path : source UDF/BDF
    output_path : destination path
    result : ClassificationResult from classifier.classify()
    base_mol_name : prefix for renamed groups (e.g. "IMC" -> "IMC_dual")
    color_attrs : dict of role -> DrawAttribute. Default = DEFAULT_COLORS.

    Returns
    -------
    str : output_path
    """
    if color_attrs is None:
        color_attrs = DEFAULT_COLORS
    role_to_name = _make_role_to_mol_name(base_mol_name, color_attrs)

    if input_path != output_path:
        shutil.copy(input_path, output_path)

    from UDFManager import UDFManager
    u = UDFManager(output_path)
    n_records = u.totalRecord()

    # Rewrite Mol_Name per molecule in each record (Mol_Name lives in
    # Set_of_Molecules which may be record-bound depending on UDF schema).
    # We update record by record AND the common (record=-1) for safety.
    for rec in [-1] + list(range(n_records)):
        try:
            u.jump(rec)
        except Exception:
            continue
        n_mol = u.size("Set_of_Molecules.molecule[]")
        if n_mol is None or n_mol <= 0:
            continue
        for role_info in result.roles:
            new_name = role_to_name[role_info.role]
            try:
                u.put(new_name,
                      f"Set_of_Molecules.molecule[{role_info.mol_index}].Mol_Name")
            except Exception:
                pass

    # Draw_Attributes.Molecule[] -- add 3 entries
    u.jump(-1)
    existing = u.size("Draw_Attributes.Molecule[]") or 0
    for idx, role in enumerate(("dual", "single", "free")):
        attr = color_attrs[role]
        new_name = role_to_name[role]
        slot = existing + idx
        try:
            u.put(new_name, f"Draw_Attributes.Molecule[{slot}].Mol_Name")
            u.put(attr.color, f"Draw_Attributes.Molecule[{slot}].color")
            u.put(attr.transparency,
                  f"Draw_Attributes.Molecule[{slot}].transparency")
            u.put(attr.radius, f"Draw_Attributes.Molecule[{slot}].radius")
        except Exception as e:
            print(f"Warning: failed to set Draw_Attributes[{slot}]: {e}")

    u.write(output_path)
    return output_path
