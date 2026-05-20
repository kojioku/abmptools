"""
colorizer.py
------------
Write H-bond classification results back to a COGNAC UDF/BDF so that
gourmet renders the participating functional groups with role-based colors.

Two strategies:

(A) ``colorize_udf`` — **Mol_Name rename** + ``Draw_Attributes.Molecule[]``.
    Renames ``Set_of_Molecules.molecule[i].Mol_Name`` per molecule role
    (``IMC_DUAL`` / ``IMC_SINGLE`` / ``IMC_FREE``) and adds 3 entries to
    ``Draw_Attributes.Molecule[]``. Works in gourmet show + OCTA post-render
    but breaks J-OCTA pre-render (Mol_Name lookup fails on renamed names).

(B) ``colorize_udf_action`` — **Python action (.act)** + Mol_Name kept.
    Emits a `_show.act` file that the gourmet ``show`` action loads to
    paint each functional group (carboxyl / amide) individually with
    sphere overlays on its atoms. The BDF text header's ``Action:`` line
    is patched to reference the new .act file (binary section untouched).
    Compatible with J-OCTA pre-render (no Mol_Name rename) and lets one
    molecule contribute to multiple roles via different functional groups.

Schema notes (Draw_Attributes.Molecule[] for strategy A):
- ``color`` is a ``select`` type with 9 named colors only:
    Red, Green, Blue, Magenta, Cyan, Yellow, White, Black, Gray.
- ``transparency`` is a float [0.0, 1.0]; **1.0 = NO transparency (opaque)**.
- ``radius`` is a float (sphere radius scale).
"""
from __future__ import annotations

import os
import re
import shutil
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

from .classifier import (
    AmideRole, CarboxylRole, ClassificationResult,
    FunctionalGroupClassification, MolRole,
)
from .functional_groups import AmideGroup, CarboxylGroup


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


# ---------------------------------------------------------------------------
# Strategy (B): Python action (.act) coloring — Mol_Name preserved
# ---------------------------------------------------------------------------

# RGBA float palette for action-mode coloring (gourmet draw functions accept
# `[r, g, b, a]` or `[r, g, b, a, radius]`).
DEFAULT_ACTION_COLORS: Dict[str, List[float]] = {
    "dual":   [1.0, 0.0, 0.0, 1.0],     # red
    "single": [0.0, 0.0, 1.0, 1.0],     # blue
    "free":   [0.7, 0.7, 0.7, 0.3],     # translucent gray
    "accept": [0.0, 0.8, 0.8, 1.0],     # cyan (amide acceptor)
}


def _patch_action_header(bdf_path: str, act_basename: str) -> None:
    """Append ``act_basename`` to the BDF text header's ``Action:`` field.

    Only the leading text section (before ``\\end{header}``) is rewritten;
    the binary section is left byte-for-byte untouched. Idempotent if
    ``act_basename`` is already listed.
    """
    with open(bdf_path, "rb") as f:
        data = f.read()
    end_marker = b"\\end{header}"
    idx = data.find(end_marker)
    if idx < 0:
        raise RuntimeError(
            f"Cannot find \\end{{header}} in {bdf_path}; "
            "is this a valid UDF/BDF file?"
        )
    header = data[:idx].decode("ascii")
    rest = data[idx:]

    m = re.search(r'^Action:"([^"]*)"', header, re.MULTILINE)
    if m:
        existing = m.group(1)
        parts = [p for p in existing.split(";") if p]
        if act_basename in parts:
            return  # already wired
        parts.append(act_basename)
        new_line = f'Action:"{";".join(parts)}"'
        header = header[:m.start()] + new_line + header[m.end():]
    else:
        # No Action line — insert right before \end{data} of the header block.
        header = header.replace(
            "\\end{data}", f'Action:"{act_basename}"\n\\end{{data}}', 1
        )

    with open(bdf_path, "wb") as f:
        f.write(header.encode("ascii") + rest)


def _render_show_body_lines(
    classification: FunctionalGroupClassification,
    carboxyls: List[CarboxylGroup],
    amides: List[AmideGroup],
    colors: Dict[str, List[float]],
) -> List[str]:
    """Common body lines (data + draw loop) shared by .act and .py outputs.

    The body uses GOURMET draw functions (``sphere`` / ``line``) and the
    UDF accessors ``size(path, [indices])`` / ``get(path, [indices])`` —
    confirmed compatible with both OCTA gourmet and J-OCTA Viewer Python
    panel.
    """
    carb_entries = []
    for cg, role in zip(carboxyls, classification.carboxyl_roles):
        carb_entries.append(
            f"    ({cg.mol_index}, {cg.c_atom}, {cg.o_atom}, "
            f"{cg.oh_atom}, {cg.ho_atom}, '{role.role}'),"
        )
    amide_entries = []
    for ag, role in zip(amides, classification.amide_roles):
        amide_entries.append(
            f"    ({ag.mol_index}, {ag.c_atom}, {ag.o_atom}, "
            f"{ag.n_atom}, '{role.role}'),"
        )

    color_lines = []
    for role, c in colors.items():
        rgba = ", ".join(f"{v:.3f}" for v in c)
        color_lines.append(f"    '{role}': [{rgba}],")

    return [
        "carboxyls = [",
        *carb_entries,
        "]",
        "amides = [",
        *amide_entries,
        "]",
        "role_color = {",
        *color_lines,
        "}",
        "sphere_radius = 0.55",
        "",
        "# Helper: get atom position by global indices",
        "def _pos(mi, ai):",
        "    return get('Structure.Position.mol[].atom[]', [mi, ai])",
        "",
        "# Backbone bonds in light gray",
        "n_mol = size('Set_of_Molecules.molecule[]')",
        "for mi in range(n_mol):",
        "    n_bonds = size('Set_of_Molecules.molecule[].bond[]', [mi])",
        "    for bi in range(n_bonds):",
        "        a1 = get('Set_of_Molecules.molecule[].bond[].atom1', [mi, bi])",
        "        a2 = get('Set_of_Molecules.molecule[].bond[].atom2', [mi, bi])",
        "        line(_pos(mi, a1), _pos(mi, a2), [0.5, 0.5, 0.5, 0.4])",
        "",
        "# Carboxyl atoms (c, o, oh, ho) coloured by role",
        "for entry in carboxyls:",
        "    mi, c_at, o_at, oh_at, ho_at, role = entry",
        "    col = role_color.get(role, role_color['free']) + [sphere_radius]",
        "    for ai in (c_at, o_at, oh_at, ho_at):",
        "        sphere(_pos(mi, ai), col)",
        "",
        "# Amide atoms (c, o, n) coloured by accept/free",
        "for entry in amides:",
        "    mi, c_at, o_at, n_at, role = entry",
        "    if role == 'accept':",
        "        key = 'accept'",
        "    else:",
        "        key = 'free'",
        "    col = role_color[key] + [sphere_radius]",
        "    for ai in (c_at, o_at, n_at):",
        "        sphere(_pos(mi, ai), col)",
    ]


def _render_show_action(
    classification: FunctionalGroupClassification,
    carboxyls: List[CarboxylGroup],
    amides: List[AmideGroup],
    colors: Dict[str, List[float]],
) -> str:
    """Build the .act content: body wrapped in ``autorun: showHbond()``.

    OCTA gourmet only — J-OCTA Viewer has been observed to crash on this
    autorun form; use :func:`_render_show_python_panel` instead for J-OCTA.
    """
    body = _render_show_body_lines(classification, carboxyls, amides, colors)
    lines = [
        "# Auto-generated by abmptools.hbond.colorize_udf_action.",
        "# Renders functional groups with role-based color overlay.",
        "# Carboxyl: dual=red, single=blue, free=gray.",
        "# Amide:   accept=cyan, free=gray.",
        "# autorun fires on file open and on each record change, so the",
        "# overlay refreshes automatically when navigating frames.",
        "# NOTE: OCTA gourmet only. J-OCTA Viewer may crash on this form;",
        "# use the companion <prefix>_show.py via the Python panel instead.",
        "",
        "autorun : showHbond() : \\begin",
        *body,
        "\\end",
        "",
    ]
    return "\n".join(lines)


def _render_show_python_panel(
    classification: FunctionalGroupClassification,
    carboxyls: List[CarboxylGroup],
    amides: List[AmideGroup],
    colors: Dict[str, List[float]],
    target_bdf_basename: str,
) -> str:
    """Build a .py script for direct execution in the Viewer's Python panel.

    Workflow (J-OCTA Viewer):
      1. File → Open ``<target_bdf_basename>`` (Mol_Name preserved)
      2. Python panel → ``Load…`` this .py (or copy-paste the body)
      3. Click ``Run``

    No ``autorun``/``\\begin`` wrapper — just a plain Python script body,
    which the J-OCTA Viewer accepts via the Python panel even when the
    autorun action form crashes.
    """
    body = _render_show_body_lines(classification, carboxyls, amides, colors)
    lines = [
        "# Auto-generated by abmptools.hbond.write_show_python_script.",
        "# Renders functional groups with role-based color overlay.",
        "# Carboxyl: dual=red, single=blue, free=gray.",
        "# Amide:   accept=cyan, free=gray.",
        "#",
        "# How to use (J-OCTA Viewer / OCTA gourmet):",
        f"#   1. File -> Open {target_bdf_basename} (Mol_Name preserved copy)",
        "#   2. Python panel -> Load this file (or paste its contents)",
        "#   3. Click Run",
        "#",
        "# Companion to the .act autorun overlay; use this .py path when",
        "# the .act form is incompatible with the host viewer.",
        "",
        *body,
        "",
    ]
    return "\n".join(lines)


def colorize_udf_action(
    input_path: str,
    output_bdf_path: str,
    output_act_path: str,
    classification: FunctionalGroupClassification,
    carboxyls: List[CarboxylGroup],
    amides: List[AmideGroup],
    action_colors: Optional[Dict[str, List[float]]] = None,
) -> Tuple[str, str]:
    """Emit a Mol_Name-preserved BDF + autorun action that paints each
    functional group with a role-based color overlay.

    For J-OCTA Viewer (which can crash on the autorun form), pair this with
    :func:`write_show_python_script` to emit a Python-panel-loadable script.

    Returns
    -------
    (bdf_path, act_path)
    """
    if action_colors is None:
        action_colors = DEFAULT_ACTION_COLORS

    if input_path != output_bdf_path:
        shutil.copy(input_path, output_bdf_path)

    act_text = _render_show_action(
        classification, carboxyls, amides, action_colors,
    )
    with open(output_act_path, "w", encoding="ascii") as f:
        f.write(act_text)

    act_basename = os.path.basename(output_act_path)
    _patch_action_header(output_bdf_path, act_basename)

    return output_bdf_path, output_act_path


def write_hbond_attributes(
    bdf_path: str,
    classification: FunctionalGroupClassification,
    carboxyls: List[CarboxylGroup],
    amides: List[AmideGroup],
    attribute_name: str = "hbond",
    value_map: Optional[Dict[str, str]] = None,
) -> str:
    """Append per-atom ``Attributes[]`` entries tagging each functional-group
    atom with its role (carboxyl: Dual/Single/Free, amide: Accept/Free).

    The new entry is **appended to the tail** of each atom's ``Attributes[]``;
    existing entries (the ``Name='1' Value='molecular:<id>'`` mol identifier
    that J-OCTA writes by default, plus any other pre-existing fields) are
    left untouched.

    Atoms tagged:
    - carboxyl ``c`` / ``o`` / ``oh`` / ``ho`` — value from COOH role
    - amide   ``c`` / ``o`` / ``n``        — value from amide role
    Other atoms are not modified. This lets J-OCTA filter / highlight by the
    ``hbond`` Attribute (e.g. show only ``Dual`` atoms).

    Default value_map keeps Title-case strings consistent across roles:
    ``{'dual':'Dual','single':'Single','free':'Free','accept':'Accept'}``.
    """
    if value_map is None:
        value_map = {
            "dual": "Dual", "single": "Single",
            "free": "Free", "accept": "Accept",
        }

    atom_value: Dict[Tuple[int, int], str] = {}
    for cg, role in zip(carboxyls, classification.carboxyl_roles):
        v = value_map.get(role.role, role.role)
        for ai in (cg.c_atom, cg.o_atom, cg.oh_atom, cg.ho_atom):
            atom_value[(cg.mol_index, ai)] = v
    for ag, role in zip(amides, classification.amide_roles):
        v = value_map.get(role.role, role.role)
        for ai in (ag.c_atom, ag.o_atom, ag.n_atom):
            # carboxyl > amide priority on rare overlap (shouldn't happen)
            atom_value.setdefault((ag.mol_index, ai), v)

    from UDFManager import UDFManager
    u = UDFManager(bdf_path)
    n_records = u.totalRecord()
    for rec in [-1] + list(range(n_records)):
        try:
            u.jump(rec)
        except Exception:
            continue
        n_mol = u.size("Set_of_Molecules.molecule[]")
        if not n_mol:
            continue
        for mi in range(n_mol):
            n_at = u.size(f"Set_of_Molecules.molecule[{mi}].atom[]") or 0
            for ai in range(n_at):
                key = (mi, ai)
                if key not in atom_value:
                    continue
                base = f"Set_of_Molecules.molecule[{mi}].atom[{ai}].Attributes"
                cur_n = u.size(f"{base}[]") or 0
                # Skip if we've already appended (idempotent re-runs)
                already = False
                for j in range(cur_n):
                    try:
                        if u.get(f"{base}[{j}].Name") == attribute_name:
                            already = True
                            break
                    except Exception:
                        pass
                if already:
                    continue
                slot = cur_n
                try:
                    u.put(attribute_name, f"{base}[{slot}].Name")
                    u.put(atom_value[key], f"{base}[{slot}].Value")
                except Exception as e:
                    print(f"Warn: failed Attributes[{slot}] at "
                          f"mol[{mi}].atom[{ai}]: {e}")
    u.write(bdf_path)
    return bdf_path


def write_show_python_script(
    output_py_path: str,
    classification: FunctionalGroupClassification,
    carboxyls: List[CarboxylGroup],
    amides: List[AmideGroup],
    target_bdf_basename: str,
    action_colors: Optional[Dict[str, List[float]]] = None,
) -> str:
    """Emit a plain Python script for direct Python-panel execution.

    Unlike the .act autorun form, this script is meant to be loaded into the
    Viewer's Python panel (Load… or paste) and executed manually via Run.
    It contains no UDF Action wiring, so the host BDF is untouched.

    Use it on the Mol_Name-preserved copy (``<prefix>.bdf``); the .act
    autorun form often crashes the J-OCTA Viewer but the Python-panel form
    runs fine.
    """
    if action_colors is None:
        action_colors = DEFAULT_ACTION_COLORS
    text = _render_show_python_panel(
        classification, carboxyls, amides, action_colors, target_bdf_basename,
    )
    with open(output_py_path, "w", encoding="ascii") as f:
        f.write(text)
    return output_py_path
