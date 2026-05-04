# -*- coding: utf-8 -*-
"""
abmptools.cg.membrane.topology_composer
----------------------------------------
Post-processing helpers for insane's outputs.

Insane writes:
  - ``bilayer.gro`` with ions named ``NA+`` / ``CL-`` (resname & atomname).
  - ``insane_topol.top`` with a generic ``#include "martini.itp"``, a
    placeholder ``Protein`` moleculetype name, and ``NA+`` / ``CL-`` in the
    ``[ molecules ]`` table.

Martini 3 ITPs that abmptools requires (``martini_v3.0.0_ions_v1.itp``) use
``NA`` / ``CL`` (no charge sign) as moleculetype/atomname. ``gmx grompp``
flags the mismatch as a warning that we'd rather avoid; we also need to
substitute the M3 ITP includes and the actual peptide moleculetype name.

This module provides:

- :func:`compose_topology` -- rewrites insane's ``topol.top`` into a
  GROMACS-ready Martini 3 topology with full ITP includes.
- :func:`normalize_ion_atom_names_gro` -- substitutes ``NA+``/``CL-`` with
  ``NA``/``CL`` in a ``.gro``, preserving column alignment.
- :func:`get_moleculetype_name_from_itp` -- reads the first
  ``[ moleculetype ]`` directive in a peptide ITP (cg.peptide writes
  ``molecule_0`` by default).
"""
from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import List, Optional

logger = logging.getLogger(__name__)


# ITPs to include in the composed topology. Order matters: protein /
# solvents / ions / phospholipids -- but martini_v3.0.0.itp must be first
# because it defines atom types referenced by the others.
_M3_ITP_INCLUDES = [
    "martini_v3.0.0.itp",
    "martini_v3.0.0_solvents_v1.itp",
    "martini_v3.0.0_ions_v1.itp",
    "martini_v3.0.0_phospholipids_v1.itp",
]


# ---------------------------------------------------------------------------
# Peptide moleculetype name extraction (re-used from cg.peptide convention)
# ---------------------------------------------------------------------------

def get_moleculetype_name_from_itp(itp_path: Path) -> str:
    """Return the first moleculetype name in *itp_path*.

    Parses the first ``[ moleculetype ]`` block in the ITP and returns the
    first non-comment, non-blank token of the next non-comment, non-blank
    line. cg.peptide writes ``molecule_0`` here.
    """
    text = Path(itp_path).read_text()
    in_block = False
    for raw in text.splitlines():
        line = raw.split(";", 1)[0].strip()
        if not line:
            continue
        if line.startswith("[") and line.endswith("]"):
            section = line[1:-1].strip().lower()
            in_block = section == "moleculetype"
            continue
        if in_block:
            return line.split()[0]
    raise ValueError(f"No [ moleculetype ] block in {itp_path}")


# ---------------------------------------------------------------------------
# topol.top rewriter
# ---------------------------------------------------------------------------

def compose_topology(
    insane_top: Path,
    output_top: Path,
    *,
    peptide_itp_path: Path,
    peptide_count: int = 1,
    peptide_itp_relpath: Optional[str] = None,
    system_title: str = "Martini 3 peptide-membrane system",
) -> Path:
    """Rewrite insane's topology with Martini 3 4-ITP includes + peptide ITP.

    Parameters
    ----------
    insane_top
        The raw ``insane_topol.top`` produced by :func:`insane_runner.run_insane`.
    output_top
        Path to write the composed topology.
    peptide_itp_path
        Path to the peptide ITP (used to extract the actual moleculetype
        name; cg.peptide writes ``molecule_0``).
    peptide_count
        How many peptide molecules to declare in ``[ molecules ]``. Default 1.
    peptide_itp_relpath
        ``#include`` path written into the topology for the peptide ITP. If
        None, the absolute path of ``peptide_itp_path`` is used. For the
        builder, this is typically a relative path like
        ``"molecules/<name>/<name>.itp"``.
    system_title
        ``[ system ]`` title line.

    Returns
    -------
    Path to the written topology.
    """
    insane_top = Path(insane_top)
    output_top = Path(output_top)
    raw = insane_top.read_text()

    mol_name = get_moleculetype_name_from_itp(Path(peptide_itp_path))
    pep_inc = (
        peptide_itp_relpath
        if peptide_itp_relpath is not None
        else str(Path(peptide_itp_path).resolve())
    )

    # 1) replace the FF #include block with the Martini 3 4-ITP block + peptide ITP
    ff_block = "\n".join(f'#include "{itp}"' for itp in _M3_ITP_INCLUDES)
    ff_block += f'\n#include "{pep_inc}"'

    # Insane's typical line is `#include "martini.itp"`.
    new = re.sub(
        r'#include\s+"martini[^"]*\.itp"',
        ff_block,
        raw,
        count=1,
    )
    if new == raw:
        # Fallback: no include line found; prepend the FF block.
        logger.warning(
            "No '#include \"martini*.itp\"' line found in %s; "
            "prepending Martini 3 includes.",
            insane_top,
        )
        new = ff_block + "\n\n" + raw

    # 2) Rewrite [ system ] body to a meaningful title.
    new = _replace_system_block(new, system_title)

    # 3) In [ molecules ]: substitute ``Protein`` -> mol_name (count line),
    #    ``NA+`` -> ``NA``, ``CL-`` -> ``CL``.
    new = _rewrite_molecules_block(
        new,
        protein_to=mol_name,
        peptide_count=peptide_count,
    )

    output_top.parent.mkdir(parents=True, exist_ok=True)
    output_top.write_text(new)
    logger.info(
        "compose_topology: wrote %s (peptide moleculetype=%s)",
        output_top, mol_name,
    )
    return output_top


def _replace_system_block(text: str, title: str) -> str:
    """Substitute the body of the ``[ system ]`` directive."""
    pattern = re.compile(
        r"(\[\s*system\s*\][ \t]*\n)"
        r"((?:[^\[\n][^\n]*\n)+)",
        flags=re.IGNORECASE,
    )
    return pattern.sub(rf"\1{title}\n", text, count=1)


def _rewrite_molecules_block(
    text: str, *, protein_to: str, peptide_count: int,
) -> str:
    """Rewrite ``Protein`` / ``NA+`` / ``CL-`` in ``[ molecules ]``.

    Only the ``[ molecules ]`` directive body is touched. Outside the block,
    occurrences of ``NA+`` / ``CL-`` (e.g. inside comments) are preserved.
    """
    # Find the [ molecules ] header.
    m = re.search(r"\[\s*molecules\s*\][ \t]*\n", text, flags=re.IGNORECASE)
    if not m:
        # No molecules block; leave text as-is and let downstream tools
        # complain. Don't silently invent content here.
        logger.warning("no [ molecules ] block found; topology untouched")
        return text

    header_end = m.end()
    body = text[header_end:]
    head = text[:header_end]

    new_lines: List[str] = []
    for line in body.splitlines(keepends=True):
        # Stop at the next directive (next ``[`` at line start). Keep the
        # rest of the file untouched.
        if line.lstrip().startswith("["):
            new_lines.append(line)
            new_lines.append(body.split(line, 1)[1] if line in body else "")
            break

        replaced = _replace_molecules_line(
            line, protein_to=protein_to, peptide_count=peptide_count,
        )
        new_lines.append(replaced)
    else:
        # Loop exhausted without breaking (no further section after molecules).
        pass

    rebuilt = "".join(new_lines)
    # Re-attach any unmodified suffix that was passed through above.
    if not rebuilt.endswith("\n") and not rebuilt.endswith("\r"):
        rebuilt += ""
    return head + rebuilt if rebuilt else head + body


def _replace_molecules_line(
    line: str, *, protein_to: str, peptide_count: int,
) -> str:
    """Per-line substitutions inside ``[ molecules ]``.

    Comment-only and blank lines pass through verbatim.
    """
    stripped = line.split(";", 1)[0]
    if not stripped.strip():
        return line

    parts = stripped.split()
    if len(parts) < 2:
        return line
    mol, rest = parts[0], parts[1:]

    # Substitute moleculetype name
    if mol == "Protein":
        mol_out = protein_to
        # Force the count to peptide_count when rewriting Protein
        rest = [str(peptide_count)] + rest[1:]
    elif mol == "NA+":
        mol_out = "NA"
    elif mol == "CL-":
        mol_out = "CL"
    else:
        return line  # unchanged

    # Preserve trailing comment if any.
    comment = ""
    if ";" in line:
        comment = ";" + line.split(";", 1)[1].rstrip("\n")
    new_body = f"{mol_out:<15s} {rest[0]}"
    suffix = "\n" if line.endswith("\n") else ""
    return new_body + ("  " + comment if comment else "") + suffix


# ---------------------------------------------------------------------------
# .gro ion atom-name normalization
# ---------------------------------------------------------------------------

# Valid ion substitutions: 3-char old -> 2-char-padded new.
# In the .gro format, the resname column is 5 chars (left-justified) and the
# atomname column is 5 chars (right-justified). We replace the 3-char string
# in-place; padding stays consistent.
_ION_REWRITES = (
    ("NA+", "NA "),    # 3 chars -> 3 chars (preserves column width)
    ("CL-", "CL "),
)


def normalize_ion_atom_names_gro(gro_path: Path, out_path: Optional[Path] = None) -> Path:
    """Substitute ``NA+`` / ``CL-`` with ``NA`` / ``CL`` in a .gro file.

    Replacements happen only in the resname (cols 6-10) and atomname
    (cols 11-15) columns of atom lines. The first 2 lines (title, atom
    count) and the last (box) are passed through unchanged. Column widths
    are preserved by padding with a trailing space.
    """
    gro_path = Path(gro_path)
    if out_path is None:
        out_path = gro_path
    text = gro_path.read_text()
    lines = text.splitlines(keepends=True)
    if len(lines) < 3:
        Path(out_path).write_text(text)
        return Path(out_path)

    title = lines[0]
    natoms_line = lines[1]
    box_line = lines[-1]
    atom_lines = lines[2:-1]

    new_atom_lines: List[str] = []
    for line in atom_lines:
        if len(line) < 20:
            new_atom_lines.append(line)
            continue
        resname = line[5:10]
        atomname = line[10:15]
        for old, new in _ION_REWRITES:
            if old in resname:
                resname = resname.replace(old, new)
            if old in atomname:
                atomname = atomname.replace(old, new)
        new_line = line[:5] + resname + atomname + line[15:]
        new_atom_lines.append(new_line)

    out = title + natoms_line + "".join(new_atom_lines) + box_line
    Path(out_path).write_text(out)
    return Path(out_path)
