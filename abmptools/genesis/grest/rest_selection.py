# -*- coding: utf-8 -*-
"""
abmptools.genesis.grest.rest_selection
--------------------------------------
Resolve a :class:`RESTSelectionSpec` to a flat list of residue
numbers (1-based, AMBER/GENESIS convention).

Two modes:

- ``mode="explicit"``: parse residue specs like ``"1-138"``,
  ``"21,96"``, ``"274-275"`` to ``List[int]``. Pure Python; no
  external tools required.
- ``mode="around"``: invoke ``cpptraj`` with a generated input
  script, parse ``resinfo`` output to get the residue list. Requires
  AmberTools ``cpptraj`` on the path; raises :class:`CommandError`
  otherwise.

Convention: GENESIS ``[SELECTION] groupN = rno:1-138`` uses 1-based
residue numbering identical to AMBER prmtop. We always emit and
expect 1-based numbers. Comparison vs cpptraj output (which is also
1-based) is direct.
"""
from __future__ import annotations

import logging
import re
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional

from ._subprocess import run_command, write_text
from .models import RESTSelectionSpec

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Result container
# ---------------------------------------------------------------------------

@dataclass
class RESTSelectionResult:
    """Resolved REST residues + metadata for ``validate`` preview.

    Attributes
    ----------
    residues
        Sorted, deduplicated 1-based residue numbers.
    n_residues
        ``len(residues)``.
    n_atoms
        Total atom count in those residues. Filled in by
        :func:`atom_count_from_prmtop` if a prmtop is supplied;
        otherwise 0.
    n_heavy_atoms
        Heavy-atom (non-H) count. Same convention as ``n_atoms``.
    selection_string
        GENESIS-style ``rno:`` selection (e.g. ``"rno:1-138"``).
    """
    residues: List[int]
    n_residues: int
    n_atoms: int = 0
    n_heavy_atoms: int = 0
    selection_string: str = ""


# ---------------------------------------------------------------------------
# Explicit-mode parser
# ---------------------------------------------------------------------------

_EXPLICIT_TOKEN_RE = re.compile(r"^\s*(\d+)\s*(?:-\s*(\d+))?\s*$")


def parse_explicit_residues(specs: List[str]) -> List[int]:
    """Parse a list of explicit residue specs to a flat 1-based int list.

    Accepts:

    - Single residue:   ``"21"`` -> ``[21]``
    - Hyphen range:      ``"1-138"`` -> ``[1, 2, ..., 138]``
    - Comma list:        ``"21,96,274"`` -> ``[21, 96, 274]``
    - Mixed:             ``"1-3,5,7-9"`` -> ``[1, 2, 3, 5, 7, 8, 9]``

    The output is sorted and deduplicated.

    Raises
    ------
    ValueError
        If a token is malformed or an ascending range has start > end.
    """
    out: List[int] = []
    for raw in specs:
        for token in raw.split(","):
            token = token.strip()
            if not token:
                continue
            m = _EXPLICIT_TOKEN_RE.match(token)
            if not m:
                raise ValueError(
                    f"Cannot parse residue token {token!r}; "
                    f"expected '<int>' or '<int>-<int>'."
                )
            start = int(m.group(1))
            end_s = m.group(2)
            end = int(end_s) if end_s is not None else start
            if start > end:
                raise ValueError(
                    f"Range start > end in token {token!r}."
                )
            out.extend(range(start, end + 1))
    return sorted(set(out))


# ---------------------------------------------------------------------------
# around-mode (cpptraj-driven)
# ---------------------------------------------------------------------------

def render_cpptraj_around_script(
    prmtop: Path,
    center: str,
    radius_A: float,
    out_txt: Path,
) -> str:
    """Render a cpptraj input script that prints ``resinfo`` for the
    residues whose centre-of-mass is within ``radius_A`` of *center*.

    The script writes one residue per line with format
    ``#Res Name First Last Natom #Orig #Mol C I`` (cpptraj default
    ``resinfo`` format) to *out_txt*.

    Parameters
    ----------
    prmtop
        AMBER prmtop file (used by cpptraj as the parm).
    center
        cpptraj-style mask, e.g. ``"rno:96"`` or ``":96"``. Both
        ``rno:`` (GENESIS-flavoured) and ``:`` (cpptraj-native) are
        normalised to the bare ``:N`` form here.
    radius_A
        Distance cutoff in Angstrom.
    out_txt
        Destination path for the resinfo dump.
    """
    # Normalise GENESIS-style "rno:96" -> cpptraj "@:96" form.
    if center.startswith("rno:"):
        cpptraj_center = ":" + center[len("rno:") :]
    elif center.startswith(":"):
        cpptraj_center = center
    else:
        cpptraj_center = ":" + center.lstrip()

    # cpptraj's "byres" (around) mask: ``(<center>)<:radius`` selects
    # all residues with at least one atom within radius_A of any atom
    # in <center>. Wrap with byres to get unique residues.
    around_mask = f"({cpptraj_center})<@{radius_A}"

    return (
        f"parm {prmtop}\n"
        f"resinfo {around_mask} & !:WAT & !:HOH & !:Na+ & !:Cl- "
        f"& !:NA & !:CL out {out_txt}\n"
        f"run\n"
        f"quit\n"
    )


_RESINFO_LINE_RE = re.compile(
    r"^\s*(\d+)\s+(\S+)\s+\d+\s+\d+\s+(\d+)"
)


def parse_cpptraj_resinfo(text: str) -> List[int]:
    """Parse a cpptraj ``resinfo`` output dump to a list of resnos.

    cpptraj's ``resinfo`` output looks like:

    ::

        #Res  Name First  Last Natom #Orig #Mol C I
            1 LYS      1    24    24     1     1
            5 ARG     68    91    24     5     1
            ...

    Lines starting with ``#`` are comments. We extract the first
    numeric column (``#Res``).
    """
    residues: List[int] = []
    for line in text.splitlines():
        if not line or line.lstrip().startswith("#"):
            continue
        m = _RESINFO_LINE_RE.match(line)
        if m:
            residues.append(int(m.group(1)))
    return sorted(set(residues))


def resolve_around(
    spec: RESTSelectionSpec,
    prmtop: Path,
    workdir: Path,
    cpptraj_path: str = "cpptraj",
) -> List[int]:
    """Run cpptraj to resolve an ``around``-mode REST selection.

    Parameters
    ----------
    spec
        Validated ``RESTSelectionSpec`` with ``mode="around"``.
    prmtop
        AMBER topology file.
    workdir
        Directory for cpptraj input/output staging (auto-created).
    cpptraj_path
        Executable name or path; default ``"cpptraj"``.

    Returns
    -------
    List[int]
        Sorted, deduplicated 1-based residue numbers.
    """
    if spec.mode != "around":
        raise ValueError(
            f"resolve_around requires mode='around', got {spec.mode!r}"
        )
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)

    inp_path = workdir / "rest_around.cpptraj"
    out_txt = workdir / "rest_residues.txt"
    script = render_cpptraj_around_script(
        prmtop=Path(prmtop),
        center=spec.center,
        radius_A=spec.radius_A,
        out_txt=out_txt,
    )
    write_text(inp_path, script)

    logger.info(
        "Resolving around-mode REST selection via cpptraj "
        "(center=%s, radius=%.2f A)",
        spec.center,
        spec.radius_A,
    )
    run_command(
        [cpptraj_path, "-i", str(inp_path)],
        cwd=workdir,
        capture=True,
    )
    return parse_cpptraj_resinfo(out_txt.read_text())


# ---------------------------------------------------------------------------
# Top-level entry
# ---------------------------------------------------------------------------

def resolve_rest_selection(
    spec: RESTSelectionSpec,
    prmtop: Optional[Path] = None,
    workdir: Optional[Path] = None,
    cpptraj_path: str = "cpptraj",
) -> RESTSelectionResult:
    """Resolve a :class:`RESTSelectionSpec` to a residue list and metadata.

    For ``mode="explicit"`` no external tools are required; ``prmtop``
    and ``workdir`` may be ``None``. For ``mode="around"`` both must
    be supplied.
    """
    if spec.mode == "explicit":
        residues = parse_explicit_residues(spec.residues)
    elif spec.mode == "around":
        if prmtop is None or workdir is None:
            raise ValueError(
                "resolve_rest_selection(mode='around') requires both "
                "prmtop and workdir."
            )
        residues = resolve_around(spec, prmtop, workdir, cpptraj_path)
    else:  # pragma: no cover -- __post_init__ rejects others
        raise ValueError(f"Unknown selection mode: {spec.mode}")

    return RESTSelectionResult(
        residues=residues,
        n_residues=len(residues),
        selection_string=format_genesis_selection(residues),
    )


def format_genesis_selection(residues: List[int]) -> str:
    """Render a residue list as a GENESIS ``rno:`` selection string.

    Contiguous runs are compacted into ranges. Examples:

    - ``[1, 2, 3]`` -> ``"rno:1-3"``
    - ``[21, 96, 274, 275]`` -> ``"rno:21,96,274-275"``
    - ``[]`` -> ``""``
    """
    if not residues:
        return ""
    sorted_r = sorted(set(residues))
    runs: List[str] = []
    start = sorted_r[0]
    prev = start
    for r in sorted_r[1:]:
        if r == prev + 1:
            prev = r
            continue
        runs.append(f"{start}-{prev}" if start != prev else str(start))
        start = r
        prev = r
    runs.append(f"{start}-{prev}" if start != prev else str(start))
    return "rno:" + ",".join(runs)
