# -*- coding: utf-8 -*-
"""
abmptools.gro2udf.make_whole
----------------------------
Put every molecule back together in every frame (``--make-whole``), from the
topology alone -- no gmx, no .tpr.

Why this exists next to ``-pbc nojump``
=======================================
``gmx trjconv -pbc nojump`` (what ``--trajectory`` runs by default) keeps each
atom's path continuous: an atom that moved more than half a box since the
previous frame is assumed to have been wrapped, and is shifted back. That
needs two things to be true, and a coarse-grained production run can break
both:

* **the first frame is whole.** nojump carries the reference over; it never
  repairs a molecule that starts out split. mdrun's own ``.gro`` / ``.xtc``
  can be split in every frame.
* **nothing moves half a box between frames.** With frames written every
  10 ns, small coarse-grained molecules can move that far, so a real move is
  mistaken for a wrap and the molecule is torn by one box length.

``gmx trjconv -pbc mol`` avoids both, but needs a ``.tpr`` -- and a ``.tpr``
can only be read by a GROMACS at least as new as the one that wrote it.

Here each frame is fixed on its own, walking the bonds of every molecule
(``[ bonds ]`` and ``[ constraints ]``) and moving each atom to the periodic
image nearest the atom it is bonded to. A molecule that crosses the box edge
moves to the other side between frames, but it is never shown split.

Only rectangular boxes are handled (``.gro`` frames here carry three box
lengths).
"""
from __future__ import annotations

from collections import defaultdict, deque
from typing import List, Optional, Sequence, Tuple

import numpy as np

from .top_model import GROFrameData

__all__ = ["unwrap_steps", "make_frames_whole"]


def unwrap_steps(raw) -> List[Tuple[np.ndarray, np.ndarray]]:
    """Return the walk as ``[(parents, children), ...]``, one entry per depth.

    Indices are 0-based positions in the full system (the order coordinates
    come in). Within one depth every child's parent is already placed, so a
    whole depth can be moved at once.
    """
    links = {}
    for name, atoms, bonds in zip(raw.mol_types, raw.atomlist, raw.bondlist):
        pairs = {(min(b[0], b[1]), max(b[0], b[1])) for b in bonds}
        pairs |= {(min(a, b), max(a, b))
                  for a, b in raw.constraint_pairs.get(name, [])}
        links[name] = (len(atoms), sorted(pairs))

    levels = defaultdict(lambda: ([], []))
    offset = 0
    for name in raw.mol_instance_list:
        n, pairs = links[name]
        adj = defaultdict(list)
        for i, j in pairs:                       # 1-based within the molecule
            adj[i - 1].append(j - 1)
            adj[j - 1].append(i - 1)
        depth = [-1] * n
        for start in range(n):
            if depth[start] >= 0:
                continue
            depth[start] = 0
            queue = deque([start])
            while queue:
                i = queue.popleft()
                for j in adj[i]:
                    if depth[j] >= 0:
                        continue
                    depth[j] = depth[i] + 1
                    levels[depth[j]][0].append(offset + i)
                    levels[depth[j]][1].append(offset + j)
                    queue.append(j)
        offset += n
    return [(np.asarray(p, dtype=int), np.asarray(c, dtype=int))
            for _, (p, c) in sorted(levels.items())]


def make_frames_whole(frames: Optional[Sequence[GROFrameData]],
                      steps: List[Tuple[np.ndarray, np.ndarray]]
                      ) -> Optional[List[GROFrameData]]:
    """Return *frames* with every molecule whole. ``None`` passes through."""
    if frames is None:
        return None
    out: List[GROFrameData] = []
    for f in frames:
        if getattr(f, "triclinic", False):
            raise ValueError(
                "--make-whole handles rectangular boxes only, and the frame at "
                "t = {} ps has a triclinic box (only its diagonal reaches "
                "here, so molecules would be put together with the wrong box). "
                "Make the molecules whole first with `gmx trjconv -pbc mol` "
                "and pass that trajectory without --make-whole.".format(f.time))
        x = np.asarray(f.coord_list, dtype=float)
        box = np.asarray(f.cell[:3], dtype=float)
        if np.any(box <= 0):
            raise ValueError("--make-whole needs a box; frame at t = {} ps has "
                             "{}".format(f.time, list(f.cell)))
        for parents, children in steps:
            d = x[children] - x[parents]
            x[children] -= np.round(d / box) * box
        out.append(GROFrameData(step=f.step, time=f.time,
                                coord_list=x.tolist(), cell=f.cell,
                                triclinic=f.triclinic))
    return out
