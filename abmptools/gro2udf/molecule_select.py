# -*- coding: utf-8 -*-
"""
abmptools.gro2udf.molecule_select
---------------------------------
Keep only some molecule types when converting (``--keep-molecules``).

A solvated system is mostly solvent: in a coarse-grained run three quarters of
the particles can be water, and a UDF of every frame of every particle is slow
to open in a viewer. Dropping the solvent used to mean writing a second
``.top`` (and cutting the ``.gro`` / ``.xtc`` to match), which someone else
then has to reproduce by hand.

Here the selection is made on read instead. The ``.top``, ``.gro`` and
trajectory are the originals; the ``[ molecules ]`` order fixes which atoms
belong to which instance, so the same index list cuts every frame.
"""
from __future__ import annotations

import copy
from collections import Counter
from typing import Iterable, List, Optional, Sequence, Tuple

from .top_model import GROFrameData

__all__ = ["select_molecules", "subset_frames"]


def select_molecules(raw, keep: Iterable[str]) -> Tuple[object, List[int]]:
    """Return ``(raw', atom_indices)`` keeping only the molecule types in *keep*.

    ``raw'`` is a copy of *raw* whose ``mol_instance_list`` holds only the
    kept instances, in their original order. ``atom_indices`` are the 0-based
    positions of their atoms in the full system -- the order coordinates come
    in -- for :func:`subset_frames`.

    Raises ``ValueError`` when a name is not a molecule type of the system,
    so a typo stops the conversion instead of silently dropping a species.
    """
    keep = list(dict.fromkeys(keep))            # de-duplicate, keep order
    n_atoms = {name: len(atoms)
               for name, atoms in zip(raw.mol_types, raw.atomlist)}
    counts = Counter(raw.mol_instance_list)
    unknown = [k for k in keep if k not in counts]
    if unknown:
        listing = ", ".join("{} x{}".format(n, c) for n, c in counts.items())
        raise ValueError(
            "--keep-molecules: {} not in [ molecules ]. The system has: {}"
            .format(", ".join(unknown), listing))

    wanted = set(keep)
    indices: List[int] = []
    kept_instances: List[str] = []
    offset = 0
    for name in raw.mol_instance_list:
        n = n_atoms[name]
        if name in wanted:
            indices.extend(range(offset, offset + n))
            kept_instances.append(name)
        offset += n

    selected = copy.copy(raw)
    selected.mol_instance_list = kept_instances
    return selected, indices


def subset_frames(frames: Optional[Sequence[GROFrameData]],
                  indices: Sequence[int],
                  n_atoms_full: Optional[int] = None
                  ) -> Optional[List[GROFrameData]]:
    """Cut every frame down to *indices*. ``None`` passes through.

    *n_atoms_full*, when given, is checked against each frame first: a frame
    of a different system would otherwise be cut without complaint, and the
    UDF would hold the wrong atoms under the right names.
    """
    if frames is None:
        return None
    out: List[GROFrameData] = []
    for f in frames:
        if n_atoms_full is not None and len(f.coord_list) != n_atoms_full:
            raise ValueError(
                "--keep-molecules: a frame has {} atoms but the .top describes "
                "{}; the coordinates are not of this system"
                .format(len(f.coord_list), n_atoms_full))
        out.append(GROFrameData(step=f.step, time=f.time,
                                coord_list=[f.coord_list[i] for i in indices],
                                cell=f.cell))
    return out
