# -*- coding: utf-8 -*-
"""
abmptools.amorphous.trajectory_ingest
-------------------------------------
Helpers to read GROMACS trajectory files (``xtc`` / ``trr``) into the same
``GROFrameData`` records that :class:`abmptools.gro2udf.top_exporter.TopExporter`
consumes when writing COGNAC UDF structure records.

Used by downstream tools that want to reuse a static
``TopModel`` (built from ``system.top`` + initial ``system.gro``) and swap
the coordinate frames in for a longer MD trajectory — e.g. fcews-manybody's
``setupgetcontact_gromacs`` injects xtc-sourced frames so the resulting
UDF can be fed into legacy contact / AJF extraction unchanged.

MDAnalysis is imported lazily so this module can be loaded in environments
that don't have the package installed.
"""
from __future__ import annotations

import logging
from typing import List, Optional, Tuple

from ..gro2udf.top_model import GROFrameData

logger = logging.getLogger(__name__)


__all__ = ["frames_from_xtc"]


def _bonds_from_top(top_path: str) -> List[Tuple[int, int]]:
    """Parse a GROMACS .top into a system-wide 0-based bond list.

    For each entry in ``[ molecules ]`` we replicate the per-mol-type
    ``[ bonds ]`` with an offset equal to the running atom count.
    Returns a list of ``(atom_i, atom_j)`` tuples in 0-based indexing,
    suitable for :meth:`MDAnalysis.Universe.add_bonds`.
    """
    from ..gro2udf.top_parser import TopParser

    raw = TopParser().parse(top_path)
    type_to_idx = {name: i for i, name in enumerate(raw.mol_types)}

    bonds: List[Tuple[int, int]] = []
    seen = set()
    offset = 0
    for instance_name in raw.mol_instance_list:
        if instance_name not in type_to_idx:
            # Unknown molecule type — skip but advance offset by 0
            # (we won't know the atom count, so we leave it; downstream
            #  make_whole will simply not unwrap this segment).
            continue
        i = type_to_idx[instance_name]
        n_atoms = len(raw.atomlist[i])
        for b in raw.bondlist[i]:
            # bondlist entry: [atom1, atom2, type_index] (1-based)
            a1, a2 = int(b[0]) - 1, int(b[1]) - 1
            if a1 == a2:
                continue
            edge = (offset + min(a1, a2), offset + max(a1, a2))
            if edge in seen:
                continue
            seen.add(edge)
            bonds.append(edge)
        offset += n_atoms

    return bonds


def frames_from_xtc(
    topology_path: str,
    xtc_path: str,
    start: int = 0,
    end: Optional[int] = None,
    top_path: Optional[str] = None,
) -> List[GROFrameData]:
    """Read a GROMACS trajectory into :class:`GROFrameData` records.

    Parameters
    ----------
    topology_path : str
        Topology file MDAnalysis builds the ``Universe`` from. Most callers
        pass ``build/system.gro`` because GROMACS-2026 ``.tpr`` (tpx 138+)
        is not yet supported by MDAnalysis' TPRParser. The atom ordering
        must match ``xtc_path``.
    xtc_path : str
        Path to the compressed trajectory (``*.xtc`` or ``*.trr``).
    start, end : int, optional
        Inclusive frame range. ``start`` defaults to 0 (first frame);
        ``end`` defaults to the last frame.
    top_path : str, optional
        Path to the GROMACS ``.top``. When provided we parse the
        ``[ bonds ]`` section to populate Universe bonds and then apply
        :func:`MDAnalysis.lib.mdamath.make_whole` per fragment per frame
        — i.e. molecules that straddle a periodic-boundary in ``xtc``
        get reconstructed as single connected clusters. Without this
        the FMO snapshot extraction downstream sees half-broken
        molecules and SCC fails (see the methanol-acetone reference
        run for the diagnostic case).

    Returns
    -------
    list of GROFrameData
        One record per selected frame, with coordinates converted to nm
        (MDAnalysis returns Å) and cell as the diagonal of the unit-cell
        matrix (orthogonal box — non-orthogonal boxes are flattened to
        their ``a/b/c`` lengths, matching the GRO frame convention).

    Raises
    ------
    RuntimeError
        When MDAnalysis is not installed; the error message points to
        ``pip install MDAnalysis``.
    """
    try:
        import MDAnalysis as mda  # type: ignore
    except ImportError as e:
        raise RuntimeError(
            "MDAnalysis is required by "
            "abmptools.amorphous.trajectory_ingest.frames_from_xtc "
            "(`pip install MDAnalysis`)."
        ) from e

    universe = mda.Universe(topology_path, xtc_path)

    apply_make_whole = False
    make_whole = None
    if top_path is not None:
        try:
            from MDAnalysis.lib.mdamath import make_whole as _make_whole
            bonds = _bonds_from_top(top_path)
            if bonds:
                universe.add_bonds(bonds)
                apply_make_whole = True
                make_whole = _make_whole
                logger.info(
                    "trajectory_ingest: added %d bonds from %s, "
                    "make_whole will run per fragment per frame",
                    len(bonds), top_path,
                )
        except Exception as e:  # noqa: BLE001
            # Non-fatal: log and continue without make_whole. Caller can
            # still pre-process xtc with `gmx trjconv -pbc mol`.
            logger.warning(
                "trajectory_ingest: failed to load bonds from %s (%s); "
                "skipping make_whole. Pre-process xtc with "
                "'gmx trjconv -pbc mol -ur compact' to ensure whole "
                "molecules.",
                top_path, e,
            )

    total = len(universe.trajectory)
    last = total - 1 if end is None else min(end, total - 1)
    first = max(0, start)

    frames: List[GROFrameData] = []
    for i, ts in enumerate(universe.trajectory):
        if i < first or i > last:
            continue
        if apply_make_whole:
            for fragment in universe.atoms.fragments:
                make_whole(fragment)
        # MDAnalysis positions are Å — convert to nm (GROMACS convention)
        coords_nm = (universe.atoms.positions / 10.0).tolist()
        # universe.dimensions: [a, b, c, alpha, beta, gamma] in Å / deg
        dims = universe.dimensions
        cell_nm = [dims[0] / 10.0, dims[1] / 10.0, dims[2] / 10.0]
        frames.append(
            GROFrameData(
                step=i,
                time=float(getattr(ts, "time", 0.0)),
                coord_list=coords_nm,
                cell=cell_nm,
            )
        )
    return frames
