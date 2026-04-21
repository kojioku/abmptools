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

from typing import List, Optional

from ..gro2udf.top_model import GROFrameData


__all__ = ["frames_from_xtc"]


def frames_from_xtc(
    tpr_path: str,
    xtc_path: str,
    start: int = 0,
    end: Optional[int] = None,
) -> List[GROFrameData]:
    """Read a GROMACS trajectory into :class:`GROFrameData` records.

    Parameters
    ----------
    tpr_path : str
        Path to the binary ``*.tpr`` (needed by MDAnalysis to get the atom
        ordering and topology consistent with the ``.gro`` used to build
        the TopModel).
    xtc_path : str
        Path to the compressed trajectory (``*.xtc`` or ``*.trr``).
    start, end : int, optional
        Inclusive frame range. ``start`` defaults to 0 (first frame);
        ``end`` defaults to the last frame.

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

    universe = mda.Universe(tpr_path, xtc_path)
    total = len(universe.trajectory)
    last = total - 1 if end is None else min(end, total - 1)
    first = max(0, start)

    frames: List[GROFrameData] = []
    for i, ts in enumerate(universe.trajectory):
        if i < first or i > last:
            continue
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
