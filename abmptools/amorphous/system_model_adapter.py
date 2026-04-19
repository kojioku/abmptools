# -*- coding: utf-8 -*-
"""
abmptools.amorphous.system_model_adapter
------------------------------------------
Convert OpenFF artefacts (Interchange) into an
``abmptools.core.system_model.SystemModel`` so downstream pipelines can
treat "amorphous build" output uniformly with "UDF-derived" output.

Current scope (Phase 1)
~~~~~~~~~~~~~~~~~~~~~~~
Populates the SystemModel fields needed to:

- share a single intermediate representation with the GAFF + FMO-charge
  path used by fcewsmb (see ``fcewsmb/gaff_adapter.py``);
- tag the object with ``ensemble_family`` so :mod:`abmptools.udf2gro`
  writers can accept or reject it;
- drive the GROMACS ``.gro`` writer.

The full bond/angle/dihedral topology (``mol_topologies``) is **not**
populated here because OpenFF Interchange already knows how to emit a
``.top`` directly. Callers that need a GROMACS ``.top`` should keep
using ``interchange.to_top(path)``. A follow-up phase will fill
``mol_topologies`` so that :class:`abmptools.udf2gro.gromacs.writers.TopWriter`
can also consume this SystemModel.

Usage
~~~~~
::

    from abmptools.amorphous.system_model_adapter import from_interchange

    inter = create_interchange(molecules, counts, box_nm, mixture_pdb)
    model = from_interchange(inter, title="run1")
    # model.ensemble_family == "gromacs_ok"
    # model.atom_positions / model.cell / model.atom_types are populated.
"""
from __future__ import annotations

import os
import tempfile
from pathlib import Path
from typing import Any, List, Optional, Tuple

from ..core.system_model import (
    AtomPosition,
    AtomType,
    CellGeometry,
    SystemModel,
    classify_ensemble,
)


def from_interchange(
    interchange: Any,
    *,
    title: str = "amorphous",
    mdp_path: Optional[str] = None,
    scratch_dir: Optional[str] = None,
) -> SystemModel:
    """Build a :class:`SystemModel` from an OpenFF ``Interchange``.

    Internally the Interchange is serialised to a temporary ``.top`` /
    ``.gro`` pair, which is then parsed into a
    :class:`~abmptools.gro2udf.top_model.TopModel`. A subset of that
    TopModel is copied into a :class:`SystemModel` (atom types, global
    positions, cell, and molecule sequence); bond/angle/dihedral records
    are **not** copied in this phase.

    Parameters
    ----------
    interchange : openff.interchange.Interchange
        Parameterised system (typically produced by
        :func:`abmptools.amorphous.parameterizer.create_interchange`).
    title : str, optional
        Assigned to :attr:`SystemModel.title` (appears on the first line
        of a GroWriter output).
    mdp_path : str, optional
        Forwarded to the TopAdapter so ``ref_t`` / ``tau_t`` / Ewald cutoff
        can be derived. Leave ``None`` if no MDP is available.
    scratch_dir : str, optional
        Directory for the transient ``.top``/``.gro`` pair. Defaults to a
        system temp dir (auto-cleaned).

    Returns
    -------
    SystemModel
        ``ensemble_family`` is set via :func:`classify_ensemble` — the
        OpenFF path only produces GROMACS-native algorithms, so this
        almost always resolves to ``"gromacs_ok"``.
    """
    from ..gro2udf.top_parser import TopParser
    from ..gro2udf.top_adapter import TopAdapter

    if scratch_dir is None:
        ctx = tempfile.TemporaryDirectory(prefix="abmp_sysmodel_")
    else:
        os.makedirs(scratch_dir, exist_ok=True)
        ctx = _NullCtx(scratch_dir)

    with ctx as work:
        gro_path = str(Path(work) / "_amorphous.gro")
        top_path = str(Path(work) / "_amorphous.top")

        interchange.to_gro(gro_path)
        interchange.to_top(top_path)

        raw = TopParser().parse(top_path)
        top_model = TopAdapter().build(raw, gro_path, mdp_path=mdp_path)

    return _topmodel_to_systemmodel(top_model, title=title)


# ---------------------------------------------------------------------------
# Internal: TopModel → SystemModel projection (minimal subset)
# ---------------------------------------------------------------------------

def _topmodel_to_systemmodel(tm: Any, *, title: str) -> SystemModel:
    """Project the subset of :class:`TopModel` that has a 1:1 mapping to
    :class:`SystemModel`.

    Fields copied:

    - atom_types    (name, mass, sigma, epsilon)
    - atom_positions (per-atom xyz from the first GRO frame)
    - cell           (from the first GRO frame)
    - mol_sequence   (collapsed run-length of ``mol_instance_list``)
    - force-field defaults (comb_rule, fudgeLJ, fudgeQQ)
    - calcQQ = 1 (OpenFF Interchange always includes electrostatics)
    - ensemble_family = "gromacs_ok"

    Fields left intentionally empty:

    - mol_topologies (bond/angle/dihedral records — not yet required;
      keep emitting .top via ``interchange.to_top`` for now)
    - sim_params / ndx_data / cluster_data / fixed_labels
    """
    atom_types = [
        AtomType(
            name=a.name,
            mass=a.mass,
            sigma=a.sigma,
            epsilon=a.epsilon,
        )
        for a in tm.atom_type_specs
    ]

    positions: List[AtomPosition] = []
    cell: Optional[CellGeometry] = None
    if tm.frames:
        frame = tm.frames[0]
        cell = CellGeometry(a=frame.cell[0], b=frame.cell[1], c=frame.cell[2])

        # Build a flat stream of AtomPosition records by walking each
        # molecule instance in order and pulling coordinates from frame.
        spec_by_name = {s.name: s for s in tm.mol_specs}
        atom_global = 0
        mol_id_running = 0
        for inst_name in tm.mol_instance_list:
            mol_id_running += 1
            mol_id_capped = min(mol_id_running, 99999)
            spec = spec_by_name.get(inst_name)
            if spec is None:
                continue
            for a in spec.atoms:
                if atom_global >= len(frame.coord_list):
                    break
                x, y, z = frame.coord_list[atom_global]
                positions.append(
                    AtomPosition(
                        mol_id=mol_id_capped,
                        mol_name_short=inst_name[:5],
                        atom_gro_name=a.atom_name[:5],
                        atom_id=min(atom_global + 1, 99999),
                        x=x, y=y, z=z,
                        vx=0.0, vy=0.0, vz=0.0,
                    )
                )
                atom_global += 1

    # Run-length collapse of mol_instance_list → [(name, count), ...]
    mol_sequence: List[Tuple[str, int]] = []
    prev_name: Optional[str] = None
    run_count = 0
    for n in tm.mol_instance_list:
        if n == prev_name:
            run_count += 1
        else:
            if prev_name is not None:
                mol_sequence.append((prev_name, run_count))
            prev_name = n
            run_count = 1
    if prev_name is not None:
        mol_sequence.append((prev_name, run_count))

    model = SystemModel(
        title=title,
        udf_path="",
        comb_rule=tm.comb_rule,
        flags14=0,
        fudgeLJ=tm.fudge_lj,
        fudgeQQ=tm.fudge_qq,
        calcQQ=1,
        atom_types=atom_types,
        mol_topologies=[],
        mol_sequence=mol_sequence,
        atom_positions=positions,
        cell=cell,
        sim_params=None,
        ndx_data=None,
    )

    # OpenFF → GROMACS native algorithms; default "gromacs_ok".
    model.ensemble_family = classify_ensemble("")
    return model


class _NullCtx:
    """Context manager that yields a pre-existing directory unchanged."""

    def __init__(self, path: str):
        self._path = path

    def __enter__(self) -> str:
        return self._path

    def __exit__(self, *exc) -> None:
        return None
