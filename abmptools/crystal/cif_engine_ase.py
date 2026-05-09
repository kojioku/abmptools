# -*- coding: utf-8 -*-
"""
abmptools.crystal.cif_engine_ase
---------------------------------
ASE-based CIF reader + PDB writer for the crystal-FMO pipeline.

Provides:

- :func:`read_cif_to_atoms`     -- ``ase.io.read`` thin wrapper
- :func:`expand_supercell`      -- ``Atoms.repeat`` with explicit shifts
- :func:`detect_molecules`      -- ``ase.neighborlist`` + connected
                                   components (Z' validation included)
- :func:`write_pdb_for_abmp`    -- ASE Atoms + molecule partition →
                                   abmptools-compatible HETATM PDB
                                   (Phase D-2)
- :func:`run_ase`               -- top-level orchestrator returning the
                                   supercell ``Atoms`` and molecule
                                   index lists
- :func:`run_ase_pipeline`      -- full Phase D-2 pipeline: CIF →
                                   layer<L>/{pdb,xyz}/<base>.{pdb,xyz}
                                   ready for :mod:`abmptools.pdb2fmo`

The ASE backend is opt-in via ``CIFEngineConfig.engine='ase'``. The
default ``legacy`` engine remains the byte-equivalent path covered by
the Phase B regression fixture; this module is intended for crystal
systems whose space group / layer count is not handled by the
self-rolled ``readcif`` parser.

**Layer convention.** ``readcif -l N`` walks ``[0, ±1, ±2, ..., ±(N//2)]``
shifts in each direction, yielding :math:`N^3` supercell copies for
odd N (``N=5`` -> 125 cells). :func:`expand_supercell` matches this by
calling ``Atoms.repeat((N, N, N))`` and shifting the origin so the
asymmetric unit sits at the centre. The atom count matches; the
absolute coordinates differ from the legacy parser by a constant
translation (the legacy layout is centred whereas
``Atoms.repeat`` starts at the cell origin). Downstream logic
(molecule detection, ``filter_molecules_by_distance``) is
translation-invariant, so this does not affect the science output.
"""
from __future__ import annotations

import os
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple


# ---------------------------------------------------------------------------
# CIF read
# ---------------------------------------------------------------------------

def read_cif_to_atoms(cif_path: str):
    """Read a CIF file and return an :class:`ase.Atoms` object.

    The returned ``Atoms`` already has symmetry-equivalent positions
    expanded (``ase.io.cif.read_cif`` parses
    ``_symmetry_equiv_pos_as_xyz``). Cell vectors and the spacegroup
    are stored in ``atoms.info``.

    Parameters
    ----------
    cif_path
        Path to the CIF file (small molecule CIF, single datablock).

    Raises
    ------
    ImportError
        If ASE is not installed (``pip install abmptools[crystal]``).
    """
    try:
        from ase.io import read
    except ImportError as exc:
        raise ImportError(
            "ASE is required for `--engine ase`; install via "
            "`mamba install ase` or `pip install abmptools[crystal]`."
        ) from exc

    # `fractional_occupancies=False` avoids a `scipy.spatial.distance.cdist`
    # TypeError that some ASE 3.28 / SciPy combinations hit on cif blocks
    # whose atom_site_occupancy column is entirely 1.0. Organic-crystal
    # FMO inputs never have partial occupancies, so disabling that branch
    # is safe.
    atoms = read(str(cif_path), format="cif", fractional_occupancies=False)
    return atoms


# ---------------------------------------------------------------------------
# Supercell expansion
# ---------------------------------------------------------------------------

def expand_supercell(atoms, layer: int):
    """Expand *atoms* into an :math:`N^3` supercell where ``N=layer``.

    Matches the legacy ``readcif -l N`` semantics for atom count: the
    output has ``layer**3`` copies of the input cell. Coordinates are
    not centred (``Atoms.repeat`` starts at the original cell origin),
    so absolute positions differ from the legacy parser by a constant
    translation; molecule detection and distance-based filtering are
    translation-invariant, so downstream stages are unaffected.

    Parameters
    ----------
    atoms
        ``ase.Atoms`` from :func:`read_cif_to_atoms`.
    layer
        Supercell layer count (positive integer; ``layer=5`` matches the
        csp7 demo).

    Returns
    -------
    ase.Atoms
        New ``Atoms`` with ``len(atoms) * layer**3`` atoms.
    """
    if layer < 1:
        raise ValueError(f"layer must be >= 1, got {layer}")
    return atoms.repeat((layer, layer, layer))


# ---------------------------------------------------------------------------
# Molecule detection
# ---------------------------------------------------------------------------

def _union_find_components(n: int, edges: Sequence[Tuple[int, int]]) -> List[List[int]]:
    """Return connected components as lists of indices.

    Plain Python union-find; sufficient for crystal supercells up to
    ~10^5 atoms (csp7 layer 5 reaches ~16k atoms, well within range).
    """
    parent = list(range(n))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for a, b in edges:
        union(int(a), int(b))

    groups: dict = {}
    for i in range(n):
        r = find(i)
        groups.setdefault(r, []).append(i)
    return list(groups.values())


def detect_molecules(
    atoms,
    atoms_in_mol: List[int],
    *,
    bond_tolerance: float = 0.4,
) -> List[List[int]]:
    """Split *atoms* into connected molecules via :mod:`ase.neighborlist`.

    Parameters
    ----------
    atoms
        Supercell ``ase.Atoms``.
    atoms_in_mol
        Expected atoms per molecule. Single value ``[N]`` means every
        molecule has ``N`` atoms (Z'=1). Two values ``[N1, N2]`` allow
        two molecule types in the asymmetric unit (Z'>=2 / heterodimer
        crystal).
    bond_tolerance
        Padding added to ``ase.neighborlist.natural_cutoffs`` (Å). The
        legacy default 0.4 Å is robust for typical organic small
        molecules.

    Returns
    -------
    List[List[int]]
        One list of atom indices per detected molecule. Total is
        ``len(atoms) // sum_per_molecule`` for homogeneous crystals.

    Raises
    ------
    ValueError
        If the detected molecule sizes do not match ``atoms_in_mol``
        within tolerance — typically indicates wrong ``atoms_in_mol``
        or a bond-tolerance mismatch.
    """
    try:
        from ase.neighborlist import natural_cutoffs, neighbor_list
    except ImportError as exc:
        raise ImportError(
            "ASE is required for molecule detection."
        ) from exc

    cutoffs = [c + bond_tolerance for c in natural_cutoffs(atoms, mult=1.0)]
    i_idx, j_idx = neighbor_list("ij", atoms, cutoffs)
    edges = list(zip(i_idx.tolist(), j_idx.tolist()))
    molecules = _union_find_components(len(atoms), edges)

    expected = set(atoms_in_mol)
    sizes_seen = sorted({len(m) for m in molecules})
    if not all(s in expected for s in sizes_seen):
        raise ValueError(
            f"Detected molecule sizes {sizes_seen} do not match "
            f"atoms_in_mol={atoms_in_mol}. Try adjusting "
            f"bond_tolerance (current {bond_tolerance} Å) or check the "
            f"asymmetric unit."
        )

    # Sort molecules by their first atom index for deterministic output.
    molecules.sort(key=lambda m: m[0])
    return molecules


def unwrap_molecules(atoms, molecules, *, bond_tolerance: float = 0.4) -> None:
    """Rewrite atom positions so each molecule sits in one piece.

    `Atoms.repeat` (the supercell expansion) carries the input cell's
    periodic boundary forward — molecules whose atoms straddle the
    asymmetric-unit boundary therefore appear split (atoms scattered to
    opposite faces of the supercell) in the raw expanded ``Atoms``.
    Both the legacy CIF parser and the historical ``readcif`` pipeline
    avoid this by always working from a centred / shifted origin; here
    we have to fix it explicitly.

    For each molecule, do a BFS along the bond graph (PBC-aware
    neighbour list) and translate every atom by the **minimum-image
    vector** to its parent. After this, the molecule is contiguous in
    real space (no boundary crossing) and ABINIT-MP can build a sane
    initial guess for the SCC procedure.

    The mutation is in-place via ``atoms.set_positions``. ``atoms.pbc``
    and ``atoms.cell`` are left untouched.

    Parameters
    ----------
    atoms
        Supercell ``ase.Atoms`` from :func:`expand_supercell`.
    molecules
        Connected-component partition from :func:`detect_molecules`.
        Each entry is a list of atom indices belonging to one molecule.
    bond_tolerance
        Same neighbour-list padding used by :func:`detect_molecules`.
    """
    try:
        from collections import deque
        from ase.geometry import find_mic
        from ase.neighborlist import natural_cutoffs, neighbor_list
    except ImportError as exc:
        raise ImportError("ASE is required for molecule unwrap.") from exc

    cutoffs = [c + bond_tolerance for c in natural_cutoffs(atoms, mult=1.0)]
    i_idx, j_idx = neighbor_list("ij", atoms, cutoffs)
    adj: dict = {}
    for a, b in zip(i_idx.tolist(), j_idx.tolist()):
        adj.setdefault(a, []).append(b)

    pos = atoms.get_positions().copy()
    cell = atoms.cell
    pbc = atoms.pbc

    for mol in molecules:
        mol_set = set(mol)
        seen = {mol[0]}
        queue = deque([mol[0]])
        while queue:
            cur = queue.popleft()
            for nb in adj.get(cur, []):
                if nb in seen or nb not in mol_set:
                    continue
                delta = pos[nb] - pos[cur]
                delta_mic, _ = find_mic(delta, cell, pbc)
                pos[nb] = pos[cur] + delta_mic
                seen.add(nb)
                queue.append(nb)

    atoms.set_positions(pos)


# ---------------------------------------------------------------------------
# Top-level orchestrator
# ---------------------------------------------------------------------------

def run_ase(
    cif: str,
    *,
    layer: int = 5,
    atoms_in_mol: Optional[List[int]] = None,
    bond_tolerance: float = 0.4,
):
    """Run the ASE-based CIF expansion pipeline.

    Returns
    -------
    Tuple[ase.Atoms, List[List[int]]]
        ``(supercell_atoms, molecules)`` where ``molecules`` is the
        connected-component partition produced by
        :func:`detect_molecules`.

    Notes
    -----
    PDB / XYZ / AJF serialisation is intentionally **not** included
    here; the orchestrator (:mod:`abmptools.crystal.builder`) drives
    that step so legacy and ASE engines share the same fragment-cut
    and AJF-emit logic.
    """
    if atoms_in_mol is None:
        atoms_in_mol = [32]

    atoms = read_cif_to_atoms(cif)
    super_atoms = expand_supercell(atoms, layer)
    molecules = detect_molecules(
        super_atoms, atoms_in_mol, bond_tolerance=bond_tolerance,
    )
    # `Atoms.repeat` carries the original cell's PBC forward, so a
    # molecule that straddles the input cell boundary lands as scattered
    # atoms in the supercell. Walk each molecule's bond graph and
    # rewrite positions to the minimum-image of the BFS parent, so
    # the molecule sits contiguously in real space.
    unwrap_molecules(super_atoms, molecules, bond_tolerance=bond_tolerance)
    return super_atoms, molecules


def _format_atom_name(element: str, index: int) -> str:
    """Format a 4-char PDB atom-name field for *element*+*index*.

    Single-char elements ('S', 'C', 'O', ...) get a leading space so the
    element symbol lands in column 13 (PDB v3.3 convention). Two-char
    elements ('Cl', 'Br') start at column 12. Indices are appended
    directly; the field truncates to 4 characters total to stay inside
    the fixed-width PDB layout.
    """
    name = f"{element}{index}"
    if len(name) > 4:
        name = name[:4]
    if len(element) == 1:
        return " " + f"{name:<3s}"
    return f"{name:<4s}"


def write_pdb_for_abmp(
    atoms,
    molecules: List[List[int]],
    *,
    molname: str = "UNK",
    output_pdb: str,
) -> str:
    """Emit an abmptools-compatible HETATM PDB from an ASE supercell.

    Layout matches the legacy ``abmptools.pdb_io.exportardpdbfull``
    output (HETATM records, blank chain id, residue sequence ==
    1-based molecule index, atom name == element + per-molecule
    serial). The resulting file feeds directly into
    :func:`abmptools.pdb2fmo.run_pdb2fmo` (with or without ``-xyz``).

    Coordinates are written verbatim (no PBC unwrap / centering).
    Callers needing the legacy "centred at origin" layout should
    translate ``atoms`` before calling this function.

    Parameters
    ----------
    atoms
        ``ase.Atoms`` (supercell). Atom positions are read from
        ``atoms.positions`` directly.
    molecules
        Per-molecule atom-index lists from :func:`detect_molecules`.
        Order determines the residue sequence numbers (1-based).
    molname
        Three-character residue name (default ``"UNK"`` matches csp7).
    output_pdb
        Output path. Existing file is overwritten.

    Returns
    -------
    str
        ``output_pdb`` (echoed for convenience in pipeline chaining).
    """
    if len(molname) > 3:
        raise ValueError(
            f"molname must be <= 3 chars (PDB residue name field), "
            f"got {molname!r}"
        )

    positions = atoms.get_positions()
    symbols = atoms.get_chemical_symbols()

    out_path = Path(output_pdb)
    lines = [
        f"COMPND    {out_path.name}",
        "AUTHOR    GENERATED BY abmptools.crystal.cif_engine_ase",
    ]
    serial = 0
    for res_seq_0, mol in enumerate(molecules):
        res_seq = res_seq_0 + 1
        elem_count: Dict[str, int] = {}
        for atom_idx in mol:
            elem = symbols[atom_idx]
            count = elem_count.get(elem, 0) + 1
            elem_count[elem] = count
            atom_name_4 = _format_atom_name(elem, count)
            x, y, z = positions[atom_idx]
            serial += 1
            lines.append(
                f"HETATM{serial:5d} {atom_name_4} {molname:<3s}  "
                f"{res_seq:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}"
                f"{1.00:6.2f}{0.00:6.2f}          {elem:>2s}  "
            )
    lines.append("END")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(lines) + "\n")
    return str(out_path)


def write_xyz_for_abmp(
    atoms,
    molecules: List[List[int]],
    *,
    output_xyz: str,
) -> str:
    """Emit a sibling XYZ file matching :func:`write_pdb_for_abmp`.

    Required when downstream :mod:`abmptools.pdb2fmo` runs in ``-xyz``
    mode: ``pdb_io.exportardxyzfull`` reads coordinates from a sibling
    ``<base>.xyz`` to populate the AJF ``&XYZ`` block at full precision.

    Format (matches the legacy ``readcif`` XYZ output):

        <atom count>
        <blank>
        <element> <x> <y> <z>
        ...

    Coordinates are written with full Python ``str(float)`` precision
    so the legacy ``-xyz`` path can extract the exact same numeric
    values that the bare ``HETATM`` (3-decimal) representation would
    truncate.
    """
    positions = atoms.get_positions()
    symbols = atoms.get_chemical_symbols()
    flat_indices = [idx for mol in molecules for idx in mol]

    lines = [str(len(flat_indices)), ""]
    for atom_idx in flat_indices:
        elem = symbols[atom_idx]
        x, y, z = positions[atom_idx]
        lines.append(f"{elem} {x} {y} {z}")
    out_path = Path(output_xyz)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(lines) + "\n")
    return str(out_path)


def run_ase_pipeline(
    cif: str,
    *,
    layer: int = 5,
    atoms_in_mol: Optional[List[int]] = None,
    bond_tolerance: float = 0.4,
    odir: str = "cifout",
    molname: str = "UNK",
    cwd: Optional[str] = None,
) -> Path:
    """Phase D-2 full pipeline: CIF → layer<L>/{pdb,xyz}/<base>.{pdb,xyz}.

    Mirrors :func:`abmptools.crystal.cif_engine_legacy.run_legacy` so
    :class:`CrystalOrchestrator` can drop in the ASE backend with the
    same interface. Output naming follows the legacy convention
    (``<cif_stem>layer<L>Zp1.{pdb,xyz}``) so downstream stages
    (pdb2fmo, postproc) keep working without changes.

    Parameters
    ----------
    cif
        Path to the CIF file (absolute or relative to *cwd*).
    layer
        Supercell layer count. ``layer=N`` -> ``N**3`` cells.
    atoms_in_mol
        Atoms per molecule (defaults to ``[32]`` matching csp7).
    bond_tolerance
        Padding added to ``natural_cutoffs`` (Å) for molecule detection.
    odir
        Output directory under cwd. Receives ``layer<L>/{pdb,xyz}/``.
    molname
        PDB residue name (3 chars).
    cwd
        Working directory (resolves *cif* and writes *odir*). Defaults
        to current process cwd.

    Returns
    -------
    Path
        ``cwd/odir/layer<L>``: contains ``pdb/`` and ``xyz/``
        subdirectories with the supercell outputs.
    """
    if atoms_in_mol is None:
        atoms_in_mol = [32]

    cif_path = Path(cif)
    base_dir = Path(cwd) if cwd is not None else Path(os.getcwd())
    if not cif_path.is_absolute():
        cif_path = (base_dir / cif).resolve()
    out_root = base_dir / odir / f"layer{layer}"
    pdb_dir = out_root / "pdb"
    xyz_dir = out_root / "xyz"

    super_atoms, molecules = run_ase(
        cif=str(cif_path),
        layer=layer,
        atoms_in_mol=atoms_in_mol,
        bond_tolerance=bond_tolerance,
    )

    base = cif_path.stem + f"layer{layer}Zp1"
    pdb_path = pdb_dir / f"{base}.pdb"
    xyz_path = xyz_dir / f"{base}.xyz"
    write_pdb_for_abmp(
        super_atoms, molecules,
        molname=molname,
        output_pdb=str(pdb_path),
    )
    write_xyz_for_abmp(
        super_atoms, molecules,
        output_xyz=str(xyz_path),
    )
    return out_root


__all__ = [
    "read_cif_to_atoms",
    "expand_supercell",
    "detect_molecules",
    "unwrap_molecules",
    "write_pdb_for_abmp",
    "write_xyz_for_abmp",
    "run_ase",
    "run_ase_pipeline",
]
