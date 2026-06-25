"""Per-molecule 2D structure diagrams marking the detected H-bond sites.

Builds an RDKit molecule from the COGNAC UDF topology (element from
``Atom_Name`` + connectivity + a single frame's 3D coordinates) and renders a
2D depiction with the **detected** donor / acceptor atoms highlighted
(red = donor, cyan = acceptor). One diagram is produced per unique molecular
species.

COGNAC UDFs (and OpenFF-interchange-exported ones in particular) usually carry
**no bond orders** (every ``order`` is 1.0), so aromaticity and C=O double
bonds are recovered from the 3D geometry via
``rdkit.Chem.rdDetermineBonds.DetermineBondOrders`` (assumes a neutral, closed
shell molecule by default — override with ``charge``).

RDKit is an optional dependency: if it is not importable, or the molecule
cannot be perceived (e.g. an ion with the wrong ``charge``), diagram generation
is skipped with a warning rather than aborting the analysis.
"""
from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple

# Two-letter elements first so "Cl"/"Br" are not truncated to "C"/"B".
_TWO_LETTER = ("Cl", "Br", "Si", "Na", "Mg", "Ca", "Fe", "Zn", "Li", "Se")
_ONE_LETTER = set("HCNOFPSKIBW")

DONOR_COLOR = (0.92, 0.30, 0.30)     # red
ACCEPTOR_COLOR = (0.20, 0.75, 0.90)  # cyan


def element_from_name(name: str) -> str:
    """Best-effort element symbol from a UDF ``Atom_Name``.

    Handles clean symbols (``"Cl"``, ``"O"``), element+index (``"C10"``,
    ``"HO"`` → ``H``), and is case-insensitive on the trailing letter so
    ``"CL"`` and ``"cl"`` both resolve to ``Cl``.
    """
    s = "".join(ch for ch in (name or "") if ch.isalpha())
    if not s:
        return ""
    cap = s[0].upper() + s[1:2].lower()
    for two in _TWO_LETTER:
        if cap.startswith(two):
            return two
    if s[0].upper() in _ONE_LETTER:
        return s[0].upper()
    return s[0].upper()


def build_mol_from_topology(topology, coords, charge: int = 0):
    """Build an RDKit ``Mol`` from a ``MoleculeTopology`` + 3D coordinates.

    Returns ``(mol, None)`` on success or ``(None, reason)`` on failure.
    ``coords`` is an ``(n_atoms, 3)`` array-like (Å). Atom order is preserved,
    so caller-supplied atom indices map 1:1 onto the returned molecule.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import rdDetermineBonds
        from rdkit.Geometry import Point3D
    except ImportError:
        return None, "rdkit not installed"

    atoms = topology.atoms
    if coords is None or len(coords) != len(atoms):
        return None, "coordinate/atom count mismatch"

    rw = Chem.RWMol()
    for a in atoms:
        sym = element_from_name(a.atom_name)
        if not sym:
            return None, f"could not resolve element for {a.atom_name!r}"
        rw.AddAtom(Chem.Atom(sym))
    for b in topology.bonds:
        rw.AddBond(int(b.atom1), int(b.atom2), Chem.BondType.SINGLE)
    mol = rw.GetMol()

    conf = Chem.Conformer(mol.GetNumAtoms())
    for i in range(mol.GetNumAtoms()):
        x, y, z = coords[i]
        conf.SetAtomPosition(i, Point3D(float(x), float(y), float(z)))
    mol.AddConformer(conf, assignId=True)

    try:
        # Perceive bond orders / aromaticity from geometry (bonds already set).
        rdDetermineBonds.DetermineBondOrders(mol, charge=int(charge))
        Chem.SanitizeMol(mol)
    except Exception as e:  # noqa: BLE001 - any perception failure → skip
        return None, f"bond-order perception failed: {type(e).__name__}: {e}"
    return mol, None


def draw_hbond_diagram(
    mol,
    donor_atoms: Sequence[int],
    acceptor_atoms: Sequence[int],
    png_path: str,
    svg_path: Optional[str] = None,
    legend: str = "",
    donor_note: str = "donor",
    acceptor_note: str = "acceptor",
) -> List[str]:
    """Render ``mol`` 2D with donor (red) / acceptor (cyan) atoms highlighted.

    ``donor_atoms`` / ``acceptor_atoms`` are atom indices into ``mol`` (the
    heavy donor atom, e.g. amide N, and the acceptor atom, e.g. carbonyl O).
    Non-polar H are removed for clarity; the highlighted heavy atoms are tracked
    across H removal via an atom property. Returns the list of written paths.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit.Chem.Draw import rdMolDraw2D

    work = Chem.Mol(mol)
    for i in donor_atoms:
        work.GetAtomWithIdx(int(i)).SetProp("_hb_role", "donor")
    for i in acceptor_atoms:
        work.GetAtomWithIdx(int(i)).SetProp("_hb_role", "acceptor")

    heavy = Chem.RemoveHs(work)
    AllChem.Compute2DCoords(heavy)  # 2D depiction, not a projection of the 3D pose

    hi_atoms: List[int] = []
    atom_cols: Dict[int, Tuple[float, float, float]] = {}
    for a in heavy.GetAtoms():
        if not a.HasProp("_hb_role"):
            continue
        role = a.GetProp("_hb_role")
        hi_atoms.append(a.GetIdx())
        atom_cols[a.GetIdx()] = DONOR_COLOR if role == "donor" else ACCEPTOR_COLOR
        a.SetProp("atomNote", donor_note if role == "donor" else acceptor_note)

    written: List[str] = []
    targets = [("png", png_path)] + ([("svg", svg_path)] if svg_path else [])
    for kind, path in targets:
        drawer = (rdMolDraw2D.MolDraw2DCairo if kind == "png"
                  else rdMolDraw2D.MolDraw2DSVG)(1100, 820)
        op = drawer.drawOptions()
        op.legendFontSize = 20
        op.bondLineWidth = 2
        rdMolDraw2D.PrepareAndDrawMolecule(
            drawer, heavy, legend=legend,
            highlightAtoms=hi_atoms, highlightAtomColors=atom_cols,
        )
        drawer.FinishDrawing()
        data = drawer.GetDrawingText()
        with open(path, "wb" if kind == "png" else "w") as f:
            f.write(data)
        written.append(path)
    return written
