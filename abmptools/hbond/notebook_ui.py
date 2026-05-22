"""
notebook_ui.py
--------------
Jupyter ipywidgets-based UI for hbond analyzer.

Use::

    from abmptools.hbond import open_panel
    open_panel("/path/to/IMC.bdf")

Provides:
- functional group detection summary (textual + 2D RDKit structure of mol[0])
- criteria selector (luzar-chandler / strict / custom sliders)
- record range selector
- output prefix
- Run button → executes Analyzer, displays stats table + count plot
- path to colored.bdf with hint "open in gourmet to visualize"

Dependencies:
- ipywidgets >= 8
- rdkit (optional; falls back to text-only summary)
- matplotlib (for inline plot)
"""
from __future__ import annotations

import io
import os
import tempfile
from typing import Optional


def _build_rdkit_mol_from_topology(topo, positions):
    """Build an RDKit Mol from MoleculeTopology + 3D positions.

    Uses Atom_Name as the element symbol (works for IMC GAFF2 output:
    Cl, O, N, C, H).
    """
    try:
        from rdkit import Chem
        from rdkit.Geometry import Point3D
    except ImportError:
        return None

    rwmol = Chem.RWMol()
    for atom in topo.atoms:
        elem = atom.atom_name.strip()
        if elem and elem[0].isupper():
            sym = elem
            if len(sym) >= 2 and sym[1].isupper():
                sym = sym[0]
            elif len(sym) >= 2 and sym[1].isdigit():
                sym = sym[0]
        else:
            sym = "C"
        try:
            ra = Chem.Atom(sym)
        except RuntimeError:
            ra = Chem.Atom("C")
        rwmol.AddAtom(ra)

    bond_type_map = {
        1.0: Chem.BondType.SINGLE,
        1.5: Chem.BondType.AROMATIC,
        2.0: Chem.BondType.DOUBLE,
        3.0: Chem.BondType.TRIPLE,
    }
    for bond in topo.bonds:
        bt = bond_type_map.get(bond.order, Chem.BondType.SINGLE)
        try:
            rwmol.AddBond(bond.atom1, bond.atom2, bt)
        except RuntimeError:
            pass

    mol = rwmol.GetMol()
    conf = Chem.Conformer(mol.GetNumAtoms())
    for i in range(mol.GetNumAtoms()):
        x, y, z = positions[i]
        conf.SetAtomPosition(i, Point3D(float(x), float(y), float(z)))
    mol.AddConformer(conf, assignId=True)
    try:
        Chem.SanitizeMol(mol, sanitizeOps=Chem.SANITIZE_ALL ^ Chem.SANITIZE_PROPERTIES)
    except Exception:
        pass
    return mol


def _render_mol_svg(mol, highlight_atoms_by_color=None, size=(400, 300)) -> str:
    """Render an RDKit Mol to SVG with optional atom highlights.

    highlight_atoms_by_color: dict like
        {'red': [3, 4, 18, 40], 'blue': [2, 5, 14]}
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem, Draw
        from rdkit.Chem.Draw import rdMolDraw2D
    except ImportError:
        return ""

    heavy_mol = Chem.RemoveHs(mol) if mol.GetNumAtoms() > 0 else mol
    try:
        AllChem.Compute2DCoords(heavy_mol)
    except Exception:
        pass

    highlight_atoms = []
    atom_colors = {}
    color_map = {
        "red": (1.0, 0.3, 0.3),
        "blue": (0.3, 0.3, 1.0),
        "green": (0.3, 0.8, 0.3),
    }
    if highlight_atoms_by_color:
        heavy_idx_of_orig = {}
        orig_to_heavy = {}
        for new_idx, atom in enumerate(heavy_mol.GetAtoms()):
            orig = atom.GetPropsAsDict().get("origIdx")
            if orig is not None:
                orig_to_heavy[orig] = new_idx
        for color_name, atom_ids in highlight_atoms_by_color.items():
            rgb = color_map.get(color_name, (0.7, 0.7, 0.7))
            for a_idx in atom_ids:
                heavy_idx = orig_to_heavy.get(a_idx, a_idx)
                if 0 <= heavy_idx < heavy_mol.GetNumAtoms():
                    highlight_atoms.append(heavy_idx)
                    atom_colors[heavy_idx] = rgb

    drawer = rdMolDraw2D.MolDraw2DSVG(size[0], size[1])
    drawer.drawOptions().addAtomIndices = True
    drawer.DrawMolecule(
        heavy_mol,
        highlightAtoms=highlight_atoms,
        highlightAtomColors=atom_colors,
    )
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    # Strip the XML preamble so the inline <svg> embeds cleanly inside an
    # ipywidgets.HTML widget. JupyterLab sanitises documents that start with
    # `<?xml ...?>` and falls back to an "image" placeholder.
    if svg.startswith("<?xml"):
        svg = svg.split("?>", 1)[1].lstrip()
    return svg


def open_panel(bdf_path: str):
    """Open the hbond analyzer UI for the given BDF file.

    Must be called from a Jupyter notebook cell.
    """
    import ipywidgets as widgets
    from IPython.display import HTML, Image, clear_output, display

    from .analyzer import Analyzer, AnalyzerConfig
    from .functional_groups import (
        detect_amides, detect_carboxyls, summarize_groups
    )
    from .hbond_detector import HBondCriteria
    from .bdf_reader import BDFTrajectory

    # Pre-load topology
    traj = BDFTrajectory(bdf_path)
    traj.load_topology()
    summary = summarize_groups(traj.molecules)

    n_sec_amide = summary['n_amides'] - summary['n_amides_tert']
    info_html = widgets.HTML(value=f"""
    <h3>abmptools.hbond panel</h3>
    <p><b>Input:</b> <code>{bdf_path}</code></p>
    <ul>
      <li>n_molecules: <b>{len(traj.molecules)}</b></li>
      <li>n_records: <b>{traj.n_records}</b></li>
      <li>force field: <b>{summary['force_field']}</b></li>
      <li>detected carboxyls: <b>{summary['n_carboxyls']}</b>
          (in {summary['n_mols_with_carboxyl']} mols)</li>
      <li>detected amides: <b>{summary['n_amides']}</b>
          (tert: {summary['n_amides_tert']}, sec/prim: {n_sec_amide})</li>
      <li>detected hydroxyls (non-carboxyl): <b>{summary['n_hydroxyls']}</b></li>
      <li>detected N-H donors: <b>{summary['n_amine_donors']}</b></li>
    </ul>
    """)

    # Try to render mol[0] 2D structure with functional groups highlighted
    svg_view = widgets.HTML(value="<i>(RDKit not available — 2D structure skipped)</i>")
    try:
        frame0 = traj.get_frame(0)
        mol0 = _build_rdkit_mol_from_topology(
            traj.molecules[0], frame0.positions[0]
        )
        if mol0 is not None:
            carb0 = [g for g in detect_carboxyls(traj.molecules) if g.mol_index == 0]
            amid0 = [g for g in detect_amides(traj.molecules) if g.mol_index == 0]
            highlights = {}
            if carb0:
                cg = carb0[0]
                highlights["red"] = [cg.c_atom, cg.o_atom, cg.oh_atom, cg.ho_atom]
            if amid0:
                ag = amid0[0]
                highlights["blue"] = [ag.c_atom, ag.o_atom, ag.n_atom]
            svg = _render_mol_svg(mol0, highlight_atoms_by_color=highlights)
            if svg:
                svg_view.value = (
                    "<h4>mol[0] structure (carboxyl=red, amide=blue)</h4>"
                    + svg
                )
    except Exception as e:
        svg_view.value = f"<i>2D render failed: {e}</i>"

    # Config widgets
    criteria_dd = widgets.Dropdown(
        options=["luzar-chandler", "strict", "custom"],
        value="luzar-chandler",
        description="Criteria:",
    )
    d_da_slider = widgets.FloatSlider(
        value=3.5, min=2.5, max=4.5, step=0.1,
        description="d(D-A) max:", disabled=True,
    )
    angle_slider = widgets.FloatSlider(
        value=120.0, min=90.0, max=180.0, step=5.0,
        description="∠ min:", disabled=True,
    )
    d_ha_slider = widgets.FloatSlider(
        value=2.5, min=1.5, max=3.5, step=0.1,
        description="d(H-A) max:", disabled=True,
    )

    def _on_criteria_change(change):
        is_custom = change.new == "custom"
        d_da_slider.disabled = not is_custom
        angle_slider.disabled = not is_custom
        d_ha_slider.disabled = not is_custom
    criteria_dd.observe(_on_criteria_change, names="value")

    rec_start = widgets.IntText(value=0, description="rec start:")
    rec_end = widgets.IntText(value=-1, description="rec end:")
    mol_name = widgets.Text(value="IMC", description="mol prefix:")
    out_prefix = widgets.Text(
        value=os.path.join(tempfile.gettempdir(), "hbond_result"),
        description="out prefix:",
        layout=widgets.Layout(width="500px"),
    )

    # Classification mode (v1.28+): IMC (4-species COOH-centric) vs generic
    # (donor-type x acceptor-type pair stats, for PVA / peptide / alcohol etc.)
    classify_mode_dd = widgets.Dropdown(
        options=[
            ("imc (COOH 4-species: dual/chain/single/free)", "imc"),
            ("generic (donor-type x acceptor-type pair stats)", "generic"),
        ],
        value="imc",
        description="Mode:",
        layout=widgets.Layout(width="400px"),
    )

    # Functional group pair selection (Phase 3 enhancement)
    donor_cb = {
        "carboxyl":     widgets.Checkbox(value=True,  description="COOH (carboxyl OH)"),
        "amide_donor":  widgets.Checkbox(value=False, description="N-H (sec/prim amide)"),
        "amine_donor":  widgets.Checkbox(value=False, description="N-H (amine)"),
        "hydroxyl":     widgets.Checkbox(value=False, description="OH (alcohol)"),
    }
    acceptor_cb = {
        "carboxyl_O":   widgets.Checkbox(value=True,  description="C=O (carboxyl)"),
        "amide_O":      widgets.Checkbox(value=True,  description="C=O (amide)"),
        "hydroxyl_O":   widgets.Checkbox(value=False, description="O-H (as acceptor)"),
        "ether_O":      widgets.Checkbox(value=False, description="O (ether/ester)"),
    }
    donor_box = widgets.VBox(
        [widgets.HTML("<b>Donor groups (D-H)</b>")] + list(donor_cb.values()),
        layout=widgets.Layout(border="1px solid lightgray", padding="6px"),
    )
    acceptor_box = widgets.VBox(
        [widgets.HTML("<b>Acceptor groups (A)</b>")] + list(acceptor_cb.values()),
        layout=widgets.Layout(border="1px solid lightgray", padding="6px"),
    )

    # Lifetime options
    do_lifetime = widgets.Checkbox(value=True, description="compute lifetime")
    gap_tol = widgets.IntText(value=0, description="gap tol:")
    dt_input = widgets.FloatText(value=1.0, description="dt:")
    lifetime_box = widgets.HBox([do_lifetime, gap_tol, dt_input])

    run_button = widgets.Button(description="Run", button_style="primary")
    output = widgets.Output()

    def _on_run(_):
        with output:
            clear_output()
            custom = None
            if criteria_dd.value == "custom":
                custom = HBondCriteria(
                    d_da_max=d_da_slider.value,
                    d_ha_max=d_ha_slider.value,
                    angle_min=angle_slider.value,
                )
            selected_donors = [k for k, w in donor_cb.items() if w.value]
            selected_acceptors = [k for k, w in acceptor_cb.items() if w.value]
            if not selected_donors or not selected_acceptors:
                print("ERROR: select at least one donor and one acceptor group")
                return
            cfg = AnalyzerConfig(
                bdf_path=bdf_path,
                out_prefix=out_prefix.value,
                criteria_mode=criteria_dd.value,
                custom_criteria=custom,
                record_start=rec_start.value,
                record_end=rec_end.value,
                base_mol_name=mol_name.value,
                verbose=True,
                classify_mode=classify_mode_dd.value,
                donor_groups=selected_donors,
                acceptor_groups=selected_acceptors,
                compute_lifetime=do_lifetime.value,
                gap_tolerance=gap_tol.value,
                dt=dt_input.value,
            )
            analyzer = Analyzer(cfg)
            analyzer.load()
            analyzer.run()
            paths = analyzer.write_outputs()

            # stats table
            fr = analyzer.frame_results[-1]
            cls = fr.classification
            if classify_mode_dd.value == "generic":
                gcls = fr.generic_classification
                print("\n--- Last frame summary (generic, donor x acceptor) ---")
                if gcls is not None:
                    for (dt, at), st in sorted(gcls.pair_stats.items()):
                        denom_d = gcls.n_donor_candidates.get(dt, 0) or 1
                        denom_a = gcls.n_acceptor_candidates.get(at, 0) or 1
                        print(
                            "  {:>13s} -> {:<13s}  n={:4d}  "
                            "donors={}/{}  acceptors={}/{}".format(
                                dt, at, st.n_hbonds,
                                st.n_uniq_donors, denom_d,
                                st.n_uniq_acceptors, denom_a,
                            )
                        )
                print("\nOutputs:")
                for kind, path in paths.items():
                    print(f"  {kind}: {path}")
                return
            print("\n--- Last frame summary (per functional group) ---")
            if cls.n_carboxyls > 0:
                print(f"  COOH dual:   {cls.n_carboxyls_dual:4d} / "
                      f"{cls.n_carboxyls}  "
                      f"({100*cls.ratio_carboxyl_dual:.1f}%)")
                print(f"  COOH chain:  {cls.n_carboxyls_chain:4d} / "
                      f"{cls.n_carboxyls}  "
                      f"({100*cls.ratio_carboxyl_chain:.1f}%)")
                print(f"  COOH single: {cls.n_carboxyls_single:4d} / "
                      f"{cls.n_carboxyls}  "
                      f"({100*cls.ratio_carboxyl_single:.1f}%)")
                print(f"  COOH free:   {cls.n_carboxyls_free:4d} / "
                      f"{cls.n_carboxyls}  "
                      f"({100*cls.ratio_carboxyl_free:.1f}%)")
            if cls.n_amides > 0:
                print(f"  amide accept: {cls.n_amides_accept:4d} / "
                      f"{cls.n_amides}  "
                      f"({100*cls.ratio_amide_accept:.1f}%)")
                print(f"  amide free:   {cls.n_amides_free:4d} / "
                      f"{cls.n_amides}  "
                      f"({100*cls.ratio_amide_free:.1f}%)")
            print(f"  H-bonds: cc={fr.n_hbonds_cc}, ca={fr.n_hbonds_ca}")
            print("\n--- Mol-level representative role (for coloring) ---")
            print(f"  dual mols:   {fr.n_dual_mols:4d}")
            print(f"  chain mols:  {fr.n_chain_mols:4d}")
            print(f"  single mols: {fr.n_single_mols:4d}")
            print(f"  free mols:   {fr.n_free_mols:4d}")
            print("\nOutputs:")
            for kind, path in paths.items():
                print(f"  {kind}: {path}")
            if "uncolored" in paths:
                print(f"\n→ {paths['uncolored']!r}: Mol_Name preserved "
                      "(open in J-OCTA pre-render).")
            if "colored" in paths:
                print(f"→ {paths['colored']!r}: Mol_Name renamed + "
                      "Draw_Attributes set (open in gourmet / OCTA "
                      "post-render; in gourmet Python panel set "
                      "show.all('line','mol','molname',...)).")

            if "plot" in paths and os.path.exists(paths["plot"]):
                # Embed via IPython.display.Image so the PNG is base64-encoded
                # into the notebook. Passing a local filesystem path to
                # `<img src="...">` does not work in Jupyter because the
                # browser cannot fetch arbitrary local files from the
                # notebook server.
                display(Image(filename=paths["plot"], width=600))

    run_button.on_click(_on_run)

    panel = widgets.VBox([
        info_html,
        svg_view,
        classify_mode_dd,
        widgets.HBox([donor_box, acceptor_box]),
        widgets.HBox([criteria_dd]),
        widgets.HBox([d_da_slider, d_ha_slider, angle_slider]),
        widgets.HBox([rec_start, rec_end, mol_name]),
        lifetime_box,
        out_prefix,
        run_button,
        output,
    ])
    display(panel)
    return panel
