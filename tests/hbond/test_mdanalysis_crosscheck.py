"""
Cross-validation of the geometric H-bond detector against MDAnalysis.

The abmptools detector (`hbond_detector.detect_hbonds`) implements the
Luzar-Chandler criteria (``d(D...A) <= 3.5 Å`` and ``∠(D-H...A) >= 120°``)
with an orthogonal PBC minimum-image convention in plain numpy. This test
feeds the *same* donor/acceptor candidate set, the *same* criteria and the
*same* box to MDAnalysis' independent ``HydrogenBondAnalysis`` and asserts
that both implementations find exactly the same H-bonds.

Only the distance/angle *implementation* differs between the two paths, so a
match rules out a bug in the geometry. The tiny residual (< 0.01 Å, < 0.1°)
is the 3-decimal rounding incurred when writing coordinates to PDB for
MDAnalysis, not an algorithmic difference.

Skipped automatically when the IMC reference BDF or MDAnalysis is absent
(MDAnalysis is an optional, validation-only dependency).
"""
import os

import numpy as np
import pytest

IMC_BDF = "/home/okuwaki/llm-project/SI/IMC_result450.0_out_rec900.bdf"

pytestmark = pytest.mark.skipif(
    not os.path.exists(IMC_BDF), reason="IMC BDF not available"
)


def _element_from_type(atom_type):
    a = atom_type.lower()
    for pre, e in (("cl", "Cl"), ("br", "Br"), ("h", "H"), ("c", "C"),
                   ("o", "O"), ("n", "N"), ("s", "S"), ("p", "P"), ("f", "F")):
        if a.startswith(pre):
            return e
    return "C"


def _write_pdb(path, molecules, positions, cell):
    with open(path, "w") as f:
        f.write("CRYST1%9.3f%9.3f%9.3f  90.00  90.00  90.00 P 1           1\n"
                % (cell.a, cell.b, cell.c))
        serial = 0
        for m, mol in enumerate(molecules):
            for j, atom in enumerate(mol.atoms):
                el = _element_from_type(atom.atom_type)
                xyz = positions[m][j]
                serial += 1
                f.write(
                    "ATOM  %5d %-4s %3s A%4d    %8.3f%8.3f%8.3f"
                    "  1.00  0.00          %2s\n"
                    % (serial % 100000, el, "MOL", 1,
                       xyz[0], xyz[1], xyz[2], el)
                )
        f.write("END\n")


def test_detector_matches_mdanalysis(tmp_path):
    mda = pytest.importorskip(
        "MDAnalysis", reason="MDAnalysis (validation-only) not installed"
    )
    from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import (
        HydrogenBondAnalysis as HBA,
    )
    from abmptools.hbond import Analyzer, AnalyzerConfig
    from abmptools.hbond.hbond_detector import (
        DonorSite, AcceptorSite, detect_hbonds, HBondCriteria,
    )

    # --- (1) analyzer official detection (the COOH→COOH + COOH→amide set) ---
    cfg = AnalyzerConfig(
        bdf_path=IMC_BDF, out_prefix=str(tmp_path / "imc"),
        classify_mode="imc", do_colorize=False, do_copy_uncolored=False,
        do_plot=False, do_distance_plots=False, verbose=False,
    )
    an = Analyzer(cfg)
    an.load()
    fr = an.run()[0]
    n_official = fr.n_hbonds_cc + fr.n_hbonds_ca

    # --- candidate donor/acceptor sites rebuilt from the detected groups ---
    donors = [DonorSite(g.mol_index, g.oh_atom, g.ho_atom)
              for g in an.carboxyls]
    acceptors = ([AcceptorSite(g.mol_index, g.o_atom) for g in an.carboxyls]
                 + [AcceptorSite(g.mol_index, g.o_atom) for g in an.amides])

    frame = an.traj.get_frame(0)
    positions, cell = frame.positions, frame.cell

    # --- (2) abmptools detector on that candidate set ---
    mine = detect_hbonds(donors, acceptors, positions, cell,
                         HBondCriteria.luzar_chandler())
    # candidate reconstruction must reproduce the analyzer's own count
    assert len(mine) == n_official

    # global atom indices (MDAnalysis is 0-based, contiguous)
    natoms = [len(m.atoms) for m in an.traj.molecules]
    offset = np.concatenate([[0], np.cumsum(natoms)])[:-1].astype(int)
    mine_keys = {
        (int(offset[hb.donor_mol] + hb.donor_h),
         int(offset[hb.acceptor_mol] + hb.acceptor_a)): (hb.d_da, hb.angle)
        for hb in mine
    }

    # --- (3) MDAnalysis independent detection, same candidates/criteria ---
    pdb = str(tmp_path / "imc_rec0.pdb")
    _write_pdb(pdb, an.traj.molecules, positions, cell)
    u = mda.Universe(pdb)
    u.dimensions = [cell.a, cell.b, cell.c, 90.0, 90.0, 90.0]

    donor_heavy_g = sorted({int(offset[d.mol_index] + d.d_local)
                            for d in donors})
    hyd_g = sorted({int(offset[d.mol_index] + d.h_local) for d in donors})
    acc_g = sorted({int(offset[a.mol_index] + a.a_local) for a in acceptors})

    hba = HBA(
        universe=u,
        donors_sel="index " + " ".join(map(str, donor_heavy_g)),
        hydrogens_sel="index " + " ".join(map(str, hyd_g)),
        acceptors_sel="index " + " ".join(map(str, acc_g)),
        d_a_cutoff=3.5, d_h_a_angle_cutoff=120.0, d_h_cutoff=1.2,
        update_selections=False,
    )
    hba.run()
    mda_keys = {
        (int(r[2]), int(r[3])): (float(r[4]), float(r[5]))
        for r in hba.results.hbonds
    }

    # --- assertions -------------------------------------------------------
    assert len(mine_keys) > 50, "sanity: IMC should have many H-bonds"
    # identical count
    assert len(mine_keys) == len(mda_keys)
    # identical set of (hydrogen, acceptor) pairs — no false +/- either way
    assert set(mine_keys) == set(mda_keys)
    # numerical agreement; residual is PDB 3-decimal rounding only
    for k in mine_keys:
        assert abs(mine_keys[k][0] - mda_keys[k][0]) < 0.01   # d_DA, Å
        assert abs(mine_keys[k][1] - mda_keys[k][1]) < 0.1    # angle, deg
