# -*- coding: utf-8 -*-
"""tests/cg_dpd/test_helpers.py — aij skeleton helper + monomer hand-craft builder."""
from __future__ import annotations

import pytest

from .conftest import requires_cognac


# --- create_empty_aij ------------------------------------------------------

def test_create_empty_aij_default():
    """3 segments → 6 pairs (3 同種 + 3 異種)."""
    from abmptools.cg.dpd import create_empty_aij
    aij = create_empty_aij(["A", "B", "W"])
    assert aij.mode == "a"
    assert aij.aii == 25.0
    assert len(aij.pairs) == 6  # 3 + 3
    same = [(i, j, v) for (i, j, v) in aij.pairs if i == j]
    diff = [(i, j, v) for (i, j, v) in aij.pairs if i != j]
    assert len(same) == 3
    assert all(v == 25.0 for (_, _, v) in same)
    assert len(diff) == 3
    assert all(v == 25.0 for (_, _, v) in diff)  # default off_diagonal=25


def test_create_empty_aij_chi_mode():
    """chi mode では同種は chi=0、 異種は off_diagonal."""
    from abmptools.cg.dpd import create_empty_aij
    aij = create_empty_aij(["A", "B"], mode="chi", off_diagonal=1.5)
    assert aij.mode == "chi"
    same = next((v for (i, j, v) in aij.pairs if i == j == "A"))
    diff = next((v for (i, j, v) in aij.pairs if (i, j) == ("A", "B")))
    assert same == 0.0
    assert diff == 1.5


def test_create_empty_aij_off_diagonal():
    """off_diagonal で異種 pair の値を変えられる."""
    from abmptools.cg.dpd import create_empty_aij
    aij = create_empty_aij(["A", "B"], aii=25.0, off_diagonal=35.0)
    same = [v for (i, j, v) in aij.pairs if i == j]
    diff = [v for (i, j, v) in aij.pairs if i != j]
    assert same == [25.0, 25.0]
    assert diff == [35.0]


def test_create_empty_aij_invalid_mode():
    from abmptools.cg.dpd import create_empty_aij
    with pytest.raises(ValueError, match="mode must be"):
        create_empty_aij(["A"], mode="invalid")


def test_create_empty_aij_writable(tmp_path):
    """生成した skeleton aij を write_aij で書き出して読み戻せる."""
    from abmptools.cg.dpd import create_empty_aij, write_aij, read_aij
    aij = create_empty_aij(["A", "B", "W"], aii=25.0, off_diagonal=30.0)
    p = tmp_path / "skel.dat"
    write_aij(p, aij)
    re = read_aij(p)
    assert re.mode == "a"
    assert len(re.pairs) == 6


# --- build_monomer ----------------------------------------------------------

def test_build_monomer_solvent_1_particle():
    """1-particle solvent (bond/angle なし)."""
    from abmptools.cg.dpd import build_monomer
    wat = build_monomer("wat", particle_names=["W"])
    assert wat.name == "wat"
    assert wat.n_particles == 1
    assert wat.particle_names == ["W"]
    assert wat.bond12 == []
    assert wat.angle13 == []


def test_build_monomer_dimer():
    """2-particle dimer (1 bond)."""
    from abmptools.cg.dpd import build_monomer
    dim = build_monomer("dimer", particle_names=["A", "B"], bond12=[(0, 1)])
    assert dim.n_particles == 2
    assert dim.bond12 == [(0, 1)]
    assert len(dim.bond12h) == 1
    # bond12h: [type, i, j, dist, stiff, 0]
    assert dim.bond12h[0] == [1, 0, 1, 0.86, 50.0, 0]


def test_build_monomer_trimer_with_angle():
    """3-particle linear chain (2 bond + 1 angle)."""
    from abmptools.cg.dpd import build_monomer
    tri = build_monomer(
        "trimer", particle_names=["A", "B", "C"],
        bond12=[(0, 1), (1, 2)], angle13=[(0, 1, 2)],
    )
    assert tri.n_particles == 3
    assert len(tri.bond12) == 2
    assert len(tri.angle13) == 1
    # angle13data: [a, b, c, eq, stiff]、 default eq=0 (180°), stiff=5
    assert tri.angle13data[0] == [0, 1, 2, 0.0, 5.0]


def test_build_monomer_custom_params():
    """bond12_dist / bond12_stiff / angle13_eq / angle13_stiff を上書き."""
    from abmptools.cg.dpd import build_monomer
    m = build_monomer(
        "x", particle_names=["A", "B", "C"],
        bond12=[(0, 1)], bond12_dist=1.0, bond12_stiff=100.0,
        angle13=[(0, 1, 2)], angle13_eq=60.0, angle13_stiff=10.0,
    )
    assert m.bond12h[0] == [1, 0, 1, 1.0, 100.0, 0]
    assert m.angle13data[0] == [0, 1, 2, 60.0, 10.0]


@requires_cognac
def test_build_monomer_with_dpd_builder(tmp_path):
    """hand-craft monomer を CGDpdBuilder に突っ込んで R1 UDF 生成まで通る."""
    from abmptools.cg.dpd import (
        build_monomer, create_empty_aij, write_aij,
        DpdSystemSpec, CGDpdBuilder, AijMatrix,
    )
    # 1-particle water monomer + 2-segment aij (W + ?)、 simple solvent system
    wat = build_monomer("wat", particle_names=["W"])
    aij = create_empty_aij(["W"], aii=25.0)
    aij_path = tmp_path / "wat_aij.dat"
    write_aij(aij_path, aij)

    spec = DpdSystemSpec(monomers=[wat], aij=aij, project_name="water-only")
    builder = CGDpdBuilder(spec=spec)
    udf = builder.build_udf(tmp_path / "wat_uin.udf")
    assert udf.exists()
    text = udf.read_text()
    assert '"W"' in text  # Atom_Type
    assert '"wat"' in text  # Set_of_Molecules
