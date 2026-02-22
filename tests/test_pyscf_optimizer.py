# -*- coding: utf-8 -*-
"""
Tests for abmptools.mlopt.pyscf_optimizer.

Parser/writer tests run without any optional dependencies.
Optimisation tests require pyscf + geometric (or berny) and are
automatically skipped when those are not installed.
"""
import math
import textwrap
from pathlib import Path

import pytest

from abmptools.mlopt.pyscf_optimizer import (
    QMOptimizerPySCF,
    _infer_element_from_atom_name,
    _parse_pdb,
    _parse_xyz,
    _write_xyz,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Small water molecule – slightly distorted so an optimiser can converge
WATER_XYZ = textwrap.dedent("""\
    3
    water
    O   0.000000   0.000000   0.119748
    H   0.000000   0.756950  -0.478993
    H   0.000000  -0.756950  -0.478993
""")

# Methane molecule
METHANE_XYZ = textwrap.dedent("""\
    5
    methane
    C   0.000000   0.000000   0.000000
    H   0.629118   0.629118   0.629118
    H  -0.629118  -0.629118   0.629118
    H  -0.629118   0.629118  -0.629118
    H   0.629118  -0.629118  -0.629118
""")

# Minimal PDB for water (element column populated)
WATER_PDB = textwrap.dedent("""\
    ATOM      1  O   HOH A   1       0.000   0.000   0.120  1.00  0.00           O
    ATOM      2  H1  HOH A   1       0.000   0.757  -0.479  1.00  0.00           H
    ATOM      3  H2  HOH A   1       0.000  -0.757  -0.479  1.00  0.00           H
    END
""")

# PDB with blank element column – element must be inferred from atom name
WATER_PDB_NO_ELEM = textwrap.dedent("""\
    ATOM      1  O   HOH A   1       0.000   0.000   0.120  1.00  0.00
    ATOM      2  H1  HOH A   1       0.000   0.757  -0.479  1.00  0.00
    ATOM      3  H2  HOH A   1       0.000  -0.757  -0.479  1.00  0.00
    END
""")

_pyscf_available = pytest.importorskip  # used below as a marker factory


def _skip_if_pyscf_missing():
    try:
        import pyscf  # noqa: F401
    except ImportError:
        pytest.skip("pyscf not installed")


def _skip_if_geomopt_missing(solver="geometric"):
    _skip_if_pyscf_missing()
    if solver == "geometric":
        try:
            import geometric  # noqa: F401
        except ImportError:
            pytest.skip("geometric not installed (pip install geometric)")
    else:
        try:
            import berny  # noqa: F401
        except ImportError:
            pytest.skip("pyberny not installed (pip install pyberny)")


# ---------------------------------------------------------------------------
# Parser tests (no optional deps)
# ---------------------------------------------------------------------------


def test_parse_xyz_water(tmp_path):
    p = tmp_path / "water.xyz"
    p.write_text(WATER_XYZ)
    atoms = _parse_xyz(p)
    assert len(atoms) == 3
    assert atoms[0][0] == "O"
    assert atoms[1][0] == "H"
    assert atoms[2][0] == "H"
    # Check approximate coordinates
    assert abs(atoms[0][3] - 0.119748) < 1e-5  # O z-coord


def test_parse_xyz_methane(tmp_path):
    p = tmp_path / "methane.xyz"
    p.write_text(METHANE_XYZ)
    atoms = _parse_xyz(p)
    assert len(atoms) == 5
    assert atoms[0][0] == "C"
    assert all(a[0] == "H" for a in atoms[1:])


def test_parse_xyz_bad_count(tmp_path):
    p = tmp_path / "bad.xyz"
    p.write_text("not_a_number\ncomment\nO 0 0 0\n")
    with pytest.raises(ValueError, match="atom count"):
        _parse_xyz(p)


def test_parse_xyz_too_short(tmp_path):
    p = tmp_path / "short.xyz"
    p.write_text("3\nonly one line follows\nO 0 0 0\n")
    with pytest.raises(ValueError):
        _parse_xyz(p)


def test_parse_pdb_water_with_element(tmp_path):
    p = tmp_path / "water.pdb"
    p.write_text(WATER_PDB)
    atoms = _parse_pdb(p)
    assert len(atoms) == 3
    assert atoms[0][0] == "O"
    assert atoms[1][0] == "H"
    assert atoms[2][0] == "H"
    assert abs(atoms[0][3] - 0.120) < 1e-3  # z-coord of O


def test_parse_pdb_water_inferred_element(tmp_path):
    p = tmp_path / "water_noelem.pdb"
    p.write_text(WATER_PDB_NO_ELEM)
    atoms = _parse_pdb(p)
    assert len(atoms) == 3
    assert atoms[0][0] == "O"
    assert atoms[1][0] == "H"


def test_parse_pdb_empty(tmp_path):
    p = tmp_path / "empty.pdb"
    p.write_text("REMARK nothing here\nEND\n")
    with pytest.raises(ValueError, match="No ATOM/HETATM"):
        _parse_pdb(p)


def test_infer_element_simple():
    assert _infer_element_from_atom_name(" O  ") == "O"
    assert _infer_element_from_atom_name(" H  ") == "H"
    assert _infer_element_from_atom_name(" N  ") == "N"
    assert _infer_element_from_atom_name(" C  ") == "C"
    # FE (iron) is a known 2-letter element
    assert _infer_element_from_atom_name("FE  ") == "Fe"
    # MG (magnesium) is a known 2-letter element
    assert _infer_element_from_atom_name("MG  ") == "Mg"
    # ZN (zinc) is a known 2-letter element
    assert _infer_element_from_atom_name("ZN  ") == "Zn"
    # CA: calcium is a known 2-letter element (Ca); ambiguous with alpha-C
    # but the lookup table takes precedence for known elements
    assert _infer_element_from_atom_name(" CA ") == "Ca"


def test_infer_element_leading_digit():
    # '1HB' -> strip '1' -> 'HB'; 'Hb' is NOT a known 2-letter element
    # -> falls back to first letter 'H'
    assert _infer_element_from_atom_name("1HB ") == "H"
    assert _infer_element_from_atom_name("2H  ") == "H"
    # '1CG' -> 'CG' -> 'Cg' not a known element -> 'C'
    assert _infer_element_from_atom_name("1CG ") == "C"


def test_infer_element_raises_on_blank():
    with pytest.raises(ValueError):
        _infer_element_from_atom_name("    ")
    with pytest.raises(ValueError):
        _infer_element_from_atom_name("123 ")


def test_write_xyz(tmp_path):
    atoms = [("O", 0.0, 0.0, 0.0), ("H", 0.0, 0.757, -0.479)]
    p = tmp_path / "out.xyz"
    _write_xyz(p, atoms, comment="test")
    lines = p.read_text().splitlines()
    assert lines[0].strip() == "2"
    assert lines[1].strip() == "test"
    assert lines[2].startswith("O   ")


# ---------------------------------------------------------------------------
# QMOptimizerPySCF constructor / validation tests (no optional deps)
# ---------------------------------------------------------------------------


def test_constructor_defaults():
    opt = QMOptimizerPySCF()
    assert opt.functional == "B3LYP"
    assert opt.basis == "def2-SVP"
    assert opt.dispersion == "d3bj"
    assert opt.charge == 0
    assert opt.spin == 0
    assert opt.max_steps == 100
    assert opt.solver == "geometric"


def test_constructor_bad_solver():
    with pytest.raises(ValueError, match="solver"):
        QMOptimizerPySCF(solver="badopt")


def test_constructor_bad_dispersion():
    with pytest.raises(ValueError, match="dispersion"):
        QMOptimizerPySCF(dispersion="d4")


def test_optimize_missing_file(tmp_path):
    _skip_if_pyscf_missing()
    opt = QMOptimizerPySCF()
    with pytest.raises(FileNotFoundError):
        opt.optimize(tmp_path / "no_such_file.xyz", tmp_path / "out.xyz")


def test_optimize_same_path(tmp_path):
    _skip_if_pyscf_missing()
    p = tmp_path / "mol.xyz"
    p.write_text(WATER_XYZ)
    opt = QMOptimizerPySCF()
    with pytest.raises(ValueError, match="different paths"):
        opt.optimize(p, p)


def test_optimize_unsupported_format(tmp_path):
    _skip_if_pyscf_missing()
    p = tmp_path / "mol.sdf"
    p.write_text("dummy")
    opt = QMOptimizerPySCF()
    with pytest.raises(ValueError, match="Unsupported input format"):
        opt.optimize(p, tmp_path / "out.xyz")


# ---------------------------------------------------------------------------
# Full optimisation tests (require pyscf + geometric)
# ---------------------------------------------------------------------------


def test_optimize_water_xyz(tmp_path):
    """Full geometry optimisation of water from an xyz file."""
    _skip_if_geomopt_missing("geometric")

    in_xyz = tmp_path / "water.xyz"
    out_xyz = tmp_path / "water_opt.xyz"
    in_xyz.write_text(WATER_XYZ)

    # Use a small basis for speed; dispersion off to avoid dep
    opt = QMOptimizerPySCF(
        functional="B3LYP",
        basis="sto-3g",
        dispersion="none",
        max_steps=50,
        verbose=0,
    )
    result = opt.optimize(str(in_xyz), str(out_xyz))

    assert out_xyz.exists(), "Output xyz not created"
    assert result["converged"] is True, "Optimisation did not converge"
    assert isinstance(result["energy"], float)
    assert not math.isnan(result["energy"])
    assert result["energy_hartree"] < 0  # DFT energy should be negative
    assert result["steps"] >= 0
    assert result["out_xyz"] == str(out_xyz)

    # Check output xyz is parseable and has correct atom count
    opt_atoms = _parse_xyz(out_xyz)
    assert len(opt_atoms) == 3
    assert opt_atoms[0][0] == "O"


def test_optimize_methane_xyz(tmp_path):
    """Full geometry optimisation of methane from an xyz file."""
    _skip_if_geomopt_missing("geometric")

    in_xyz = tmp_path / "methane.xyz"
    out_xyz = tmp_path / "methane_opt.xyz"
    in_xyz.write_text(METHANE_XYZ)

    opt = QMOptimizerPySCF(
        functional="B3LYP",
        basis="sto-3g",
        dispersion="none",
        max_steps=50,
        verbose=0,
    )
    result = opt.optimize(str(in_xyz), str(out_xyz))

    assert out_xyz.exists()
    assert isinstance(result["energy"], float)
    assert not math.isnan(result["energy"])
    opt_atoms = _parse_xyz(out_xyz)
    assert len(opt_atoms) == 5


def test_optimize_water_pdb(tmp_path):
    """Full geometry optimisation of water from a pdb file; output is xyz."""
    _skip_if_geomopt_missing("geometric")

    in_pdb = tmp_path / "water.pdb"
    out_xyz = tmp_path / "water_opt.xyz"
    in_pdb.write_text(WATER_PDB)

    opt = QMOptimizerPySCF(
        functional="B3LYP",
        basis="sto-3g",
        dispersion="none",
        max_steps=50,
        verbose=0,
    )
    result = opt.optimize(str(in_pdb), str(out_xyz))

    assert out_xyz.exists()
    assert result["converged"] is True
    opt_atoms = _parse_xyz(out_xyz)
    assert len(opt_atoms) == 3


def test_optimize_water_berny(tmp_path):
    """Full geometry optimisation using the berny solver."""
    _skip_if_geomopt_missing("berny")

    in_xyz = tmp_path / "water.xyz"
    out_xyz = tmp_path / "water_berny_opt.xyz"
    in_xyz.write_text(WATER_XYZ)

    opt = QMOptimizerPySCF(
        functional="B3LYP",
        basis="sto-3g",
        dispersion="none",
        solver="berny",
        max_steps=50,
        verbose=0,
    )
    result = opt.optimize(str(in_xyz), str(out_xyz))

    assert out_xyz.exists()
    assert isinstance(result["energy"], float)
    opt_atoms = _parse_xyz(out_xyz)
    assert len(opt_atoms) == 3


def test_result_dict_keys(tmp_path):
    """Verify all expected keys are present in the result dict."""
    _skip_if_geomopt_missing("geometric")

    in_xyz = tmp_path / "water.xyz"
    out_xyz = tmp_path / "water_opt.xyz"
    in_xyz.write_text(WATER_XYZ)

    opt = QMOptimizerPySCF(
        functional="B3LYP",
        basis="sto-3g",
        dispersion="none",
        max_steps=50,
        verbose=0,
    )
    result = opt.optimize(str(in_xyz), str(out_xyz))

    required_keys = {"energy", "energy_hartree", "steps", "converged", "out_xyz"}
    assert required_keys.issubset(result.keys()), (
        f"Missing keys: {required_keys - set(result.keys())}"
    )
    assert isinstance(result["converged"], bool)
    assert isinstance(result["steps"], int)
    assert isinstance(result["out_xyz"], str)
