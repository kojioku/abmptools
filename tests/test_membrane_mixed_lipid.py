# -*- coding: utf-8 -*-
"""
tests/test_membrane_mixed_lipid.py
----------------------------------
Unit tests for the per-lipid APL handling introduced in v1.17.1
(mixed-lipid support).

These tests cover the **pure-function** parts:
    - LipidSpec.apl_angstrom2 default 0.0 → auto-lookup
    - DEFAULT_LIPID_APL coverage of common species
    - _resolve_apl precedence: explicit > table > fallback
    - estimate_distxy_angstrom for pure / binary / ternary mixtures
    - assemble_packmol_memgen_cmd ratio reduction with multi-lipid

End-to-end packmol-memgen / tleap calls are exercised in the
integration smoke (`tests/integration/run_membrane_us_smoke.sh`)
and can be re-run there with a mixed-lipid config.
"""
from __future__ import annotations

from math import sqrt

import pytest

from abmptools.membrane import (
    LipidSpec, PeptideSpec, MembraneConfig, IonSpec,
    USProtocol, EquilibrationProtocol, PullingProtocol,
)
from abmptools.membrane.bilayer import (
    DEFAULT_LIPID_APL,
    _classify_lipid_head,
    _resolve_apl,
    assemble_packmol_memgen_cmd,
    estimate_distxy_angstrom,
    list_known_lipids,
)


def _basic_cfg(lipids):
    """Build a minimal valid MembraneConfig with the given lipid list."""
    return MembraneConfig(
        backend="amber",
        lipids=lipids,
        peptide=PeptideSpec(name="aa5", sequence="AAAAA"),
        output_dir="/tmp/membrane_mixed_lipid_test",
    )


# ---------------------------------------------------------------------------
# DEFAULT_LIPID_APL coverage
# ---------------------------------------------------------------------------

class TestDefaultLipidAPL:
    @pytest.mark.parametrize("resname,expected_min,expected_max", [
        ("POPC", 65.0, 70.0),    # canonical phospholipid, ~67
        ("DOPC", 70.0, 75.0),    # di-unsaturated, larger
        ("DPPC", 60.0, 66.0),    # di-saturated, smaller (Lα)
        ("POPE", 50.0, 60.0),    # PE head, smaller than PC
        ("CHL1", 35.0, 42.0),    # cholesterol, much smaller
        ("CHOL", 35.0, 42.0),    # alias
    ])
    def test_known_lipid_in_range(self, resname, expected_min, expected_max):
        assert resname in DEFAULT_LIPID_APL
        assert expected_min <= DEFAULT_LIPID_APL[resname] <= expected_max

    def test_table_has_chl_aliases(self):
        # CHL / CHL1 / CHOL all map to the same value (cholesterol)
        assert DEFAULT_LIPID_APL["CHL"] == DEFAULT_LIPID_APL["CHL1"]
        assert DEFAULT_LIPID_APL["CHOL"] == DEFAULT_LIPID_APL["CHL1"]


# ---------------------------------------------------------------------------
# _resolve_apl precedence
# ---------------------------------------------------------------------------

class TestResolveAPL:
    def test_explicit_overrides_table(self):
        # User-supplied APL takes precedence over DEFAULT_LIPID_APL
        l = LipidSpec(resname="POPC", n_per_leaflet=10, apl_angstrom2=80.0)
        assert _resolve_apl(l) == 80.0

    def test_table_lookup_when_default(self):
        # apl_angstrom2 default 0.0 → auto-lookup
        l = LipidSpec(resname="POPC", n_per_leaflet=10)
        assert _resolve_apl(l) == DEFAULT_LIPID_APL["POPC"]

    def test_unknown_residue_uses_fallback(self):
        l = LipidSpec(resname="UNKNOWN_LIPID", n_per_leaflet=10)
        assert _resolve_apl(l, fallback=72.0) == 72.0

    def test_fallback_default_is_65(self):
        l = LipidSpec(resname="UNKNOWN_LIPID", n_per_leaflet=10)
        assert _resolve_apl(l) == 65.0


# ---------------------------------------------------------------------------
# estimate_distxy_angstrom — pure / binary / ternary mixtures
# ---------------------------------------------------------------------------

class TestEstimateDistxy:
    def test_pure_popc(self):
        cfg = _basic_cfg([LipidSpec("POPC", n_per_leaflet=64)])
        # area = 64 * 67 = 4288, sqrt ~ 65.48
        assert abs(estimate_distxy_angstrom(cfg) - sqrt(64 * 67)) < 1e-6

    def test_binary_popc_chl1(self):
        # 4:1 mole ratio (POPC:CHL1) — typical "raft" model
        cfg = _basic_cfg([
            LipidSpec("POPC", n_per_leaflet=80),
            LipidSpec("CHL1", n_per_leaflet=20),
        ])
        expected_area = 80 * 67 + 20 * 38   # = 5360 + 760 = 6120
        assert abs(estimate_distxy_angstrom(cfg) - sqrt(expected_area)) < 1e-6

    def test_ternary_popc_pope_chl1(self):
        cfg = _basic_cfg([
            LipidSpec("POPC", n_per_leaflet=40),
            LipidSpec("POPE", n_per_leaflet=20),
            LipidSpec("CHL1", n_per_leaflet=10),
        ])
        expected_area = 40 * 67 + 20 * 56 + 10 * 38
        assert abs(estimate_distxy_angstrom(cfg) - sqrt(expected_area)) < 1e-6

    def test_explicit_apl_overrides_table(self):
        # Force a non-default APL on POPC (e.g. for low-T gel-phase studies)
        cfg = _basic_cfg([
            LipidSpec("POPC", n_per_leaflet=64, apl_angstrom2=50.0),
        ])
        assert abs(estimate_distxy_angstrom(cfg) - sqrt(64 * 50)) < 1e-6

    def test_unknown_lipid_uses_fallback(self):
        cfg = _basic_cfg([
            LipidSpec("MYNEWLIPID", n_per_leaflet=64),
        ])
        # default fallback = 65
        assert abs(estimate_distxy_angstrom(cfg) - sqrt(64 * 65)) < 1e-6
        # custom fallback = 70
        assert abs(estimate_distxy_angstrom(cfg, apl_angstrom2=70.0)
                   - sqrt(64 * 70)) < 1e-6

    def test_mixed_known_unknown_lipid(self):
        # Known POPC uses table (67); unknown falls back to apl_angstrom2 param
        cfg = _basic_cfg([
            LipidSpec("POPC", n_per_leaflet=50),
            LipidSpec("MYNEW", n_per_leaflet=10),
        ])
        # POPC: table → 67; MYNEW: fallback 65
        expected = 50 * 67 + 10 * 65
        assert abs(estimate_distxy_angstrom(cfg) - sqrt(expected)) < 1e-6


# ---------------------------------------------------------------------------
# assemble_packmol_memgen_cmd — multi-lipid ratio + distxy_fix
# ---------------------------------------------------------------------------

class TestAssembleCmdMultiLipid:
    def test_pure_popc_ratio_one(self):
        cfg = _basic_cfg([LipidSpec("POPC", n_per_leaflet=64)])
        cmd = assemble_packmol_memgen_cmd(
            config=cfg, peptide_pdb="p.pdb", output_pdb="out.pdb",
        )
        assert "--lipids" in cmd
        assert "POPC" in cmd[cmd.index("--lipids") + 1]
        assert "--ratio" in cmd
        assert cmd[cmd.index("--ratio") + 1] == "1"

    def test_binary_popc_chl1_4to1(self):
        cfg = _basic_cfg([
            LipidSpec("POPC", n_per_leaflet=80),
            LipidSpec("CHL1", n_per_leaflet=20),
        ])
        cmd = assemble_packmol_memgen_cmd(
            config=cfg, peptide_pdb="p.pdb", output_pdb="out.pdb",
        )
        # lipid string colon-separated
        assert cmd[cmd.index("--lipids") + 1] == "POPC:CHL1"
        # ratio gcd-reduced: 80:20 → 4:1
        assert cmd[cmd.index("--ratio") + 1] == "4:1"

    def test_ternary_2to1to1(self):
        cfg = _basic_cfg([
            LipidSpec("POPC", n_per_leaflet=40),
            LipidSpec("POPE", n_per_leaflet=20),
            LipidSpec("CHL1", n_per_leaflet=20),
        ])
        cmd = assemble_packmol_memgen_cmd(
            config=cfg, peptide_pdb="p.pdb", output_pdb="out.pdb",
        )
        assert cmd[cmd.index("--lipids") + 1] == "POPC:POPE:CHL1"
        assert cmd[cmd.index("--ratio") + 1] == "2:1:1"

    def test_distxy_fix_uses_per_lipid_apl(self):
        cfg = _basic_cfg([
            LipidSpec("POPC", n_per_leaflet=80),
            LipidSpec("CHL1", n_per_leaflet=20),
        ])
        cmd = assemble_packmol_memgen_cmd(
            config=cfg, peptide_pdb="p.pdb", output_pdb="out.pdb",
        )
        # find --distxy_fix and parse the value (Å, two decimals)
        idx = cmd.index("--distxy_fix")
        distxy = float(cmd[idx + 1])
        # expected: sqrt(80*67 + 20*38) ≈ sqrt(6120) ≈ 78.23 Å
        expected = sqrt(80 * 67 + 20 * 38)
        assert abs(distxy - expected) < 0.5  # within 0.5 Å of the formula

    def test_table_size_v1172(self):
        """v1.17.2 expanded the table to ~60 entries (Lipid21 PE/PG/PS/PA + SM)."""
        # Sanity: at least 50 entries (PC ~11, PE ~10, PG ~10, PS ~10, PA ~9, SM ~7, sterol ~3)
        assert len(DEFAULT_LIPID_APL) >= 50
        # All known head groups represented
        heads = {_classify_lipid_head(r) for r in DEFAULT_LIPID_APL}
        assert {"PC", "PE", "PG", "PS", "PA", "SM", "sterol"}.issubset(heads)


class TestClassifyLipidHead:
    @pytest.mark.parametrize("resname,expected", [
        ("POPC", "PC"), ("DOPC", "PC"), ("AHPC", "PC"),
        ("POPE", "PE"), ("DPPE", "PE"),
        ("POPG", "PG"), ("DAPG", "PG"),
        ("POPS", "PS"), ("DSPS", "PS"),
        ("POPA", "PA"), ("DPPA", "PA"),
        ("PSM", "SM"), ("OSM", "SM"), ("HSM", "SM"),
        ("CHL1", "sterol"), ("CHOL", "sterol"), ("CHL", "sterol"),
        ("ZZZZZ", "other"),
    ])
    def test_known_classification(self, resname, expected):
        assert _classify_lipid_head(resname) == expected


class TestListKnownLipids:
    def test_no_filter_returns_all(self):
        rows = list_known_lipids()
        assert len(rows) == len(DEFAULT_LIPID_APL)

    def test_pc_filter(self):
        rows = list_known_lipids(head_group="PC")
        assert all(head == "PC" for _name, _apl, head in rows)
        # Includes well-known PC species
        names = {r[0] for r in rows}
        assert {"POPC", "DOPC", "DPPC"}.issubset(names)

    def test_sm_filter(self):
        rows = list_known_lipids(head_group="SM")
        # 7 SM variants in the table (LSM..HSM)
        assert len(rows) == 7
        # Sorted alphabetically
        names = [r[0] for r in rows]
        assert names == sorted(names)

    def test_sterol_filter(self):
        rows = list_known_lipids(head_group="sterol")
        names = {r[0] for r in rows}
        assert {"CHL", "CHL1", "CHOL"}.issubset(names)

    def test_apl_values_present(self):
        rows = list_known_lipids()
        for resname, apl, _head in rows:
            assert apl > 0
            assert apl == DEFAULT_LIPID_APL[resname]


    def test_charmm_backend_adds_charmm_flag_with_mix(self):
        cfg = MembraneConfig(
            backend="charmm36",
            charmm_ff_dir="/dev/null/dummy",
            lipids=[
                LipidSpec("POPC", n_per_leaflet=80),
                LipidSpec("POPE", n_per_leaflet=20),
            ],
            peptide=PeptideSpec(name="aa5", sequence="AAAAA"),
            output_dir="/tmp/x",
        )
        cmd = assemble_packmol_memgen_cmd(
            config=cfg, peptide_pdb="p.pdb", output_pdb="out.pdb",
        )
        assert "--charmm" in cmd
        assert cmd[cmd.index("--lipids") + 1] == "POPC:POPE"
        assert cmd[cmd.index("--ratio") + 1] == "4:1"
