# -*- coding: utf-8 -*-
"""tests/cg_dpd/conftest.py — 共有 fixture (cholesterol を CG segment 化した結果)."""
from __future__ import annotations

import subprocess
import sys
import tempfile
from pathlib import Path

import pytest


def _cognac_schema_available() -> bool:
    """UDFManager + cognac112.udf スキーマが解決でき named-path put が可能か。

    build-udf (UDFManager ベース writer) はスキーマ解決が必須なので、OCTA が
    無い / cognac112.udf が見つからない環境では writer 系テストを skip する。
    """
    try:
        from UDFManager import UDFManager  # noqa
    except Exception:
        return False
    d = tempfile.mkdtemp()
    p = Path(d) / "probe.udf"
    p.write_text('\\include{"cognac112.udf"}\n\\begin{data}\n\\end{data}\n')
    try:
        u = UDFManager(str(p))
        u.put("x", "Interactions.Pair_Interaction[0].Name")
        return True
    except Exception:
        return False


# 全 cg_dpd テストで共有する skip マーカー (UDFManager ベース writer 用)
requires_cognac = pytest.mark.skipif(
    not _cognac_schema_available(),
    reason="UDFManager + cognac112.udf schema unavailable (OCTA 非導入環境)",
)


@pytest.fixture(scope="session")
def repo_root() -> Path:
    """abmptools repo root path."""
    return Path(__file__).resolve().parents[2]


@pytest.fixture(scope="session")
def cholesterol_pdb(repo_root) -> Path:
    p = repo_root / "sample" / "cg_segmenter" / "cholesterol_rdkit.pdb"
    if not p.exists():
        pytest.skip(f"cholesterol PDB sample not found: {p}")
    return p


@pytest.fixture(scope="session")
def cholesterol_cg(tmp_path_factory, cholesterol_pdb) -> dict:
    """cg_segmenter で cholesterol を 5 segment + monomer/calc_sett 出力する。

    Returns
    -------
    dict with keys:
        ``td`` : tmp output dir
        ``monomer`` : path to chol_monomer
        ``calc_sett`` : path to chol_calc_sett
    """
    td = tmp_path_factory.mktemp("chol_cg")
    subprocess.run(
        [
            sys.executable, "-m", "abmptools.fragmenter.cg_segmenter", "dpdgen",
            "--pdb", str(cholesterol_pdb),
            "--output-dir", str(td),
            "--monomer-name", "chol", "--box", "12",
        ],
        check=True, capture_output=True,
    )
    return {
        "td": td,
        "monomer": td / "chol_monomer",
        "calc_sett": td / "chol_calc_sett",
    }


@pytest.fixture(scope="session")
def sample_aij_a_mode(tmp_path_factory):
    """5-segment 'a' mode aij.dat (cholesterol P0..P4 用 dummy)."""
    from abmptools.cg.dpd import AijMatrix, write_aij
    td = tmp_path_factory.mktemp("aij_a")
    aij = AijMatrix(
        segments=[f"P{i}" for i in range(5)],
        pairs=[
            ("P0","P0",25.0), ("P1","P1",25.0), ("P2","P2",25.0),
            ("P3","P3",25.0), ("P4","P4",25.0),
            ("P0","P1",28.0), ("P0","P2",30.0), ("P0","P3",32.0), ("P0","P4",40.0),
            ("P1","P2",28.0), ("P1","P3",30.0), ("P1","P4",40.0),
            ("P2","P3",28.0), ("P2","P4",40.0), ("P3","P4",40.0),
        ],
        mode="a", aii=25.0,
    )
    p = td / "aij_a.dat"
    write_aij(p, aij)
    return p


@pytest.fixture(scope="session")
def sample_aij_chi_mode(tmp_path_factory):
    """3-segment 'chi' mode aij.dat (chi → a 変換テスト用)."""
    from abmptools.cg.dpd import AijMatrix, write_aij
    td = tmp_path_factory.mktemp("aij_chi")
    aij = AijMatrix(
        segments=["A", "B", "W"],
        pairs=[("A","B", -4.108), ("A","W", 2.0), ("B","W", 3.5)],
        mode="chi", aii=25.0,
    )
    p = td / "aij_chi.dat"
    write_aij(p, aij)
    return p


@pytest.fixture(scope="session")
def dpm_template_path() -> Path:
    """user-provided 空 dpm template (man/octa/dpdfile-test/dpm-sample.dpm)。

    abmptools repo 外なので存在チェックして無ければ skip。
    """
    p = Path("/home/okuwaki/llm-project/fcews-workspace/man/octa/dpdfile-test/dpm-sample.dpm")
    if not p.exists():
        pytest.skip(f"dpm template not found (B 案 user-provided): {p}")
    return p


@pytest.fixture(scope="session")
def virtual_mom_template_path() -> Path:
    p = Path(
        "/home/okuwaki/llm-project/fcews-workspace/man/octa/dpdfile-test"
        "/dpm-sample/monomer-lib/segA/Virtual.mom"
    )
    if not p.exists():
        pytest.skip(f"Virtual.mom template not found: {p}")
    return p
