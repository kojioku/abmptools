"""Numeric regression tests for the csp7 R00001 layer3 HF/6-31G reference.

The fixture under tests/regression/reference/main/crystal_csp7/R00001/
contains the canonical getifiepieda outputs frozen on 2026-05-09:

    expected_layer3_hf_631g_ifiesum.csv     (1 row, target-frag-1 sums)
    expected_layer3_hf_631g_ifiedt.csv      (16 dimer-es=False pairs around frag 1)

Both CSVs are produced by `python -m abmptools.getifiepieda --multi 1
-dimeres -imd -zp 5 -t 1 1 1 -nof90 -i '[..,..]'` against the abinitmp
v2r8 log of the 9 h reference run. Re-running the extract script on
the same log must reproduce the same bytes.

Two test layers:

1. **Shape sanity** (always runs, ~0.01 s): committed CSVs have the
   expected columns / row counts / numeric ranges.

2. **Live numeric regression** (`@pytest.mark.slow`, requires
   abinitmp + ENABLE_FMO_LIVE_REGRESSION=1): run the full
   `abmp-crystal pipeline --run-local`, then run the extract script
   (= getifiepieda subprocess), and compare the produced CSVs to the
   committed ones byte-for-byte. Skipped by default; ~9 hours wall
   time on 1 core.

The earlier excerpt-log roundtrip path was removed because the new
extractor delegates to getifiepieda, which needs the *full* log to
parse all of `## MONOMER ENERGY`, `## HF-IFIE`, `## PIEDA` consistently.
The legacy `excerpt_layer3_hf_631g.log` and `expected_layer3_hf_631g.json`
are kept on disk for archival / cross-reference but are no longer
read by the test suite.
"""
from __future__ import annotations

import csv
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
REF_DIR = REPO_ROOT / "tests/regression/reference/main/crystal_csp7/R00001"
EXTRACT_SCRIPT = REF_DIR / "extract_layer3_hf_631g.py"
EXPECTED_SUM_CSV = REF_DIR / "expected_layer3_hf_631g_ifiesum.csv"
EXPECTED_DT_CSV = REF_DIR / "expected_layer3_hf_631g_ifiedt.csv"


def _read_csv(path: Path) -> list[list[str]]:
    with path.open() as fh:
        return list(csv.reader(fh))


def test_expected_csvs_have_committed_shape() -> None:
    """Sanity-check the committed reference CSVs."""
    sum_rows = _read_csv(EXPECTED_SUM_CSV)
    assert len(sum_rows) == 2  # header + 1 data row
    sum_header = sum_rows[0]
    for col in ("HF-IFIE", "ES", "EX", "CT-mix", "MonomerEnergy(1)"):
        assert col in sum_header, f"column {col} missing from ifiesum header"

    sum_data = sum_rows[1]
    sum_map = dict(zip(sum_header, sum_data))
    # MP2 columns must be 0 (HF log defensive padding in getmomenedf/getfiltifpifd).
    assert float(sum_map["MP2-IFIE"]) == 0.0
    assert float(sum_map["DI(MP2)"]) == 0.0
    # MonomerEnergy(1) is fragment 1's HF SCF energy in hartree (~-1435.43).
    monomer = float(sum_map["MonomerEnergy(1)"])
    assert -1436.0 < monomer < -1434.0
    # Sum of HF-IFIE around frag 1 (in kcal/mol). The reference value
    # is -23.74; allow slack so the test is robust to formatting.
    hf_sum = float(sum_map["HF-IFIE"])
    assert -30.0 < hf_sum < -10.0

    dt_rows = _read_csv(EXPECTED_DT_CSV)
    assert len(dt_rows) >= 2  # header + at least 1 row
    dt_header = dt_rows[0]
    for col in ("I", "J", "DIST", "DIMER-ES", "HF-IFIE", "ES", "EX", "CT-mix"):
        assert col in dt_header, f"column {col} missing from ifiedt header"
    # All rows are dimer-es=False (the `-dimeres` filter).
    es_idx = dt_header.index("DIMER-ES")
    for row in dt_rows[1:]:
        assert row[es_idx] == "F", f"unexpected DIMER-ES value: {row[es_idx]!r}"


@pytest.mark.slow
@pytest.mark.skipif(
    os.environ.get("ENABLE_FMO_LIVE_REGRESSION", "0") != "1",
    reason="set ENABLE_FMO_LIVE_REGRESSION=1 to run the ~9 h abinitmp regression",
)
def test_live_layer3_hf_631g_against_reference(tmp_path: Path) -> None:
    """Run the full pipeline + abinitmp + getifiepieda and compare to
    the committed CSVs byte-for-byte.

    Wall time: ~9 hours on 1 core (abinitmp v2r8). Run manually only.
    Inputs (cif + UNK.ajf) live in abmptools-sample at
    ``$ABMPTOOLS_SAMPLE_DIR/sample/csp7_ciftest/crystal_reference/R00001_layer3_hf_631g/``.
    """
    sample_dir = Path(os.environ.get(
        "ABMPTOOLS_SAMPLE_DIR",
        os.path.expanduser("~/repos/abmptools-sample"),
    ))
    inputs = sample_dir / "sample/csp7_ciftest/crystal_reference/R00001_layer3_hf_631g"
    if not inputs.is_dir():
        pytest.skip(
            f"crystal_reference inputs not found at {inputs}. "
            "Set ABMPTOOLS_SAMPLE_DIR to your abmptools-sample checkout."
        )

    work = tmp_path / "csp7_r00001_l3_hf_631g"
    work.mkdir()
    shutil.copy(inputs / "XXXI-MMFF-R00001.cif", work)
    shutil.copy(inputs / "UNK.ajf", work)

    abinit_dir = os.environ.get("ABINIT_DIR", "/home/okuwaki/bin/abinitmp/v2r8")
    binary_name = os.environ.get("ABINIT_BIN", "abinitmp")
    if not (Path(abinit_dir) / binary_name).is_file() and shutil.which(binary_name) is None:
        pytest.skip(f"abinitmp not found: {abinit_dir}/{binary_name}")

    yaml_text = f"""
project_name: csp7_r00001_l3_hf_631g
output_dir: ./out
inputs:
  - cif: XXXI-MMFF-R00001.cif
    layer: 3
    atoms_in_mol: [32]
cif_engine:
  engine: legacy
fragment:
  cutmode: around
  solutes: [0]
  criteria: 6.0
  molname: [UNK]
  template_ajf: ./UNK.ajf
fmo:
  method: HF
  basis_set: 6-31G
  memory: 4000
  abinit_ver: rev23
  npro: 1
  is_xyz: true
hpc:
  scheduler: local
  abinit_dir: '{abinit_dir}'
  binary_name: '{binary_name}'
  nodes: 1
  proc_per_node: 1
  omp_threads: 1
  elapse: '12:00:00'
postproc:
  enable: false
""".strip()
    (work / "crystal.yaml").write_text(yaml_text + "\n")

    subprocess.run(
        [sys.executable, "-m", "abmptools.crystal", "pipeline",
         "--config", "crystal.yaml", "--run-local"],
        cwd=work,
        check=True,
    )

    log = next((work / "out").rglob("*-around_ar6.0.log"))
    actual_dir = work / "actual"
    subprocess.run(
        [sys.executable, str(EXTRACT_SCRIPT),
         str(log), "--out-dir", str(actual_dir)],
        check=True,
    )

    actual_sum = actual_dir / "expected_layer3_hf_631g_ifiesum.csv"
    actual_dt = actual_dir / "expected_layer3_hf_631g_ifiedt.csv"
    assert actual_sum.read_bytes() == EXPECTED_SUM_CSV.read_bytes(), (
        "ifiesum.csv drift vs committed reference"
    )
    assert actual_dt.read_bytes() == EXPECTED_DT_CSV.read_bytes(), (
        "ifiedt.csv drift vs committed reference"
    )
