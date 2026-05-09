"""Generate the numeric reference for csp7 R00001 layer3 HF/6-31G.

Phase D-3 originally extracted Total energy / monomer / IFIE with a
local Python regex parser. Phase D-3 (revised, 2026-05-09) replaces
that with a thin wrapper around `python -m abmptools.getifiepieda`
to avoid duplicating the in-tree post-processor: a single source of
truth, shared with production csp7 1500-structure runs.

Usage (from repo root, after a successful `--run-local` run):

    python tests/regression/reference/main/crystal_csp7/R00001/extract_layer3_hf_631g.py \\
        .smoke_artifacts/csp7_r00001_l3_hf_631g/out/XXXI-MMFF-R00001/cifout/layer3/pdb/for_abmp/XXXI-MMFF-R00001layer3Zp1-around_ar6.0.log \\
        --out-dir tests/regression/reference/main/crystal_csp7/R00001

Two CSV files land in --out-dir:

    expected_layer3_hf_631g_ifiesum.csv    (1 row, target-frag-1 sums + MonomerEnergy(1))
    expected_layer3_hf_631g_ifiedt.csv     (16 rows, dimer-es=False pairs around frag 1)

These are the canonical getifiepieda outputs renamed for repo storage.
"""
from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
from pathlib import Path


def split_log_for_getifiepieda(log: Path) -> tuple[str, str, int]:
    """Split a log filename like `XXXI-MMFF-R00001layer3Zp1-around_ar6.0.log`
    into (prefix_with_dir, suffix, structure_number) for getifiepieda's
    `-i '["prefix", "suffix"]'` form.
    """
    m = re.search(r"^(.*)(\d{5})(.*)$", log.name)
    if m is None:
        raise ValueError(
            f"log basename does not contain a 5-digit zero-padded structure id: {log.name!r}"
        )
    prefix = str(log.parent / m.group(1))
    structure_id = int(m.group(2))
    suffix = m.group(3)
    return prefix, suffix, structure_id


def run_getifiepieda(log: Path, work_dir: Path) -> Path:
    """Invoke getifiepieda with the SI-doc-canonical flags. Returns the
    `csv/` directory that getifiepieda emitted under work_dir.
    """
    prefix, suffix, sid = split_log_for_getifiepieda(log)
    arg_i = f'["{prefix}", "{suffix}"]'
    cmd = [
        sys.executable, "-m", "abmptools.getifiepieda",
        "--multi", "1",
        "-dimeres",
        "-imd",
        "-zp", "5",
        "-t", str(sid), str(sid), "1",
        "-nof90",
        "-i", arg_i,
    ]
    subprocess.run(cmd, check=True, cwd=str(work_dir))
    csv_dir = work_dir / "csv"
    if not csv_dir.is_dir():
        raise RuntimeError(f"getifiepieda did not produce {csv_dir}")
    return csv_dir


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("log", type=Path,
                    help="Full ABINIT-MP log (must have a 5-digit zero-padded id "
                         "in the basename, e.g. XXXI-MMFF-R00001layer3Zp1-...log).")
    ap.add_argument("--out-dir", type=Path, required=True,
                    help="Destination directory for the renamed reference CSVs.")
    args = ap.parse_args()

    # `.absolute()` instead of `.resolve()` so we don't dereference
    # symlinks that exist only to give getifiepieda its 5-digit
    # zero-padded structure id (used by the public-molecule samples
    # added in Phase D-5).
    log = args.log.absolute()
    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    # getifiepieda expects to write csv/ to the cwd. Run it from the
    # log's parent dir's parent so cwd has space for csv/ without
    # contaminating the for_abmp/ directory.
    work_dir = log.parent.parent
    csv_dir = run_getifiepieda(log, work_dir)

    sum_src = csv_dir / "frag1-dimer-es-false-ifiesum.csv"
    dt_src = csv_dir / "frag1-dimer-es-false-ifiedt.csv"
    if not sum_src.is_file() or not dt_src.is_file():
        raise RuntimeError(
            f"expected getifiepieda outputs not found in {csv_dir}: "
            f"sum={sum_src.is_file()}, dt={dt_src.is_file()}"
        )

    sum_dst = out_dir / "expected_layer3_hf_631g_ifiesum.csv"
    dt_dst = out_dir / "expected_layer3_hf_631g_ifiedt.csv"
    shutil.copy(sum_src, sum_dst)
    shutil.copy(dt_src, dt_dst)

    print(f"wrote {sum_dst} ({sum_dst.stat().st_size} bytes)")
    print(f"wrote {dt_dst} ({dt_dst.stat().st_size} bytes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
