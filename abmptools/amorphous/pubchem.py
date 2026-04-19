# -*- coding: utf-8 -*-
"""
abmptools.amorphous.pubchem
----------------------------
Thin wrapper around the PubChem PUG REST API to fetch 3D SDF conformers
or canonical SMILES without requiring an extra Python dependency.

Only the standard library (`urllib`) is used, so this module is cheap to
import and safe to use in environments where `requests` / `pubchempy`
may not be installed.

Usage (Python)::

    from abmptools.amorphous.pubchem import fetch_3d_sdf, fetch_smiles

    sdf_text = fetch_3d_sdf(3825)                      # CID
    sdf_text = fetch_3d_sdf("ketoprofen", by="name")   # compound name
    smi      = fetch_smiles(3825)

Usage (CLI)::

    python -m abmptools.amorphous.pubchem --cid 3825 --output ketoprofen.sdf
    python -m abmptools.amorphous.pubchem --name ketoprofen --smiles

If the requested compound has no 3D conformer on PubChem (common for
polymers, metal complexes, or unusual salts), :class:`PubChemNo3DError`
is raised with an actionable message.
"""
from __future__ import annotations

import logging
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)

PUBCHEM_REST = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
_DEFAULT_UA = "abmptools/amorphous-pubchem"
_VALID_BY = ("cid", "name", "smiles", "inchi", "inchikey")


class PubChemError(RuntimeError):
    """Base class for PubChem-related failures."""


class PubChemNotFoundError(PubChemError):
    """Raised when PubChem returns 404 (no such compound / no such record)."""


class PubChemNo3DError(PubChemNotFoundError):
    """Raised when the requested compound exists but has no 3D conformer."""


def _validate_by(by: str) -> str:
    by = by.lower()
    if by not in _VALID_BY:
        raise ValueError(
            f"'by' must be one of {_VALID_BY}, got {by!r}."
        )
    return by


def _fetch(url: str, timeout: float = 30.0) -> str:
    """GET ``url`` and return the response body decoded as text.

    Raises :class:`PubChemNotFoundError` on 404 and :class:`PubChemError`
    on any other HTTP / network error.
    """
    logger.debug("PubChem GET %s", url)
    req = urllib.request.Request(url, headers={"User-Agent": _DEFAULT_UA})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            return resp.read().decode("utf-8")
    except urllib.error.HTTPError as e:
        if e.code == 404:
            raise PubChemNotFoundError(
                f"PubChem returned 404 for {url}"
            ) from e
        raise PubChemError(
            f"PubChem API error (HTTP {e.code}): {url}"
        ) from e
    except urllib.error.URLError as e:
        raise PubChemError(
            f"PubChem connection error for {url}: {e.reason}"
        ) from e


def _query_segment(query: object, by: str) -> str:
    """URL-encode the query term for the PUG REST path."""
    return urllib.parse.quote(str(query), safe="")


def fetch_3d_sdf(query: object, by: str = "cid",
                 timeout: float = 30.0) -> str:
    """Fetch a 3D SDF (with hydrogens, MMFF94-optimized) from PubChem.

    Parameters
    ----------
    query : object
        CID (int/str) or compound name / SMILES / InChI / InChIKey.
    by : str
        Lookup type. One of ``"cid"`` (default), ``"name"``, ``"smiles"``,
        ``"inchi"``, ``"inchikey"``.
    timeout : float
        HTTP timeout in seconds.

    Returns
    -------
    str
        The full SDF text.

    Raises
    ------
    PubChemNo3DError
        If the compound exists but PubChem has no 3D conformer for it.
    PubChemNotFoundError
        If the compound itself is not found.
    PubChemError
        For other API or network errors.
    """
    by = _validate_by(by)
    url = (
        f"{PUBCHEM_REST}/compound/{by}/{_query_segment(query, by)}"
        f"/SDF?record_type=3d"
    )
    try:
        return _fetch(url, timeout=timeout)
    except PubChemNotFoundError as e:
        raise PubChemNo3DError(
            f"No 3D conformer available on PubChem for {by}={query!r} "
            "(the compound may not exist, or may exist only as a 2D "
            "record — common for polymers, salts, and metal complexes). "
            "Consider using fetch_smiles() and letting OpenFF generate "
            "a conformer instead."
        ) from e


def fetch_smiles(query: object, by: str = "cid",
                 timeout: float = 30.0) -> str:
    """Fetch the canonical SMILES for a compound from PubChem.

    Parameters
    ----------
    query, by, timeout : see :func:`fetch_3d_sdf`.

    Returns
    -------
    str
        Canonical SMILES string (whitespace stripped).

    Raises
    ------
    PubChemNotFoundError
        If the compound is not found.
    """
    by = _validate_by(by)
    url = (
        f"{PUBCHEM_REST}/compound/{by}/{_query_segment(query, by)}"
        f"/property/CanonicalSMILES/TXT"
    )
    return _fetch(url, timeout=timeout).strip()


def download_3d_sdf(query: object, output_path: str,
                    by: str = "cid", timeout: float = 30.0) -> str:
    """Fetch a 3D SDF and write it to ``output_path``.

    The parent directory is created if necessary. Returns the absolute
    path of the written file.
    """
    sdf_text = fetch_3d_sdf(query, by=by, timeout=timeout)
    out = Path(output_path)
    if out.parent and str(out.parent) != "":
        out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(sdf_text)
    logger.info("Wrote 3D SDF (%d chars) to %s", len(sdf_text), out)
    return str(out.resolve())


# ---- CLI ----

def _build_cli_parser():
    import argparse

    p = argparse.ArgumentParser(
        prog="abmptools.amorphous.pubchem",
        description=(
            "Download a 3D SDF or canonical SMILES from PubChem. "
            "Raises an explicit error if no 3D conformer is available."
        ),
    )
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument("--cid", help="PubChem Compound ID (integer)")
    group.add_argument("--name", help="Compound name (e.g. 'ketoprofen')")
    group.add_argument("--smiles", dest="smiles_query",
                       help="SMILES to look up")
    group.add_argument("--inchi", help="InChI to look up")
    group.add_argument("--inchikey", help="InChIKey to look up")

    p.add_argument("--output", "-o",
                   help="Output SDF path. Default: <query>.sdf in cwd.")
    p.add_argument("--smiles-only", action="store_true",
                   help="Print canonical SMILES to stdout instead of "
                        "downloading the SDF.")
    p.add_argument("--timeout", type=float, default=30.0,
                   help="HTTP timeout in seconds (default: 30)")
    p.add_argument("-v", "--verbose", action="store_true",
                   help="Enable debug logging")
    return p


def main(argv=None) -> int:
    """CLI entry point."""
    p = _build_cli_parser()
    args = p.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    # Determine query value and lookup key
    if args.cid is not None:
        query, by = args.cid, "cid"
    elif args.name is not None:
        query, by = args.name, "name"
    elif args.smiles_query is not None:
        query, by = args.smiles_query, "smiles"
    elif args.inchi is not None:
        query, by = args.inchi, "inchi"
    else:
        query, by = args.inchikey, "inchikey"

    try:
        if args.smiles_only:
            smi = fetch_smiles(query, by=by, timeout=args.timeout)
            print(smi)
            return 0

        # Default: download 3D SDF
        output = args.output
        if not output:
            safe = urllib.parse.quote(str(query), safe="").replace("%", "_")
            output = f"{safe}.sdf"
        path = download_3d_sdf(query, output, by=by, timeout=args.timeout)
        print(path)
        return 0

    except PubChemNo3DError as e:
        logger.error("No 3D conformer: %s", e)
        return 3
    except PubChemNotFoundError as e:
        logger.error("Not found: %s", e)
        return 2
    except PubChemError as e:
        logger.error("PubChem error: %s", e)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
