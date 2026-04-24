# -*- coding: utf-8 -*-
"""
abmptools.amorphous.cli
-------------------------
Command-line interface for the amorphous structure builder.
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import List, Optional

from .models import BuildConfig, ComponentSpec


def _parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="amorphous-builder",
        description="Build amorphous multi-component systems for GROMACS MD.",
    )

    # --- Input ---
    inp = p.add_argument_group("Input molecules")
    inp.add_argument("--mol", nargs="+", metavar="FILE",
                     help="SDF/MOL files for each component")
    inp.add_argument("--smiles", nargs="+", metavar="SMILES",
                     help="SMILES strings for each component")
    inp.add_argument("--name", nargs="+", metavar="NAME",
                     help="Names for each component (optional)")
    inp.add_argument("--pubchem_cid", nargs="+", metavar="CID",
                     help="PubChem CIDs; each CID is resolved to a 3D SDF "
                          "and treated as --mol input for that component. "
                          "Raises an error if PubChem has no 3D conformer.")
    inp.add_argument("--pubchem_name", nargs="+", metavar="NAME",
                     help="PubChem compound names; resolved to 3D SDF like "
                          "--pubchem_cid. Useful when the CID is unknown.")
    inp.add_argument("--pubchem_cache_dir", type=str, default=None,
                     help="Directory to cache downloaded PubChem SDFs "
                          "(default: <output_dir>/input/)")

    # --- Composition ---
    comp = p.add_argument_group("Composition")
    comp.add_argument("--n_mol", nargs="+", type=int, metavar="N",
                      help="Number of molecules for each component")
    comp.add_argument("--weight_fraction", nargs="+", type=float, metavar="W",
                      help="Weight fraction for each component (sum ~1.0)")
    comp.add_argument("--total_molecules", type=int, default=0,
                      help="Total molecules (with --weight_fraction)")

    # --- Box / density ---
    box = p.add_argument_group("Box / density")
    box.add_argument("--box", type=float, default=0.0,
                     help="Box edge length [nm] (0 = auto from density)")
    box.add_argument("--density", type=float, default=0.8,
                     help="Target density [g/cm^3] (default: 0.8)")

    # --- MD parameters ---
    md = p.add_argument_group("MD parameters")
    md.add_argument("--temperature", type=float, default=300.0,
                    help="Final temperature [K] (default: 300)")
    md.add_argument("--T_high", type=float, default=600.0,
                    help="High temperature [K] (default: 600)")
    md.add_argument("--seed", type=int, default=None,
                    help="Random seed")
    md.add_argument("--forcefield", type=str,
                    default="openff_unconstrained-2.1.0.offxml",
                    help="OpenFF force field (default: openff_unconstrained-2.1.0.offxml)")
    md.add_argument("--charge_method", type=str, default="",
                    choices=["", "am1bcc", "nagl", "gasteiger"],
                    help="Partial charge backend: '' / 'am1bcc' (default; "
                         "Interchange's AM1-BCC via AmberTools sqm — requires "
                         "Linux/macOS), 'nagl' (ML AM1-BCC via openff-nagl; "
                         "Windows-compatible), or 'gasteiger' (fast fallback).")
    md.add_argument("--nagl_model", type=str,
                    default="openff-gnn-am1bcc-0.1.0-rc.3.pt",
                    help="openff-nagl model file (used only with "
                         "--charge_method nagl). Check available models via "
                         "openff.nagl_models.list_available_nagl_models().")

    # --- Packmol ---
    pk = p.add_argument_group("Packmol")
    pk.add_argument("--packmol_tolerance", type=float, default=2.0,
                    help="Packmol tolerance [A] (default: 2.0)")
    pk.add_argument("--packmol_path", type=str, default="packmol",
                    help="Path to packmol binary")

    # --- Output ---
    out = p.add_argument_group("Output")
    out.add_argument("--output_dir", type=str, default=".",
                     help="Output directory (default: .)")

    # --- JSON config ---
    p.add_argument("--config", type=str, default=None,
                   help="JSON config file (overrides all other args)")

    # --- Verbosity ---
    p.add_argument("-v", "--verbose", action="store_true",
                   help="Enable verbose logging")

    return p.parse_args(argv)


def _resolve_pubchem_inputs(args: argparse.Namespace) -> List[str]:
    """Download PubChem 3D SDFs requested via --pubchem_cid / --pubchem_name.

    Returns a list of local SDF paths, in the order they should be
    appended to the mol/smiles input lists.
    """
    cids = args.pubchem_cid or []
    names = args.pubchem_name or []
    if not cids and not names:
        return []

    from .pubchem import download_3d_sdf  # lazy import

    cache_dir = args.pubchem_cache_dir
    if cache_dir is None:
        # Default: <output_dir>/input/
        cache_dir = str(Path(args.output_dir) / "input")
    Path(cache_dir).mkdir(parents=True, exist_ok=True)

    resolved: List[str] = []
    for cid in cids:
        out = str(Path(cache_dir) / f"pubchem_cid{cid}.sdf")
        path = download_3d_sdf(cid, out, by="cid")
        resolved.append(path)
    for name in names:
        safe = "".join(c if c.isalnum() else "_" for c in name)
        out = str(Path(cache_dir) / f"pubchem_name_{safe}.sdf")
        path = download_3d_sdf(name, out, by="name")
        resolved.append(path)
    return resolved


def _build_config_from_args(args: argparse.Namespace) -> BuildConfig:
    """Build a BuildConfig from parsed CLI args."""
    if args.config:
        return BuildConfig.from_json(args.config)

    # Resolve PubChem queries (if any) to local SDF paths and prepend /
    # append them to --mol inputs.
    pubchem_sdfs = _resolve_pubchem_inputs(args)

    # Determine number of components
    n_mol_inputs = list(args.mol or []) + pubchem_sdfs
    n_smiles_inputs = list(args.smiles or [])
    n_components = max(len(n_mol_inputs), len(n_smiles_inputs))

    if n_components == 0:
        raise ValueError(
            "Provide --mol, --smiles, --pubchem_cid, or --pubchem_name for "
            "at least one component."
        )

    names = args.name or []
    n_mol_list = args.n_mol or []
    wf_list = args.weight_fraction or []

    components = []
    for i in range(n_components):
        smiles = n_smiles_inputs[i] if i < len(n_smiles_inputs) else ""
        sdf_path = n_mol_inputs[i] if i < len(n_mol_inputs) else ""
        name = names[i] if i < len(names) else ""
        n_mol = n_mol_list[i] if i < len(n_mol_list) else 0
        wf = wf_list[i] if i < len(wf_list) else 0.0

        components.append(ComponentSpec(
            name=name,
            smiles=smiles,
            sdf_path=sdf_path,
            n_mol=n_mol,
            weight_fraction=wf,
        ))

    return BuildConfig(
        components=components,
        total_molecules=args.total_molecules,
        box_size_nm=args.box,
        density_g_cm3=args.density,
        temperature=args.temperature,
        T_high=args.T_high,
        seed=args.seed,
        forcefield=args.forcefield,
        packmol_tolerance=args.packmol_tolerance,
        packmol_path=args.packmol_path,
        output_dir=args.output_dir,
        charge_method=args.charge_method,
        nagl_model=args.nagl_model,
    )


def main(argv: Optional[List[str]] = None) -> None:
    """CLI entry point."""
    args = _parse_args(argv)

    level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    config = _build_config_from_args(args)

    from .builder import AmorphousBuilder
    builder = AmorphousBuilder(config)
    result = builder.build()

    print(f"\nBuild complete!")
    print(f"  Output dir : {result['output_dir']}")
    print(f"  GRO        : {result['gro']}")
    print(f"  TOP        : {result['top']}")
    print(f"  NDX        : {result['ndx']}")
    print(f"  Box size   : {result['box_nm']:.4f} nm")
    print(f"  Counts     : {result['counts']}")
    print(f"  Run script : {result['run_script']}")
    print(f"  Wrap script: {result['wrap_script']}")
