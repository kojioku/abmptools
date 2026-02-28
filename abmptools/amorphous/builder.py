# -*- coding: utf-8 -*-
"""
abmptools.amorphous.builder
-----------------------------
Orchestrator: integrates all stages of amorphous structure building.

Data flow:
  SMILES/SDF → molecule_prep → density → packing → parameterizer → mdp/ndx
"""
from __future__ import annotations

import json
import logging
import os
from pathlib import Path
from typing import Any, Dict, List, Optional

from ..core.system_model import AnnealProtocol
from .models import BuildConfig, ComponentSpec
from .density import estimate_box_size_nm, weight_fractions_to_counts
from .molecule_prep import (
    prepare_molecule,
    get_molecular_weight,
    write_single_mol_pdb,
)
from .packing import run_packmol
from .parameterizer import create_interchange, export_gromacs
from .mdp_protocol import write_all_mdp, write_run_script
from .ndx_writer import write_ndx

logger = logging.getLogger(__name__)


class AmorphousBuilder:
    """Build an amorphous multi-component system for GROMACS MD.

    Usage::

        from abmptools.amorphous import AmorphousBuilder, BuildConfig, ComponentSpec

        config = BuildConfig(
            components=[
                ComponentSpec(smiles="CCCCC", name="pentane", n_mol=200),
                ComponentSpec(smiles="c1ccccc1", name="benzene", n_mol=50),
            ],
            density_g_cm3=0.8,
            temperature=300,
            output_dir="./output",
        )
        builder = AmorphousBuilder(config)
        result = builder.build()
    """

    def __init__(self, config: BuildConfig) -> None:
        self.config = config
        self._molecules: List[Any] = []
        self._mol_weights: List[float] = []
        self._counts: List[int] = []
        self._box_nm: float = 0.0
        self._pdb_paths: List[str] = []

    def build(self) -> Dict[str, Any]:
        """Run the full build pipeline.

        Returns
        -------
        dict
            Keys: "output_dir", "gro", "top", "ndx", "mdp_files",
                  "run_script", "config_json", "box_nm", "counts".
        """
        cfg = self.config
        out = Path(cfg.output_dir)
        input_dir = out / "input"
        build_dir = out / "build"
        md_dir = out / "md"
        for d in [input_dir, build_dir, md_dir]:
            d.mkdir(parents=True, exist_ok=True)

        # Save config for reproducibility
        config_json = str(input_dir / "config.json")
        cfg.to_json(config_json)

        # Stage 1: Prepare molecules
        logger.info("=== Stage 1: Molecule preparation ===")
        self._prepare_molecules(input_dir)

        # Stage 2: Determine counts and box size
        logger.info("=== Stage 2: Density / box estimation ===")
        self._compute_counts_and_box()

        # Stage 3: Packmol packing
        logger.info("=== Stage 3: Packmol packing ===")
        mixture_pdb = str(build_dir / "mixture.pdb")
        self._run_packing(mixture_pdb, str(build_dir))

        # Stage 4: Parameterize with OpenFF
        logger.info("=== Stage 4: OpenFF parameterization ===")
        gro_path = str(build_dir / "system.gro")
        top_path = str(build_dir / "system.top")
        gromacs_files = self._parameterize(mixture_pdb, gro_path, top_path)

        # Stage 5: Write NDX
        logger.info("=== Stage 5: NDX index groups ===")
        ndx_path = str(build_dir / "system.ndx")
        self._write_ndx(ndx_path)

        # Stage 6: Write MDP files and run script
        logger.info("=== Stage 6: MDP protocol ===")
        protocol = self._make_protocol()
        tc_grps = self._tc_grps_string()
        mdp_files = write_all_mdp(protocol, str(md_dir), tc_grps=tc_grps)
        run_script = write_run_script(str(md_dir))

        logger.info("=== Build complete: %s ===", cfg.output_dir)
        return {
            "output_dir": str(out.resolve()),
            "gro": gromacs_files["gro"],
            "top": gromacs_files["top"],
            "ndx": str(Path(ndx_path).resolve()),
            "mdp_files": mdp_files,
            "run_script": run_script,
            "config_json": config_json,
            "box_nm": self._box_nm,
            "counts": list(self._counts),
        }

    # ---- internal stages ----

    def _prepare_molecules(self, input_dir: Path) -> None:
        """Prepare OpenFF Molecules and write single-molecule PDBs."""
        for i, comp in enumerate(self.config.components):
            mol = prepare_molecule(
                smiles=comp.smiles,
                sdf_path=comp.sdf_path,
                name=comp.name or f"comp_{i}",
            )
            mw = get_molecular_weight(mol)
            comp.molecular_weight = mw
            if not comp.name:
                comp.name = mol.name

            pdb_path = str(input_dir / f"component_{i}.pdb")
            write_single_mol_pdb(mol, pdb_path)

            self._molecules.append(mol)
            self._mol_weights.append(mw)
            self._pdb_paths.append(pdb_path)

    def _compute_counts_and_box(self) -> None:
        """Determine molecule counts and box size."""
        cfg = self.config

        # Determine counts
        has_n_mol = any(c.n_mol > 0 for c in cfg.components)
        has_wf = any(c.weight_fraction > 0 for c in cfg.components)

        if has_n_mol:
            self._counts = [c.n_mol for c in cfg.components]
        elif has_wf:
            wfs = [c.weight_fraction for c in cfg.components]
            total = cfg.total_molecules if cfg.total_molecules > 0 else 200
            self._counts = weight_fractions_to_counts(
                self._mol_weights, wfs, total
            )
        else:
            raise ValueError(
                "Specify either n_mol or weight_fraction for each component."
            )

        # Determine box size
        if cfg.box_size_nm > 0:
            self._box_nm = cfg.box_size_nm
        else:
            self._box_nm = estimate_box_size_nm(
                self._mol_weights, self._counts, cfg.density_g_cm3
            )
        logger.info("Box size: %.4f nm, Counts: %s", self._box_nm, self._counts)

    def _run_packing(self, output_pdb: str, build_dir: str) -> str:
        return run_packmol(
            pdb_paths=self._pdb_paths,
            counts=self._counts,
            box_size_nm=self._box_nm,
            output_pdb=output_pdb,
            build_dir=build_dir,
            tolerance=self.config.packmol_tolerance,
            seed=self.config.seed,
            packmol_path=self.config.packmol_path,
        )

    def _parameterize(self, mixture_pdb: str,
                      gro_path: str, top_path: str) -> Dict[str, str]:
        interchange = create_interchange(
            molecules=self._molecules,
            counts=self._counts,
            box_size_nm=self._box_nm,
            mixture_pdb=mixture_pdb,
            forcefield_name=self.config.forcefield,
        )
        return export_gromacs(interchange, gro_path, top_path)

    def _write_ndx(self, ndx_path: str) -> str:
        comp_names = [c.name or f"comp_{i}"
                      for i, c in enumerate(self.config.components)]
        atom_counts = [mol.n_atoms for mol in self._molecules]
        return write_ndx(comp_names, atom_counts, self._counts, ndx_path)

    def _make_protocol(self) -> AnnealProtocol:
        cfg = self.config
        return AnnealProtocol(
            T_high=cfg.T_high,
            T_low=cfg.temperature,
            P_ref=cfg.pressure,
            em_steps=cfg.em_steps,
            em_tol=cfg.em_tol,
            nvt_high_nsteps=cfg.nvt_high_nsteps,
            npt_high_nsteps=cfg.npt_high_nsteps,
            anneal_nsteps=cfg.anneal_nsteps,
            npt_low_nsteps=cfg.npt_low_nsteps,
            dt=cfg.dt,
            tau_t=cfg.tau_t,
            tau_p=cfg.tau_p,
            nstxout_compressed=cfg.nstxout_compressed,
            nstenergy=cfg.nstenergy,
            gen_seed=cfg.seed,
        )

    def _tc_grps_string(self) -> str:
        """Build the tc-grps string for MDP (e.g. 'API Polymer')."""
        names = [c.name or f"comp_{i}"
                 for i, c in enumerate(self.config.components)]
        return " ".join(names)
