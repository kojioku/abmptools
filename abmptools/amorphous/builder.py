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
from .mdp_protocol import (
    write_all_mdp,
    write_run_script,
    write_udf_export_script,
    write_wrap_script,
)
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

        # Stage 5b: position_restraints for trimer cluster center (if frozen
        # atoms specified). freezegrps was tried first but failed Fugaku MPI:
        # LINCS/SETTLE matrix inversion failed (`determinant = -inf`) because
        # hard freeze + DD + rigid water constraints are incompatible. Switch
        # to harmonic position_restraints (k=10000 kJ/mol/nm² default) which
        # plays nice with all constraints and keeps atoms within ~0.01 Å of
        # initial position — functionally equivalent for trimer center fix.
        frozen = list(getattr(self.config, "frozen_atom_indices", []) or [])
        define_posres = None
        if frozen:
            self._add_trimer_posres_to_top(top_path, frozen)
            define_posres = "POSRES_TRIMER"

        # Stage 6: Write MDP files and run script
        logger.info("=== Stage 6: MDP protocol ===")
        protocol = self._make_protocol()
        tc_grps = self._tc_grps_string()
        mdp_files = write_all_mdp(protocol, str(md_dir), tc_grps=tc_grps,
                                  define_posres=define_posres)
        run_script = write_run_script(str(md_dir))
        wrap_script = write_wrap_script(str(md_dir))
        udf_script = write_udf_export_script(str(md_dir))

        logger.info("=== Build complete: %s ===", cfg.output_dir)
        return {
            "output_dir": str(out.resolve()),
            "gro": gromacs_files["gro"],
            "top": gromacs_files["top"],
            "ndx": str(Path(ndx_path).resolve()),
            "mdp_files": mdp_files,
            "run_script": run_script,
            "wrap_script": wrap_script,
            "udf_script": udf_script,
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
                pdb_path=comp.pdb_path,
                name=comp.name or f"comp_{i}",
                charge_method=self.config.charge_method,
                nagl_model=self.config.nagl_model,
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
        # Optional pre-built cluster (UDF cluster_file equivalent).
        # When cfg.cluster_pdb_path is set, packmol places the cluster
        # rigidly at box center first. We reduce the first component's
        # count by the cluster's mol count (n_trimer derived from
        # frozen_atom_indices) so total mol count stays the same.
        cluster_pdb = str(getattr(self.config, 'cluster_pdb_path', '') or '')
        cluster_pdb = cluster_pdb if cluster_pdb else None
        adjusted_counts = list(self._counts)
        if cluster_pdb and self._molecules:
            atoms_per_first = int(self._molecules[0].n_atoms)
            frozen = list(getattr(self.config, 'frozen_atom_indices', []) or [])
            if frozen:
                n_cluster_mol = max((int(a) - 1) // atoms_per_first
                                    for a in frozen) + 1
                adjusted_counts[0] = max(0, adjusted_counts[0] - n_cluster_mol)
                logger.info(
                    "cluster_pdb=%s placed as rigid block; first-component "
                    "count reduced by %d → %d", cluster_pdb, n_cluster_mol,
                    adjusted_counts[0],
                )
        return run_packmol(
            pdb_paths=self._pdb_paths,
            counts=adjusted_counts,
            box_size_nm=self._box_nm,
            output_pdb=output_pdb,
            build_dir=build_dir,
            tolerance=self.config.packmol_tolerance,
            seed=self.config.seed,
            packmol_path=self.config.packmol_path,
            cluster_pdb=cluster_pdb,
        )

    def _parameterize(self, mixture_pdb: str,
                      gro_path: str, top_path: str) -> Dict[str, str]:
        # When the user opted into nagl / gasteiger we pre-assigned charges
        # on each molecule; tell Interchange to reuse them instead of
        # invoking AM1-BCC via sqm (which has no Windows build).
        use_precomputed = self.config.charge_method in ("nagl", "gasteiger")
        interchange = create_interchange(
            molecules=self._molecules,
            counts=self._counts,
            box_size_nm=self._box_nm,
            mixture_pdb=mixture_pdb,
            forcefield_name=self.config.forcefield,
            use_precomputed_charges=use_precomputed,
        )
        return export_gromacs(interchange, gro_path, top_path)

    def _write_ndx(self, ndx_path: str) -> str:
        comp_names = [c.name or f"comp_{i}"
                      for i, c in enumerate(self.config.components)]
        atom_counts = [mol.n_atoms for mol in self._molecules]
        frozen = list(getattr(self.config, "frozen_atom_indices", []) or [])
        return write_ndx(comp_names, atom_counts, self._counts, ndx_path,
                         frozen_atom_indices=frozen)

    def _add_trimer_posres_to_top(self, top_path: str,
                                  frozen_atoms_1based: list) -> None:
        """Add [ position_restraints ] for trimer cluster center to system.top.

        Derives local atom indices and trimer molecule count from the global
        1-based ``frozen_atoms_1based`` and the first component's atom count.
        If the first moleculetype's molecule count > trimer count, splits the
        moletype: a ``<name>_TRIMER`` copy (trimer count instances, with
        posres) is inserted BEFORE the original moleculetype (now with
        reduced count). Otherwise, the posres block is appended to the
        existing first moletype.

        The posres block is wrapped in ``#ifdef POSRES_TRIMER`` so anneal /
        sampling mdp can toggle via ``define = -DPOSRES_TRIMER``.

        Why posres (not freezegrps): hard freeze on rigid water O atoms
        breaks LINCS / SETTLE matrix inversion under MPI Domain
        Decomposition (`determinant = -inf`). Harmonic posres
        (k = posres_force_constant, default 10000 kJ/mol/nm²) keeps the
        atoms within ~0.01 Å of initial position without any constraint
        conflict — functionally equivalent for trimer cluster fix.
        """
        if not self._molecules:
            return
        atoms_per_first_mol = int(self._molecules[0].n_atoms)
        # Local 1-based atom indices within the first moletype
        local_atoms = sorted({(int(a) - 1) % atoms_per_first_mol + 1
                              for a in frozen_atoms_1based})
        # Trimer mol count = max mol-index used + 1
        n_trimer = max((int(a) - 1) // atoms_per_first_mol
                       for a in frozen_atoms_1based) + 1
        fc = float(getattr(self.config, 'posres_force_constant', 10000.0))

        text = Path(top_path).read_text()

        # Find first [ moleculetype ] block. Block extends until next
        # [ moleculetype ] OR [ system ] marker.
        import re
        starts = [m.start() for m in re.finditer(r'^\[ moleculetype \]',
                                                  text, flags=re.MULTILINE)]
        if not starts:
            return
        end_match = re.search(r'^\[ (?:moleculetype|system) \]', text[starts[0] + 1:],
                              flags=re.MULTILINE)
        if end_match is None:
            block_end = len(text)
        else:
            block_end = starts[0] + 1 + end_match.start()
        first_block = text[starts[0]:block_end]

        # First non-comment line after [ moleculetype ] is "<name>\t<nrexcl>"
        moltype_name = None
        for line in first_block.split('\n')[1:]:
            s = line.strip()
            if not s or s.startswith(';'):
                continue
            moltype_name = s.split()[0]
            break
        if not moltype_name:
            return

        # Find first entry in [ molecules ] section
        mol_section_match = re.search(r'^\[ molecules \]', text, flags=re.MULTILINE)
        if mol_section_match is None:
            return
        mol_start = mol_section_match.end()
        first_count = None
        first_line_text = None
        for line in text[mol_start:].split('\n'):
            s = line.strip()
            if not s or s.startswith(';') or s.startswith('['):
                continue
            parts = s.split()
            if len(parts) >= 2 and parts[0] == moltype_name:
                try:
                    first_count = int(parts[1])
                    first_line_text = line
                    break
                except ValueError:
                    pass
        if first_count is None:
            return

        # Build posres block
        posres_lines = ['', '#ifdef POSRES_TRIMER',
                        '[ position_restraints ]',
                        '; ai funct kx ky kz']
        for ai in local_atoms:
            posres_lines.append(f'  {ai}  1  {fc:g}  {fc:g}  {fc:g}')
        posres_lines.append('#endif')
        posres_str = '\n'.join(posres_lines) + '\n'

        needs_split = (first_count > n_trimer)
        if needs_split:
            # Insert TRIMER copy BEFORE the original. Rename only the
            # moltype declaration line (first non-comment line after the
            # [ moleculetype ] marker).
            trimer_block_lines = first_block.split('\n')
            for i in range(1, len(trimer_block_lines)):
                s = trimer_block_lines[i].strip()
                if not s or s.startswith(';'):
                    continue
                if s.split()[0] == moltype_name:
                    trimer_block_lines[i] = trimer_block_lines[i].replace(
                        moltype_name, moltype_name + '_TRIMER', 1)
                    break
            trimer_block_text = '\n'.join(trimer_block_lines)
            # GROMACS forbids 2 moltypes both having [ settles ]
            # (`Only one such is allowed`). Convert TRIMER copy's settles
            # to equivalent [ constraints ] (LINCS-based) so the original
            # moltype keeps its faster SETTLE and TRIMER becomes
            # constraint-based rigid water. For TIP3P/SPC 3-site water:
            #   [ settles ]
            #     1  1  dOH  dHH
            # → [ constraints ]
            #     1 2 1 dOH    ; O-H1
            #     1 3 1 dOH    ; O-H2
            #     2 3 1 dHH    ; H1-H2
            settles_block_pattern = re.compile(
                r'(\[ settles \].*?)(?=\n\[ |\Z)', re.DOTALL,
            )
            settles_match = settles_block_pattern.search(trimer_block_text)
            if settles_match:
                settles_text = settles_match.group(1)
                # Parse the first non-comment row: "i funct dOH dHH"
                dOH = dHH = None
                for line in settles_text.split('\n')[1:]:
                    s = line.strip()
                    if not s or s.startswith(';') or s.startswith('['):
                        continue
                    parts = s.split()
                    if len(parts) >= 4:
                        try:
                            dOH = float(parts[2])
                            dHH = float(parts[3])
                            break
                        except ValueError:
                            continue
                if dOH is not None:
                    constraints_text = (
                        f'[ constraints ]\n'
                        f'; ai aj funct value (TIP3P/SPC rigid water, '
                        f'converted from [settles])\n'
                        f'  1 2 1 {dOH}\n'
                        f'  1 3 1 {dOH}\n'
                        f'  2 3 1 {dHH}\n'
                    )
                    trimer_block_text = trimer_block_text.replace(
                        settles_text, constraints_text, 1)
            trimer_block = trimer_block_text.rstrip() + '\n' + posres_str + '\n'
            text = text[:starts[0]] + trimer_block + text[starts[0]:]
            # Update the [ molecules ] line
            new_first_line = (
                f'{moltype_name}_TRIMER\t{n_trimer}\n'
                f'{moltype_name}\t{first_count - n_trimer}'
            )
            text = text.replace(first_line_text, new_first_line, 1)
        else:
            # No split — append posres at end of first moltype block
            new_block = first_block.rstrip() + '\n' + posres_str + '\n'
            text = text.replace(first_block, new_block, 1)

        Path(top_path).write_text(text)
        logger.info(
            "added [ position_restraints ] (k=%g) for moltype %s (n_trimer=%d, "
            "local_atoms=%s, split=%s)",
            fc, moltype_name, n_trimer, local_atoms, needs_split,
        )

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
        """Build the tc-grps string for MDP (e.g. 'API Polymer').

        Names are deduplicated while preserving first-seen order so that
        pure-component pairs (e.g. ``ComponentSpec(name='B_water', ...)``
        × 2) emit a single group, matching the single ``ref-t`` /
        ``tau-t`` value emitted by :func:`mdp_protocol._thermostat_block`.
        Without dedup, grompp rejects the mdp with
        ``Invalid T coupling input: 2 groups, 1 ref-t values``.
        """
        seen: dict[str, None] = {}
        for i, c in enumerate(self.config.components):
            name = c.name or f"comp_{i}"
            if name not in seen:
                seen[name] = None
        return " ".join(seen.keys())
