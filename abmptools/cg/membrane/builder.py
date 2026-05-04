# -*- coding: utf-8 -*-
"""
abmptools.cg.membrane.builder
------------------------------
End-to-end Martini 3 peptide-membrane US system builder.

7-stage pipeline (:meth:`MembraneCGBuilder.build`):

    0. ``_copy_ff_files``     -- 4 Martini 3 ITPs -> ``output_dir/``
    1. ``_stage1_cg_peptide`` -- :class:`PeptideCGBuilder` sub-call
                                  (``solvent_enabled=False, mdp_*=False``)
                                  produces ``molecules/<name>/<name>_cg.pdb``
                                  + ``<name>.itp``.
    2. ``_stage2_insane``     -- ``insane`` embeds peptide in POPC bilayer.
                                  Outputs raw ``insane_topol.top`` and
                                  ``bilayer.gro``.
    3. ``_stage3_topology``   -- :func:`topology_composer.compose_topology`
                                  rewrites ITP includes / molecule name /
                                  ion names. Output: ``topol.top``.
    4. ``_stage4_ions``       -- normalize ``NA+/CL-`` in .gro, then
                                  ``gmx grompp`` + ``gmx genion`` (re-using
                                  cg.peptide ``add_ions``). Output:
                                  ``system_ions.gro``.
    5. ``_stage5_index``      -- :func:`system_packer.write_ndx_from_gro_cg`
                                  -> ``index.ndx`` with named groups.
    6. ``_stage6_mdps``       -- em / nvt / npt / pull / window MDPs.
    7. ``_stage7_run_script`` -- ``run.sh``.

Output layout::

    output_dir/
      input/config.json           # config dump (re-loaded by run.sh stages)
      martini_v3.0.0*.itp (4 files)
      molecules/<name>/{name}_atomistic.pdb, {name}_cg.pdb, {name}.itp
      bilayer.gro                 # raw insane output
      insane_topol.top            # raw insane topology (kept for audit)
      system_ions.gro             # final coordinates
      topol.top                   # composed topology
      index.ndx
      ions.{mdp,tpr}              # genion bookkeeping
      mdp/{em,nvt,npt}.mdp        # equilibration MDPs
      pull/pull.mdp               # pulling MDP
      windows/win_NNN/window.mdp  # 13 per-window MDPs (NN = 000..012)
      run.sh                      # top-level orchestrator
"""
from __future__ import annotations

import logging
import shutil
from pathlib import Path
from typing import Any, Dict, Optional

from abmptools.cg.peptide.builder import PeptideCGBuilder
from abmptools.cg.peptide.models import (
    PeptideBuildConfig,
    PeptideSpec,
)

from . import (
    forcefield_check,
    insane_runner,
    mdp_templates,
    pulling,
    system_packer,
    topology_composer,
    umbrella,
)
from ._subprocess import ensure_dir, write_text
from .models import MembraneCGBuildConfig

logger = logging.getLogger(__name__)


class MembraneCGBuilder:
    """Build a Martini 3 peptide-membrane US system from a
    :class:`MembraneCGBuildConfig`."""

    def __init__(self, config: MembraneCGBuildConfig):
        self.config = config
        self.output_dir = ensure_dir(Path(config.output_dir).resolve())
        self.input_dir = ensure_dir(self.output_dir / "input")
        self.molecules_dir = ensure_dir(self.output_dir / "molecules")
        self.mdp_dir = ensure_dir(self.output_dir / "mdp")
        self.pull_dir = ensure_dir(self.output_dir / "pull")
        self.windows_dir = ensure_dir(self.output_dir / "windows")
        self.peptide_cg_pdb: Optional[Path] = None
        self.peptide_itp: Optional[Path] = None
        self.peptide_itp_relpath: Optional[str] = None

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def build(self) -> Dict[str, Any]:
        """Run the full 7-stage pipeline. Returns a result dictionary.

        The pipeline produces all GROMACS inputs but does not run any MD
        (``gmx mdrun``). The user invokes ``bash run.sh`` separately;
        ``run.sh`` orchestrates em -> nvt -> npt -> pull -> per-window mdrun
        -> wham.
        """
        if self.config.peptide is None:
            raise ValueError("MembraneCGBuildConfig.peptide is required")

        ff_dir = self._resolve_ff_dir()

        logger.info("=== Stage 0/7: copy Martini 3 ITPs to output_dir ===")
        self._copy_ff_files(ff_dir)

        logger.info("=== Stage 1/7: cg.peptide sub-call (CG peptide) ===")
        self._stage1_cg_peptide(ff_dir)

        logger.info("=== Stage 2/7: insane (peptide in bilayer) ===")
        bilayer_gro, insane_top = self._stage2_insane()

        logger.info("=== Stage 3/7: compose topology ===")
        topol_top = self._stage3_topology(insane_top)

        logger.info("=== Stage 4/7: normalize ion names + add_ions ===")
        ions_gro = self._stage4_ions(bilayer_gro, topol_top)

        logger.info("=== Stage 5/7: index file ===")
        ndx_path = self._stage5_index(ions_gro)

        logger.info("=== Stage 6/7: MDP files (em/nvt/npt/pull/windows) ===")
        equil_mdps, pull_mdp, window_mdps = self._stage6_mdps()

        logger.info("=== Stage 7/7: run.sh ===")
        run_script = self._stage7_run_script(
            top=topol_top, gro=ions_gro, ndx=ndx_path,
            equil_mdps=equil_mdps,
            pull_mdp=pull_mdp,
            window_mdps=window_mdps,
        )

        config_json = self.input_dir / "config.json"
        self.config.to_json(str(config_json))

        result: Dict[str, Any] = {
            "output_dir": self.output_dir,
            "gro": ions_gro,
            "top": topol_top,
            "ndx": ndx_path,
            "equil_mdps": equil_mdps,
            "pull_mdp": pull_mdp,
            "window_mdps": window_mdps,
            "run_script": run_script,
            "config_json": config_json,
            "n_windows": self.config.umbrella.n_windows,
            "lipid_resname": self.config.lipids[0].resname,
            "n_lipids_total": (
                2 * sum(lp.n_per_leaflet for lp in self.config.lipids)
            ),
        }
        return result

    # ------------------------------------------------------------------
    # Stage helpers
    # ------------------------------------------------------------------

    def _stage1_cg_peptide(self, ff_dir: Path) -> None:
        """Sub-call :class:`PeptideCGBuilder` (stages 1+2 only) to get CG PDB+ITP.

        cg.peptide's builder writes into ``output_dir/molecules/<name>/`` and
        also drops bookkeeping files (``packed.gro`` / ``topol.top`` /
        ``index.ndx`` / ``run/run.sh``) at the sub-call's ``output_dir`` root.
        We use a dedicated subdirectory so those bookkeeping files don't
        clobber the membrane builder's outputs, then **copy** just the
        CG PDB + ITP into the final ``molecules/<name>/`` location.
        """
        spec = self.config.peptide
        # If user supplied pre-built CG files, skip the sub-call.
        if spec.cg_pdb_path and spec.cg_itp_path:
            self.peptide_cg_pdb = Path(spec.cg_pdb_path).resolve()
            self.peptide_itp = Path(spec.cg_itp_path).resolve()
            self.peptide_itp_relpath = str(self.peptide_itp)
            logger.info(
                "Using pre-built CG peptide files: %s + %s",
                self.peptide_cg_pdb, self.peptide_itp,
            )
            return

        # Build via cg.peptide in a dedicated scratch subdirectory.
        scratch = self.output_dir / "_cg_peptide_subbuild"
        scratch.mkdir(parents=True, exist_ok=True)
        sub_cfg = PeptideBuildConfig(
            peptides=[
                PeptideSpec(name=spec.name, sequence=spec.sequence, count=1),
            ],
            box_size_nm=4.0,                      # placeholder, not used
            martini_itp_dir=str(ff_dir),
            output_dir=str(scratch),
            solvent_enabled=False,                # we don't want gmx solvate
            mdp_em=False, mdp_nvt=False,
            mdp_npt=False, mdp_md=False,
            martinize2_path=self.config.martinize2_path,
            gmx_path=self.config.gmx_path,
            tleap_path=self.config.tleap_path,
        )
        PeptideCGBuilder(sub_cfg).build()

        # cg.peptide wrote molecules/<name>/{name}_cg.pdb + {name}.itp under
        # scratch. Copy the two artefacts to the final molecules/<name>/.
        src_dir = scratch / "molecules" / spec.name
        dst_dir = self.molecules_dir / spec.name
        dst_dir.mkdir(parents=True, exist_ok=True)
        cg_pdb_src = src_dir / f"{spec.name}_cg.pdb"
        itp_src = src_dir / f"{spec.name}.itp"
        cg_pdb_dst = dst_dir / f"{spec.name}_cg.pdb"
        itp_dst = dst_dir / f"{spec.name}.itp"
        shutil.copy2(cg_pdb_src, cg_pdb_dst)
        shutil.copy2(itp_src, itp_dst)
        # Optionally also copy the atomistic PDB for audit (not strictly needed).
        atom_src = src_dir / f"{spec.name}_atomistic.pdb"
        if atom_src.exists():
            shutil.copy2(atom_src, dst_dir / f"{spec.name}_atomistic.pdb")

        self.peptide_cg_pdb = cg_pdb_dst
        self.peptide_itp = itp_dst
        # Relative path so the final topol.top is portable.
        self.peptide_itp_relpath = str(
            self.peptide_itp.relative_to(self.output_dir)
        )

    def _stage2_insane(self) -> tuple:
        """Run insane: embed peptide in POPC bilayer + solvent + ions."""
        assert self.peptide_cg_pdb is not None, (
            "_stage1_cg_peptide must run before _stage2_insane"
        )
        lipid = self.config.lipids[0]
        bilayer_gro, insane_top = insane_runner.run_insane(
            output_dir=self.output_dir,
            lipid_resname=lipid.resname,
            n_per_leaflet=lipid.n_per_leaflet,
            box_d_nm=self.config.insane_d_nm,
            box_z_nm=self.config.box_z_nm,
            peptide_pdb=self.peptide_cg_pdb,
            peptide_z_offset_nm=self.config.peptide.initial_z_offset_nm,
            salt_concentration_M=self.config.nacl_molar,
            solvent=self.config.solvent_type,
            pbc=self.config.insane_pbc,
            extra_args=self.config.insane_extra_args,
            insane_path=self.config.insane_path,
        )
        return bilayer_gro, insane_top

    def _stage3_topology(self, insane_top: Path) -> Path:
        """Rewrite insane's topology with M3 ITP includes + peptide ITP."""
        topol_top = self.output_dir / "topol.top"
        assert self.peptide_itp is not None
        return topology_composer.compose_topology(
            insane_top, topol_top,
            peptide_itp_path=self.peptide_itp,
            peptide_count=1,
            peptide_itp_relpath=self.peptide_itp_relpath,
        )

    def _stage4_ions(self, bilayer_gro: Path, topol_top: Path) -> Path:
        """Normalize NA+/CL- in .gro, then add_ions if neutralization needed.

        Insane already adds the salt at the requested molar concentration
        (``-salt 0.15``) and balances charge (``-charge auto``). We only need
        to (a) normalize the atom names so ``gmx grompp`` doesn't warn, and
        (b) optionally re-run ``gmx genion`` for additional neutralisation
        if the user changed ``nacl_molar`` post-insane (rare).

        For v1 we trust insane's salt placement and only normalize names.
        Returns the path to the final coordinates.
        """
        # Normalize NA+/CL- -> NA / CL in-place (preserving column widths).
        normalized = self.output_dir / "system_ions.gro"
        topology_composer.normalize_ion_atom_names_gro(
            bilayer_gro, out_path=normalized,
        )
        return normalized

    def _stage5_index(self, gro_path: Path) -> Path:
        """Write index.ndx with Bilayer / Peptide / W / NA / CL / Non_Bilayer."""
        ndx_path = self.output_dir / "index.ndx"
        system_packer.write_ndx_from_gro_cg(
            gro_path=str(gro_path), ndx_path=str(ndx_path),
        )
        return ndx_path

    def _stage6_mdps(self) -> tuple:
        """Render em / nvt / npt MDPs + pull.mdp + per-window MDPs.

        ``pull_init_nm`` is initialized to the peptide's prescribed z offset
        (config.peptide.initial_z_offset_nm). After running NPT
        equilibration, the user (or run.sh) may want to recompute it from
        the equilibrated .gro using
        :func:`pulling.estimate_initial_pull_coord`; we record the initial
        guess in pull.mdp so ``gmx grompp`` accepts the file at build time.
        """
        equil = mdp_templates.write_equilibration_mdps(
            config=self.config, equil_dir=self.mdp_dir,
        )
        pull_mdp = pulling.write_pulling_mdp_cg(
            config=self.config, pull_dir=self.pull_dir,
            pull_init_nm=self.config.peptide.initial_z_offset_nm,
        )
        window_mdps = umbrella.write_window_mdps(
            config=self.config, windows_dir=self.windows_dir,
        )
        return equil, pull_mdp, window_mdps

    def _stage7_run_script(
        self, *, top: Path, gro: Path, ndx: Path,
        equil_mdps: Dict[str, Path],
        pull_mdp: Path,
        window_mdps: Dict[int, Path],
    ) -> Path:
        return umbrella.write_run_script(
            config=self.config,
            top=top, gro=gro, ndx=ndx,
            equil_mdps=equil_mdps,
            pull_mdp=pull_mdp,
            window_mdps=window_mdps,
            output_dir=self.output_dir,
        )

    # ------------------------------------------------------------------
    # ITP / ff_dir handling
    # ------------------------------------------------------------------

    def _resolve_ff_dir(self) -> Path:
        if self.config.martini_itp_dir:
            return Path(self.config.martini_itp_dir).resolve()
        return self.output_dir / "ff"

    def _copy_ff_files(self, ff_dir: Path) -> None:
        """Copy required Martini 3 ITP files into output_dir.

        Mirrors the cg.peptide builder's behaviour: insane writes
        ``topol.top`` with bare-name ``#include`` directives, and
        ``gmx grompp`` resolves them relative to the topology's directory.
        """
        for fname in forcefield_check.REQUIRED_MARTINI_FILES:
            src = ff_dir / fname
            dst = self.output_dir / fname
            if src.exists():
                if dst.exists():
                    continue
                shutil.copy2(src, dst)
                logger.debug("Copied %s -> %s", src, dst)
            else:
                logger.warning(
                    "Required Martini 3 file %s not found in %s; "
                    "gmx grompp will fail without it",
                    fname, ff_dir,
                )
