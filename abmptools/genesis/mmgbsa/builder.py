# -*- coding: utf-8 -*-
"""
abmptools.genesis.mmgbsa.builder
--------------------------------
End-to-end orchestrator for the GENESIS MM/GBSA pipeline.

4-stage pipeline (:meth:`MMGBSAOrchestrator.run`):

1. ``divide``        -- :mod:`pdb_splitter` per target
                        (-> ``<basename>/<basename>_{ligand,receptor}_<resno>.pdb``)
2. ``parameterize``  -- :mod:`ligand_parameterize` + :mod:`system_builder`
                        per target (-> ``<basename>/{complex,ligand,receptor}.prmtop+inpcrd``)
3. ``run_gbsa``      -- :mod:`gbsa_runner` per target
                        (-> ``<basename>/{complex,ligand,receptor}.{inp,log,dcd,rst}``)
4. ``analyze``       -- :mod:`analysis` over all targets
                        (-> ``analysis_results.csv``, ``dg_bind_plot.png``)

Stages can be invoked individually for re-runs and debugging:

    orch = MMGBSAOrchestrator(cfg)
    orch.divide()
    orch.parameterize()
    orch.run_gbsa()
    result = orch.analyze()

or all at once via :meth:`run`. Per-target failures are recorded in
``result["failures"]`` and (depending on ``cfg.fail_fast``) either halt
the pipeline or are skipped.
"""
from __future__ import annotations

import json
import logging
from dataclasses import asdict
from pathlib import Path
from typing import Any, Dict, List

from . import (
    analysis,
    gbsa_runner,
    inp_writer,
    ligand_parameterize,
    pdb_splitter,
    system_builder,
)
from ._subprocess import CommandError, ensure_dir
from .models import MMGBSABuildConfig, TargetSpec

logger = logging.getLogger(__name__)


class MMGBSAOrchestrator:
    """Run the 4-stage MM/GBSA pipeline from a :class:`MMGBSABuildConfig`."""

    def __init__(self, config: MMGBSABuildConfig):
        self.config = config
        self.output_dir = ensure_dir(Path(config.output_dir).resolve())
        self.input_dump_dir = ensure_dir(self.output_dir / "input")

        # Per-target intermediate state.
        self._splits: Dict[str, Dict[str, Path]] = {}
        self._params: Dict[str, system_builder.SystemBuildResult] = {}
        self._runs: Dict[str, gbsa_runner.GBSARunResult] = {}
        self._failures: Dict[str, str] = {}

        # If folder mode (input_dir set, targets empty), enumerate now.
        self._resolve_folder_mode_targets()

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def run(self) -> Dict[str, Any]:
        """Run all 4 stages in sequence. Returns the result dict."""
        logger.info(
            "=== MMGBSAOrchestrator.run start (project=%s, %d targets) ===",
            self.config.project_name, len(self.config.targets),
        )
        self._dump_config()
        self.divide()
        self.parameterize()
        self.run_gbsa()
        result = self.analyze()
        logger.info("=== MMGBSAOrchestrator.run complete ===")
        return result

    def divide(self) -> Dict[str, Dict[str, Path]]:
        """Stage 1: split each target's PDB into receptor + ligand."""
        logger.info("=== Stage 1/4: PDB splitter ===")
        input_dir = Path(self.config.input_dir) if self.config.input_dir else None
        for target in self.config.targets:
            tname = self._target_name(target)
            try:
                receptor, ligand = pdb_splitter.split_target(
                    target=target,
                    output_dir=self.output_dir,
                    input_dir=input_dir,
                )
                self._splits[tname] = {
                    "receptor": receptor,
                    "ligand": ligand,
                    "target_dir": receptor.parent,
                }
                logger.info(
                    "  %s -> receptor=%s ligand=%s",
                    tname, receptor.name, ligand.name,
                )
            except Exception as exc:
                self._record_failure(tname, "divide", exc)
                if self.config.fail_fast:
                    raise
        # Persist dirnames.in for POC compatibility.
        self._write_dirnames_in()
        return self._splits

    def parameterize(self) -> Dict[str, system_builder.SystemBuildResult]:
        """Stage 2: acpype + tleap for the 3 systems per target."""
        logger.info("=== Stage 2/4: parameterize (acpype + tleap) ===")
        for target in self.config.targets:
            tname = self._target_name(target)
            if tname not in self._splits:
                # divide stage failed previously; honour fail_fast policy.
                continue
            split = self._splits[tname]
            target_dir = split["target_dir"]
            try:
                acp = ligand_parameterize.run_acpype(
                    ligand_pdb=split["ligand"],
                    workdir=target_dir,
                    config=self.config.ligand,
                    acpype_path=self.config.acpype_path,
                )
                params = system_builder.build_three_systems(
                    receptor_pdb=split["receptor"],
                    frcmod=acp.frcmod,
                    mol2=acp.mol2,
                    workdir=target_dir,
                    ff=self.config.force_field,
                    tleap_path=self.config.tleap_path,
                )
                self._params[tname] = params
                logger.info("  %s -> 3 systems parameterized", tname)
            except Exception as exc:
                self._record_failure(tname, "parameterize", exc)
                if self.config.fail_fast:
                    raise
        return self._params

    def run_gbsa(self) -> Dict[str, gbsa_runner.GBSARunResult]:
        """Stage 3: atdyn GBSA single-point on the 3 systems per target."""
        logger.info("=== Stage 3/4: atdyn GBSA single-point ===")
        inps = inp_writer.render_three_inps(
            energy=self.config.energy,
            minimize=self.config.minimize,
        )
        for target in self.config.targets:
            tname = self._target_name(target)
            if tname not in self._params:
                continue
            target_dir = self._splits[tname]["target_dir"]
            try:
                result = gbsa_runner.run_gbsa_one_target(
                    target_dir=target_dir,
                    inps=inps,
                    atdyn_path=self.config.atdyn_path,
                    mpirun_path=self.config.mpirun_path,
                    mpi_processes=self.config.mpi_processes,
                )
                self._runs[tname] = result
                logger.info("  %s -> 3 GBSA logs", tname)
            except Exception as exc:
                self._record_failure(tname, "run_gbsa", exc)
                if self.config.fail_fast:
                    raise
        return self._runs

    def analyze(self) -> Dict[str, Any]:
        """Stage 4: aggregate logs + ΔG_bind CSV + bar plot."""
        logger.info("=== Stage 4/4: analysis ===")
        analysis_dir = ensure_dir(self.output_dir / "analysis")
        target_dirs = []
        for target in self.config.targets:
            tname = self._target_name(target)
            if tname not in self._runs:
                continue
            target_dirs.append(self._splits[tname]["target_dir"])

        if not target_dirs:
            logger.warning("No successful targets to analyze.")
            return self._collect_result(analyze=None)

        analyze_result = analysis.analyze(
            target_dirs=target_dirs,
            out_csv=analysis_dir / "analysis_results.csv",
            out_png=analysis_dir / "dg_bind_plot.png",
        )
        return self._collect_result(analyze=analyze_result)

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _target_name(self, target: TargetSpec) -> str:
        """Stable display name (POC convention: PDB stem)."""
        if target.name:
            return target.name
        return Path(target.pdb).stem

    def _resolve_folder_mode_targets(self) -> None:
        """If folder mode, populate ``targets`` from ``input_dir/*.pdb``."""
        if self.config.targets or not self.config.input_dir:
            return
        input_dir = Path(self.config.input_dir)
        if not input_dir.is_dir():
            raise FileNotFoundError(
                f"input_dir not found: {input_dir}"
            )
        # All TargetSpecs share the same ligand_resno/chain in folder mode;
        # those need to come from a separate config knob -- for v1.22.0 we
        # require the user to set this via shortcut helper or to populate
        # targets explicitly. CLI's "pipeline -i ... -r 201" synthesises a
        # config with explicit targets before instantiating the orchestrator,
        # so this branch is reached only for accidental folder-mode JSON
        # configs without targets, which we treat as an error.
        raise ValueError(
            f"input_dir={self.config.input_dir} is set but targets is empty. "
            "Folder-mode requires the CLI shortcut "
            "(pipeline -i <dir> -r <ligand_resno>) which synthesises "
            "TargetSpecs at parse time."
        )

    def _dump_config(self) -> None:
        """Persist the resolved config for reproducibility."""
        cfg_path = self.input_dump_dir / "config.json"
        cfg_path.write_text(
            json.dumps(asdict(self.config), indent=2, ensure_ascii=False)
        )

    def _write_dirnames_in(self) -> None:
        """POC-compatibility: ``dirnames.in`` lists per-target dirs."""
        path = self.output_dir / "dirnames.in"
        names = [self._target_name(t) for t in self.config.targets]
        path.write_text("\n".join(names) + "\n")

    def _record_failure(self, tname: str, stage: str, exc: Exception) -> None:
        msg = f"{stage}: {type(exc).__name__}: {exc}"
        self._failures[tname] = msg
        logger.error("FAILURE [%s] %s", tname, msg)

    def _collect_result(
        self,
        analyze: Any,
    ) -> Dict[str, Any]:
        """Build the result dict returned by :meth:`run` / :meth:`analyze`."""
        return {
            "output_dir": self.output_dir,
            "config_json": self.input_dump_dir / "config.json",
            "n_targets": len(self.config.targets),
            "n_succeeded": len(self._runs),
            "n_failed": len(self._failures),
            "failures": dict(self._failures),
            "splits": {
                k: {kk: str(vv) for kk, vv in v.items()}
                for k, v in self._splits.items()
            },
            "csv_path": (
                analyze.csv_path if analyze is not None else None
            ),
            "png_path": (
                analyze.png_path if analyze is not None else None
            ),
            "targets": (
                [
                    {
                        "name": t.name,
                        "complex_e": t.complex_e,
                        "complex_s": t.complex_s,
                        "ligand_e": t.ligand_e,
                        "ligand_s": t.ligand_s,
                        "receptor_e": t.receptor_e,
                        "receptor_s": t.receptor_s,
                        "dg_bind": t.dg_bind,
                    }
                    for t in analyze.targets
                ]
                if analyze is not None
                else []
            ),
        }


def synthesize_targets_from_folder(
    input_dir: Path, ligand_resno: int, chain: str = None,
) -> List[TargetSpec]:
    """Enumerate ``*.pdb`` in *input_dir* into a list of TargetSpecs.

    Helper for the CLI ``pipeline -i ... -r 201`` folder-mode shortcut.
    """
    input_dir = Path(input_dir)
    if not input_dir.is_dir():
        raise FileNotFoundError(f"input_dir not found: {input_dir}")
    pdbs = sorted(input_dir.glob("*.pdb"))
    if not pdbs:
        raise FileNotFoundError(f"No *.pdb files in {input_dir}")
    return [
        TargetSpec(pdb=p.name, ligand_resno=ligand_resno, chain=chain)
        for p in pdbs
    ]


__all__ = [
    "MMGBSAOrchestrator",
    "synthesize_targets_from_folder",
]
