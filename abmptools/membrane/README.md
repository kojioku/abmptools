# abmptools.membrane

Lipid-bilayer umbrella-sampling (US) system builder.

## Scope

Build GROMACS-ready systems for peptide–membrane PMF calculations:

```
peptide (in water)
    ↓ pulling along z
peptide (in bilayer)
    ↓ window-by-window US
PMF(z)
```

Sister package to `abmptools.amorphous` (3D random boxes).
Both share `abmptools.core.system_model` for the simulation-parameter IR.

## Backends

Two parameterisation routes, selected by `MembraneConfig.backend`:

| Backend | Force field | License |
|---|---|---|
| `"amber"`    | ff19SB + Lipid21 + TIP3P (Joung-Cheatham ions) | AmberTools, commercial OK |
| `"charmm36"` | CHARMM36m + CHARMM36 lipids + CHARMM-modified TIP3P | parameters commercial OK (see below) |

## Commercial-use rule (strict)

This subpackage is intended for **company use without separate licensing**.
Therefore the following are **forbidden** by design:

- ❌ **CGenFF web server** output (https://cgenff.umaryland.edu) —
  commercial use requires a Silcsbio subscription. The pipeline does
  **not** call `cgenff` / `match` / Silcsbio binaries.
- ❌ **CHARMM-GUI** automated system / topology generation —
  academic-only without paid subscription. The pipeline does **not**
  produce or consume CHARMM-GUI output.

The following are explicitly **allowed**:

- ✅ **AmberTools** (tleap, antechamber, parmchk2, packmol-memgen) —
  free, redistributable, commercial OK.
- ✅ **CHARMM36 force-field parameter values** when sourced from
  MacKerell lab's free distribution or Klauda lab's GROMACS port.
  MacKerell lab's stated license: "free of charge to academic and
  industrial researchers."
- ✅ **Packmol** (Martínez et al.) — free.
- ✅ **GROMACS** — LGPL.
- ✅ **parmed** (Shirts et al.) — LGPL.

If a future feature would require a forbidden tool, it must use an
allowed alternative (e.g., parametrise novel small molecules with
GAFF2/antechamber instead of CGenFF).

## Module map

| Module | Role |
|---|---|
| `models.py`              | Dataclasses (MembraneConfig, LipidSpec, …) |
| `builder.py`             | Orchestrator (MembraneUSBuilder) |
| `bilayer.py`             | packmol-memgen wrapper for bilayer + peptide |
| `parameterize_amber.py`  | tleap → parmed → GROMACS top |
| `parameterize_charmm.py` | CHARMM36 ff staging → GROMACS top |
| `mdp_us_protocol.py`     | semiisotropic NPT + pull-code MDP writers |
| `pulling.py`             | reaction-coordinate generation |
| `umbrella.py`            | per-window MDP + run script |
| `pmf.py`                 | gmx wham / pymbar wrapper |

## Status

- Phase A (skeleton) — done (commit `a0d2bb9`)
- Phase B (AMBER end-to-end) — in progress
  - B-1 (bilayer.py) — done (smoke verified)
  - B-2 (parameterize_amber.py) — pending
  - B-3..B-7 — pending
- Phase C (CHARMM36 backend) — pending
- Phase D (tutorial doc + integration test) — pending

## Environment requirements

- `AMBERHOME` must be resolvable (auto-detected from
  `packmol-memgen` binary location for conda installs; otherwise set
  the env var manually).
- `PATH` must include `tleap`, `packmol-memgen`, `gmx` (or pass
  full paths via `MembraneConfig.tleap_path` / `.packmol_memgen_path` /
  `.gmx_path`).

### packmol-memgen + NumPy ≥ 1.24 patch

`packmol-memgen` 2023.2.24 ships an internal copy of `pdbremix`
(`packmol_memgen/lib/pdbremix/v3numpy.py`) that uses the deprecated
`np.float` alias, which was removed in NumPy 1.24. If you see

    AttributeError: module 'numpy' has no attribute 'float'

apply this 2-line patch (idempotent — safe to re-run):

```bash
sed -i 's/np\.zeros(3, dtype=np\.float)/np.zeros(3, dtype=float)/' \
       <env>/lib/python3.10/site-packages/packmol_memgen/lib/pdbremix/v3numpy.py
sed -i 's/np\.array(args, dtype=np\.float, copy=True)/np.array(args, dtype=float, copy=True)/' \
       <env>/lib/python3.10/site-packages/packmol_memgen/lib/pdbremix/v3numpy.py
```

`<env>` is your conda env root (e.g.
`~/.local/share/mamba/envs/abmptoolsenv`). The fix only changes
`np.float` → `float` (Python builtin), which is the recommended
replacement per NumPy's deprecation note.
