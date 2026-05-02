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

- Phase A (skeleton) — in progress
- Phase B (AMBER end-to-end) — pending
- Phase C (CHARMM36 backend) — pending
- Phase D (tutorial doc + integration test) — pending
