# `insulin_smoke` — 1 ns insulin × caprate × taurocholate smoke

Recreation of Hossain et al. 2023 (*Nanoscale* 15, 19180) **system 4**
(2 × insulin + 32 × caprate + 2 × taurocholate, 10 nm cubic, fasted
state) using commercial-permissive force fields:

- Insulin (PDB **2G4M**): A-chain (21 aa) + B-chain (30 aa), 3 disulfide
  bridges. Free download from RCSB; not redistributed by abmptools.
- AMBER ff14SB protein force field
- GAFF2 + AM1-BCC for sodium caprate and taurocholate
- TIP3P water + Joung-Cheatham Na/Cl

## Composition (Hossain 2023 system 4)

| Species | Count | Notes |
|---|---|---|
| Insulin | 2 (× 2 chains each = 4) | PDB 2G4M, 3 disulfides per copy |
| Na-caprate (neutral) | 16 | GAFF2 / AM1-BCC |
| Na-caprate (charged) | 16 | GAFF2 / AM1-BCC |
| Taurocholate | 2 | charge -1, GAFF2 / AM1-BCC |
| Box | 10 nm cubic | |
| Salt | 0.15 M NaCl | Joung-Cheatham |

Expected size: ~99,000 atoms.

## Setup

```bash
# 1. Download insulin PDB 2G4M (license-free per RCSB terms)
cd sample/formulation/insulin_smoke
curl -O https://files.rcsb.org/download/2G4M.pdb
mv 2G4M.pdb insulin_2G4M.pdb

# 2. Validate
micromamba run -n abmptoolsenv python -m abmptools.formulation validate \
    --config config.json

# 3. Build (~10-20 min depending on acpype + tleap solvatebox)
micromamba run -n abmptoolsenv python -m abmptools.formulation build \
    --config config.json --output-dir /tmp/insulin_smoke

# 4. Run MD (em + nvt + npt + 1 ns prod). On a GPU this is ~10-20 min;
#    on CPU only it is ~6-8 h. Set GMX/NT to taste.
cd /tmp/insulin_smoke
NT=8 MDRUN_OPTS="-nb gpu" bash run.sh

# 5. Analyze
micromamba run -n abmptoolsenv python -m pip install MDAnalysis networkx
micromamba run -n abmptoolsenv python -m abmptools.formulation analyze \
    --traj /tmp/insulin_smoke/prod/prod.xtc \
    --top /tmp/insulin_smoke/system.top \
    --out /tmp/insulin_smoke/analysis \
    --enhancer-resnames CPRN,CPRC \
    --bile-salt-resnames TCH
```

## Force-field caveat

This build is **not** a direct reproduction of Hossain 2023 absolute
PMF values. The paper used CHARMM36m + CGenFF 1.0.0 (academic
license); we use ff14SB + GAFF2 + AM1-BCC (commercial-permissive).
Qualitative trends (aggregate formation, peptide release order
mediated by enhancers vs. bile salt) are expected to reproduce within
~5-15 kJ/mol; absolute PMF values are not directly comparable.

To run a longer, statistically meaningful production simulation,
bump `production.nsteps` from `500000` (1 ns) to `250000000` (500 ns)
in `config.json`.
