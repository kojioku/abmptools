# FAQ

## How do I install ABMPTools?

```bash
pip install --user .
```

This installs the `abmptools` package and attempts to compile the Fortran shared library. If `gfortran` is not available, the install still succeeds — the Fortran acceleration is optional.

Ensure `numpy` and `pandas` are installed in your environment (they are not auto-installed; see [dependencies.md](dependencies.md)).

## What Python version is required?

Python 3.6+ is assumed based on code features. Python 3.8+ is recommended. No explicit version constraint is declared in `setup.py`.

## How do I run a specific tool?

All tools are invoked as Python modules:

```bash
python -m abmptools.<module_name> [options]
```

Examples:
```bash
python -m abmptools.generateajf -i protein.pdb --method MP2 -basis 6-31G*
python -m abmptools.getifiepieda --frag 10 -d 8.0 -i calc.log
python -m abmptools.log2cpf -i calc.log -o output.cpf
python -m abmptools.convertcpf -i large.cpf -f 1-100 -v 23
```

Use `-h` with any module to see available options:
```bash
python -m abmptools.getifiepieda -h
```

## What is CPF and which versions are supported?

CPF (Coordinate Property File) is ABINIT-MP's output format containing atomic coordinates, fragment definitions, monomer/dimer energies, and PIEDA components. Supported versions:

| Version | Notes |
|---------|-------|
| 4.201 | Legacy format |
| 7.0 | MIZUHO variant |
| 10 | Default since v1.8.0 |
| 23 | Latest (Open 1.0 Rev.23); used for DIFIE output |

`CPFManager.parse()` auto-detects the version. Use `convertcpf` to convert between versions:
```bash
python -m abmptools.convertcpf -i input.cpf -v 23
```

## What is DIFIE?

DIFIE (Dynamic IFIE) is a time-averaged CPF file generated from multiple MD trajectory snapshots. It computes the mean and standard deviation of IFIE/PIEDA values across snapshots and stores them in a single CPF file with `M-` (mean) and `S-` (std) prefixed column headers.

```bash
python -m abmptools.generate_difie -i traj-xxx.cpf -t 1 10 1 -f 1-100 -np 4
```

The output CPF can be visualized in Biostation Viewer.

## Do I need gfortran?

No. gfortran is optional and only needed to compile the Fortran shared library (`readifiepiedalib.so`) that accelerates IFIE/PIEDA reading from log files.

- **With gfortran**: Run `make` or install via `pip install --user .` to compile the library.
- **Without gfortran**: All functionality works via pure Python. Use the `-nof90` flag with `getifiepieda` to explicitly use the Python path.

Note that some features are only available via specific paths:
- MP3/MP4 extraction requires the Fortran module.
- PB-IFIE, BSSE-IFIE, monomer/dimer energies require the pure Python mode (`-nof90`).

## How do I use the Fortran-accelerated IFIE reader?

If the Fortran library is compiled (present at `abmptools/f90/bin/readifiepiedalib.so`), it is loaded automatically by `getifiepieda.py`. To disable it:

```bash
python -m abmptools.getifiepieda --frag 10 -d 8.0 -i calc.log -nof90
```

To compile manually:
```bash
make    # in repository root
```

## What are the supported ABINIT-MP versions?

| ABINIT-MP Version | Revisions | AJF Version Flag |
|--------------------|-----------|------------------|
| v1 | Rev.10–23 | `rev22`, `rev23`, etc. |
| v2 | Rev.4–8 | `v2rev4`, `v2rev8` |

Specify the version when generating AJF files:
```bash
python -m abmptools.generateajf -i protein.pdb -ajfv v2rev8
```

## How do I handle nucleic acid systems?

Nucleic acid support was added in v1.14.6. The `LOGManager` (`logmanager.py`) and `log2config` module handle:
- Pure nucleic acid systems.
- Protein-nucleic acid complexes.
- CYS disulfide bridge handling in mixed systems.

Use the standard `log2config` or `log2cpf` workflows — nucleic acid detection is automatic.

## Where are sample workflows?

The `sample/` directory contains working examples:

| Directory | Workflow | Run |
|-----------|----------|-----|
| `sample/convertcpf/` | CPF conversion/filtering | `bash run.sh` |
| `sample/generate_difie/` | DIFIE averaging (TrpCage, CS4) | `bash run.sh` |
| `sample/gbsa/` | GBSA solvation setup | `bash run.sh` |
| `sample/generateajf/` | AJF generation | `bash run.sh` |
| `sample/log2config/` | Log to config extraction | `bash run.sh` |
| `sample/log2cpf/` | Log to CPF conversion | `bash run.sh` |
| `sample/rmsd/` | RMSD analysis | `bash run.sh` |

## How do I use CPFManager programmatically?

```python
import abmptools

cpf = abmptools.CPFManager()
cpf.parse("input.cpf")

# Access data as pandas DataFrames
print(cpf.atominfo)     # Per-atom data
print(cpf.fraginfo)     # Fragment definitions
print(cpf.diminfo)      # Dimer interaction data (IFIE/PIEDA)
print(cpf.mominfo)      # Monomer energies
print(cpf.static_data)  # Summary statistics

# Write to new CPF
cpf.write("Title", "output.cpf")
```

## How do I speed up multi-sample analysis?

Use the `-np` flag to enable parallel processing:

```bash
python -m abmptools.getifiepieda --multi 10 -d 8.0 -t 1 100 1 -i template.log -np 8
python -m abmptools.generate_difie -i traj-xxx.cpf -t 1 50 1 -np 8
```

This uses Python's `multiprocessing.Pool` to read and process log/CPF files in parallel. See [parallelization.md](parallelization.md) for details.

## Amorphous builder: can I fetch 3D SDFs automatically from PubChem?

Yes, from `abmptools 1.15.3+`. Use `--pubchem_cid` / `--pubchem_name` on `build_amorphous.py`, or import the helpers directly:

```python
from abmptools.amorphous import fetch_3d_sdf, fetch_smiles, PubChemNo3DError
sdf = fetch_3d_sdf(3825)                    # by CID
sdf = fetch_3d_sdf("ketoprofen", by="name") # by name
```

If the compound has no 3D conformer on PubChem (common for polymers, salts, metal complexes), `PubChemNo3DError` is raised. In that case, fall back to `fetch_smiles(...)` and let OpenFF generate a conformer from SMILES, or prepare a 3D SDF locally (RDKit/OpenBabel + AddHs/Embed).

Network access to `https://pubchem.ncbi.nlm.nih.gov` is required. No extra Python dependency is needed — `urllib` is used.

## Amorphous builder: Packmol fails with `Fortran runtime error: Illegal seek`?

This happens with `packmol 21.2.1` from conda-forge when its input is piped through a non-seekable `stdin`. Fixed in `abmptools/amorphous/packing.py` (1.15.1+) by passing the input file via a real file descriptor and resolving PDB / output paths to absolute paths. If you are on an older abmptools, upgrade (`pip install -U abmptools`) or manually patch `run_packmol()` to use `stdin=open(inp_path, "rb")`.

## Amorphous builder: `ModuleNotFoundError: No module named 'pkg_resources'`?

Triggered by `openff.amber_ff_ports` (used indirectly via `openff-toolkit`) on setuptools 82+, where `pkg_resources` has been spun out. Pin setuptools below 81 until the upstream migration lands:

```bash
pip install "setuptools<81"
```

This applies to any environment that mixes OpenFF (Amber-style ports) with a recent setuptools.

## Amorphous MD on WSL2: GROMACS does not detect my NVIDIA GPU?

conda-forge's `gromacs` is typically built with OpenCL for GPU support, but WSL2 does not ship a NVIDIA OpenCL ICD (`libcuda.so` is present, but `/etc/OpenCL/vendors/nvidia.icd` is not). The mdrun log will say `GPU detection failed: No valid OpenCL driver found`.

Options:

1. **Run on CPU** — for small systems (~a few thousand atoms) 8 CPU cores are often enough (our ketoprofen ×50 / 1.3 ns test finished in ~10 minutes).
2. **Install a CUDA-enabled GROMACS** via conda-forge:
   ```bash
   micromamba create -n gmxcuda -c conda-forge gromacs=2025.4=mpi_openmpi_cuda_*
   ```
   Requires a working CUDA driver inside WSL2 (already the case with recent `libcuda.so`).

## CG (Martini 3): How do I install the Martini 3 force field files?

`abmptools.cg.peptide` (1.18.0+) and `abmptools.cg.membrane` (1.19.0+) do **not** bundle the Martini 3 ITPs because cgmartini.nl distributes them under unspecified terms. Users download `martini_v300.zip` and unzip the required ITPs into `./ff/`:

```bash
mkdir ff
curl -L -o /tmp/m300.zip \
    https://cgmartini-library.s3.ca-central-1.amazonaws.com/1_Downloads/ff_parameters/martini3/martini_v300.zip
unzip -o -j /tmp/m300.zip \
    "martini_v300/martini_v3.0.0.itp" \
    "martini_v300/martini_v3.0.0_solvents_v1.itp" \
    "martini_v300/martini_v3.0.0_ions_v1.itp" \
    "martini_v300/martini_v3.0.0_phospholipids_v1.itp" \
    -d ff/
```

`cg.peptide` needs the first 3 ITPs; `cg.membrane` adds the phospholipid ITP. The `validate` subcommand (`python -m abmptools.cg.peptide validate ...` / `python -m abmptools.cg.membrane validate ...`) checks all required files are present. Please cite Souza et al. 2021 *Nat. Methods* in any publication using the Martini 3 force field.

## CG: My cg.membrane PMF has 30-80 kJ/mol left/right asymmetry — bug?

No. This is a known limitation of **single-direction umbrella sampling**: the peptide is pulled from +z above the bilayer to -z below in one trajectory, so each window inherits a sampling history biased toward that direction. The two water-side minima should be physically equivalent for a homogeneous bilayer, but with finite per-window MD time the WHAM integration anchors PMF=0 at one end and the other end drifts.

Mitigations (not implemented in 1.19.0):

1. **Bidirectional pulling + Bennett acceptance** — pull `+z → -z` and `-z → +z` separately, combine forward/reverse forces.
2. **Replica-exchange umbrella sampling** — use `plumed` to swap windows during MD.
3. **Post-process symmetrize** — fit and force-symmetrize about z=0 (only valid for a symmetric system).

The same artifact appears in AA membrane Phase D (CHARMM36) and is documented in `docs/membrane.md` Phase C+. See `docs/tutorial_cg_membrane_us.md` §5 for details.

## CG: My pull stage hangs at step 0 with 99% CPU and no progress?

Symptom: `pull/pull.log` shows only `Step 0` and `pull.xtc` does not grow.

Root cause: pull MDP has both `pull-coord1-geometry = direction-periodic` AND a `pull-group1-pbcatom = N` line. With Martini 3 dt=20 fs, this combination puts the integrator in a deadlock loop. AA membrane (dt=2 fs) does not see this pathology.

Fix in 1.19.0+: `MembraneCGBuilder._stage6_mdps` only injects pbcatom into **window MDPs** (which use `direction` geometry, dynamic-box compatible). The pull MDP relies on `direction-periodic`'s internal periodic image handling and does not need pbcatom. If you build with an older builder version, manually remove `pull-group1-pbcatom` from `pull/pull.mdp` and re-grompp.

## CG: `gmx grompp` on a window MDP fails with "Pull group 1 is larger than half the box"?

Symptom (typically at the bilayer-center window, e.g. `win_018/window.mdp` for 31-window protocol):

```
ERROR 1 [file windows/win_018/window.mdp]:
  When the maximum distance from a pull group reference atom to other
  atoms in the group is larger than 0.5 times half the box size a
  centrally placed atom should be chosen as pbcatom. Pull group 1 is
  larger than that and does not have a specific atom selected as
  reference atom.
```

Root cause: the bilayer xy-extent exceeds half the post-NPT compressed shortest box vector, but the window MDP lacks an explicit `pull-group1-pbcatom` directive.

Fix in 1.19.0+: `MembraneCGBuilder._stage6_mdps` calls `find_pbc_center_atom` on `system_ions.gro` + `index.ndx` to locate the atom closest to the bilayer COM, then injects `pull-group1-pbcatom = <atom>` into every window MDP. If you build with an older builder version, manually compute it (`python -c "from abmptools.cg.membrane.pulling import find_pbc_center_atom; print(find_pbc_center_atom(gro_path='system_ions.gro', ndx_path='index.ndx', group_name='Bilayer'))"`) and append the line to each `windows/win_NNN/window.mdp`.

## What is the `tips/` directory?

The `tips/` directory contains auxiliary scripts for MD post-processing workflows with AMBER, GROMACS, and NAMD. These are standalone helper scripts, not part of the core `abmptools` package.
