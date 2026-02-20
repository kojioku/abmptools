# I/O Specification

## Overview

ABMPTools reads and writes multiple file formats used in Fragment Molecular Orbital (FMO) calculations. This document specifies each format, the modules that handle it, and relevant details.

## Input Formats

### PDB (Protein Data Bank)

- **Module**: `pdb_io.py` (`readpdb()`)
- **Read modes**: `TER` (split on TER records), `resnum` (split on residue number), `rfile` (read from file metadata)
- **Usage**: Primary structural input for protein and molecular systems.
- **Notes**: Supports extended atom numbering (5+ digits). Fixed-column format parsing.

### ABINIT-MP LOG

- **Module**: `logmanager.py` (`LOGManager.parse()`), `anlfmo.py` (`readlog()`)
- **Content**: Calculation output including energies, fragment info, IFIE/PIEDA tables.
- **Versions supported**: ABINIT-MP v1 (Rev.10–23), v2 (Rev.4–8).
- **Compressed**: `.gz` files supported transparently.
- **Special cases**: `Nprint=0` logs, nucleic acid systems, CYS disulfide bridges.

### CPF (Coordinate Property File)

- **Module**: `cpfmanager.py` (`CPFManager.parse()`)
- **Format**: Mixed binary/text format with sections for atoms, fragments, dimers, monomers, conditions.
- **Versions**: See [CPF Format Versions](#cpf-format-versions) below.
- **Compressed**: `.gz` files supported transparently.

### UDF (OCTA COGNAC)

- **Module**: `udf_io.py`, `udfrm_io.py`
- **Content**: Binary molecular dynamics trajectory data from OCTA COGNAC simulations.
- **Dependency**: Requires optional `UDFManager` library.
- **Usage**: Extract frames for FMO calculation setup.

### CIF (Crystallographic Information File)

- **Module**: `readcif.py`
- **Functions**: `getcartesiancellvec()`, `getcartesianmol()`, `intocell()`
- **Usage**: Convert crystal structures to Cartesian coordinates for FMO input.

### XYZ

- **Module**: `mol_io.py` (`read_xyz()`)
- **Format**: Standard XYZ molecular coordinate format (atom count, comment, atom lines).

### Config Dictionary (`segment_data.dat`)

- **Module**: `ajf2config.py` (writer), `pdb2fmo.py` (reader)
- **Format**: Python dictionary serialized to file, containing fragment definitions for non-protein systems.
- **Content**: Molecule name → list of fragment atom ranges.

### Parameter File (`input_param`)

- **Module**: `pdb2fmo.py`
- **Format**: Python-evaluable dictionary with FMO setup parameters.
- **Key fields**:

```python
{
    'cutmode': 'around',        # sphere | cube | around | neutral | none
    'solutes': [0, 1],          # solute molecule indices
    'criteria': 8.0,            # distance cutoff (Å)
    'ionmode': 'remain',        # remain | remove
    'ionname': ['NA'],          # ion residue names
    'molname': ['000', '001'],  # molecule type names
    'getmode': 'rfile',         # PDB read mode
}
```

## Output Formats

### AJF (ABINIT-MP Input File)

- **Module**: `generateajf.py`, `addsolvfrag.py`
- **Content**: Complete ABINIT-MP input specification — coordinates, fragment definitions, method, basis set, options.
- **Version-aware**: Output format varies by ABINIT-MP version (`rev22`, `v2rev4`, `v2rev8`, etc.).

### CPF (Output)

- **Module**: `cpfmanager.py` (`CPFManager.write()`), `log2cpf.py`, `generate_difie.py`, `convertcpf.py`
- **Format**: Same mixed binary/text format as input, with version specification.
- **DIFIE variant**: Columns prefixed with `M-` (mean) and `S-` (standard deviation).

### CSV

- **Module**: `getifiepieda.py`, `cpf2ifielist.py`, `getcharge.py`
- **Content**: IFIE/PIEDA interaction data, distance matrices, charge data.
- **See**: [CSV Output Columns](#csv-output-columns) below.

### PDB (Output)

- **Module**: `pdbmodify.py`, `udfrm_io.py`, `pdb2fmo.py`
- **Usage**: Edited/converted PDB structures.

### XYZ (Output)

- **Module**: `mol_io.py`, `udfrm_io.py`
- **Usage**: Coordinate export from UDF or internal structures.

## CPF Format Versions

| Version | Origin | Notes |
|---------|--------|-------|
| 4.201 | Legacy | Older ABINIT-MP format |
| 7.0 | MIZUHO | MIZUHO-specific variant (added v1.13.2) |
| 10 | Standard | Default version since v1.8.0 |
| 23 | Open 1.0 Rev.23 | Latest format; used for DIFIE output |

Each version differs in:
- Header format and metadata fields.
- Available charge types (MUL-HF, MUL-MP2, NPA-HF, NPA-MP2, ESP-HF, ESP-MP2).
- Dimer data columns (PIEDA components, BSSE, SCS-MP2).
- Monomer data columns (NR, HF, MP2, MP3).

`CPFManager.parse()` auto-detects the version from the file header.

## CSV Output Columns

### IFIE/PIEDA Per-Pair Output

Standard columns in single-point and multi-sample CSV files:

| Column | Description | Unit |
|--------|-------------|------|
| `I` | Reference fragment (name + ID) | — |
| `J` | Partner fragment (name + ID) | — |
| `DIST` | Minimum inter-fragment distance | Å |
| `DIMER-ES` | Dimer-ES approximation flag (0 = full, 1 = approximate) | — |
| `HF-IFIE` | Hartree-Fock IFIE | kcal/mol |
| `MP2-IFIE` | MP2 IFIE (HF + MP2 correlation) | kcal/mol |
| `PR-TYPE1` | PR-MP2 corrected IFIE | kcal/mol |
| `GRIMME` | SCS-MP2 (Grimme) | kcal/mol |
| `JUNG` | SCS-MP2 (Jung) | kcal/mol |
| `HILL` | SCS-MP2 (Hill) | kcal/mol |
| `ES` | Electrostatic component | kcal/mol |
| `EX` | Exchange repulsion component | kcal/mol |
| `CT-mix` | Charge transfer + mix component | kcal/mol |
| `DI(MP2)` | Dispersion-like (MP2 correlation) component | kcal/mol |
| `q(I=>J)` | Charge transfer amount | e |
| `TIMES` | Timestep (multi-sample only) | — |

### ffmatrix Output Files

Fragment-fragment matrix mode produces one CSV per PIEDA component:

| File | Content |
|------|---------|
| `Distance-ffmatrix.csv` | Minimum distance (Å) |
| `ES-ffmatrix.csv` | Electrostatic |
| `EX-ffmatrix.csv` | Exchange |
| `CT-ffmatrix.csv` | Charge transfer |
| `HF-ffmatrix.csv` | HF-IFIE |
| `MP2corr-ffmatrix.csv` | MP2 correlation only |
| `MP2total-ffmatrix.csv` | HF + MP2 correlation |
| `PRMP2corr-ffmatrix.csv` | PR-MP2 correction |
| `PRMP2total-ffmatrix.csv` | HF + PR-MP2 total |

Matrix format: rows and columns labeled as `ResidueName(FragID)`.

### cpf2ifielist Output

| Column | Description | Unit |
|--------|-------------|------|
| `fragi` | Fragment I index | — |
| `fragj` | Fragment J index | — |
| `min-dist` | Minimum distance | Å |
| `ES` | Electrostatic | kcal/mol |
| `DI` | Dispersion | kcal/mol |
| `PR-MP2` | PR-MP2 | kcal/mol |
| `SCS-MP2(Grimme)` | SCS-MP2 | kcal/mol |
| `EX` | Exchange | kcal/mol |
| `CT` | Charge transfer | kcal/mol |
| `DQ` | Higher-order term | kcal/mol |
| `MP2-Total` | Total MP2 IFIE | kcal/mol |

## Unit Conventions

| Quantity | Unit | Conversion |
|----------|------|------------|
| Energy (IFIE/PIEDA) | kcal/mol | × 627.5095 from Hartree |
| Distance | Å (Angstrom) | × 0.529177 from Bohr |
| Charge transfer | e (electron) | — |
| Coordinates (internal) | Å | — |

## Compressed File Support

Both `CPFManager` and `LOGManager` transparently handle `.gz` compressed files:
- Detection: filename ends with `.gz`.
- Implementation: `gzip.open()` in read mode.
- All CLI tools that accept LOG or CPF files support `.gz` inputs.
