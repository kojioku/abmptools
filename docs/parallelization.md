# Parallelization

## Overview

ABMPTools uses Python's `multiprocessing.Pool` for embarrassingly parallel workloads — primarily reading multiple log/CPF files and processing independent trajectory frames. Seven modules implement parallel execution.

## Modules Using `multiprocessing.Pool`

| Module | Parallelized Operation | Typical Use Case |
|--------|----------------------|------------------|
| `getifiepieda.py` | Reading multiple log files in multi-sample mode | Time-series IFIE extraction across trajectory |
| `anlfmo.py` | Log parsing and fragment interaction analysis | Multi-sample aggregation |
| `generate_difie.py` | Loading CPF files for DIFIE averaging | DIFIE from many trajectory snapshots |
| `molcalc.py` | Distance matrix / particle interaction calculations | Large-system distance computation |
| `abinit_io.py` | Fragment setup operations | Fragment assignment for large systems |
| `mol_io.py` | Coordinate transformations | Batch coordinate processing |
| `udfrm_io.py` | UDF frame conversion (`run_convert()`) | Converting multiple MD frames to PDB |

**Assumption**: `udfcreate.py` also appears to use multiprocessing based on imports, but the primary parallelization targets are the modules listed above.

## `-np` CLI Flag

Several CLI tools expose the `-np` / `--pynp` flag to control the number of parallel workers:

```bash
# IFIE extraction with 8 parallel workers
python -m abmptools.getifiepieda --multi 10 -d 8.0 -t 1 100 1 -i template.log -np 8

# DIFIE generation with 5 workers
python -m abmptools.generate_difie -i traj-xxx.cpf -t 1 50 1 -np 5

# Log to config conversion with 4 workers
python -m abmptools.log2config -i calculation.log -np 4
```

If `-np` is not specified, modules typically default to serial execution (1 process).

## Parallelization Patterns

### Pattern: Parallel File Reading

The most common pattern reads N independent files in parallel, then aggregates results:

```python
from multiprocessing import Pool

def read_one_file(args):
    """Worker function: parse a single log/CPF file."""
    filename, options = args
    # ... parse and return results ...
    return result

# Main process
with Pool(processes=np) as pool:
    results = pool.map(read_one_file, file_args_list)

# Aggregate results (e.g., concatenate DataFrames, compute averages)
```

This pattern is used in:
- `getifiepieda.py` — Parallel log reading for `--multi` and `--tfmatrix` modes.
- `generate_difie.py` — Parallel CPF loading for DIFIE computation.
- `anlfmo.py` — Parallel log analysis.

### Pattern: Parallel Frame Conversion

Used for MD trajectory processing:

```python
with Pool(processes=np) as pool:
    pool.map(convert_single_frame, frame_indices)
```

This pattern is used in:
- `udfrm_io.py` — Converting UDF records to PDB files in parallel.

### Pattern: Parallel Computation

Used for distance matrix and coordinate operations:

```python
with Pool(processes=np) as pool:
    distances = pool.map(calc_distance_chunk, atom_chunks)
```

This pattern is used in:
- `molcalc.py` — Distance matrix calculation for large systems.

## Performance Considerations

1. **I/O bound workloads**: Multi-sample log/CPF reading is typically I/O bound. The speedup from `-np` depends on storage throughput (spinning disk vs. SSD vs. network filesystem).

2. **Process startup overhead**: `multiprocessing.Pool` creates worker processes via `fork()`. For small numbers of files (< 5), the overhead may negate the benefit.

3. **Memory usage**: Each worker process loads its own copy of the data. For large CPF files, memory usage scales linearly with `-np`. Monitor memory when processing many large files.

4. **GIL independence**: Since `multiprocessing` uses separate processes (not threads), the Python GIL does not limit parallelism. Each worker has its own interpreter.

5. **Recommended `-np` values**:
   - For file-heavy workloads (many log/CPF files): match the number of CPU cores.
   - For memory-heavy workloads (large CPF files): keep `-np` modest (4–8) to avoid memory pressure.
   - On shared systems: respect other users' resource allocations.

6. **Fortran library interaction**: The Fortran shared library (`readifiepiedalib.so`) is loaded independently in each worker process, so parallel Fortran-accelerated reading works correctly.
