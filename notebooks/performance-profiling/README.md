# Performance Profiling

Scripts for benchmarking CAPRICHO's data curation pipeline. These were used to generate the timing and memory claims reported in the manuscript.

## Scripts

- **`benchmark_case_studies.sh`** — Runs each case study command 3 times and reports mean/std wall-clock time and peak memory (RSS). Uses `measure_command.py` internally.
- **`measure_command.py`** — Wraps any command to measure wall time and peak resident set size via `resource.getrusage(RUSAGE_CHILDREN)`.

## Usage

```bash
chmod +x benchmark_case_studies.sh measure_command.py
./benchmark_case_studies.sh
```

Requires:
- CAPRICHO installed in a `capricho` micromamba environment
- ChEMBL database downloaded locally (`capricho download`)

Results are written to `benchmark_results.csv` and summarized at the end.

## Results

Measured on a MacBook Air (Apple M4, 24 GB RAM) using CAPRICHO 1.0.0 with ChEMBL 36 (local SQL backend).

| Case Study | Targets/Assays | Bioactivity Type | Mean Time | Std Time | Peak Memory |
|---|---|---|---|---|---|
| Max-curation (Ki) | 26 targets | Ki | 3.3 min | 0.2 min | 5.1 GB |
| Max-curation (IC50) | 39 targets | IC50 | 8.1 min | 0.2 min | 5.4 GB |
| CYP inhibition | 14 targets | Mixed | 3.3 min | 0.02 min | 2.8 GB |
| Caco-2 (A→B + B→A) | 141 assays | Permeability | < 15 s | — | < 300 MB |

All runs used `n = 3` repetitions.
