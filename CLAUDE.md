# CLAUDE.md

## Project Overview

**dynamic_multiplex** is a dual-language (Python and R) package for multiplex network analysis with customizable interlayer ties. It is a community detection algorithm that identifies temporal communities by only connecting nodes in adjacent periods, with explicit control over which layers influence which layers via a `layer_links` argument.

## Repository Structure

- `python_code/` — Python package (`src/dynamic_multiplex/`)
- `r_code/` — R package (`R/`)
- `scripts/` — Benchmark scripts comparing against competing approaches
- `manuscript/` — Manuscript pipeline (synthetic experiments + empirical application)
- `.github/workflows/` — CI/CD (R package testing on macOS)

Both R and Python implementations are functionally equivalent with near 1:1 correspondence.

## Build & Test

### Python
```bash
pip install -e "./python_code[dev,louvain]"
pytest python_code/tests
```

### R
```bash
Rscript -e "remotes::install_local('r_code')"
Rscript -e "testthat::test_dir('r_code/tests/testthat')"
```

## Benchmarks

```bash
bash scripts/run_benchmarks.sh          # both R and Python
bash scripts/run_benchmarks.sh python   # Python only
bash scripts/run_benchmarks.sh r        # R only
```

## Manuscript Pipeline

```bash
bash manuscript/run_pipeline.sh             # full pipeline
bash manuscript/run_pipeline.sh synthetic   # synthetic experiments only
bash manuscript/run_pipeline.sh empirical   # empirical (peacesciencer) only
```

This legacy synthetic pipeline writes to `manuscript/output/`. **Do not treat
`manuscript/output/` as the paper's data** — see the Data & Analysis Guardrails
below. The authoritative results for the paper live in
`replication/extended/output/`.

## Key Dependencies

- **Python**: networkx, pandas, numpy; optional: python-louvain, python-igraph, leidenalg
- **R**: igraph; optional: ggplot2, gganimate, ggalluvial, RColorBrewer, peacesciencer

## Coding Conventions

- Snake_case for all function and variable names in both languages
- Python uses type annotations and dataclasses
- R uses roxygen2 documentation and S3 classes
- Both languages support Louvain and Leiden community detection algorithms
- All fitting functions return a dict/list with: algorithm, layer_communities, layer_links, interlayer_ties, directed
- Node and community indices are 1-indexed in both languages

## Data & Analysis Guardrails

These rules exist because prior work repeatedly drifted onto stale or wrong data.
They are mandatory for any figure, table, or statistic in the paper.

### A. Canonical data source
1. The paper's authoritative data is ROAR `replication/extended/output/`, most
   recent run only.
2. Off-limits as inputs: `manuscript/output/`, `paper_assets/`, any pulled cloud
   snapshot, any other vintage.
3. Never select a data file by name match — resolve to the last-run directory
   explicitly.

### B. Pre-run provenance check (enforcement)
4. Before generating any figure, table, or statistic, print the date/mtime of
   every input file and confirm it is from the last run. Show the evidence; do
   not assert it.
5. If more than one vintage exists for an input, stop and ask. Do not guess.

### C. Metric & comparison correctness
6. The simulation metric is `nmi_joint` (temporal), never per-layer NMI. Per-layer
   NMI flatters untracked methods and is not the comparison metric.
7. The only no-coupling baseline is Cross-sectional + Hungarian. Raw independent /
   cross-sectional detection is not a valid comparator on temporal tasks.

### D. Cluster safety (ROAR)
8. Never `scancel` anything except a specific job ID the user named; never operate
   on all of jfe4's jobs; confirm before any irreversible cluster action.
9. On ROAR, archive (`mv`) rather than delete, and confirm scope first.
