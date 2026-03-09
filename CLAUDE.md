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

Outputs go to `manuscript/output/`.

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
