# CLAUDE.md

## Project Overview

**dynamic_multiplex** is a dual-language (Python and R) package for multiplex network analysis with customizable interlayer ties. It enables community detection across multiple network layers with explicit control over which layers influence each other via a `layer_links` argument.

## Repository Structure

- `python_code/` — Python package (`src/dynamic_multiplex/`)
- `r_code/` — R package (`R/`)
- `.github/workflows/` — CI/CD (R package testing on macOS)

Both implementations are functionally equivalent with near 1:1 correspondence.

## Build & Test

### Python
```bash
cd python_code
pip install -e ".[dev,louvain]"
pytest
```

### R
```bash
cd r_code
Rscript -e "devtools::test()"
# or full check:
R CMD check .
```

## Key Dependencies

- **Python**: networkx, pandas, numpy; optional: python-louvain, python-igraph, leidenalg
- **R**: igraph; optional: ggplot2, gganimate, ggalluvial, RColorBrewer

## Coding Conventions

- Snake_case for all function and variable names in both languages
- Python uses type annotations and dataclasses
- R uses roxygen2 documentation and S3 classes
- Both languages support Louvain and Leiden community detection algorithms
- All fitting functions return a dict/list with: algorithm, layer_communities, layer_links, interlayer_ties, directed
