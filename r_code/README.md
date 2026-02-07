# bayesnet

`bayesnet` is being repurposed into an R package for multiplex community modeling with customizable interlayer ties.

## Why this package

Most multilayer workflows assume every network layer connects to every other layer. This package adds explicit control over *which layers influence which layers* via a `layer_links` argument (`from`, `to`, `weight`).

## Current scope

This repository now keeps only the multiplex/dynamic community-detection workflow; legacy ERGM-stage functions from prior iterations were removed to keep the package focused.

## Current prototype functions

> Note: Louvain in `igraph` is undirected. When `directed = TRUE` and `algorithm = "louvain"`, layers are collapsed to undirected weighted graphs before clustering.


- `fit_multilayer_jaccard()`
  - Fits Louvain or Leiden communities on each layer (supports directed layers).
  - Builds interlayer ties between communities using weighted Jaccard similarity across selected layer pairs.
- `fit_multilayer_overlap()`
  - Fits Louvain or Leiden communities on each layer (supports directed layers).
  - Builds interlayer ties between communities using weighted overlap coefficient across selected layer pairs.
- `fit_multilayer_identity_ties()`
  - Fits Louvain or Leiden communities on each layer (supports directed layers).
  - Builds interlayer ties only for the same node across selected adjacent layers.
- `simulate_and_fit_multilayer()`
  - Simulates multiplex layers from a planted partition process and runs one of the three fitting strategies.

## Installation (development)

```r
# install.packages("remotes")
remotes::install_local("r_code")
```

## Quick example

```r
library(bayesnet)

sim <- simulate_and_fit_multilayer(
  directed = TRUE,
  n_nodes = 50,
  n_layers = 4,
  n_communities = 3,
  fit_type = "jaccard",
  algorithm = "louvain",
  seed = 123
)

head(sim$fit$interlayer_ties)

# Custom layer influence map
custom_links <- data.frame(
  from = c(1, 2),
  to = c(2, 4),
  weight = c(1, 0.6)
)

fit_overlap <- fit_multilayer_overlap(
  sim$layers,
  algorithm = "leiden",
  layer_links = custom_links,
  min_similarity = 0.1
)
```


## Simulation benchmark for manuscripts

You can benchmark whether the package recovers *interlayer community ties* better
than multinet-style all-to-all coupling assumptions.

```r
bench <- benchmark_interlayer_tie_recovery(
  n_reps = 100,
  simulation_args = list(
    n_nodes = 200,
    n_layers = 8,
    n_communities = 5,
    persistence = 0.90,
    p_in = 0.20,
    p_out = 0.02,
    directed = TRUE
  ),
  algorithm = "leiden"
)

bench$summary  # includes precision, recall, F1, elapsed_sec, F1/sec, ARI, modularity
```

A manuscript-ready script is included at:

- `inst/scripts/manuscript_simulation_benchmark.R`


## Python companion package

A Python scaffold with matching core APIs is available in `python_code/`.

```bash
pip install -e python_code
pytest -q python_code/tests
```

## Planned next steps

- Add S3 print/summary/plot methods for fit objects.
- Add benchmarking utilities against other multilayer tooling.
- Add package tests and validation diagnostics.
