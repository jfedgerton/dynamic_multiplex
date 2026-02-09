# dynamic_multiplex

`dynamic_multiplex` is an R package for multiplex community modeling with customizable interlayer ties.

## Why this package

Most multilayer workflows assume every network layer connects to every other layer. This package adds explicit control over *which layers influence which layers* via a `layer_links` argument (`from`, `to`, `weight`).

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

## Installation

Development install:

```r
# install.packages("remotes")
remotes::install_local(".")
```

GitHub install (recommended for collaborators):

```r
# install.packages("remotes")
remotes::install_github("jarededgerton/bayesnet", subdir = "r_code")
```

## Quick example

```r
library(dynamic_multiplex)

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

## Planned next steps

- Add S3 print/summary/plot methods for fit objects.
- Add benchmarking utilities against other multilayer tooling.
- Add package tests and validation diagnostics.
