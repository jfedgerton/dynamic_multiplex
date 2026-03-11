# dynamic_multiplex

`dynamic_multiplex` is an R package for multiplex community modeling with customizable interlayer ties.

## Why this package

Standard multislice approaches (Mucha et al. 2010) connect all layers to all layers, meaning community structure at distant time periods influences assignments everywhere. When applied to temporal data, this creates a pooling problem: the community configuration at 2010 affects the structure detected at 1950. `dynamic_multiplex` provides explicit control over *which layers influence which layers* via a `layer_links` argument (`from`, `to`, `weight`), defaulting to adjacent-only temporal coupling.

## Method

This is a **two-stage** approach:

1. **Per-layer community detection**: Louvain or Leiden is run independently on each layer. Community assignments within a layer are determined solely by that layer's topology.
2. **Post-hoc interlayer ties**: After communities are detected, interlayer ties are constructed between communities (or nodes, for identity ties) in selected layer pairs using similarity metrics.

The interlayer ties describe relationships between layers but **do not retroactively influence** the per-layer community assignments. This contrasts with true coupled multislice methods (e.g., Mucha et al. 2010), where a single optimization jointly determines communities across all layers. The advantage of the two-stage approach is that each layer's communities reflect only that layer's structure, avoiding the temporal pooling problem where distant layers distort local community detection.

## Functions

> Note: Louvain in `igraph` is undirected. When `directed = TRUE` and `algorithm = "louvain"`, layers are collapsed to undirected weighted graphs before clustering.

- `fit_multilayer_jaccard()`
  - Detects communities independently on each layer (stage 1, supports directed layers).
  - Builds interlayer ties between communities using Jaccard similarity across selected layer pairs (stage 2).
- `fit_multilayer_overlap()`
  - Detects communities independently on each layer (stage 1, supports directed layers).
  - Builds interlayer ties between communities using overlap coefficient across selected layer pairs (stage 2).
- `fit_multilayer_weighted_jaccard()`
  - Detects communities independently on each layer (stage 1, supports directed layers).
  - Builds interlayer ties using node-strength weighted Jaccard similarity across selected layer pairs (stage 2).
- `fit_multilayer_weighted_overlap()`
  - Detects communities independently on each layer (stage 1, supports directed layers).
  - Builds interlayer ties using node-strength weighted overlap coefficient across selected layer pairs (stage 2).
- `fit_multilayer_identity_ties()`
  - Detects communities independently on each layer (stage 1, supports directed layers).
  - Builds node-level interlayer ties linking the same node across selected adjacent layers (stage 2).
- `simulate_and_fit_multilayer()`
  - Simulates multiplex layers from a planted partition process and runs one of the five fitting strategies.

All `fit_multilayer_*` functions accept a `seed` parameter for reproducible community detection. All layers must share the same node set.

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
  min_similarity = 0.1,
  add_self_loops = TRUE,
  self_loop_multiplier = 1
)
```

## Plotting options

```r
# 1) Static series of network panels
plot_multilayer_series(sim$layers, fit = sim$fit, directed = TRUE, palette = "Dark2")

# 2) GIF animation with colorblind-friendly community colors (`Set2`/`Dark2` from RColorBrewer)
animate_multilayer_gif(
  sim$layers,
  fit = sim$fit,
  output_file = "multilayer_communities.gif",
  directed = TRUE,
  fps = 2,
  palette = "Set2"
)

# 3) Alluvial plot for community flow across temporal layers
plot_multilayer_alluvial(sim$fit, max_nodes = 100, palette = "Dark2")
```
