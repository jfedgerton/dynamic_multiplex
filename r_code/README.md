# dynamicmultiplex

`dynamicmultiplex` is an R package for multiplex community modeling with customizable interlayer ties.

## Why this package

Standard multislice approaches (Mucha et al. 2010) connect all layers to all layers, meaning community structure at distant time periods influences assignments everywhere. When applied to temporal data, this creates a pooling problem: the community configuration at 2010 affects the structure detected at 1950. `dynamicmultiplex` provides explicit control over *which layers influence which layers* via a `layer_links` argument (`from`, `to`, `weight`), defaulting to adjacent-only temporal coupling.

## Functions

> Note: `igraph`'s Louvain and Leiden are undirected. When `directed = TRUE`, layers are collapsed to undirected weighted graphs (reciprocal arcs summed) before clustering; the Python package keeps edge direction for Leiden.

- `fit_multilayer_jaccard()`
  - Fits Leiden (default) or Louvain communities on each layer, links communities in the selected layer pairs by the Jaccard index of their members, and groups the linked communities into cross-layer meta-communities with a second community detection on the community graph.
- `fit_multilayer_overlap()`
  - Fits Leiden (default) or Louvain communities on each layer, links communities in the selected layer pairs by the overlap coefficient of their members, and groups the linked communities into cross-layer meta-communities with a second community detection on the community graph.
- `fit_multilayer_weighted_jaccard()`
  - Fits Leiden (default) or Louvain communities on each layer, links communities in the selected layer pairs by a node-strength weighted Jaccard index, and groups the linked communities into cross-layer meta-communities with a second community detection on the community graph.
- `fit_multilayer_weighted_overlap()`
  - Fits Leiden (default) or Louvain communities on each layer, links communities in the selected layer pairs by a node-strength weighted overlap coefficient, and groups the linked communities into cross-layer meta-communities with a second community detection on the community graph.
- `fit_multilayer_identity_ties()`
  - Multislice modularity (Mucha et al. 2010): ties each node to its own copy in the selected layer pairs (weight = layer-link weight times `omega`) and detects meta-communities jointly on the coupled layers.
- `fit_multilayer_hungarian()`
  - Fits communities on each layer independently, then aligns labels between consecutive layers with Hungarian matching on community overlap (a snapshot-and-match baseline).
- `simulate_and_fit_multilayer()`
  - Simulates multiplex layers from a planted partition process and runs one of the five fit types (`jaccard`, `overlap`, `weighted_jaccard`, `weighted_overlap`, `identity`).
- `bootstrap_multilayer()`, `partition_stability()`, `co_assignment_ci()`
  - Parametric network bootstrap of a fit, the stability score, and co-assignment intervals.
- `plot_multilayer_series()`, `animate_multilayer_gif()`, `plot_multilayer_alluvial()`
  - Plots colored by the tracked meta-communities.

## Input handling

- Layers must share one node set unless `allow_unequal_nodes = TRUE`. With unequal node sets, every layer needs vertex names (igraph names or adjacency-matrix dimnames); nodes are matched across layers by name.
- Within-layer self-loops and zero-weight edges are dropped (the same rule as the Python package).
- `objective` sets the first-stage (per-layer) objective only; the second stage always maximizes modularity.
- `seed = NULL` (the default) uses the caller's random number stream; pass an integer for a reproducible fit. The Python package defaults to `seed=123`.

## Installation

Development install:

```r
# install.packages("remotes")
remotes::install_local(".")
```

GitHub install (recommended for collaborators):

```r
# install.packages("remotes")
remotes::install_github("jfedgerton/dynamic_multiplex", subdir = "r_code")
```

## Quick example

```r
library(dynamicmultiplex)

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
