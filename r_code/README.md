# dynamicmultiplex

`dynamicmultiplex` is the R implementation of DynMux. DynMux is a suite of tools for temporal networks that change over time. Nodes enter and leave, groups form, split, merge and dissolve, and periods can differ in size. The temporal fitters track nodes and communities across periods. A community keeps one id for as long as it persists, and a newly formed group gets a new id. The two-stage fitters detect communities in each period and link them to communities in the linked periods. `layer_links` sets which periods are directly linked, and by default that is adjacent periods only. With the default `tracking = "smoothed"`, all periods are fitted jointly, so later periods can change earlier labels. `tracking = "filtered"` prevents that. Multislice modularity (Mucha et al. 2010) and Hungarian matching are included as comparison methods. The bootstrap stability score helps choose between DynMux and multislice.

## Why this package

Standard multislice approaches (Mucha et al. 2010) connect all layers to all layers, meaning community structure at distant time periods influences assignments everywhere. When applied to temporal data, this creates a pooling problem: the community configuration at 2010 affects the structure detected at 1950. `dynamicmultiplex` provides explicit control over *which layers influence which layers* via a `layer_links` argument (`from`, `to`, `weight`), defaulting to adjacent-only temporal coupling.

## Functions

> Note: `igraph`'s Louvain and Leiden are undirected. When `directed = TRUE`, layers are collapsed to undirected weighted graphs (reciprocal arcs summed) before per-layer clustering; the Python package keeps edge direction for Leiden. The multislice fit (`fit_multilayer_identity_ties()`) keeps edge direction in each slice when `directed = TRUE`.

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
- `layers_from_membership()`
  - Builds one weighted layer per time slice from group memberships (for example states in IGOs): shared-group counts, size-weighted, or Kumar hypergraph reweighting. No weighting is the default; each slice is built on its own.
- `bootstrap_multilayer()`, `partition_stability()`
  - Parametric network bootstrap of a fit and the stability score s, for comparing two fits of the same network. Since 1.4.0 no interval or accuracy-floor outputs are reported.
- `plot_multilayer_series()`, `animate_multilayer_gif()`, `plot_multilayer_alluvial()`
  - Plots colored by the tracked meta-communities.

## Input handling

- Layers must share one node set unless `allow_unequal_nodes = TRUE`. Nodes may then enter and leave, with or without a change in layer size; every layer needs vertex names (igraph names or adjacency-matrix dimnames), and nodes are matched across layers by name. This holds for every fit function and for `bootstrap_multilayer()`.
- Within-layer self-loops and zero-weight edges are dropped (the same rule as the Python package).
- Two-stage fits: `repair_chains = TRUE` rejoins groups that stage 2 cut along time (mutual best match with similarity at least `repair_min_similarity`, default 0.2); `tracking = "filtered"` makes the meta ids of layers 1..t depend only on layers 1..t. Both are off by default.
- `objective` sets the first-stage (per-layer) objective only; the second stage always maximizes modularity.
- `seed = 123` by default (since 1.4.0, as in the Python package); pass `seed = NULL` for unseeded detection. The caller's random number stream is restored after the fit.

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
