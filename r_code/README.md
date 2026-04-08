# dynamic_multiplex

`dynamic_multiplex` is an R package for multiplex community modeling with customizable interlayer ties.

## Why this package

Standard multislice approaches (Mucha et al. 2010) connect all layers to all layers, meaning community structure at distant time periods influences assignments everywhere. When applied to temporal data, this creates a pooling problem: the community configuration at 2010 affects the structure detected at 1950. `dynamic_multiplex` provides explicit control over *which layers influence which layers* via a `layer_links` argument (`from`, `to`, `weight`), defaulting to adjacent-only temporal coupling.

## Functions

> Note: Louvain in `igraph` is undirected. When `directed = TRUE` and `algorithm = "louvain"`, layers are collapsed to undirected weighted graphs before clustering.


- `fit_multilayer_jaccard()`
  - Fits Louvain or Leiden communities on each layer (supports directed layers).
  - Builds interlayer ties between communities using Jaccard similarity across selected layer pairs.
- `fit_multilayer_overlap()`
  - Fits Louvain or Leiden communities on each layer (supports directed layers).
  - Builds interlayer ties between communities using overlap coefficient across selected layer pairs.
- `fit_multilayer_weighted_jaccard()`
  - Fits Louvain or Leiden communities on each layer (supports directed layers).
  - Builds interlayer ties using node-strength weighted Jaccard similarity across selected layer pairs.
- `fit_multilayer_weighted_overlap()`
  - Fits Louvain or Leiden communities on each layer (supports directed layers).
  - Builds interlayer ties using node-strength weighted overlap coefficient across selected layer pairs.
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
