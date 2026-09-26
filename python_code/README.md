# dynamic_multiplex (Python)

`dynamic_multiplex` is a Python package for multiplex community modeling with customizable interlayer ties.

## Why this package

Standard multislice approaches (Mucha et al. 2010) connect all layers to all layers, meaning community structure at distant time periods influences assignments everywhere. When applied to temporal data, this creates a pooling problem: the community configuration at 2010 affects the structure detected at 1950. `dynamic_multiplex` provides explicit control over *which layers influence which layers* via a `layer_links` argument (`from`, `to`, `weight`), defaulting to adjacent-only temporal coupling.

## Functions

> Note: with `directed=True`, Leiden keeps edge direction in the per-layer detection; Louvain collapses directed layers to undirected weighted graphs (reciprocal arcs summed). The second stage (community graph, or multislice for identity ties) is undirected.

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

## Input handling

- Layers must share one node set unless `allow_unequal_nodes=True`. With unequal node sets, pass networkx graphs with node labels; nodes are matched across layers by label and keep their labels as ids.
- Within-layer self-loops and zero-weight edges are dropped (the same rule as the R package).
- `objective` sets the first-stage (per-layer) objective only; the second stage always maximizes modularity.
- `seed=123` by default; pass `seed=None` for unseeded detection. The R package defaults to `seed = NULL` (the caller's random number stream).

## Installation (development)

```bash
pip install -e ./python_code
```

Optional algorithms:

```bash
pip install -e ./python_code[louvain]
pip install -e ./python_code[leiden]
```

## Quick example

```python
from dynamic_multiplex import simulate_and_fit_multilayer, fit_multilayer_overlap

sim = simulate_and_fit_multilayer(
    directed=True,
    n_nodes=50,
    n_layers=4,
    n_communities=3,
    fit_type="jaccard",
    algorithm="louvain",
    seed=123,
)

print(sim["fit"]["interlayer_ties"].head())

custom_links = [
    {"from": 1, "to": 2, "weight": 1.0},
    {"from": 2, "to": 4, "weight": 0.6},
]

fit_overlap = fit_multilayer_overlap(
    sim["layers"],
    algorithm="leiden",
    layer_links=custom_links,
    min_similarity=0.1,
    add_self_loops=True,
    self_loop_multiplier=1.0,
)
```

## Testing

```bash
pip install -e ./python_code[dev,louvain]
pytest
```
