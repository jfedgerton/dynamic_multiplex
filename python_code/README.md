# dynamic_multiplex (Python)

`dynamic_multiplex` is the Python implementation of DynMux. DynMux is a suite of tools for temporal networks that change over time. Nodes enter and leave, groups form, split, merge and dissolve, and periods can differ in size. Every function tracks nodes and communities across periods. A community keeps one id for as long as it persists, and a newly formed group gets a new id. The two-stage fitters detect communities in each period and link them to communities in the linked periods. `layer_links` sets which periods are linked, and by default that is adjacent periods only, so distant periods are not pooled. Multislice modularity (Mucha et al. 2010) and Hungarian matching are included as comparison methods. The bootstrap stability score helps choose between DynMux and multislice.

## Why this package

Standard multislice approaches (Mucha et al. 2010) connect all layers to all layers, meaning community structure at distant time periods influences assignments everywhere. When applied to temporal data, this creates a pooling problem: the community configuration at 2010 affects the structure detected at 1950. `dynamic_multiplex` provides explicit control over *which layers influence which layers* via a `layer_links` argument (`from`, `to`, `weight`), defaulting to adjacent-only temporal coupling.

## Functions

> Note: with `directed=True`, Leiden keeps edge direction in the per-layer detection; Louvain collapses directed layers to undirected weighted graphs (reciprocal arcs summed). The community-graph second stage is undirected; the multislice stage for identity ties keeps edge direction in each slice (directed null model), with undirected identity ties.

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
- `partition_membership_bipartite()`, `partition_membership_strict()`, `fit_multilayer_from_partitions()` (Python only)
  - Per-slice partitions of membership data without projection: Barber's bipartite modularity (Leiden) and strict (all-or-nothing) hypergraph modularity (hypernetx; `pip install dynamic-multiplex[hypergraph]`, which needs pandas < 3). `fit_multilayer_from_partitions()` tracks any per-slice partitions over time with the DynMux second stage or Hungarian matching. These are comparison specifications; in pilot simulations both did worse than the projections. The R package deliberately has no counterpart.
- `bootstrap_multilayer()`, `partition_stability()`
  - Parametric network bootstrap of a fit and the stability score s, for comparing two fits of the same network. Since 1.4.0 no interval or accuracy-floor outputs are reported.

## Input handling

- Layers must share one node set unless `allow_unequal_nodes=True`. Nodes may then enter and leave, with or without a change in layer size; pass networkx graphs with node labels, and nodes are matched across layers by label and keep their labels as ids. This holds for every fit function and for `bootstrap_multilayer()`.
- Within-layer self-loops and zero-weight edges are dropped (the same rule as the R package).
- Two-stage fits: `repair_chains=True` rejoins groups that stage 2 cut along time (mutual best match with similarity at least `repair_min_similarity`, default 0.2); `tracking="filtered"` makes the meta ids of layers 1..t depend only on layers 1..t. Both are off by default.
- `objective` sets the first-stage (per-layer) objective only; the second stage always maximizes modularity.
- `seed=123` by default, as in the R package since 1.4.0; pass `seed=None` for unseeded detection.

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
