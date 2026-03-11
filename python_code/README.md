# dynamic_multiplex (Python)

`dynamic_multiplex` is a Python package for multiplex community modeling with customizable interlayer ties.

## Why this package

Standard multislice approaches (Mucha et al. 2010) connect all layers to all layers, meaning community structure at distant time periods influences assignments everywhere. When applied to temporal data, this creates a pooling problem: the community configuration at 2010 affects the structure detected at 1950. `dynamic_multiplex` provides explicit control over *which layers influence which layers* via a `layer_links` argument (`from`, `to`, `weight`), defaulting to adjacent-only temporal coupling.

## Method

This is a **two-stage** approach:

1. **Per-layer community detection**: Louvain or Leiden is run independently on each layer. Community assignments within a layer are determined solely by that layer's topology.
2. **Post-hoc interlayer ties**: After communities are detected, interlayer ties are constructed between communities (or nodes, for identity ties) in selected layer pairs using similarity metrics.

The interlayer ties describe relationships between layers but **do not retroactively influence** the per-layer community assignments. This contrasts with true coupled multislice methods (e.g., Mucha et al. 2010), where a single optimization jointly determines communities across all layers. The advantage of the two-stage approach is that each layer's communities reflect only that layer's structure, avoiding the temporal pooling problem where distant layers distort local community detection.

## Functions

- `fit_multilayer_jaccard()`
  - Detects communities independently on each layer (stage 1).
  - Builds interlayer ties between communities using Jaccard similarity across selected layer pairs (stage 2).
- `fit_multilayer_overlap()`
  - Detects communities independently on each layer (stage 1).
  - Builds interlayer ties between communities using overlap coefficient across selected layer pairs (stage 2).
- `fit_multilayer_weighted_jaccard()`
  - Detects communities independently on each layer (stage 1).
  - Builds interlayer ties using node-strength weighted Jaccard similarity across selected layer pairs (stage 2).
- `fit_multilayer_weighted_overlap()`
  - Detects communities independently on each layer (stage 1).
  - Builds interlayer ties using node-strength weighted overlap coefficient across selected layer pairs (stage 2).
- `fit_multilayer_identity_ties()`
  - Detects communities independently on each layer (stage 1).
  - Builds node-level interlayer ties linking the same node across selected adjacent layers (stage 2).
- `simulate_and_fit_multilayer()`
  - Simulates multiplex layers from a planted partition process and runs one of the five fitting strategies.

All `fit_multilayer_*` functions accept a `seed` parameter for reproducible community detection. All layers must share the same node set; graphs with arbitrary (string or non-contiguous integer) node labels are supported via automatic internal remapping.

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
