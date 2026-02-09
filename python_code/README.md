# dynamic_multiplex (Python)

`dynamic_multiplex` is a Python package for multiplex community modeling with customizable interlayer ties.

## Why this package

Most multilayer workflows assume every network layer connects to every other layer. This package adds explicit control over *which layers influence which layers* via a `layer_links` argument (`from`, `to`, `weight`).

## Current prototype functions

- `fit_multilayer_jaccard()`
  - Fits Louvain or Leiden communities on each layer.
  - Builds interlayer ties between communities using weighted Jaccard similarity across selected layer pairs.
- `fit_multilayer_overlap()`
  - Fits Louvain or Leiden communities on each layer.
  - Builds interlayer ties between communities using weighted overlap coefficient across selected layer pairs.
- `fit_multilayer_identity_ties()`
  - Fits Louvain or Leiden communities on each layer.
  - Builds interlayer ties only for the same node across selected adjacent layers.
- `simulate_and_fit_multilayer()`
  - Simulates multiplex layers from a planted partition process and runs one of the three fitting strategies.

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
)
```

## Testing

```bash
pip install -e ./python_code[dev,louvain]
pytest
```
