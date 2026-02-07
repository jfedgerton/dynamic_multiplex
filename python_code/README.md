# bayesnet-py

Python companion package for multiplex and dynamic community modeling.

## Features

- Louvain and Leiden community detection across network layers.
- Configurable interlayer tie construction:
  - weighted Jaccard
  - weighted overlap coefficient
  - node-identity ties
- Directed-layer support.
- Simulation and benchmarking tools for evolving systems.

## Install

```bash
pip install -e python_code
```

## Quickstart

```python
from bayesnet_py import simulate_and_fit_multilayer

res = simulate_and_fit_multilayer(
    n_nodes=50,
    n_layers=4,
    fit_type="jaccard",
    algorithm="leiden",
    directed=True,
    seed=123,
)

print(res["fit"]["interlayer_ties"].head())
```
