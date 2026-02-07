from __future__ import annotations

from typing import Iterable

import igraph as ig
import numpy as np
import pandas as pd


def prepare_multilayer_graphs(layers: Iterable, directed: bool = False) -> list[ig.Graph]:
    layers = list(layers)
    if len(layers) < 2:
        raise ValueError("layers must contain at least two network layers")

    graphs: list[ig.Graph] = []
    for i, layer in enumerate(layers, start=1):
        if isinstance(layer, ig.Graph):
            graphs.append(layer)
            continue

        arr = np.asarray(layer)
        if arr.ndim != 2 or arr.shape[0] != arr.shape[1]:
            raise ValueError(f"Layer {i} must be a square adjacency matrix or igraph Graph")

        mode = "directed" if directed else "undirected"
        g = ig.Graph.Weighted_Adjacency(arr.tolist(), mode=mode, loops=False)
        graphs.append(g)

    return graphs


def make_layer_links(n_layers: int, layer_links: pd.DataFrame | None = None) -> pd.DataFrame:
    if layer_links is None:
        return pd.DataFrame(
            {
                "from": np.arange(1, n_layers),
                "to": np.arange(2, n_layers + 1),
                "weight": np.ones(n_layers - 1, dtype=float),
            }
        )

    if not {"from", "to"}.issubset(layer_links.columns):
        raise ValueError("layer_links must include columns: from, to")

    links = layer_links.copy()
    if "weight" not in links.columns:
        links["weight"] = 1.0

    if ((links["from"] < 1) | (links["to"] < 1) | (links["from"] > n_layers) | (links["to"] > n_layers)).any():
        raise ValueError("layer_links indices must be between 1 and number of layers")

    return links[["from", "to", "weight"]].reset_index(drop=True)


def weighted_jaccard(a: set[int], b: set[int]) -> float:
    union = len(a.union(b))
    return 0.0 if union == 0 else len(a.intersection(b)) / union


def weighted_overlap(a: set[int], b: set[int]) -> float:
    denom = min(len(a), len(b))
    return 0.0 if denom == 0 else len(a.intersection(b)) / denom
