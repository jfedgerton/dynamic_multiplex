from __future__ import annotations

import numpy as np

from .fit import (
    fit_multilayer_identity_ties,
    fit_multilayer_jaccard,
    fit_multilayer_overlap,
)


def simulate_and_fit_multilayer(
    n_nodes: int = 100,
    n_layers: int = 4,
    n_communities: int = 4,
    p_in: float = 0.2,
    p_out: float = 0.05,
    fit_type: str = "jaccard",
    algorithm: str = "louvain",
    layer_links=None,
    min_similarity: float = 0.0,
    seed: int | None = None,
    directed: bool = False,
):
    rng = np.random.default_rng(seed)
    membership = rng.integers(1, n_communities + 1, size=n_nodes)

    layers = []
    for _ in range(n_layers):
        mat = np.zeros((n_nodes, n_nodes), dtype=float)
        if directed:
            for i in range(n_nodes):
                for j in range(n_nodes):
                    if i == j:
                        continue
                    p = p_in if membership[i] == membership[j] else p_out
                    mat[i, j] = rng.binomial(1, p)
        else:
            for i in range(n_nodes - 1):
                for j in range(i + 1, n_nodes):
                    p = p_in if membership[i] == membership[j] else p_out
                    tie = rng.binomial(1, p)
                    mat[i, j] = tie
                    mat[j, i] = tie
        layers.append(mat)

    if fit_type == "jaccard":
        fit = fit_multilayer_jaccard(layers, algorithm, layer_links, min_similarity, directed=directed)
    elif fit_type == "overlap":
        fit = fit_multilayer_overlap(layers, algorithm, layer_links, min_similarity, directed=directed)
    elif fit_type == "identity":
        fit = fit_multilayer_identity_ties(layers, algorithm, layer_links, directed=directed)
    else:
        raise ValueError("fit_type must be jaccard, overlap, or identity")

    return {"layers": layers, "true_membership": membership, "fit": fit, "directed": directed}


def simulate_evolving_multilayer(
    n_nodes: int = 150,
    n_layers: int = 6,
    n_communities: int = 4,
    persistence: float = 0.85,
    p_in: float = 0.20,
    p_out: float = 0.03,
    directed: bool = False,
    seed: int | None = None,
):
    rng = np.random.default_rng(seed)
    membership = np.zeros((n_nodes, n_layers), dtype=int)
    membership[:, 0] = rng.integers(1, n_communities + 1, size=n_nodes)

    for t in range(1, n_layers):
        keep = rng.binomial(1, persistence, size=n_nodes).astype(bool)
        membership[:, t] = membership[:, t - 1]
        membership[~keep, t] = rng.integers(1, n_communities + 1, size=np.sum(~keep))

    layers = []
    for t in range(n_layers):
        mat = np.zeros((n_nodes, n_nodes), dtype=float)
        if directed:
            for i in range(n_nodes):
                for j in range(n_nodes):
                    if i == j:
                        continue
                    p = p_in if membership[i, t] == membership[j, t] else p_out
                    mat[i, j] = rng.binomial(1, p)
        else:
            for i in range(n_nodes - 1):
                for j in range(i + 1, n_nodes):
                    p = p_in if membership[i, t] == membership[j, t] else p_out
                    tie = rng.binomial(1, p)
                    mat[i, j] = tie
                    mat[j, i] = tie
        layers.append(mat)

    return {"layers": layers, "membership": membership}
