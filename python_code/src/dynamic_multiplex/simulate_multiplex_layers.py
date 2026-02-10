from __future__ import annotations

import numpy as np

from .fit_multilayer_identity_ties import fit_multilayer_identity_ties
from .fit_multilayer_jaccard import fit_multilayer_jaccard
from .fit_multilayer_overlap import fit_multilayer_overlap


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
    if fit_type not in {"jaccard", "overlap", "identity"}:
        raise ValueError("`fit_type` must be one of {'jaccard', 'overlap', 'identity'}.")

    rng = np.random.default_rng(seed)
    memberships = rng.integers(1, n_communities + 1, size=n_nodes)

    layers = []
    for _ in range(n_layers):
        mat = np.zeros((n_nodes, n_nodes))
        if directed:
            for i in range(n_nodes):
                for j in range(n_nodes):
                    if i == j:
                        continue
                    prob = p_in if memberships[i] == memberships[j] else p_out
                    mat[i, j] = rng.binomial(1, prob)
        else:
            for i in range(n_nodes - 1):
                for j in range(i + 1, n_nodes):
                    prob = p_in if memberships[i] == memberships[j] else p_out
                    tie = rng.binomial(1, prob)
                    mat[i, j] = tie
                    mat[j, i] = tie
        layers.append(mat)

    if fit_type == "jaccard":
        fit = fit_multilayer_jaccard(
            layers,
            algorithm=algorithm,
            layer_links=layer_links,
            min_similarity=min_similarity,
            directed=directed,
        )
    elif fit_type == "overlap":
        fit = fit_multilayer_overlap(
            layers,
            algorithm=algorithm,
            layer_links=layer_links,
            min_similarity=min_similarity,
            directed=directed,
        )
    else:
        fit = fit_multilayer_identity_ties(
            layers,
            algorithm=algorithm,
            layer_links=layer_links,
            directed=directed,
        )

    return {
        "layers": layers,
        "true_membership": memberships,
        "fit": fit,
        "directed": directed,
    }
