from __future__ import annotations

import pandas as pd

from .multilayer_utils import fit_layer_communities, make_layer_links, prepare_multilayer_graphs


def fit_multilayer_identity_ties(
    layers,
    algorithm: str = "louvain",
    layer_links=None,
    resolution_parameter: float = 1.0,
    directed: bool = False,
):
    graph_layers = prepare_multilayer_graphs(layers, directed=directed)
    links = make_layer_links(len(graph_layers), layer_links)
    fit = fit_layer_communities(
        graph_layers,
        algorithm=algorithm,
        resolution_parameter=resolution_parameter,
        directed=directed,
    )

    n_nodes = graph_layers[0].number_of_nodes()
    if not all(g.number_of_nodes() == n_nodes for g in graph_layers):
        raise ValueError("All layers must have the same number of nodes for identity ties.")

    ties = []
    for _, row in links.iterrows():
        for node in range(1, n_nodes + 1):
            ties.append(
                {
                    "from_layer": int(row["from"]),
                    "to_layer": int(row["to"]),
                    "node": node,
                    "layer_weight": float(row["weight"]),
                }
            )

    return {
        "algorithm": algorithm,
        "layer_communities": fit,
        "layer_links": links,
        "interlayer_ties": pd.DataFrame(ties),
        "directed": directed,
    }
