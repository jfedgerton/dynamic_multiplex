from __future__ import annotations

import pandas as pd

from .multilayer_utils import (
    _is_zero_indexed,
    fit_layer_communities,
    make_layer_links,
    prepare_multilayer_graphs,
)


def fit_multilayer_identity_ties(
    layers,
    algorithm: str = "louvain",
    layer_links=None,
    resolution_parameter: float = 1.0,
    directed: bool = False,
    objective: str | None = None,
):
    graph_layers = prepare_multilayer_graphs(layers, directed=directed)
    links = make_layer_links(len(graph_layers), layer_links)
    fit = fit_layer_communities(
        graph_layers,
        algorithm=algorithm,
        resolution_parameter=resolution_parameter,
        directed=directed,
        objective=objective,
    )

    ties = []
    for _, row in links.iterrows():
        g_from = graph_layers[int(row["from"]) - 1]
        g_to = graph_layers[int(row["to"]) - 1]

        shared = sorted(set(g_from.nodes()) & set(g_to.nodes()))
        both_zero = _is_zero_indexed(g_from) and _is_zero_indexed(g_to)

        for node in shared:
            node_id = node + 1 if both_zero else node
            ties.append(
                {
                    "from_layer": int(row["from"]),
                    "to_layer": int(row["to"]),
                    "node": node_id,
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
