from __future__ import annotations

import pandas as pd

from .multilayer_utils import (
    _is_zero_indexed,
    detect_interlayer_communities,
    fit_layer_communities,
    make_layer_links,
    prepare_multilayer_graphs,
)


def fit_multilayer_identity_ties(
    layers,
    algorithm: str = "leiden",
    layer_links=None,
    resolution_parameter: float = 1.0,
    directed: bool = False,
    objective: str | None = None,
):
    """Fit per-layer communities with identity (node-level) interlayer ties.

    Returns
    -------
    dict
        Keys include ``layer_communities`` (per-layer detection),
        ``meta_communities`` (the cross-layer tracked partition; identity
        ties are node-level and lack community columns, so the second stage
        falls back to per-layer communities made globally distinct),
        ``meta_ids``, ``interlayer_ties``, and ``layer_links``.
    """
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

    interlayer_ties = pd.DataFrame(ties)

    # Second-stage detection. Identity ties are node-level (columns
    # from_layer, to_layer, node, layer_weight) and lack community columns, so
    # detect_interlayer_communities falls back to per-layer communities made
    # globally distinct (no cross-layer merging). Resolution has no argument
    # here, so pass 1.0.
    meta = detect_interlayer_communities(
        layer_communities=fit,
        interlayer_ties=interlayer_ties,
        algorithm=algorithm,
        resolution_parameter=1.0,
    )

    return {
        "algorithm": algorithm,
        "layer_communities": fit,
        "meta_communities": meta["membership"],
        "meta_ids": meta["meta_ids"],
        "layer_links": links,
        "interlayer_ties": interlayer_ties,
        "directed": directed,
    }
