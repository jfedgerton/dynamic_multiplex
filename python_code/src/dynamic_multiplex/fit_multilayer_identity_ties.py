from __future__ import annotations

import pandas as pd

from .multilayer_utils import (
    _is_zero_indexed,
    detect_multislice_communities,
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
        ``meta_communities`` (the node-level Mucha multislice partition: one
        supra-graph stacking each layer's adjacency plus identity interlayer
        ties, detected in a single pass, so a node's meta-community can be
        pulled across layers through the coupling), ``meta_ids`` (``None``),
        ``interlayer_ties``, and ``layer_links``.
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

    # Node-level Mucha multislice second stage: stack layers (intra = original
    # adjacency) with identity interlayer ties and run a single detection on
    # the supra-graph, so a node's meta-community can be pulled across layers
    # through the coupling. Per-layer detection (layer_communities) is unchanged.
    meta_membership = detect_multislice_communities(
        graph_layers=graph_layers,
        interlayer_ties=interlayer_ties,
        algorithm=algorithm,
    )

    return {
        "algorithm": algorithm,
        "layer_communities": fit,
        "meta_communities": meta_membership,
        "meta_ids": None,
        "layer_links": links,
        "interlayer_ties": interlayer_ties,
        "directed": directed,
    }
