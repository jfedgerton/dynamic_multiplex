from __future__ import annotations

import pandas as pd

from .multilayer_utils import (
    _id_offset,
    _sort_key,
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
    omega: float = 1.0,
    directed: bool = False,
    objective: str | None = None,
    seed: int | None = 123,
    allow_unequal_nodes: bool = False,
):
    """Fit per-layer communities with identity (node-level) interlayer ties.

    Parameters
    ----------
    omega : float
        Interlayer coupling strength for the multislice supra-graph (Mucha's
        omega). Multiplies the interlayer identity-edge weights on top of any
        ``layer_links`` weights. Larger ``omega`` couples layers more strongly
        and, past a point, collapses everything into one meta-community;
        smaller ``omega`` decouples toward independent per-layer detection.
    objective : str | None
        Objective for the first-stage (per-layer) detection stored in
        ``layer_communities``: ``"modularity"`` or ``"cpm"`` (Leiden only).
        Defaults to ``"cpm"`` for directed layers and ``"modularity"``
        otherwise. It does not affect ``meta_communities``, which always
        optimise multislice modularity.
    directed : bool
        If True, layers are directed graphs. The multislice stage keeps arc
        direction in every slice (directed null model ``k_out * k_in / m_s``;
        a fully reciprocated network gives the same ``meta_communities`` as
        its undirected twin at the same ``omega``); the interlayer identity
        ties stay undirected. Every layer must then be a directed graph or an
        adjacency matrix. Directed layers fitted with ``directed=False`` are
        symmetrised (``A + A^T``) with a warning.
    allow_unequal_nodes : bool
        If True, layers may hold different node sets; they must then be
        networkx graphs with node labels, matched across layers by label.
    resolution_parameter : float
        Also forwarded to the multislice supra-graph detection as Mucha's
        modularity resolution (larger values yield more, smaller
        meta-communities).

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
    graph_layers = prepare_multilayer_graphs(
        layers, directed=directed, require_same_nodes=not allow_unequal_nodes
    )
    links = make_layer_links(len(graph_layers), layer_links)
    fit = fit_layer_communities(
        graph_layers,
        algorithm=algorithm,
        resolution_parameter=resolution_parameter,
        directed=directed,
        objective=objective,
        seed=seed,
    )

    ties = []
    for _, row in links.iterrows():
        g_from = graph_layers[int(row["from"]) - 1]
        g_to = graph_layers[int(row["to"]) - 1]

        shared = sorted(set(g_from.nodes()) & set(g_to.nodes()), key=_sort_key)
        # every layer uses the same label-to-id rule (see _id_offset), so a
        # shared node has one id in both layers ----
        off = _id_offset(g_from)

        for node in shared:
            node_id = node + off if off else node
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
        omega=omega,
        resolution_parameter=resolution_parameter,
        seed=seed,
        directed=directed,
    )

    return {
        "algorithm": algorithm,
        "layer_communities": fit,
        "meta_communities": meta_membership,
        "meta_ids": None,
        "layer_links": links,
        "interlayer_ties": interlayer_ties,
        "directed": directed,
        "node_labels": graph_layers[0].graph.get("node_labels"),
    }
