from __future__ import annotations

from .multilayer_utils import (
    add_community_self_loops,
    community_overlap_edges,
    detect_interlayer_communities,
    fit_layer_communities,
    layer_node_strengths,
    make_layer_links,
    prepare_multilayer_graphs,
)


def fit_multilayer_weighted_overlap(
    layers,
    algorithm: str = "leiden",
    layer_links=None,
    min_similarity: float = 0.0,
    resolution_parameter: float = 1.0,
    directed: bool = False,
    add_self_loops: bool = True,
    self_loop_multiplier: float = 1.0,
    objective: str | None = None,
    seed: int | None = 123,
    allow_unequal_nodes: bool = False,
):
    """Fit per-layer communities and node-strength-weighted overlap ties.

    Parameters
    ----------
    objective : str | None
        Objective for the first-stage (per-layer) detection: ``"modularity"``
        or ``"cpm"`` (Leiden only). Defaults to ``"cpm"`` for directed layers
        and ``"modularity"`` otherwise. It does not affect the second stage,
        which always maximises modularity on the community graph.
    allow_unequal_nodes : bool
        If True, layers may hold different node sets (actors entering or
        leaving). Layers must then be networkx graphs with node labels; nodes
        are matched across layers by label and only nodes present in both
        layers of a linked pair contribute to the similarity. Default False.

    Returns
    -------
    dict
        Keys include ``layer_communities`` (per-layer detection),
        ``meta_communities`` (the cross-layer tracked partition from the
        second-stage detection, validated by ``bootstrap_multilayer``; see
        ``extract_meta_membership``), ``meta_ids``, ``interlayer_ties``, and
        ``layer_links``.
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

    node_weights = layer_node_strengths(graph_layers, directed=directed)

    interlayer_ties = community_overlap_edges(
        fit=fit,
        layer_links=links,
        metric="overlap",
        min_similarity=min_similarity,
        node_weights_by_layer=node_weights,
    )

    if add_self_loops:
        interlayer_ties = add_community_self_loops(
            edge_df=interlayer_ties,
            fit=fit,
            layer_links=links,
            self_loop_multiplier=self_loop_multiplier,
            min_similarity=min_similarity,
            directed=directed,
        )

    # Second-stage detection: group per-layer communities into cross-layer
    # meta-communities from the interlayer ties (the tracked partition).
    meta = detect_interlayer_communities(
        layer_communities=fit,
        interlayer_ties=interlayer_ties,
        algorithm=algorithm,
        resolution_parameter=resolution_parameter,
        seed=seed,
    )

    return {
        "algorithm": algorithm,
        "layer_communities": fit,
        "meta_communities": meta["membership"],
        "meta_ids": meta["meta_ids"],
        "layer_links": links,
        "interlayer_ties": interlayer_ties,
        "directed": directed,
        "weighting": "node_strength",
    }
