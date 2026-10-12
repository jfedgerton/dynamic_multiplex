from __future__ import annotations

from .multilayer_utils import (
    community_overlap_edges,
    fit_layer_communities,
    make_layer_links,
    prepare_multilayer_graphs,
    run_two_stage,
)


def fit_multilayer_overlap(
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
    tracking: str = "smoothed",
    repair_chains: bool = False,
    repair_min_similarity: float = 0.2,
):
    """Fit per-layer communities and interlayer weighted overlap ties.

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

    tracking : str
        ``"smoothed"`` (default, the 1.3.2 behaviour): stage 2 is fitted on all
        layers at once, so a layer's meta ids can depend on later layers.
        ``"filtered"``: stage 2 runs forward in time and the meta ids of layers
        1..t depend only on layers 1..t (no future information; see
        ``detect_interlayer_communities_filtered``). Filtered mode assumes the
        layers are in time order and needs a seed (not None).
    repair_chains : bool
        If True, rejoin meta-communities that stage 2 cut along time: a group
        that ends at layer s and one that starts at a linked layer t are joined
        when their communities are mutual best matches with similarity at
        least ``repair_min_similarity`` (see ``repair_meta_chains``). Default
        False (1.3.2 behaviour).
    repair_min_similarity : float
        Similarity floor for chain repair, in [0, 1]. Default 0.2.

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

    pair_ties = community_overlap_edges(
        fit=fit,
        layer_links=links,
        metric="overlap",
        min_similarity=min_similarity,
    )

    # Second-stage detection: group per-layer communities into cross-layer
    # meta-communities from the interlayer ties (the tracked partition).
    meta = run_two_stage(
        fit, pair_ties, links,
        add_self_loops=add_self_loops,
        self_loop_multiplier=self_loop_multiplier,
        min_similarity=min_similarity,
        directed=directed,
        algorithm=algorithm,
        resolution_parameter=resolution_parameter,
        seed=seed,
        tracking=tracking,
        repair_chains=repair_chains,
        repair_min_similarity=repair_min_similarity,
    )
    interlayer_ties = meta["interlayer_ties"]

    return {
        "algorithm": algorithm,
        "layer_communities": fit,
        "meta_communities": meta["membership"],
        "meta_ids": meta["meta_ids"],
        "layer_links": links,
        "interlayer_ties": interlayer_ties,
        "directed": directed,
        "tracking": meta["tracking"],
        "repair_chains": meta["repair_chains"],
        "repair_min_similarity": meta["repair_min_similarity"],
        "chain_repairs": meta["chain_repairs"],
        "node_labels": graph_layers[0].graph.get("node_labels"),
    }
