from __future__ import annotations

from .multilayer_utils import (
    add_community_self_loops,
    community_overlap_edges,
    fit_layer_communities,
    make_layer_links,
    prepare_multilayer_graphs,
)
from .multilayer_utils import community_overlap_edges, fit_layer_communities, make_layer_links, prepare_multilayer_graphs


def fit_multilayer_overlap(
    layers,
    algorithm: str = "louvain",
    layer_links=None,
    min_similarity: float = 0.0,
    resolution_parameter: float = 1.0,
    directed: bool = False,
    add_self_loops: bool = True,
    self_loop_multiplier: float = 1.0,
):
    graph_layers = prepare_multilayer_graphs(layers, directed=directed)
    links = make_layer_links(len(graph_layers), layer_links)
    fit = fit_layer_communities(
        graph_layers,
        algorithm=algorithm,
        resolution_parameter=resolution_parameter,
        directed=directed,
    )

    interlayer_ties = community_overlap_edges(
        fit=fit,
        layer_links=links,
        metric="overlap",
        min_similarity=min_similarity,
    )

    if add_self_loops:
        interlayer_ties = add_community_self_loops(
            edge_df=interlayer_ties,
            fit=fit,
            layer_links=links,
            self_loop_multiplier=self_loop_multiplier,
            min_similarity=min_similarity,
        )

    return {
        "algorithm": algorithm,
        "layer_communities": fit,
        "layer_links": links,
        "interlayer_ties": interlayer_ties,
        "directed": directed,
    }
