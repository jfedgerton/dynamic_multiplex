from __future__ import annotations

from .multilayer_utils import community_overlap_edges, fit_layer_communities, make_layer_links, prepare_multilayer_graphs


def fit_multilayer_jaccard(
    layers,
    algorithm: str = "louvain",
    layer_links=None,
    min_similarity: float = 0.0,
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

    interlayer_ties = community_overlap_edges(
        fit=fit,
        layer_links=links,
        metric="jaccard",
        min_similarity=min_similarity,
    )

    return {
        "algorithm": algorithm,
        "layer_communities": fit,
        "layer_links": links,
        "interlayer_ties": interlayer_ties,
        "directed": directed,
    }
