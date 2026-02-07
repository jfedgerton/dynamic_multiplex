from __future__ import annotations

from typing import Callable

import igraph as ig
import pandas as pd

from .multilayer_utils import (
    make_layer_links,
    prepare_multilayer_graphs,
    weighted_jaccard,
    weighted_overlap,
)


def _fit_layer_communities(
    graph_layers: list[ig.Graph], algorithm: str = "louvain", resolution_parameter: float = 1.0, directed: bool = False
) -> list[dict]:
    out: list[dict] = []
    for g in graph_layers:
        g_for_louvain = g
        if algorithm == "louvain":
            if directed and g.is_directed():
                g_for_louvain = g.as_undirected(combine_edges="sum")
            cl = g_for_louvain.community_multilevel(weights=g_for_louvain.es["weight"] if "weight" in g_for_louvain.es.attributes() else None)
            g_used = g_for_louvain
        elif algorithm == "leiden":
            objective = "CPM" if directed else "modularity"
            cl = g.community_leiden(
                objective_function=objective,
                resolution_parameter=resolution_parameter,
                weights=g.es["weight"] if "weight" in g.es.attributes() else None,
            )
            g_used = g
        else:
            raise ValueError("algorithm must be 'louvain' or 'leiden'")

        membership = cl.membership
        comms: dict[int, list[int]] = {}
        for node, cid in enumerate(membership, start=1):
            comms.setdefault(cid + 1, []).append(node)

        modularity = None
        if not g_used.is_directed():
            modularity = g_used.modularity(membership, weights=g_used.es["weight"] if "weight" in g_used.es.attributes() else None)

        out.append(
            {
                "membership": membership,
                "modularity": modularity,
                "communities": comms,
            }
        )
    return out


def _community_overlap_edges(
    fit: list[dict], layer_links: pd.DataFrame, metric: str, min_similarity: float
) -> pd.DataFrame:
    sim_fun: Callable[[set[int], set[int]], float]
    if metric == "jaccard":
        sim_fun = weighted_jaccard
    elif metric == "overlap":
        sim_fun = weighted_overlap
    else:
        raise ValueError("metric must be jaccard or overlap")

    rows: list[dict] = []
    for _, link in layer_links.iterrows():
        from_idx = int(link["from"])
        to_idx = int(link["to"])
        w = float(link["weight"])

        from_comms = fit[from_idx - 1]["communities"]
        to_comms = fit[to_idx - 1]["communities"]

        for from_cid, from_nodes in from_comms.items():
            for to_cid, to_nodes in to_comms.items():
                sim = sim_fun(set(from_nodes), set(to_nodes))
                ws = sim * w
                if ws >= min_similarity:
                    rows.append(
                        {
                            "from_layer": from_idx,
                            "to_layer": to_idx,
                            "from_community": from_cid,
                            "to_community": to_cid,
                            "similarity": sim,
                            "layer_weight": w,
                            "weighted_similarity": ws,
                        }
                    )

    return pd.DataFrame(rows)


def fit_multilayer_jaccard(
    layers,
    algorithm: str = "louvain",
    layer_links: pd.DataFrame | None = None,
    min_similarity: float = 0.0,
    resolution_parameter: float = 1.0,
    directed: bool = False,
):
    graphs = prepare_multilayer_graphs(layers, directed=directed)
    links = make_layer_links(len(graphs), layer_links)
    fit = _fit_layer_communities(graphs, algorithm, resolution_parameter, directed)
    ties = _community_overlap_edges(fit, links, "jaccard", min_similarity)
    return {"algorithm": algorithm, "directed": directed, "layer_communities": fit, "layer_links": links, "interlayer_ties": ties}


def fit_multilayer_overlap(
    layers,
    algorithm: str = "louvain",
    layer_links: pd.DataFrame | None = None,
    min_similarity: float = 0.0,
    resolution_parameter: float = 1.0,
    directed: bool = False,
):
    graphs = prepare_multilayer_graphs(layers, directed=directed)
    links = make_layer_links(len(graphs), layer_links)
    fit = _fit_layer_communities(graphs, algorithm, resolution_parameter, directed)
    ties = _community_overlap_edges(fit, links, "overlap", min_similarity)
    return {"algorithm": algorithm, "directed": directed, "layer_communities": fit, "layer_links": links, "interlayer_ties": ties}


def fit_multilayer_identity_ties(
    layers,
    algorithm: str = "louvain",
    layer_links: pd.DataFrame | None = None,
    resolution_parameter: float = 1.0,
    directed: bool = False,
):
    graphs = prepare_multilayer_graphs(layers, directed=directed)
    links = make_layer_links(len(graphs), layer_links)
    fit = _fit_layer_communities(graphs, algorithm, resolution_parameter, directed)

    n_nodes = graphs[0].vcount()
    if any(g.vcount() != n_nodes for g in graphs):
        raise ValueError("All layers must have same number of nodes for identity ties")

    rows = []
    for _, link in links.iterrows():
        for node in range(1, n_nodes + 1):
            rows.append({"from_layer": int(link["from"]), "to_layer": int(link["to"]), "node": node, "layer_weight": float(link["weight"])})

    ties = pd.DataFrame(rows)
    return {"algorithm": algorithm, "directed": directed, "layer_communities": fit, "layer_links": links, "interlayer_ties": ties}
