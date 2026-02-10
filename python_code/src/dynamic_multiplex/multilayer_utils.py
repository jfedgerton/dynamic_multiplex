from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import networkx as nx
import numpy as np
import pandas as pd


@dataclass
class LayerCommunityFit:
    membership: dict[int, int]
    modularity: float | None
    communities: dict[int, list[int]]


def _as_graph(layer, directed: bool) -> nx.Graph | nx.DiGraph:
    if isinstance(layer, (nx.Graph, nx.DiGraph)):
        return layer

    matrix = np.asarray(layer)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("Each layer must be a square adjacency matrix or a NetworkX graph.")

    graph_type = nx.DiGraph if directed else nx.Graph
    graph = nx.from_numpy_array(matrix, create_using=graph_type)
    return graph


def prepare_multilayer_graphs(layers: list, directed: bool = False) -> list[nx.Graph | nx.DiGraph]:
    if not isinstance(layers, list) or len(layers) < 2:
        raise ValueError("`layers` must be a list with at least two network layers.")

    return [_as_graph(layer, directed=directed) for layer in layers]


def make_layer_links(n_layers: int, layer_links: list[dict] | pd.DataFrame | None = None) -> pd.DataFrame:
    if layer_links is None:
        return pd.DataFrame(
            {
                "from": np.arange(1, n_layers),
                "to": np.arange(2, n_layers + 1),
                "weight": 1.0,
            }
        )

    links = pd.DataFrame(layer_links).copy()
    if not {"from", "to"}.issubset(links.columns):
        raise ValueError("`layer_links` must contain `from` and `to` columns.")

    if "weight" not in links.columns:
        links["weight"] = 1.0

    if ((links[["from", "to"]] < 1).any().any()) or ((links[["from", "to"]] > n_layers).any().any()):
        raise ValueError("`layer_links` indices must be between 1 and number of layers.")

    return links[["from", "to", "weight"]]


def _communities_to_membership(communities: Iterable[set[int]]) -> dict[int, int]:
    membership: dict[int, int] = {}
    for idx, nodes in enumerate(communities, start=1):
        for n in nodes:
            membership[n + 1] = idx
    return membership


def fit_layer_communities(
    graph_layers: list[nx.Graph | nx.DiGraph],
    algorithm: str = "louvain",
    resolution_parameter: float = 1.0,
    directed: bool = False,
) -> list[LayerCommunityFit]:
    algorithm = algorithm.lower()
    if algorithm not in {"louvain", "leiden"}:
        raise ValueError("`algorithm` must be one of {'louvain', 'leiden'}.")

    fits: list[LayerCommunityFit] = []

    for g in graph_layers:
        if algorithm == "louvain":
            try:
                import community as community_louvain
            except ImportError as exc:
                raise ImportError("Install optional dependency `python-louvain` for Louvain support.") from exc

            g_input = g.to_undirected() if directed else g
            partition = community_louvain.best_partition(g_input, weight="weight", resolution=resolution_parameter)
            communities = {}
            for node_zero, comm in partition.items():
                communities.setdefault(comm, []).append(node_zero + 1)

            mod = None
            if not directed:
                sets = [set(n - 1 for n in nodes) for nodes in communities.values()]
                mod = nx.algorithms.community.modularity(g_input, sets, weight="weight")

            fits.append(LayerCommunityFit(membership={k + 1: v + 1 for k, v in partition.items()}, modularity=mod, communities={k + 1: v for k, v in communities.items()}))
        else:
            try:
                import igraph as ig
                import leidenalg
            except ImportError as exc:
                raise ImportError("Install optional dependencies `python-igraph` and `leidenalg` for Leiden support.") from exc

            g_input = g
            edges = [(u, v) for u, v in g_input.edges()]
            weights = [g_input[u][v].get("weight", 1.0) for u, v in edges]
            ig_graph = ig.Graph(
                n=g_input.number_of_nodes(),
                edges=edges,
                directed=g_input.is_directed(),
            )
            if weights:
                ig_graph.es["weight"] = weights

            objective = leidenalg.CPMVertexPartition if directed else leidenalg.ModularityVertexPartition
            partition = leidenalg.find_partition(
                ig_graph,
                objective,
                weights=weights if weights else None,
                resolution_parameter=resolution_parameter,
            )
            membership = {idx + 1: comm + 1 for idx, comm in enumerate(partition.membership)}
            comms = {i + 1: [node + 1 for node in sorted(nodes)] for i, nodes in enumerate(partition)}

            mod = None if directed else ig_graph.modularity(partition.membership, weights=weights if weights else None)
            fits.append(LayerCommunityFit(membership=membership, modularity=mod, communities=comms))

    return fits


def weighted_jaccard(a: list[int], b: list[int]) -> float:
    inter = len(set(a).intersection(b))
    union = len(set(a).union(b))
    return 0.0 if union == 0 else inter / union


def weighted_overlap(a: list[int], b: list[int]) -> float:
    inter = len(set(a).intersection(b))
    min_size = min(len(a), len(b))
    return 0.0 if min_size == 0 else inter / min_size


def community_overlap_edges(
    fit: list[LayerCommunityFit],
    layer_links: pd.DataFrame,
    metric: str = "jaccard",
    min_similarity: float = 0.0,
) -> pd.DataFrame:
    if metric not in {"jaccard", "overlap"}:
        raise ValueError("`metric` must be one of {'jaccard', 'overlap'}.")

    sim_fun = weighted_jaccard if metric == "jaccard" else weighted_overlap
    rows: list[dict] = []

    for _, row in layer_links.iterrows():
        from_idx = int(row["from"])
        to_idx = int(row["to"])
        layer_weight = float(row["weight"])

        from_comms = fit[from_idx - 1].communities
        to_comms = fit[to_idx - 1].communities

        for from_c, from_nodes in from_comms.items():
            for to_c, to_nodes in to_comms.items():
                sim = sim_fun(from_nodes, to_nodes)
                weighted_sim = sim * layer_weight
                if weighted_sim >= min_similarity:
                    rows.append(
                        {
                            "from_layer": from_idx,
                            "to_layer": to_idx,
                            "from_community": from_c,
                            "to_community": to_c,
                            "similarity": sim,
                            "layer_weight": layer_weight,
                            "weighted_similarity": weighted_sim,
                        }
                    )

    return pd.DataFrame(rows)
