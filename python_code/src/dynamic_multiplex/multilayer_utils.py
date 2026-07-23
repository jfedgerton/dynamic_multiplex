from __future__ import annotations

from dataclasses import dataclass

import networkx as nx
import numpy as np
import pandas as pd


@dataclass
class LayerCommunityFit:
    membership: dict[int, int]
    modularity: float | None
    communities: dict[int, list[int]]


def _is_zero_indexed(graph: nx.Graph) -> bool:
    """Check if graph uses standard 0-indexed contiguous integer nodes."""
    return set(graph.nodes()) == set(range(graph.number_of_nodes()))


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


def fit_layer_communities(
    graph_layers: list[nx.Graph | nx.DiGraph],
    algorithm: str = "leiden",
    resolution_parameter: float = 1.0,
    directed: bool = False,
    objective: str | None = None,
) -> list[LayerCommunityFit]:
    algorithm = algorithm.lower()
    if algorithm not in {"louvain", "leiden"}:
        raise ValueError("`algorithm` must be one of {'louvain', 'leiden'}.")

    if objective is not None:
        objective = objective.lower()
        if objective not in {"modularity", "cpm"}:
            raise ValueError("`objective` must be one of {'modularity', 'cpm'}.")
        if algorithm == "louvain" and objective == "cpm":
            raise ValueError("Louvain does not support the CPM objective. Use algorithm='leiden'.")

    # Resolve effective objective: explicit choice, or default based on direction
    effective_objective = objective if objective is not None else ("cpm" if directed else "modularity")

    fits: list[LayerCommunityFit] = []

    if effective_objective == "modularity":
        for i, g in enumerate(graph_layers):
            weights = [d.get("weight", 1.0) for _, _, d in g.edges(data=True)]
            if any(w < 0 for w in weights):
                raise ValueError(
                    f"Layer {i + 1} contains negative edge weights. "
                    "Modularity-based methods do not support negative weights. "
                    'Use objective="cpm" to select the CPM objective, '
                    "which handles negative weights correctly."
                )

    for g in graph_layers:
        zero_indexed = _is_zero_indexed(g)

        if algorithm == "louvain":
            try:
                import community as community_louvain
            except ImportError as exc:
                raise ImportError("Install optional dependency `python-louvain` for Louvain support.") from exc

            g_input = g.to_undirected() if directed else g
            # Fixed detection seed so the partition is a deterministic function
            # of the graph: bootstrap variability then comes from the resampled
            # data, not from solver tie-breaking. Without this, best_partition
            # runs off unseeded global RNG and the pipeline is not reproducible
            # at a fixed bootstrap seed.
            partition = community_louvain.best_partition(
                g_input, weight="weight", resolution=resolution_parameter,
                random_state=123,
            )
            communities = {}
            for node, comm in partition.items():
                node_id = node + 1 if zero_indexed else node
                communities.setdefault(comm, []).append(node_id)

            mod = None
            if not directed and effective_objective != "cpm" and g_input.number_of_edges() > 0:
                sets = [set(n - 1 if zero_indexed else n for n in nodes) for nodes in communities.values()]
                mod = nx.algorithms.community.modularity(g_input, sets, weight="weight")

            membership = {
                (node + 1 if zero_indexed else node): comm + 1
                for node, comm in partition.items()
            }
            fits.append(LayerCommunityFit(membership=membership, modularity=mod, communities={k + 1: v for k, v in communities.items()}))
        else:
            try:
                import igraph as ig
                import leidenalg
            except ImportError as exc:
                raise ImportError("Install optional dependencies `python-igraph` and `leidenalg` for Leiden support.") from exc

            g_input = g
            node_list = sorted(g_input.nodes())
            node_to_idx = {n: i for i, n in enumerate(node_list)}

            edges = [(node_to_idx[u], node_to_idx[v]) for u, v in g_input.edges()]
            weights = [g_input[u][v].get("weight", 1.0) for u, v in g_input.edges()]
            ig_graph = ig.Graph(
                n=g_input.number_of_nodes(),
                edges=edges,
                directed=g_input.is_directed(),
            )
            if weights:
                ig_graph.es["weight"] = weights

            partition_type = leidenalg.CPMVertexPartition if effective_objective == "cpm" else leidenalg.RBConfigurationVertexPartition
            # Fixed detection seed for reproducibility (see the Louvain branch).
            partition = leidenalg.find_partition(
                ig_graph,
                partition_type,
                weights=weights if weights else None,
                resolution_parameter=resolution_parameter,
                seed=123,
            )

            membership = {}
            comms = {}
            for idx, comm in enumerate(partition.membership):
                node = node_list[idx]
                node_id = node + 1 if zero_indexed else node
                membership[node_id] = comm + 1
                comms.setdefault(comm + 1, []).append(node_id)
            for k in comms:
                comms[k] = sorted(comms[k])

            mod = None if directed or effective_objective == "cpm" or ig_graph.ecount() == 0 else ig_graph.modularity(partition.membership, weights=weights if weights else None)
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


def weighted_jaccard_similarity(a: list[int], b: list[int], weights_a: dict[int, float], weights_b: dict[int, float]) -> float:
    nodes = set(a).union(b)
    if not nodes:
        return 0.0

    inter = set(a).intersection(b)
    inter_weight = sum(min(weights_a.get(node, 0.0), weights_b.get(node, 0.0)) for node in inter)
    union_weight = sum(max(weights_a.get(node, 0.0), weights_b.get(node, 0.0)) for node in nodes)
    return 0.0 if union_weight == 0 else inter_weight / union_weight


def weighted_overlap_similarity(a: list[int], b: list[int], weights_a: dict[int, float], weights_b: dict[int, float]) -> float:
    inter = set(a).intersection(b)
    inter_weight = sum(min(weights_a.get(node, 0.0), weights_b.get(node, 0.0)) for node in inter)
    a_weight = sum(weights_a.get(node, 0.0) for node in a)
    b_weight = sum(weights_b.get(node, 0.0) for node in b)
    min_weight = min(a_weight, b_weight)
    return 0.0 if min_weight == 0 else inter_weight / min_weight


def layer_node_strengths(graph_layers: list[nx.Graph | nx.DiGraph], directed: bool = False) -> list[dict]:
    strengths: list[dict] = []
    for graph in graph_layers:
        zero_indexed = _is_zero_indexed(graph)
        if directed and isinstance(graph, nx.DiGraph):
            node_strength = {
                (node + 1 if zero_indexed else node): float(graph.in_degree(node, weight="weight") + graph.out_degree(node, weight="weight"))
                for node in graph.nodes()
            }
        else:
            node_strength = {
                (node + 1 if zero_indexed else node): float(graph.degree(node, weight="weight"))
                for node in graph.nodes()
            }
        strengths.append(node_strength)

    return strengths


def community_overlap_edges(
    fit: list[LayerCommunityFit],
    layer_links: pd.DataFrame,
    metric: str = "jaccard",
    min_similarity: float = 0.0,
    node_weights_by_layer: list[dict[int, float]] | None = None,
) -> pd.DataFrame:
    if metric not in {"jaccard", "overlap"}:
        raise ValueError("`metric` must be one of {'jaccard', 'overlap'}.")

    sim_fun = weighted_jaccard if metric == "jaccard" else weighted_overlap
    weighted_sim_fun = weighted_jaccard_similarity if metric == "jaccard" else weighted_overlap_similarity
    rows: list[dict] = []

    for _, row in layer_links.iterrows():
        from_idx = int(row["from"])
        to_idx = int(row["to"])
        layer_weight = float(row["weight"])

        from_comms = fit[from_idx - 1].communities
        to_comms = fit[to_idx - 1].communities

        for from_c, from_nodes in from_comms.items():
            for to_c, to_nodes in to_comms.items():
                if node_weights_by_layer is None:
                    sim = sim_fun(from_nodes, to_nodes)
                else:
                    sim = weighted_sim_fun(
                        from_nodes,
                        to_nodes,
                        node_weights_by_layer[from_idx - 1],
                        node_weights_by_layer[to_idx - 1],
                    )
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


def add_community_self_loops(
    edge_df: pd.DataFrame,
    fit: list[LayerCommunityFit],
    layer_links: pd.DataFrame,
    self_loop_multiplier: float = 1.0,
    min_similarity: float = 0.0,
    directed: bool = False,
) -> pd.DataFrame:
    rows: list[dict] = []

    # Undirected self-loops count each internal edge twice (A[i,j] + A[j,i])
    self_sim = 1.0 if directed else 2.0

    # Compute max layer_weight for each unique layer across all links
    layer_weights: dict[int, float] = {}
    for _, row in layer_links.iterrows():
        from_idx = int(row["from"])
        to_idx = int(row["to"])
        w = float(row["weight"])
        layer_weights[from_idx] = max(layer_weights.get(from_idx, 0.0), w)
        layer_weights[to_idx] = max(layer_weights.get(to_idx, 0.0), w)

    for layer_idx in sorted(layer_weights.keys()):
        layer_weight = layer_weights[layer_idx]
        comms = fit[layer_idx - 1].communities

        for comm_idx in comms.keys():
            weighted_sim = self_sim * layer_weight * self_loop_multiplier
            if weighted_sim >= min_similarity:
                rows.append(
                    {
                        "from_layer": layer_idx,
                        "to_layer": layer_idx,
                        "from_community": int(comm_idx),
                        "to_community": int(comm_idx),
                        "similarity": self_sim,
                        "layer_weight": layer_weight,
                        "weighted_similarity": weighted_sim,
                    }
                )

    if not rows:
        return edge_df

    return pd.concat([edge_df, pd.DataFrame(rows)], ignore_index=True)


def detect_interlayer_communities(
    layer_communities: list[LayerCommunityFit],
    interlayer_ties: pd.DataFrame,
    algorithm: str = "leiden",
    resolution_parameter: float = 1.0,
) -> dict:
    """Detect cross-layer (meta) communities from interlayer ties.

    Second-stage community detection. Treats each per-layer community as a
    node in a "community graph" whose edges are the interlayer similarity
    ties (plus the community self-loops), then runs community detection on
    that graph to group per-layer communities into cross-layer
    *meta-communities*. This is the step that makes custom ``layer_links``
    and the interlayer coupling actually affect the returned membership.

    Parameters
    ----------
    layer_communities : list[LayerCommunityFit]
        Per-layer detection results (each with ``membership`` and
        ``communities``).
    interlayer_ties : pandas.DataFrame
        Interlayer ties with columns ``from_layer``, ``to_layer``,
        ``from_community``, ``to_community``, ``weighted_similarity``
        (including self-loops).
    algorithm : str
        Second-stage algorithm: ``"louvain"`` or ``"leiden"`` (match the
        per-layer algorithm).
    resolution_parameter : float
        Resolution for the second-stage detection.

    Returns
    -------
    dict
        Dictionary with keys:
        - ``meta_ids``: dict mapping each supernode key ``"L<layer>C<community>"``
          to its meta-community id (or ``None`` in the fallback path).
        - ``membership``: list with one ``np.ndarray`` per layer giving each
          node's meta-community assignment (node order).
    """
    algorithm = algorithm.lower()
    n_layers = len(layer_communities)

    def key(layer: int, comm: int) -> str:
        return f"L{layer}C{comm}"

    # The second stage needs community-level ties. Node-level identity ties
    # have different columns (from_layer, to_layer, node, layer_weight) and are
    # not handled here; fall back to per-layer communities made globally
    # distinct (offset each layer's labels), i.e. no cross-layer merging.
    has_comm_ties = (
        interlayer_ties is not None
        and len(interlayer_ties) > 0
        and {"from_community", "to_community", "weighted_similarity"}.issubset(
            interlayer_ties.columns
        )
    )
    if not has_comm_ties:
        membership = []
        offset = 0
        for t in range(n_layers):
            mem_dict = layer_communities[t].membership
            node_ids = sorted(mem_dict.keys())
            mem = np.array([mem_dict[nid] for nid in node_ids], dtype=int)
            membership.append(mem + offset)
            offset += int(mem.max()) if mem.size else 0
        return {"meta_ids": None, "membership": membership}

    # Enumerate every per-layer community as a super-node (1-indexed layers).
    supernodes: list[str] = []
    seen = set()
    for t in range(n_layers):
        for comm_id in layer_communities[t].communities.keys():
            k = key(t + 1, int(comm_id))
            if k not in seen:
                seen.add(k)
                supernodes.append(k)

    # Build the community graph from interlayer ties (incl. self-loops).
    edges = []  # (from_key, to_key, weight)
    for _, row in interlayer_ties.iterrows():
        edges.append(
            (
                key(int(row["from_layer"]), int(row["from_community"])),
                key(int(row["to_layer"]), int(row["to_community"])),
                float(row["weighted_similarity"]),
            )
        )

    if len(edges) == 0:
        # No ties: every per-layer community is its own meta-community.
        meta_ids = {k: i + 1 for i, k in enumerate(supernodes)}
    else:
        if algorithm == "louvain":
            try:
                import community as community_louvain
                import networkx as nx
            except ImportError as exc:
                raise ImportError(
                    "Install optional dependency `python-louvain` for Louvain support."
                ) from exc

            cg = nx.Graph()
            cg.add_nodes_from(supernodes)
            for u, v, w in edges:
                if cg.has_edge(u, v):
                    cg[u][v]["weight"] += w
                else:
                    cg.add_edge(u, v, weight=w)
            partition = community_louvain.best_partition(
                cg, weight="weight", resolution=resolution_parameter,
                random_state=123,
            )
            meta_ids = {node: comm + 1 for node, comm in partition.items()}
        else:
            try:
                import igraph as ig
                import leidenalg
            except ImportError as exc:
                raise ImportError(
                    "Install optional dependencies `python-igraph` and `leidenalg` "
                    "for Leiden support."
                ) from exc

            node_to_idx = {k: i for i, k in enumerate(supernodes)}
            ig_edges = [(node_to_idx[u], node_to_idx[v]) for u, v, _ in edges]
            weights = [w for _, _, w in edges]
            ig_graph = ig.Graph(n=len(supernodes), edges=ig_edges, directed=False)
            ig_graph.es["weight"] = weights
            partition = leidenalg.find_partition(
                ig_graph,
                leidenalg.RBConfigurationVertexPartition,
                weights=weights,
                resolution_parameter=resolution_parameter,
                seed=123,
            )
            meta_ids = {
                supernodes[idx]: comm + 1
                for idx, comm in enumerate(partition.membership)
            }

    # Map each node to its meta-community, per layer (node order).
    membership = []
    for t in range(n_layers):
        mem_dict = layer_communities[t].membership
        node_ids = sorted(mem_dict.keys())
        mem = np.array(
            [meta_ids[key(t + 1, int(mem_dict[nid]))] for nid in node_ids],
            dtype=int,
        )
        membership.append(mem)

    return {"meta_ids": meta_ids, "membership": membership}


def detect_multislice_communities(
    graph_layers: list[nx.Graph | nx.DiGraph],
    interlayer_ties: pd.DataFrame,
    algorithm: str = "leiden",
) -> list[np.ndarray]:
    """Multislice (Mucha) meta-communities for node-identity coupling.

    Node-level second stage for the identity specification. Builds a single
    supra-graph by stacking the layers (intra-layer edges are each layer's own
    adjacency) and adding interlayer identity edges (each node tied to its own
    copies in the coupled layers, weighted by the layer-link weight), then runs
    one community detection on the whole supra-graph. This is Mucha et al.
    (2010) multislice modularity with the coupling given by ``layer_links``: a
    node's meta-community can be pulled across layers through the identity ties.

    Parameters
    ----------
    graph_layers : list[networkx.Graph | networkx.DiGraph]
        Per-layer graphs, as produced by ``prepare_multilayer_graphs``.
    interlayer_ties : pandas.DataFrame
        Node-level identity ties (columns ``from_layer``, ``to_layer``,
        ``node``, ``layer_weight``).
    algorithm : str
        ``"louvain"`` or ``"leiden"`` for the supra-graph detection.

    Returns
    -------
    list[np.ndarray]
        One array per layer giving each node's meta-community assignment
        (node order).
    """
    algorithm = algorithm.lower()
    n_layers = len(graph_layers)

    def vkey(layer: int, node) -> str:
        return f"L{layer}N{node}"

    # Per-layer canonical node ids: match the ids used in the per-layer
    # membership dict (node + 1 when the layer is zero-indexed, else the node
    # name). The identity ties' `node` values are built the same way in
    # fit_multilayer_identity_ties, so intra and inter keys line up.
    layer_nodes = []  # per layer: list of (graph_node, canonical_id) in node order
    for g in graph_layers:
        zero = _is_zero_indexed(g)
        nodes_sorted = sorted(g.nodes())
        layer_nodes.append(
            [(n, (n + 1 if zero else n)) for n in nodes_sorted]
        )

    # Supra vertices: one per (layer, node), layers 1-indexed.
    supra: list[str] = []
    seen = set()
    for t in range(n_layers):
        for _, cid in layer_nodes[t]:
            k = vkey(t + 1, cid)
            if k not in seen:
                seen.add(k)
                supra.append(k)

    edges = []  # (from_key, to_key, weight)

    # Intra-layer edges: each layer's own adjacency (original network).
    for t in range(n_layers):
        g = graph_layers[t]
        cid_map = {n: cid for n, cid in layer_nodes[t]}
        for u, v, d in g.edges(data=True):
            edges.append(
                (
                    vkey(t + 1, cid_map[u]),
                    vkey(t + 1, cid_map[v]),
                    float(d.get("weight", 1.0)),
                )
            )

    # Interlayer edges: identity ties (node to its own copy), weighted by omega.
    # Access columns directly (not iterrows, which upcasts a mixed-dtype row to
    # a single dtype and would turn integer node ids into floats).
    if interlayer_ties is not None and len(interlayer_ties) > 0:
        from_layer_col = interlayer_ties["from_layer"].tolist()
        to_layer_col = interlayer_ties["to_layer"].tolist()
        node_col = interlayer_ties["node"].tolist()
        if "layer_weight" in interlayer_ties.columns:
            weight_col = interlayer_ties["layer_weight"].tolist()
        else:
            weight_col = [1.0] * len(node_col)
        for fl, tl, nd, w in zip(from_layer_col, to_layer_col, node_col, weight_col):
            edges.append(
                (vkey(int(fl), nd), vkey(int(tl), nd), float(w))
            )

    # Single detection on the supra-graph.
    if len(edges) == 0:
        meta = {k: i + 1 for i, k in enumerate(supra)}
    elif algorithm == "louvain":
        try:
            import community as community_louvain
        except ImportError as exc:
            raise ImportError(
                "Install optional dependency `python-louvain` for Louvain support."
            ) from exc

        cg = nx.Graph()
        cg.add_nodes_from(supra)
        for u, v, w in edges:
            if cg.has_edge(u, v):
                cg[u][v]["weight"] += w
            else:
                cg.add_edge(u, v, weight=w)
        partition = community_louvain.best_partition(
            cg, weight="weight", random_state=123
        )
        meta = {node: comm + 1 for node, comm in partition.items()}
    else:
        try:
            import igraph as ig
            import leidenalg
        except ImportError as exc:
            raise ImportError(
                "Install optional dependencies `python-igraph` and `leidenalg` "
                "for Leiden support."
            ) from exc

        node_to_idx = {k: i for i, k in enumerate(supra)}
        ig_edges = [(node_to_idx[u], node_to_idx[v]) for u, v, _ in edges]
        weights = [w for _, _, w in edges]
        ig_graph = ig.Graph(n=len(supra), edges=ig_edges, directed=False)
        ig_graph.es["weight"] = weights
        # n_iterations=3 mirrors the R cluster_leiden call and gives the
        # supra-graph optimizer enough refinement passes to avoid collapsing a
        # whole slice into a single community.
        partition = leidenalg.find_partition(
            ig_graph,
            leidenalg.RBConfigurationVertexPartition,
            weights=weights,
            seed=123,
            n_iterations=3,
        )
        meta = {
            supra[idx]: comm + 1 for idx, comm in enumerate(partition.membership)
        }

    # Map back to per-layer node order.
    membership = []
    for t in range(n_layers):
        arr = np.array(
            [meta[vkey(t + 1, cid)] for _, cid in layer_nodes[t]],
            dtype=int,
        )
        membership.append(arr)
    return membership
