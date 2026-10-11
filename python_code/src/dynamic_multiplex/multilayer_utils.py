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


def _id_offset(graph: nx.Graph) -> int:
    """Offset from a graph's node labels to the node ids used in memberships.

    Layers that share one node universe are relabelled to 0..n-1 by
    ``prepare_multilayer_graphs`` and reported as 1..n (offset 1). Layers with
    different node sets keep their own labels (offset 0 for every layer), so a
    node has the same id in every layer and ids are never computed from label
    arithmetic.
    """
    if graph.graph.get("dm_raw_ids", False):
        return 0
    return 1 if _is_zero_indexed(graph) else 0


def _sort_key(node):
    """Deterministic sort key for node labels of mixed types."""
    return (type(node).__name__, node)


def _clean_graph(graph: nx.Graph | nx.DiGraph) -> nx.Graph | nx.DiGraph:
    """Copy of ``graph`` without self-loops or zero-weight edges.

    Same input policy as the R package: within-layer self-loops (including the
    diagonal of an adjacency matrix) and zero-weight edges are dropped.
    """
    out = graph.copy()
    out.remove_edges_from(list(nx.selfloop_edges(out)))
    zero = [(u, v) for u, v, d in out.edges(data=True) if d.get("weight", 1.0) == 0]
    out.remove_edges_from(zero)
    return out


def _as_graph(layer, directed: bool) -> nx.Graph | nx.DiGraph:
    if isinstance(layer, (nx.Graph, nx.DiGraph)):
        return layer

    matrix = np.asarray(layer)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("Each layer must be a square adjacency matrix or a NetworkX graph.")

    graph_type = nx.DiGraph if directed else nx.Graph
    graph = nx.from_numpy_array(matrix, create_using=graph_type)
    return graph


def _collapse_directed(graph: nx.DiGraph) -> nx.Graph:
    """Undirected copy of a directed graph with reciprocal arc weights summed.

    Matches R's ``as_undirected(mode = "collapse", weight = "sum")``;
    ``nx.DiGraph.to_undirected()`` would keep only one of the two weights.
    """
    out = nx.Graph()
    out.add_nodes_from(graph.nodes())
    for u, v, d in graph.edges(data=True):
        w = d.get("weight", 1.0)
        if out.has_edge(u, v):
            out[u][v]["weight"] += w
        else:
            out.add_edge(u, v, weight=w)
    return out


def prepare_multilayer_graphs(
    layers: list,
    directed: bool = False,
    require_same_nodes: bool = True,
) -> list[nx.Graph | nx.DiGraph]:
    if not isinstance(layers, list) or len(layers) < 2:
        raise ValueError("`layers` must be a list with at least two network layers.")

    is_matrix = [not isinstance(layer, (nx.Graph, nx.DiGraph)) for layer in layers]
    graphs = [_clean_graph(_as_graph(layer, directed=directed)) for layer in layers]

    for i, g in enumerate(graphs):
        if g.number_of_nodes() == 0:
            raise ValueError(f"Layer {i + 1} has no vertices. Drop empty layers before fitting.")

    # Validate a shared node universe across all layers. Layers with
    # different node sets (or adjacency matrices of different sizes) cannot
    # be coupled unless the caller opts in: interlayer ties and
    # meta-communities then match nodes across layers by label.
    node_sets = [set(g.nodes()) for g in graphs]
    aligned = all(ns == node_sets[0] for ns in node_sets[1:])
    # Nodes can enter and leave while layer sizes stay equal, so with
    # allow_unequal_nodes (require_same_nodes=False) every layer must carry
    # node labels; positions in an adjacency matrix carry no identity (1.4.0:
    # checked whatever the layer sizes).
    if not require_same_nodes and any(is_matrix):
        raise ValueError(
            "allow_unequal_nodes=True requires networkx graphs with node "
            "labels on every layer so that nodes can be matched across "
            "layers; adjacency matrices have no labels."
        )
    if not aligned:
        if require_same_nodes:
            raise ValueError(
                "All layers must share the same node set. "
                "Found layers with different nodes; align the node universe "
                "(adding isolates where needed) before fitting, or pass "
                "allow_unequal_nodes=True with labelled networkx graphs."
            )
        # Positions in adjacency matrices of different sizes carry no
        # identity, so they cannot be matched across layers (as in R).
        if any(is_matrix):
            raise ValueError(
                "allow_unequal_nodes=True requires networkx graphs with node "
                "labels on every layer so that nodes can be matched across "
                "layers; adjacency matrices of different sizes have no labels."
            )
        # Keep raw labels as node ids in every layer (no +1 offset), so a node
        # has one id across all layers.
        for g in graphs:
            g.graph["dm_raw_ids"] = True
        return graphs

    # Normalize arbitrary node labels (strings, non-contiguous integers) to
    # canonical 0..n-1 integers for internal processing, eliminating numeric
    # +1 assumptions on raw labels. The original labels are preserved on each
    # graph as ``graph.graph["node_labels"]`` (index i -> original label of
    # internal node i, i.e. of node id i+1 in the returned membership).
    reference_nodes = sorted(node_sets[0], key=_sort_key)
    if reference_nodes != list(range(len(reference_nodes))):
        mapping = {node: i for i, node in enumerate(reference_nodes)}
        graphs = [nx.relabel_nodes(g, mapping, copy=True) for g in graphs]
        for g in graphs:
            g.graph["node_labels"] = list(reference_nodes)

    return graphs


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
    seed: int | None = 123,
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
        off = _id_offset(g)

        if algorithm == "louvain":
            try:
                import community as community_louvain
            except ImportError as exc:
                raise ImportError("Install optional dependency `python-louvain` for Louvain support.") from exc

            # python-louvain is undirected: collapse directed layers
            g_input = _collapse_directed(g) if g.is_directed() else g
            # Seeded detection (default 123) so the partition is a deterministic
            # function of the graph: bootstrap variability then comes from the
            # resampled data, not from solver tie-breaking. Pass seed=None for
            # unseeded (non-reproducible) detection.
            partition = community_louvain.best_partition(
                g_input, weight="weight", resolution=resolution_parameter,
                random_state=seed,
            )
            communities = {}
            for node, comm in partition.items():
                node_id = node + off if off else node
                communities.setdefault(comm, []).append(node_id)

            mod = None
            if not directed and effective_objective != "cpm" and g_input.number_of_edges() > 0:
                sets = [set(n - off if off else n for n in nodes) for nodes in communities.values()]
                mod = nx.algorithms.community.modularity(g_input, sets, weight="weight")

            membership = {(node + off if off else node): comm + 1 for node, comm in partition.items()}
            fits.append(LayerCommunityFit(membership=membership, modularity=mod, communities={k + 1: v for k, v in communities.items()}))
        else:
            try:
                import igraph as ig
                import leidenalg
            except ImportError as exc:
                raise ImportError("Install optional dependencies `python-igraph` and `leidenalg` for Leiden support.") from exc

            g_input = g
            node_list = sorted(g_input.nodes(), key=_sort_key)
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
            # Seeded detection for reproducibility (see the Louvain branch).
            partition = leidenalg.find_partition(
                ig_graph,
                partition_type,
                weights=weights if weights else None,
                resolution_parameter=resolution_parameter,
                seed=seed,
            )

            membership = {}
            comms = {}
            for idx, comm in enumerate(partition.membership):
                node = node_list[idx]
                node_id = node + off if off else node
                membership[node_id] = comm + 1
                comms.setdefault(comm + 1, []).append(node_id)
            for k in comms:
                comms[k] = sorted(comms[k], key=_sort_key)

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
    # weights_a / weights_b are LAYER-WIDE node strengths: a node's strength in
    # layer a counts toward community a only when the node is a MEMBER of a
    # (likewise for b). Nodes outside a community contribute 0 on that side.
    set_a = set(a)
    set_b = set(b)
    nodes = set_a | set_b
    if not nodes:
        return 0.0

    def wa(node):
        return weights_a.get(node, 0.0) if node in set_a else 0.0

    def wb(node):
        return weights_b.get(node, 0.0) if node in set_b else 0.0

    inter_weight = sum(min(wa(node), wb(node)) for node in set_a & set_b)
    union_weight = sum(max(wa(node), wb(node)) for node in nodes)
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
        off = _id_offset(graph)
        if directed and isinstance(graph, nx.DiGraph):
            node_strength = {
                (node + off if off else node): float(graph.in_degree(node, weight="weight") + graph.out_degree(node, weight="weight"))
                for node in graph.nodes()
            }
        else:
            node_strength = {
                (node + off if off else node): float(graph.degree(node, weight="weight"))
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
            from_set = set(from_nodes)
            for to_c, to_nodes in to_comms.items():
                # communities with no shared node have zero similarity under
                # every metric; they are not stored (as in R) ----
                if from_set.isdisjoint(to_nodes):
                    continue
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
    seed: int | None = 123,
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
            node_ids = sorted(mem_dict.keys(), key=_sort_key)
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
                random_state=seed,
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
                seed=seed,
            )
            meta_ids = {
                supernodes[idx]: comm + 1
                for idx, comm in enumerate(partition.membership)
            }

    # Map each node to its meta-community, per layer (node order).
    membership = []
    for t in range(n_layers):
        mem_dict = layer_communities[t].membership
        node_ids = sorted(mem_dict.keys(), key=_sort_key)
        mem = np.array(
            [meta_ids[key(t + 1, int(mem_dict[nid]))] for nid in node_ids],
            dtype=int,
        )
        membership.append(mem)

    return {"meta_ids": meta_ids, "membership": membership}


# ---------------------------------------------------------------------------
# 1.4.0: chain repair and filtered (past-only) tracking for the two-stage fits
# ---------------------------------------------------------------------------

def _parse_key(k: str) -> tuple[int, int]:
    layer, comm = k[1:].split("C")
    return int(layer), int(comm)


def _membership_from_meta(layer_communities: list[LayerCommunityFit], meta_ids: dict) -> list[np.ndarray]:
    """Per-layer node membership (sorted node-id order) from supernode meta ids."""
    membership = []
    for t, fit in enumerate(layer_communities, start=1):
        mem_dict = fit.membership
        node_ids = sorted(mem_dict.keys(), key=_sort_key)
        membership.append(np.array([meta_ids[f"L{t}C{int(mem_dict[n])}"] for n in node_ids], dtype=int))
    return membership


def _renumber_first_appearance(meta_ids: dict) -> dict:
    """Relabel meta ids 1..K in order of first appearance (layer, then community id)."""
    new: dict = {}
    out = {}
    for k in sorted(meta_ids, key=_parse_key):
        new.setdefault(meta_ids[k], len(new) + 1)
        out[k] = new[meta_ids[k]]
    return out


def _best_matches(pair_ties: pd.DataFrame):
    """Best match of every supernode into every linked layer, by similarity.

    Returns ``best[(layer, comm, other_layer)] = (similarity, (other_layer, other_comm))``.
    Ties are broken by the lower community id.
    """
    best: dict = {}
    for r in pair_ties.itertuples(index=False):
        a = (int(r.from_layer), int(r.from_community))
        b = (int(r.to_layer), int(r.to_community))
        if a[0] == b[0]:
            continue
        sim = float(r.similarity)
        for x, y in ((a, b), (b, a)):
            cur = best.get((x[0], x[1], y[0]))
            if cur is None or sim > cur[0] or (sim == cur[0] and y[1] < cur[1][1]):
                best[(x[0], x[1], y[0])] = (sim, y)
    return best


def _check_repair_args(tracking: str, repair_chains, repair_min_similarity) -> None:
    if tracking not in ("smoothed", "filtered"):
        raise ValueError('`tracking` must be "smoothed" or "filtered".')
    if not isinstance(repair_chains, (bool, np.bool_)):
        raise ValueError("`repair_chains` must be True or False.")
    if not (0.0 <= float(repair_min_similarity) <= 1.0):
        raise ValueError("`repair_min_similarity` must be between 0 and 1.")


_REPAIR_COLUMNS = ["from_layer", "from_community", "to_layer", "to_community", "similarity"]


def repair_meta_chains(
    meta_ids: dict,
    pair_ties: pd.DataFrame,
    layer_links: pd.DataFrame,
    min_similarity: float = 0.2,
) -> tuple[dict, pd.DataFrame]:
    """Rejoin meta-communities that stage 2 cut along time (chain repair).

    Stage-2 modularity tends to cut a long-lived community into time segments
    (a resolution-limit effect on chains of communities). For each layer link
    (s, t) with s < t, in increasing order, and each community a in layer s
    (by community id), the meta-community of community b in layer t is joined
    to the meta-community of a when all of the following hold:

    * a and b are each other's best match across (s, t) by the fit's own
      similarity (mutual best match, as in Kawadia and Sreenivasan 2012);
    * their similarity is at least ``min_similarity`` (a matching floor, as in
      Greene et al. 2010);
    * a's meta-community has no community after layer s and b's has none
      before layer t, so two groups that coexist are never merged.

    Returns the repaired meta ids, renumbered 1..K by first appearance, and a
    log of the joins.
    """
    root = {m: m for m in set(meta_ids.values())}
    layers_of: dict = {}
    for k, m in meta_ids.items():
        layers_of.setdefault(m, set()).add(_parse_key(k)[0])

    def find(m):
        while root[m] != m:
            root[m] = root[root[m]]
            m = root[m]
        return m

    best = _best_matches(pair_ties)
    by_layer: dict = {}
    for k in meta_ids:
        layer, comm = _parse_key(k)
        by_layer.setdefault(layer, []).append(comm)
    link_pairs = sorted({(min(int(f), int(t)), max(int(f), int(t)))
                         for f, t in zip(layer_links["from"], layer_links["to"]) if int(f) != int(t)})
    log = []
    for s, t in link_pairs:
        for ca in sorted(by_layer.get(s, [])):
            fwd = best.get((s, ca, t))
            if fwd is None:
                continue
            sim, (_, cb) = fwd
            bwd = best.get((t, cb, s))
            if bwd is None or bwd[1] != (s, ca) or sim < min_similarity:
                continue
            ra, rb = find(meta_ids[f"L{s}C{ca}"]), find(meta_ids[f"L{t}C{cb}"])
            if ra == rb:
                continue
            if max(layers_of[ra]) == s and min(layers_of[rb]) == t:
                root[rb] = ra
                layers_of[ra] |= layers_of.pop(rb)
                log.append({"from_layer": s, "from_community": ca, "to_layer": t,
                            "to_community": cb, "similarity": sim})
    repaired = {k: find(m) for k, m in meta_ids.items()}
    return _renumber_first_appearance(repaired), pd.DataFrame(log, columns=_REPAIR_COLUMNS)


def _cluster_weighted(nodes: list, edges: dict, algorithm: str, resolution_parameter: float, seed):
    """Stage-2 solver call (same library calls as ``detect_interlayer_communities``)."""
    if algorithm == "louvain":
        import community as community_louvain
        g = nx.Graph()
        g.add_nodes_from(nodes)
        for (u, v), w in edges.items():
            g.add_edge(u, v, weight=w)
        part = community_louvain.best_partition(g, weight="weight", resolution=resolution_parameter,
                                                random_state=seed)
        return {n: part[n] for n in nodes}
    import igraph as ig
    import leidenalg
    idx = {n: i for i, n in enumerate(nodes)}
    w = list(edges.values())
    g = ig.Graph(n=len(nodes), edges=[(idx[u], idx[v]) for (u, v) in edges], directed=False)
    g.es["weight"] = w
    part = leidenalg.find_partition(g, leidenalg.RBConfigurationVertexPartition, weights=w,
                                    resolution_parameter=resolution_parameter, seed=seed)
    return {n: part.membership[idx[n]] for n in nodes}


def detect_interlayer_communities_filtered(
    layer_communities: list[LayerCommunityFit],
    pair_ties: pd.DataFrame,
    layer_links: pd.DataFrame,
    algorithm: str = "leiden",
    resolution_parameter: float = 1.0,
    seed: int | None = 123,
    add_self_loops: bool = True,
    self_loop_multiplier: float = 1.0,
    min_similarity: float = 0.0,
    directed: bool = False,
    repair_chains: bool = False,
    repair_min_similarity: float = 0.2,
) -> dict:
    """Past-only (filtered) second stage: the meta id of a layer-t community
    depends only on layers 1..t.

    At each layer t, the package's community graph is built on layers 1..t from
    the links whose ends are both at or before t. Every meta-community already
    assigned in layers before t is contracted into one read-only anchor node
    (contraction keeps degrees and internal weights, so the modularity of any
    partition that keeps those groups together is unchanged) and the graph is
    clustered with the same Leiden or Louvain call as the smoothed stage 2. A
    layer-t community takes the id of the anchor in its cluster with the
    strongest direct tie to it (then the strongest tie to the cluster's other
    layer-t communities, then the lower id); a cluster without an anchor
    starts a new id. Anchors are never relabelled. With ``repair_chains``, a
    layer-t community that started a new id is then joined to a group that
    ended at a linked earlier layer s under the rule of ``repair_meta_chains``,
    using layers s and t only.
    """
    algorithm = algorithm.lower()
    n_layers = len(layer_communities)
    meta: dict = {}
    for c in sorted(int(c) for c in layer_communities[0].communities):
        meta[f"L1C{c}"] = len(meta) + 1
    next_id = len(meta) + 1
    last_layer = {m: 1 for m in meta.values()}
    log = []
    best = _best_matches(pair_ties) if repair_chains and len(pair_ties) else {}

    for t in range(2, n_layers + 1):
        lk = layer_links[(layer_links["from"] <= t) & (layer_links["to"] <= t)].reset_index(drop=True)
        if len(pair_ties):
            pt = pair_ties[(pair_ties["from_layer"] <= t) & (pair_ties["to_layer"] <= t)]
        else:
            pt = pair_ties
        ties = pt
        if add_self_loops and len(lk):
            ties = add_community_self_loops(edge_df=pt, fit=layer_communities[:t], layer_links=lk,
                                            self_loop_multiplier=self_loop_multiplier,
                                            min_similarity=min_similarity, directed=directed)
        cur = [f"L{t}C{c}" for c in sorted(int(c) for c in layer_communities[t - 1].communities)]
        start_id = next_id

        def node_of(k):
            return ("A", meta[k]) if _parse_key(k)[0] < t else ("C", k)

        anchors = sorted({("A", m) for m in meta.values()})
        nodes = anchors + [("C", k) for k in cur]
        order = {n: i for i, n in enumerate(nodes)}
        edges: dict = {}
        for r in ties.itertuples(index=False) if len(ties) else []:
            u = node_of(f"L{int(r.from_layer)}C{int(r.from_community)}")
            v = node_of(f"L{int(r.to_layer)}C{int(r.to_community)}")
            key_uv = (u, v) if order[u] <= order[v] else (v, u)
            edges[key_uv] = edges.get(key_uv, 0.0) + float(r.weighted_similarity)

        if not edges:
            for k in cur:
                meta[k] = next_id
                next_id += 1
        else:
            mem = _cluster_weighted(nodes, edges, algorithm, resolution_parameter, seed)

            def w_between(a, b):
                return edges.get((a, b) if order[a] <= order[b] else (b, a), 0.0)

            clusters: dict = {}
            for n in nodes:
                clusters.setdefault(mem[n], []).append(n)
            groups = []
            for members in clusters.values():
                frees = [n for n in members if n[0] == "C"]
                if frees:
                    groups.append((min(_parse_key(n[1]) for n in frees), [n for n in members if n[0] == "A"], frees))
            for _, anc, frees in sorted(groups, key=lambda g: g[0]):
                if not anc:
                    for n in frees:
                        meta[n[1]] = next_id
                    next_id += 1
                    continue
                for n in frees:
                    best_anchor = max(anc, key=lambda a: (w_between(a, n), sum(w_between(a, f) for f in frees), -a[1]))
                    meta[n[1]] = best_anchor[1]

        if repair_chains:
            link_pairs = sorted({(min(int(f), int(tt)), max(int(f), int(tt)))
                                 for f, tt in zip(lk["from"], lk["to"]) if int(f) != int(tt) and max(int(f), int(tt)) == t})
            for s, _ in link_pairs:
                for k in cur:
                    if meta[k] < start_id:                     # only groups that start at t
                        continue
                    cb = _parse_key(k)[1]
                    bwd = best.get((t, cb, s))
                    if bwd is None:
                        continue
                    sim, (_, ca) = bwd
                    fwd = best.get((s, ca, t))
                    if fwd is None or fwd[1] != (t, cb) or sim < repair_min_similarity:
                        continue
                    ma = meta[f"L{s}C{ca}"]
                    # a's group must have ended at s: last seen at s and absent from t
                    if last_layer.get(ma) != s or any(meta[kk] == ma for kk in cur):
                        continue
                    mb = meta[k]
                    for kk in cur:
                        if meta[kk] == mb:
                            meta[kk] = ma
                    log.append({"from_layer": s, "from_community": ca, "to_layer": t,
                                "to_community": cb, "similarity": sim})
        for k in cur:
            last_layer[meta[k]] = t

    meta = _renumber_first_appearance(meta)
    return {"meta_ids": meta, "membership": _membership_from_meta(layer_communities, meta),
            "chain_repairs": pd.DataFrame(log, columns=_REPAIR_COLUMNS)}


def run_two_stage(
    fit: list[LayerCommunityFit],
    pair_ties: pd.DataFrame,
    links: pd.DataFrame,
    *,
    add_self_loops: bool,
    self_loop_multiplier: float,
    min_similarity: float,
    directed: bool,
    algorithm: str,
    resolution_parameter: float,
    seed,
    tracking: str,
    repair_chains: bool,
    repair_min_similarity: float,
) -> dict:
    """Second stage shared by the four two-stage fitters.

    ``tracking="smoothed"`` with ``repair_chains=False`` is the 1.3.2 code path,
    unchanged.
    """
    _check_repair_args(tracking, repair_chains, repair_min_similarity)
    interlayer_ties = pair_ties
    if add_self_loops:
        interlayer_ties = add_community_self_loops(
            edge_df=pair_ties, fit=fit, layer_links=links, self_loop_multiplier=self_loop_multiplier,
            min_similarity=min_similarity, directed=directed,
        )
    if tracking == "smoothed":
        meta = detect_interlayer_communities(
            layer_communities=fit, interlayer_ties=interlayer_ties, algorithm=algorithm,
            resolution_parameter=resolution_parameter, seed=seed,
        )
        repairs = pd.DataFrame(columns=_REPAIR_COLUMNS)
        if repair_chains and meta["meta_ids"] is not None:
            ids, repairs = repair_meta_chains(meta["meta_ids"], pair_ties, links, repair_min_similarity)
            meta = {"meta_ids": ids, "membership": _membership_from_meta(fit, ids)}
    else:
        meta = detect_interlayer_communities_filtered(
            fit, pair_ties, links, algorithm=algorithm, resolution_parameter=resolution_parameter, seed=seed,
            add_self_loops=add_self_loops, self_loop_multiplier=self_loop_multiplier,
            min_similarity=min_similarity, directed=directed, repair_chains=repair_chains,
            repair_min_similarity=repair_min_similarity,
        )
        repairs = meta["chain_repairs"]
    return {"meta_ids": meta["meta_ids"], "membership": meta["membership"],
            "interlayer_ties": interlayer_ties, "chain_repairs": repairs,
            "tracking": tracking, "repair_chains": bool(repair_chains),
            "repair_min_similarity": float(repair_min_similarity)}


def detect_multislice_communities(
    graph_layers: list[nx.Graph | nx.DiGraph],
    interlayer_ties: pd.DataFrame,
    algorithm: str = "leiden",
    omega: float = 1.0,
    resolution_parameter: float = 1.0,
    seed: int | None = 123,
    directed: bool = False,
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
        Kept for API compatibility; the multislice objective is always
        optimised with leidenalg's multiplex optimiser (Leiden moves).
    omega : float
        Interlayer coupling strength (Mucha's omega). Multiplies the
        interlayer identity-edge weights on top of any ``layer_links`` weights.
        Larger ``omega`` couples layers more strongly and, past a point,
        collapses everything into one meta-community; smaller ``omega``
        decouples toward independent per-layer detection.
    resolution_parameter : float
        Resolution for the supra-graph detection (Mucha's modularity
        resolution). Larger values yield more, smaller meta-communities.
    directed : bool
        If ``True``, every layer must be a ``networkx.DiGraph`` and each slice
        keeps its arc directions: the slice term uses Leicht and Newman's
        directed null model ``k_out_i * k_in_j / m_s``. A directed network in
        which every arc is reciprocated then gives the same fit as its
        undirected twin at the same ``omega``. The interlayer identity ties
        stay undirected (Mucha's coupling depends only on ``C_isr + C_irs``).
        If ``False`` (default), directed layers are symmetrised (``A + A^T``,
        so a reciprocated pair counts twice) with a warning, as in 1.3.1.

    Returns
    -------
    list[np.ndarray]
        One array per layer giving each node's meta-community assignment
        (node order).
    """
    algorithm = algorithm.lower()
    n_layers = len(graph_layers)

    is_dir = [g.is_directed() for g in graph_layers]
    if directed and not all(is_dir):
        raise ValueError(
            "directed=True requires every layer to be a directed graph "
            f"(networkx.DiGraph); layers {[i + 1 for i, d in enumerate(is_dir) if not d]} "
            "are undirected. Pass directed adjacency matrices or DiGraphs, or "
            "use directed=False."
        )
    if not directed and any(is_dir):
        import warnings

        warnings.warn(
            "Directed layers were passed with directed=False: the multislice "
            "fit symmetrised them (A + A^T, a reciprocated pair counts twice). "
            "Use directed=True for the directed null model.",
            UserWarning,
            stacklevel=2,
        )

    def vkey(layer: int, node) -> tuple:
        # tuple, not a formatted string: labels 1 and "1" must stay distinct
        # (networkx itself treats them as different nodes)
        return (layer, node)

    # Per-layer canonical node ids: match the ids used in the per-layer
    # membership dict (see ``_id_offset``). The identity ties' `node` values are built the same way in
    # fit_multilayer_identity_ties, so intra and inter keys line up.
    layer_nodes = []  # per layer: list of (graph_node, canonical_id) in node order
    for g in graph_layers:
        off = _id_offset(g)
        nodes_sorted = sorted(g.nodes(), key=_sort_key)
        layer_nodes.append([(n, n + off if off else n) for n in nodes_sorted])

    # Supra vertices: one per (layer, node), layers 1-indexed.
    supra: list[tuple] = []
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
                (vkey(int(fl), nd), vkey(int(tl), nd), float(w) * omega)
            )

    # Mucha et al. (2010) multislice modularity, optimised with leidenalg's
    # multiplex optimiser: one RBConfiguration (modularity) partition per slice
    # -- each slice keeps its OWN null model k_is k_js / 2m_s -- plus a CPM
    # layer holding the interlayer identity ties with resolution 0 (no null
    # term). Plain modularity on the stacked supra-graph is NOT the same
    # objective: its pooled null model is ~T times too small and penalises
    # cross-slice co-membership, so at omega <= 1 it returns each slice as one
    # community (the behaviour of this function before version 1.2.1).
    if len(edges) == 0:
        meta = {k: i + 1 for i, k in enumerate(supra)}
    else:
        try:
            import igraph as ig
            import leidenalg
        except ImportError as exc:
            raise ImportError(
                "Install optional dependencies `python-igraph` and `leidenalg` "
                "for multislice detection."
            ) from exc

        node_to_idx = {k: i for i, k in enumerate(supra)}
        supra_slice = np.array(
            [t for t in range(n_layers) for _ in layer_nodes[t]], dtype=int
        )
        intra_edges: list[tuple[int, int]] = []
        intra_w: list[float] = []
        inter_edges: list[tuple[int, int]] = []
        inter_w: list[float] = []
        for u, v, w in edges:
            iu, iv = node_to_idx[u], node_to_idx[v]
            if iu == iv:
                continue
            if supra_slice[iu] == supra_slice[iv]:
                intra_edges.append((iu, iv))
                intra_w.append(w)
            else:
                inter_edges.append((iu, iv))
                inter_w.append(w)

        # slice layers: every supra vertex is present in every layer graph, but
        # only its own slice's edges carry weight, so each slice has its own
        # degree sequence and 2m_s in its RBConfiguration null model ----
        n_supra = len(supra)
        partitions = []
        for t in range(n_layers):
            keep = [i for i, (u, v) in enumerate(intra_edges) if supra_slice[u] == t]
            # directed slices keep arc direction, so RBConfiguration uses
            # the k_out k_in / m_s null model (Leicht & Newman 2008) ----
            g_t = ig.Graph(
                n=n_supra, edges=[intra_edges[i] for i in keep], directed=directed
            )
            g_t.es["weight"] = [intra_w[i] for i in keep]
            partitions.append(
                leidenalg.RBConfigurationVertexPartition(
                    g_t, weights="weight", resolution_parameter=resolution_parameter
                )
            )
        # identity ties stay undirected for directed and undirected slices ----
        g_inter = ig.Graph(n=n_supra, edges=inter_edges, directed=False)
        assert not g_inter.is_directed()
        g_inter.es["weight"] = inter_w
        partitions.append(
            leidenalg.CPMVertexPartition(
                g_inter, weights="weight", resolution_parameter=0.0
            )
        )
        optimiser = leidenalg.Optimiser()
        if seed is not None:
            optimiser.set_rng_seed(int(seed))
        optimiser.optimise_partition_multiplex(partitions, n_iterations=-1)
        memb = partitions[0].membership
        meta = {supra[idx]: comm + 1 for idx, comm in enumerate(memb)}

    # Map back to per-layer node order.
    membership = []
    for t in range(n_layers):
        arr = np.array(
            [meta[vkey(t + 1, cid)] for _, cid in layer_nodes[t]],
            dtype=int,
        )
        membership.append(arr)
    return membership
