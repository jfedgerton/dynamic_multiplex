"""Bipartite and all-or-nothing hypergraph partitions of membership slices (Python only).

These two specifications treat a membership slice (nodes in groups) without
projecting it to a node-by-node network:

* ``partition_membership_bipartite``: nodes and groups are the two sides of a
  bipartite graph; communities maximise Barber's (2007) bipartite modularity,
  optimised with leidenalg's bipartite construction (Leiden only).
* ``partition_membership_strict``: groups are hyperedges; communities maximise
  the strict (all-or-nothing) hypergraph modularity of Kaminski et al., in
  which a group counts only when all its members share a community. Uses the
  reference implementation in hypernetx (optional extra; it requires
  pandas < 3).

Both work slice by slice, so slice t never depends on another slice.
``fit_multilayer_from_partitions`` then tracks any per-slice partitions over
time with the DynMux second stage or Hungarian matching, so these
specifications can be compared with ``layers_from_membership`` on equal terms.

The R package deliberately has no counterpart: the libraries these rely on
exist in Python and not in R.

In small pilot simulations both specifications did worse than the projections
in ``layers_from_membership`` (bipartite over-split once noise groups were
present; for the strict objective, the planted partition scored below the
partitions found). Treat them as comparison specifications, not defaults.
"""
from __future__ import annotations

import random
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from .layers_from_membership import _prepare_slice
from .multilayer_utils import (
    LayerCommunityFit,
    _check_repair_args,
    _sort_key,
    community_overlap_edges,
    make_layer_links,
    run_two_stage,
)


@dataclass
class MembershipPartition:
    """Partition of one membership slice.

    nodes : node label -> community (1..K)
    groups : group label -> community (bipartite only; None otherwise)
    quality : Barber modularity (bipartite) or strict hypergraph modularity
    initial_quality : strict modularity of the Kumar starting partition (strict only)
    groups_dropped : groups removed for having fewer than ``min_group_size`` members
    """

    nodes: dict
    groups: dict | None
    quality: float
    initial_quality: float | None = None
    groups_dropped: int = 0
    specification: str = field(default="")


def _check_common(memberships, resolution_parameter, min_group_size):
    if not isinstance(memberships, list):
        raise ValueError("`memberships` must be a list with one entry per time slice.")
    if len(memberships) == 0:
        raise ValueError("`memberships` must contain at least one slice.")
    if isinstance(resolution_parameter, bool) or not isinstance(resolution_parameter, (int, float)) \
            or not resolution_parameter > 0:
        raise ValueError("`resolution_parameter` must be a number > 0.")
    if isinstance(min_group_size, bool) or not isinstance(min_group_size, (int, np.integer)) or min_group_size < 2:
        raise ValueError("`min_group_size` must be an integer >= 2.")


def _relabel(raw) -> np.ndarray:
    """Community labels 1..K in order of first appearance."""
    _, first, inv = np.unique(np.asarray(raw), return_index=True, return_inverse=True)
    order = np.argsort(np.argsort(first))
    return order[inv] + 1


def _barber_modularity(B: np.ndarray, node_comm, group_comm, resolution_parameter: float = 1.0) -> float:
    """Barber's bipartite modularity: (1/m) sum_ij (B_ij - gamma k_i d_j / m) [c_i == c_j]."""
    B = np.asarray(B, dtype=float)
    node_comm = np.asarray(node_comm); group_comm = np.asarray(group_comm)
    m = B.sum()
    k = B.sum(axis=1); d = B.sum(axis=0)
    same = node_comm[:, None] == group_comm[None, :]
    return float(((B - resolution_parameter * np.outer(k, d) / m) * same).sum() / m)


def partition_membership_bipartite(
    memberships: list,
    resolution_parameter: float = 1.0,
    seed: int | None = 123,
    min_group_size: int = 2,
    algorithm: str = "leiden",
) -> list[MembershipPartition]:
    """Partition each membership slice by Barber's bipartite modularity.

    Nodes and groups form a bipartite graph (one edge per membership).
    Communities hold both nodes and groups and maximise

        Q_B = (1/m) sum_{i,j} (B_ij - gamma k_i d_j / m) [c_i == c_j]

    over node-group pairs (Barber 2007), with k_i the number of groups of node
    i, d_j the size of group j and m the number of memberships. It is
    optimised with leidenalg's bipartite construction: three CPM layers with
    degrees as node sizes and inter-class resolution gamma / m, combined with
    layer weights [1, -1, -1]; leidenalg's combined quality is then 2 m Q_B
    (checked at run time). Leiden only.

    Parameters
    ----------
    memberships : list
        One entry per slice, as in ``layers_from_membership``.
    resolution_parameter : float
        gamma (> 0); 1 is Barber's modularity.
    seed : int or None
        Seed for leidenalg, reused for every slice (default 123).
    min_group_size : int
        Groups with fewer members (default 2) are dropped, with a warning.
    algorithm : str
        Must be ``"leiden"``; Louvain has no bipartite construction here.

    Returns
    -------
    list[MembershipPartition]
        One per slice; ``quality`` is Barber's modularity of the partition.
    """
    if algorithm != "leiden":
        raise ValueError("The bipartite specification is available with Leiden only (algorithm='leiden').")
    _check_common(memberships, resolution_parameter, min_group_size)
    import igraph as ig
    import leidenalg

    out = []
    for t, entry in enumerate(memberships, start=1):
        labels, groups, B, _, n_drop = _prepare_slice(entry, t, min_group_size)
        n, M = B.shape
        rows, cols = np.nonzero(B)
        g = ig.Graph(n=n + M, edges=list(zip(rows.tolist(), (cols + n).tolist())), directed=False)
        types = [0] * n + [1] * M
        m = float(len(rows))
        p01, p0, p1 = leidenalg.CPMVertexPartition.Bipartite(
            g, resolution_parameter_01=resolution_parameter / m,
            degree_as_node_size=True, types=types)
        opt = leidenalg.Optimiser()
        if seed is not None:
            opt.set_rng_seed(seed)
        opt.optimise_partition_multiplex([p01, p0, p1], layer_weights=[1, -1, -1], n_iterations=-1)
        memb = _relabel(p01.membership)
        node_c, group_c = memb[:n], memb[n:]
        q = _barber_modularity(B, node_c, group_c, resolution_parameter)
        # leidenalg's combined quality (CPM with binomial node-size terms, each
        # internal edge counted twice) equals 2 m Q_B for this construction
        combined = p01.quality() - p0.quality() - p1.quality()
        assert abs(combined - 2.0 * m * q) <= 1e-8 * max(1.0, abs(2.0 * m * q)), "bipartite quality mismatch"
        out.append(MembershipPartition(
            nodes={labels[i]: int(node_c[i]) for i in range(n)},
            groups={groups[j]: int(group_c[j]) for j in range(M)},
            quality=q, groups_dropped=n_drop, specification="bipartite"))
    return out


def partition_membership_strict(
    memberships: list,
    seed: int | None = 123,
    min_group_size: int = 2,
    delta: float = 0.01,
) -> list[MembershipPartition]:
    """Partition each membership slice by strict (all-or-nothing) hypergraph modularity.

    Groups are hyperedges. Under the strict weighting a group contributes to
    modularity only when all its members are in one community (Kaminski et
    al.). Uses the hypernetx 2.4.3 reference pipeline: a Kumar starting
    partition (``hypergraph_modularity.kumar``) improved by node moves
    (``last_step`` with the ``strict`` weighting). Python's ``random`` (used by
    igraph) and NumPy's global generator (used by ``last_step``) are seeded
    with ``seed`` for every slice and restored afterwards.

    Requires the optional extra ``pip install dynamic-multiplex[hypergraph]``
    (hypernetx, which needs pandas < 3).

    Returns
    -------
    list[MembershipPartition]
        One per slice; ``quality`` is the strict modularity of the partition
        and ``initial_quality`` that of the Kumar start.
    """
    _check_common(memberships, 1.0, min_group_size)
    if isinstance(delta, bool) or not isinstance(delta, (int, float)) or not delta > 0:
        raise ValueError("`delta` must be a number > 0.")
    try:
        import hypernetx as hnx
        from hypernetx.algorithms.clustering import hypergraph_modularity as hmod
    except ImportError as exc:
        raise ImportError(
            "partition_membership_strict needs hypernetx: pip install 'dynamic-multiplex[hypergraph]' "
            "(hypernetx requires pandas < 3)."
        ) from exc

    out = []
    for t, entry in enumerate(memberships, start=1):
        labels, groups, B, _, n_drop = _prepare_slice(entry, t, min_group_size)
        edges = {g: [labels[i] for i in np.flatnonzero(B[:, j])] for j, g in enumerate(groups)}
        py_state, np_state = random.getstate(), np.random.get_state()
        try:
            if seed is not None:
                random.seed(seed)
                np.random.seed(seed)
            H = hnx.Hypergraph(edges)
            start = hmod.kumar(H, delta=delta)
            q0 = float(hmod.modularity(H, start, wdc=hmod.strict))
            final = hmod.last_step(H, start, wdc=hmod.strict, delta=delta)
            q = float(hmod.modularity(H, final, wdc=hmod.strict))
        finally:
            random.setstate(py_state)
            np.random.set_state(np_state)
        assign = {v: k for k, block in enumerate(final) for v in block}
        memb = _relabel([assign[v] for v in labels])
        out.append(MembershipPartition(
            nodes={labels[i]: int(memb[i]) for i in range(len(labels))}, groups=None,
            quality=q, initial_quality=q0, groups_dropped=n_drop, specification="strict"))
    return out


def fit_multilayer_from_partitions(
    partitions: list,
    method: str = "jaccard",
    algorithm: str = "leiden",
    layer_links=None,
    min_similarity: float = 0.0,
    resolution_parameter: float = 1.0,
    add_self_loops: bool = True,
    self_loop_multiplier: float = 1.0,
    seed: int | None = 123,
    tracking: str = "smoothed",
    repair_chains: bool = False,
    repair_min_similarity: float = 0.2,
):
    """Track supplied per-slice partitions over time.

    Runs the DynMux second stage (``method="jaccard"`` or ``"overlap"``, with
    the same options as ``fit_multilayer_jaccard``) or Hungarian matching of
    consecutive slices (``method="hungarian"``) on partitions computed
    elsewhere, for example by ``partition_membership_bipartite`` or
    ``partition_membership_strict``. Given the per-layer partitions a fitter
    computes itself, it returns that fitter's tracked partition exactly.

    Parameters
    ----------
    partitions : list
        One entry per slice: a dict node label -> community, or a
        ``MembershipPartition``. Node sets may differ across slices; nodes are
        matched by label.
    method : str
        ``"jaccard"``, ``"overlap"`` or ``"hungarian"``.
    algorithm, layer_links, min_similarity, resolution_parameter, add_self_loops,
    self_loop_multiplier, seed, tracking, repair_chains, repair_min_similarity
        As in ``fit_multilayer_jaccard`` (stage 2 only; ignored by
        ``"hungarian"`` except ``seed``, which it does not use).

    Returns
    -------
    dict
        ``meta_communities`` (one array per slice in ``nodes`` order),
        ``nodes`` (per-slice sorted node labels), ``layer_communities``,
        ``meta_ids``, ``interlayer_ties``, ``layer_links``, ``method``,
        ``tracking``, ``repair_chains``, ``repair_min_similarity``,
        ``chain_repairs``.
    """
    from .fit_multilayer_hungarian import _match_communities_hungarian

    if method not in ("jaccard", "overlap", "hungarian"):
        raise ValueError('`method` must be "jaccard", "overlap" or "hungarian".')
    if not isinstance(partitions, list) or len(partitions) < 2:
        raise ValueError("`partitions` must be a list with at least two slices.")
    mems = []
    for t, p in enumerate(partitions, start=1):
        if isinstance(p, MembershipPartition):
            p = p.nodes
        if not isinstance(p, dict) or not p:
            raise ValueError(f"Slice {t}: each partition must be a non-empty dict (node -> community) "
                             "or a MembershipPartition.")
        mems.append(dict(p))
    nodes = [sorted(m, key=_sort_key) for m in mems]

    if method == "hungarian":
        meta = _match_communities_hungarian(mems)
        return {"method": "hungarian", "meta_communities": meta, "nodes": nodes, "meta_ids": None,
                "layer_communities": None, "interlayer_ties": None, "layer_links": None,
                "tracking": "smoothed", "repair_chains": False, "repair_min_similarity": None,
                "chain_repairs": None}

    _check_repair_args(tracking, repair_chains, repair_min_similarity)
    algorithm = algorithm.lower()
    if algorithm not in ("leiden", "louvain"):
        raise ValueError("`algorithm` must be one of {'louvain', 'leiden'}.")
    fit = []
    for m in mems:
        comms: dict = {}
        for node, c in m.items():                      # insertion order, as in the fitters
            comms.setdefault(int(c), []).append(node)
        comms = {k: sorted(v, key=_sort_key) for k, v in comms.items()}
        fit.append(LayerCommunityFit(membership={k: int(v) for k, v in m.items()}, modularity=None,
                                     communities=comms))
    links = make_layer_links(len(fit), layer_links)
    pair_ties = community_overlap_edges(fit=fit, layer_links=links, metric=method, min_similarity=min_similarity)
    meta = run_two_stage(fit, pair_ties, links, add_self_loops=add_self_loops,
                         self_loop_multiplier=self_loop_multiplier, min_similarity=min_similarity, directed=False,
                         algorithm=algorithm, resolution_parameter=resolution_parameter, seed=seed,
                         tracking=tracking, repair_chains=repair_chains,
                         repair_min_similarity=repair_min_similarity)
    return {"method": method, "algorithm": algorithm, "meta_communities": meta["membership"], "nodes": nodes,
            "meta_ids": meta["meta_ids"], "layer_communities": fit, "interlayer_ties": meta["interlayer_ties"],
            "layer_links": links, "tracking": meta["tracking"], "repair_chains": meta["repair_chains"],
            "repair_min_similarity": meta["repair_min_similarity"], "chain_repairs": meta["chain_repairs"]}
