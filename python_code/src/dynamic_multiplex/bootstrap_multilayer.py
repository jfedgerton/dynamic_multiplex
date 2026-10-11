"""Parametric network bootstrap for multilayer community detection.

Redraws the full edge set B times from the fitted partition, refits the same
method and records how well each replicate's tracked partition agrees with the
original. Since 1.4.0 the package uses the bootstrap only for the stability
score s (``partition_stability``), which compares fits of the same network.
The interval and accuracy-floor outputs of earlier versions
(``co_assignment_ci``, ``community_est``, the calibrated floor) were removed:
their coverage was not reliable.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from .fit_multilayer_identity_ties import fit_multilayer_identity_ties
from .multilayer_utils import _sort_key
from .fit_multilayer_jaccard import fit_multilayer_jaccard
from .fit_multilayer_overlap import fit_multilayer_overlap
from .fit_multilayer_weighted_jaccard import fit_multilayer_weighted_jaccard
from .fit_multilayer_weighted_overlap import fit_multilayer_weighted_overlap


_FIT_FNS = {
    "jaccard": fit_multilayer_jaccard,
    "overlap": fit_multilayer_overlap,
    "weighted_jaccard": fit_multilayer_weighted_jaccard,
    "weighted_overlap": fit_multilayer_weighted_overlap,
    "identity": fit_multilayer_identity_ties,
}


@dataclass
class BootstrapResult:
    """Results from bootstrap_multilayer.

    Attributes
    ----------
    n_boot : int
        Number of bootstrap replicates completed.
    co_assignment : list[np.ndarray]
        Per-layer n_nodes x n_nodes matrices of co-assignment probabilities
        on the cross-layer meta-communities (probability two nodes share a
        persistent community).
    node_stability : list[np.ndarray]
        Per-layer array of length n_nodes giving the fraction of bootstrap
        replicates in which each node was assigned to its modal community.
    modularity_samples : list[np.ndarray]
        Per-layer array of length n_boot with bootstrap modularity values.
    community_count_reproducibility : list[float]
        Per-layer share of completed replicates whose community count equals
        the observed-network (point-estimate) count. A descriptive stability
        measure. The raw per-replicate community counts are intentionally not
        retained.
    point_estimate : dict
        The fit result from the original (unperturbed) data.
    stability_samples : dict
        ``{"nmi": array, "ari": array}`` of shape (completed replicates,
        layers): agreement of each replicate's meta-partition with the
        point-estimate partition. Summarised by ``partition_stability``.
    node_jaccard_stability : list[np.ndarray]
        Per-layer arrays: for each node, the mean Jaccard overlap between its
        replicate community and its point-estimate community.
    """

    n_boot: int
    co_assignment: list[np.ndarray]
    node_stability: list[np.ndarray]
    modularity_samples: list[np.ndarray]
    community_count_reproducibility: list[float]
    point_estimate: dict
    stability_samples: dict | None = None
    node_jaccard_stability: list[np.ndarray] | None = None


def bootstrap_multilayer(
    layers,
    fit_type: str = "jaccard",
    algorithm: str = "leiden",
    n_boot: int = 100,
    layer_links=None,
    min_similarity: float = 0.0,
    resolution_parameter: float = 1.0,
    directed: bool = False,
    seed: int | None = None,
    objective: str | None = None,
    tracking: str = "smoothed",
    repair_chains: bool = False,
    repair_min_similarity: float = 0.2,
    allow_unequal_nodes: bool = False,
) -> BootstrapResult:
    """Bootstrap confidence intervals for multilayer community detection.

    Refits communities on ``n_boot`` resampled networks using a parametric
    network bootstrap: within- and between-community edge probabilities
    (and, for weighted networks, edge-weight pools) are estimated from the
    observed network using the point-estimate partition; each replicate
    redraws the full edge set from those estimates, reproducing the
    variability of fresh data, including which edges exist.

    Uncertainty is quantified on the cross-layer *meta-communities* (the
    tracked partition from the second-stage detection), not the
    independently-detected per-layer communities. Co-assignment therefore
    answers "do these two nodes belong to the same persistent community,"
    and the community count is the number of meta-communities per layer.

    Versions before 1.1.0 instead used a Bayesian bootstrap on edge
    weights (Exponential(1) multipliers on a fixed topology). That scheme
    was removed: it understates the variability of fresh data, and in
    simulation studies intervals built from it undercovered substantially
    (~45-48% at nominal 95% for pairwise co-assignment).

    Parameters
    ----------
    layers : list
        List of adjacency matrices (numpy arrays) or NetworkX graphs.
    fit_type : str
        One of 'jaccard', 'overlap', 'weighted_jaccard', 'weighted_overlap',
        'identity'.
    algorithm : str
        Community detection algorithm ('louvain' or 'leiden').
    n_boot : int
        Number of bootstrap replicates.
    layer_links : optional
        Custom layer connectivity (DataFrame or list of dicts).
    min_similarity : float
        Minimum weighted similarity for interlayer ties.
    resolution_parameter : float
        Resolution parameter for community detection.
    directed : bool
        Whether networks are directed.
    seed : int or None
        Random seed for reproducibility.
    tracking, repair_chains, repair_min_similarity
        Passed to the two-stage fitters (see ``fit_multilayer_jaccard``) so
        that every refit uses the same settings as the point estimate. Not
        available for ``fit_type="identity"``.
    allow_unequal_nodes : bool
        If True, layers may hold different node sets (nodes entering and
        leaving, with or without a change in layer size). Layers must then be
        networkx graphs with node labels; each layer is redrawn on its own
        nodes and every refit matches nodes across layers by label. Default
        False: all layers must share one node set.

    Returns
    -------
    BootstrapResult
        Dataclass with co-assignment matrices, node stability, modularity
        samples, per-layer community-count reproducibility, and the point
        estimate.
    """
    if fit_type not in _FIT_FNS:
        raise ValueError(
            f"`fit_type` must be one of {set(_FIT_FNS.keys())}."
        )

    fit_fn = _FIT_FNS[fit_type]

    import networkx as nx

    is_graph = [hasattr(layer, "nodes") for layer in layers]
    graph_sets = [set(layer.nodes()) for layer, g in zip(layers, is_graph) if g]
    same_graph_nodes = all(ns == graph_sets[0] for ns in graph_sets[1:]) if graph_sets else True
    node_lists = None            # per-layer node labels (unequal mode only)
    if allow_unequal_nodes:
        if not all(is_graph):
            raise ValueError(
                "allow_unequal_nodes=True requires networkx graphs with node labels on every layer."
            )
        if not same_graph_nodes:
            node_lists = [sorted(layer.nodes(), key=_sort_key) for layer in layers]
    elif not same_graph_nodes:
        # Positions would silently pair different nodes (e.g. same-size layers
        # whose members changed), so fail loudly.
        raise ValueError(
            "Layers have different node sets; pass allow_unequal_nodes=True "
            "(with labelled networkx graphs) to match nodes by label."
        )

    # Convert layers to numpy arrays for resampling
    np_layers = []
    for i, layer in enumerate(layers):
        if node_lists is not None:
            np_layers.append(nx.to_numpy_array(layer, nodelist=node_lists[i], weight="weight"))
        elif is_graph[i]:
            np_layers.append(nx.to_numpy_array(layer, weight="weight"))
        else:
            np_layers.append(np.asarray(layer, dtype=float))

    n_layers = len(np_layers)
    sizes = [m.shape[0] for m in np_layers]

    def as_fit_input(mats):
        if node_lists is None:
            return mats
        out = []
        for i, m in enumerate(mats):
            g = nx.from_numpy_array(m, create_using=nx.DiGraph if directed else nx.Graph)
            out.append(nx.relabel_nodes(g, dict(enumerate(node_lists[i]))))
        return out

    # Build common kwargs
    fit_kwargs = {
        "algorithm": algorithm,
        "layer_links": layer_links,
        "directed": directed,
        "objective": objective,
    }
    if fit_type != "identity":
        fit_kwargs["min_similarity"] = min_similarity
        fit_kwargs["resolution_parameter"] = resolution_parameter
        # 1.4.0 options are passed only when set, so the default refits are
        # the 1.3.2 calls exactly.
        if tracking != "smoothed":
            fit_kwargs["tracking"] = tracking
        if repair_chains:
            fit_kwargs["repair_chains"] = repair_chains
            fit_kwargs["repair_min_similarity"] = repair_min_similarity
    elif tracking != "smoothed" or repair_chains:
        raise ValueError("`tracking` and `repair_chains` apply to the two-stage fit types, not 'identity'.")
    if node_lists is not None:
        fit_kwargs["allow_unequal_nodes"] = True

    # Point estimate on original data
    point_estimate = fit_fn(as_fit_input(np_layers), **fit_kwargs)

    # Precompute per-layer edge models for the parametric network bootstrap
    # (estimated once from the observed network + point-estimate partition)
    edge_models = []
    for layer_idx in range(n_layers):
        A = np_layers[layer_idx]
        n_nodes = sizes[layer_idx]
        memd = point_estimate["layer_communities"][layer_idx].membership
        if node_lists is None:
            mem_vec = np.array([memd[i + 1] for i in range(n_nodes)])
        else:
            mem_vec = np.array([memd[v] for v in node_lists[layer_idx]])
        same = mem_vec[:, None] == mem_vec[None, :]
        if directed:
            sel = ~np.eye(n_nodes, dtype=bool)
        else:
            sel = np.triu(np.ones((n_nodes, n_nodes), dtype=bool), k=1)
        edge_present = A > 0
        in_dyads = sel & same
        out_dyads = sel & ~same
        p_all = float(edge_present[sel].mean()) if sel.any() else 0.0
        p_in = float(edge_present[in_dyads].mean()) if in_dyads.any() else p_all
        p_out = float(edge_present[out_dyads].mean()) if out_dyads.any() else p_all
        w_all = A[sel & edge_present]
        if w_all.size == 0:
            w_all = np.array([1.0])
        w_in = A[in_dyads & edge_present]
        w_out = A[out_dyads & edge_present]
        if w_in.size == 0:
            w_in = w_all
        if w_out.size == 0:
            w_out = w_all
        edge_models.append(
            {"same": same, "sel": sel, "p_in": p_in, "p_out": p_out,
             "w_in": w_in, "w_out": w_out}
        )

    # Accumulators
    co_assign_accum = [np.zeros((n, n)) for n in sizes]
    membership_records = [[[] for _ in range(n)] for n in sizes]
    mod_samples = [[] for _ in range(n_layers)]
    count_samples = [[] for _ in range(n_layers)]

    # Agreement of every replicate with the point estimate (partition level:
    # NMI and ARI per layer; node level: Jaccard of each node's community).
    point_meta = [np.asarray(m) for m in point_estimate["meta_communities"]]
    nmi_rows: list[list[float]] = []
    ari_rows: list[list[float]] = []
    node_jaccard_accum = [np.zeros(n) for n in sizes]

    rng = np.random.default_rng(seed)

    for _b in range(n_boot):
        perturbed = []
        # Parametric network bootstrap: redraw the full edge set
        for layer_idx, em in enumerate(edge_models):
            n_nodes = sizes[layer_idx]
            probs = np.where(em["same"], em["p_in"], em["p_out"])
            draw = (rng.random((n_nodes, n_nodes)) < probs) & em["sel"]
            mat_new = np.zeros((n_nodes, n_nodes))
            on = np.where(draw)
            k = on[0].size
            if k > 0:
                same_on = em["same"][on]
                w = np.empty(k)
                n_in = int(same_on.sum())
                if n_in > 0:
                    w[same_on] = rng.choice(em["w_in"], size=n_in,
                                            replace=True)
                if k - n_in > 0:
                    w[~same_on] = rng.choice(em["w_out"], size=k - n_in,
                                             replace=True)
                mat_new[on] = w
            if not directed:
                mat_new = mat_new + mat_new.T
            np.fill_diagonal(mat_new, 0.0)
            perturbed.append(mat_new)

        try:
            boot_fit = fit_fn(as_fit_input(perturbed), **fit_kwargs)
        except Exception:
            continue

        boot_meta = [np.asarray(m) for m in boot_fit["meta_communities"]]
        nmi_rows.append([_nmi(boot_meta[t], point_meta[t]) for t in range(n_layers)])
        ari_rows.append([_ari(boot_meta[t], point_meta[t]) for t in range(n_layers)])
        for t in range(n_layers):
            node_jaccard_accum[t] += _node_jaccard(boot_meta[t], point_meta[t])

        for layer_idx in range(n_layers):
            lc = boot_fit["layer_communities"][layer_idx]

            # Meta (cross-layer) membership is the validated partition; the
            # co-assignment and community counts below are computed on it.
            # mem is an array in node order (index i -> node i+1).
            mem = np.asarray(boot_fit["meta_communities"][layer_idx])

            # Co-assignment matrix (on meta labels)
            comms: dict[int, list[int]] = {}
            for node_pos, meta_id in enumerate(mem):
                comms.setdefault(int(meta_id), []).append(node_pos)
            for comm_nodes in comms.values():
                for i_pos in range(len(comm_nodes)):
                    for j_pos in range(i_pos + 1, len(comm_nodes)):
                        ni = comm_nodes[i_pos]
                        nj = comm_nodes[j_pos]
                        co_assign_accum[layer_idx][ni, nj] += 1
                        co_assign_accum[layer_idx][nj, ni] += 1

            # Record meta membership for each node
            for node_pos, meta_id in enumerate(mem):
                membership_records[layer_idx][node_pos].append(int(meta_id))

            # Modularity (per-layer detection value, descriptive)
            mod_val = lc.modularity
            mod_samples[layer_idx].append(
                mod_val if mod_val is not None else np.nan
            )

            # Community count = number of distinct meta communities in layer
            count_samples[layer_idx].append(len(comms))

    # Completed bootstrap count (some may have failed)
    n_completed = len(mod_samples[0]) if mod_samples[0] else 0

    # Normalize co-assignment by completed replicates
    if n_completed > 0:
        co_assignment = [m / n_completed for m in co_assign_accum]
        # Diagonal = 1.0 (node always co-assigned with itself)
        for m in co_assignment:
            np.fill_diagonal(m, 1.0)
    else:
        co_assignment = co_assign_accum

    # Node stability: fraction of times in modal community
    node_stability = []
    for layer_idx in range(n_layers):
        n_nodes = sizes[layer_idx]
        stab = np.zeros(n_nodes)
        for node_idx in range(n_nodes):
            records = membership_records[layer_idx][node_idx]
            if records:
                from collections import Counter

                counts = Counter(records)
                stab[node_idx] = counts.most_common(1)[0][1] / len(records)
        node_stability.append(stab)

    # Per-layer bootstrap reproducibility of the community count: the share of
    # completed replicates whose community count equals the observed-network
    # (point-estimate) count. Computed here, from the raw per-replicate counts,
    # so those raw counts are never retained on the returned object.
    community_count_reproducibility = []
    for layer_idx in range(n_layers):
        est = len(np.unique(point_estimate["meta_communities"][layer_idx]))
        s = np.asarray(count_samples[layer_idx])
        community_count_reproducibility.append(
            float(np.mean(s == est)) if s.size else float("nan")
        )

    return BootstrapResult(
        n_boot=n_completed,
        co_assignment=co_assignment,
        node_stability=node_stability,
        modularity_samples=[np.array(s) for s in mod_samples],
        community_count_reproducibility=community_count_reproducibility,
        point_estimate=point_estimate,
        stability_samples={"nmi": np.array(nmi_rows), "ari": np.array(ari_rows)},
        node_jaccard_stability=[
            (v / n_completed) if n_completed > 0 else v for v in node_jaccard_accum
        ],
    )


def _contingency(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    _, ia = np.unique(a, return_inverse=True)
    _, ib = np.unique(b, return_inverse=True)
    ct = np.zeros((ia.max() + 1, ib.max() + 1))
    np.add.at(ct, (ia, ib), 1)
    return ct


def _nmi(a: np.ndarray, b: np.ndarray) -> float:
    """Normalized mutual information, 2I/(H(a)+H(b)), matching igraph::compare(method="nmi")."""
    ct = _contingency(a, b); n = ct.sum()
    pa = ct.sum(1) / n; pb = ct.sum(0) / n; pij = ct / n
    nz = pij > 0
    mi = float((pij[nz] * np.log(pij[nz] / np.outer(pa, pb)[nz])).sum())
    ha = float(-(pa[pa > 0] * np.log(pa[pa > 0])).sum())
    hb = float(-(pb[pb > 0] * np.log(pb[pb > 0])).sum())
    if ha + hb == 0:
        return 1.0
    return 2 * mi / (ha + hb)


def _ari(a: np.ndarray, b: np.ndarray) -> float:
    """Adjusted Rand index (Hubert and Arabie), matching igraph::compare(method="adjusted.rand")."""
    ct = _contingency(a, b); n = ct.sum()
    comb = lambda x: x * (x - 1) / 2.0
    sum_ij = comb(ct).sum(); sum_a = comb(ct.sum(1)).sum(); sum_b = comb(ct.sum(0)).sum()
    expected = sum_a * sum_b / comb(n) if n > 1 else 0.0
    max_index = (sum_a + sum_b) / 2.0
    if max_index == expected:
        return 1.0
    return float((sum_ij - expected) / (max_index - expected))


def _node_jaccard(mem_a: np.ndarray, mem_b: np.ndarray) -> np.ndarray:
    same_a = mem_a[:, None] == mem_a[None, :]
    same_b = mem_b[:, None] == mem_b[None, :]
    return (same_a & same_b).sum(1) / (same_a | same_b).sum(1)


def partition_stability(
    boot_result: BootstrapResult,
    metric: str = "nmi",
) -> dict:
    """Bootstrap stability score s of the tracked partition.

    ``stability`` is the mean agreement (NMI by default, or ARI) between each
    bootstrap replicate's meta-partition and the point-estimate partition,
    averaged over layers and replicates. It is a point estimate with Monte
    Carlo standard error ``stability_mc_se``. s has no validated mapping to
    accuracy: use it only to compare two fits of the same network (for
    example DynMux against multislice), never as an estimate of how accurate
    either fit is. Since 1.4.0 the package reports no interval or accuracy
    floor.

    Parameters
    ----------
    boot_result : BootstrapResult
        From ``bootstrap_multilayer`` (1.3.0 or later).
    metric : str
        ``"nmi"`` (default) or ``"ari"``.

    Returns
    -------
    dict
        Keys ``stability``, ``stability_by_layer``, ``stability_mc_se``,
        ``metric``, ``report``.
    """
    if metric not in ("nmi", "ari"):
        raise ValueError('metric must be "nmi" or "ari"')
    if boot_result.stability_samples is None:
        raise ValueError("boot_result has no stability_samples; rerun bootstrap_multilayer (>= 1.3.0).")
    S = np.asarray(boot_result.stability_samples[metric])
    if S.shape[0] < 2:
        raise ValueError("At least two completed bootstrap replicates are needed.")
    n_layers = S.shape[1]
    by_layer = S.mean(0)
    per_rep = S.mean(1)
    s = float(per_rep.mean())
    mc_se = float(per_rep.std(ddof=1) / np.sqrt(len(per_rep)))
    report = (
        f"Partition stability {s:.2f} ({metric.upper()} over {boot_result.n_boot} bootstrap "
        f"replicates and {n_layers} layers; MC s.e. {mc_se:.3f}). Use it only to compare fits "
        "of the same network."
    )
    return {"stability": s, "stability_by_layer": by_layer, "stability_mc_se": mc_se,
            "metric": metric, "report": report}
