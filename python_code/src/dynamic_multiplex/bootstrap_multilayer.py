"""Bootstrap confidence intervals for multilayer community detection.

Provides nonparametric bootstrap uncertainty quantification by resampling
edge weights (Bayesian bootstrap), re-running community detection B times,
and computing co-assignment probabilities, community count CIs, node-pair
co-assignment intervals, and node-level stability measures.

.. warning::
   A large simulation study found that the nominal 95% community count
   interval covers the truth in only ~40% of simulations on small networks
   (n = 50 nodes), recovering to nominal for n >= 100. A modularity CI was
   removed in version 1.1.0 because its empirical coverage is never close
   to nominal at any network size (community detection maximizes
   modularity, so the bootstrap interval concentrates around an upwardly
   biased value). Raw ``modularity_samples`` remain available for
   descriptive use.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from .fit_multilayer_identity_ties import fit_multilayer_identity_ties
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
        Per-layer n_nodes x n_nodes matrices of co-assignment probabilities.
    node_stability : list[np.ndarray]
        Per-layer array of length n_nodes giving the fraction of bootstrap
        replicates in which each node was assigned to its modal community.
    modularity_samples : list[np.ndarray]
        Per-layer array of length n_boot with bootstrap modularity values.
    community_count_samples : list[np.ndarray]
        Per-layer array of length n_boot with bootstrap community counts.
    point_estimate : dict
        The fit result from the original (unperturbed) data.
    """

    n_boot: int
    co_assignment: list[np.ndarray]
    node_stability: list[np.ndarray]
    modularity_samples: list[np.ndarray]
    community_count_samples: list[np.ndarray]
    point_estimate: dict


def bootstrap_multilayer(
    layers,
    fit_type: str = "jaccard",
    algorithm: str = "louvain",
    n_boot: int = 100,
    layer_links=None,
    min_similarity: float = 0.0,
    resolution_parameter: float = 1.0,
    directed: bool = False,
    seed: int | None = None,
    objective: str | None = None,
) -> BootstrapResult:
    """Bootstrap confidence intervals for multilayer community detection.

    Refits communities on ``n_boot`` resampled networks using a parametric
    network bootstrap: within- and between-community edge probabilities
    (and, for weighted networks, edge-weight pools) are estimated from the
    observed network using the point-estimate partition; each replicate
    redraws the full edge set from those estimates, reproducing the
    variability of fresh data, including which edges exist.

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

    Returns
    -------
    BootstrapResult
        Dataclass with co-assignment matrices, node stability, modularity
        samples, community count samples, and the point estimate.
    """
    if fit_type not in _FIT_FNS:
        raise ValueError(
            f"`fit_type` must be one of {set(_FIT_FNS.keys())}."
        )

    fit_fn = _FIT_FNS[fit_type]

    # Convert layers to numpy arrays for resampling
    np_layers = []
    for layer in layers:
        if hasattr(layer, "nodes"):  # NetworkX graph
            import networkx as nx

            np_layers.append(nx.to_numpy_array(layer, weight="weight"))
        else:
            np_layers.append(np.asarray(layer, dtype=float))

    n_layers = len(np_layers)
    n_nodes = np_layers[0].shape[0]

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

    # Point estimate on original data
    point_estimate = fit_fn(np_layers, **fit_kwargs)

    # Precompute per-layer edge models for the parametric network bootstrap
    # (estimated once from the observed network + point-estimate partition)
    edge_models = []
    for layer_idx in range(n_layers):
        A = np_layers[layer_idx]
        memd = point_estimate["layer_communities"][layer_idx].membership
        mem_vec = np.array([memd[i + 1] for i in range(n_nodes)])
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
    co_assign_accum = [np.zeros((n_nodes, n_nodes)) for _ in range(n_layers)]
    membership_records = [[[] for _ in range(n_nodes)] for _ in range(n_layers)]
    mod_samples = [[] for _ in range(n_layers)]
    count_samples = [[] for _ in range(n_layers)]

    rng = np.random.default_rng(seed)

    for _b in range(n_boot):
        perturbed = []
        # Parametric network bootstrap: redraw the full edge set
        for em in edge_models:
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
            boot_fit = fit_fn(perturbed, **fit_kwargs)
        except Exception:
            continue

        for layer_idx in range(n_layers):
            lc = boot_fit["layer_communities"][layer_idx]
            mem = lc.membership  # dict: node_id (1-indexed) -> community

            # Co-assignment matrix
            for comm_nodes in lc.communities.values():
                for i_pos in range(len(comm_nodes)):
                    for j_pos in range(i_pos + 1, len(comm_nodes)):
                        ni = comm_nodes[i_pos] - 1  # 0-indexed
                        nj = comm_nodes[j_pos] - 1
                        co_assign_accum[layer_idx][ni, nj] += 1
                        co_assign_accum[layer_idx][nj, ni] += 1

            # Record membership for each node
            for node_id, comm_id in mem.items():
                membership_records[layer_idx][node_id - 1].append(comm_id)

            # Modularity
            mod_val = lc.modularity
            mod_samples[layer_idx].append(
                mod_val if mod_val is not None else np.nan
            )

            # Community count
            count_samples[layer_idx].append(len(lc.communities))

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
        stab = np.zeros(n_nodes)
        for node_idx in range(n_nodes):
            records = membership_records[layer_idx][node_idx]
            if records:
                from collections import Counter

                counts = Counter(records)
                stab[node_idx] = counts.most_common(1)[0][1] / len(records)
        node_stability.append(stab)

    return BootstrapResult(
        n_boot=n_completed,
        co_assignment=co_assignment,
        node_stability=node_stability,
        modularity_samples=[np.array(s) for s in mod_samples],
        community_count_samples=[np.array(s) for s in count_samples],
        point_estimate=point_estimate,
    )


def community_ci(
    boot_result: BootstrapResult,
    alpha: float = 0.05,
) -> dict:
    """Summarize bootstrap results into confidence intervals.

    .. warning::
       **Community count intervals undercover badly on small networks.**
       In a simulation study on planted-partition multilayer networks
       (n = 50-400 nodes, 3-10 communities, 5-15 layers), the nominal 95%
       ``community_count_ci`` contained the true community count in only
       ~40% of simulations at n = 50 nodes. Coverage recovers to at or
       above nominal for n >= 100. A ``UserWarning`` is raised when layers
       have fewer than 100 nodes; on such networks treat the intervals as
       descriptive stability summaries, not calibrated CIs.

       Earlier versions also returned ``modularity_ci``. It was removed in
       1.1.0: its empirical coverage is never close to nominal at any
       network size (~0.40 at n = 50, 0.00 at n = 100, vacuously 1.00 at
       n = 200) because community detection maximizes modularity and the
       bootstrap interval concentrates around that optimized, upwardly
       biased value. Raw draws remain in
       ``BootstrapResult.modularity_samples``.

    Parameters
    ----------
    boot_result : BootstrapResult
        Output from ``bootstrap_multilayer``.
    alpha : float
        Significance level (default 0.05 for 95% CIs).

    Returns
    -------
    dict
        Dictionary with keys:
        - ``community_count_ci``: DataFrame with layer, estimate, lower, upper
        - ``mean_node_stability``: DataFrame with layer, mean_stability
        - ``node_stability``: list of per-layer stability arrays
        - ``co_assignment``: list of per-layer co-assignment matrices

    See Also
    --------
    co_assignment_ci : Wilson intervals for node-pair co-assignment.
    """
    if boot_result.n_boot == 0:
        raise ValueError("No completed bootstrap replicates.")

    lower_q = alpha / 2
    upper_q = 1 - alpha / 2

    n_layers = len(boot_result.modularity_samples)
    point = boot_result.point_estimate

    # Warn loudly on small networks: coverage study showed ~40% empirical
    # coverage for the nominal 95% community count CI at n = 50 nodes,
    # recovering to nominal at n >= 100.
    n_nodes = len(boot_result.node_stability[0])
    if n_nodes < 100:
        warnings.warn(
            f"Layers have {n_nodes} nodes (< 100). In simulation studies the "
            "nominal 95% community_count_ci covered the true community count "
            "in only ~40% of small-network simulations (n = 50). Treat these "
            "intervals as descriptive stability summaries, not calibrated "
            "confidence intervals.",
            UserWarning,
            stacklevel=2,
        )

    # Community count CIs
    count_rows = []
    for i in range(n_layers):
        lc = point["layer_communities"][i]
        est = len(lc.communities)
        samples = boot_result.community_count_samples[i]
        lo, hi = np.quantile(samples, [lower_q, upper_q])
        count_rows.append(
            {"layer": i + 1, "estimate": est, "lower": lo, "upper": hi}
        )

    # Mean node stability per layer
    stab_rows = []
    for i in range(n_layers):
        stab_rows.append(
            {
                "layer": i + 1,
                "mean_stability": float(np.mean(boot_result.node_stability[i])),
            }
        )

    return {
        "community_count_ci": pd.DataFrame(count_rows),
        "mean_node_stability": pd.DataFrame(stab_rows),
        "node_stability": boot_result.node_stability,
        "co_assignment": boot_result.co_assignment,
    }


def co_assignment_ci(
    boot_result: BootstrapResult,
    alpha: float = 0.05,
) -> list[dict]:
    """Wilson confidence intervals for node-pair co-assignment.

    For every pair of nodes in every layer, computes a Wilson score
    interval for the co-clustering propensity: the probability that the
    fitted community detection procedure places the two nodes in the same
    community when the data are perturbed. The point estimate is the
    co-assignment probability from ``bootstrap_multilayer`` (the share of
    bootstrap replicates in which the pair was co-assigned), and the
    interval treats the ``n_boot`` replicates as binomial draws.

    Because co-assignment is label-invariant (it never compares community
    labels across replicates, only whether two nodes sit together), it
    avoids the label-switching problem that makes per-node membership
    intervals ill-defined.

    .. warning::
       These intervals quantify the stability of the detection procedure,
       not the probability that two nodes truly share a community.
       Interpret cautiously on networks with fewer than 100 nodes, where
       community detection itself is unstable.

    Parameters
    ----------
    boot_result : BootstrapResult
        Output from ``bootstrap_multilayer``.
    alpha : float
        Significance level (default 0.05 for 95% intervals).

    Returns
    -------
    list[dict]
        One dict per layer with keys ``estimate``, ``lower``, ``upper``,
        each an ``n_nodes x n_nodes`` ndarray. Diagonals are 1 by
        construction.

    See Also
    --------
    community_ci : Community count intervals and node stability summaries.
    """
    if boot_result.n_boot == 0:
        raise ValueError("No completed bootstrap replicates.")

    from statistics import NormalDist  # stdlib; avoids a scipy dependency

    b = boot_result.n_boot
    z = NormalDist().inv_cdf(1 - alpha / 2)
    z2 = z**2

    layer_cis = []
    for phat in boot_result.co_assignment:
        denom = 1 + z2 / b
        center = (phat + z2 / (2 * b)) / denom
        half = z * np.sqrt(phat * (1 - phat) / b + z2 / (4 * b**2)) / denom
        lower = np.clip(center - half, 0.0, 1.0)
        upper = np.clip(center + half, 0.0, 1.0)
        np.fill_diagonal(lower, 1.0)
        np.fill_diagonal(upper, 1.0)
        layer_cis.append(
            {"estimate": phat, "lower": lower, "upper": upper}
        )

    return layer_cis
