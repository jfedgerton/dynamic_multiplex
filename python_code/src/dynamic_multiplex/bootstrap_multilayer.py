"""Bootstrap confidence intervals for multilayer community detection.

Provides nonparametric bootstrap uncertainty quantification by resampling
edge weights (Bayesian bootstrap), re-running community detection B times,
and computing co-assignment probabilities, modularity CIs, and node-level
stability measures.
"""

from __future__ import annotations

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
) -> BootstrapResult:
    """Bootstrap confidence intervals for multilayer community detection.

    Uses a Bayesian bootstrap on edge weights: each bootstrap replicate
    multiplies every edge weight by an independent Exponential(1) draw,
    preserving graph topology while perturbing the weighting that drives
    community detection.

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
    }
    if fit_type != "identity":
        fit_kwargs["min_similarity"] = min_similarity
        fit_kwargs["resolution_parameter"] = resolution_parameter

    # Point estimate on original data
    point_estimate = fit_fn(np_layers, **fit_kwargs)

    # Accumulators
    co_assign_accum = [np.zeros((n_nodes, n_nodes)) for _ in range(n_layers)]
    membership_records = [[[] for _ in range(n_nodes)] for _ in range(n_layers)]
    mod_samples = [[] for _ in range(n_layers)]
    count_samples = [[] for _ in range(n_layers)]

    rng = np.random.default_rng(seed)

    for _b in range(n_boot):
        # Bayesian bootstrap: multiply edge weights by Exp(1) draws
        perturbed = []
        for mat in np_layers:
            noise = rng.exponential(1.0, size=mat.shape)
            # Symmetrize noise for undirected networks
            if not directed:
                noise = (noise + noise.T) / 2.0
            perturbed_mat = mat * noise
            np.fill_diagonal(perturbed_mat, 0.0)
            perturbed.append(perturbed_mat)

        # Derive a per-replicate seed for community detection reproducibility
        boot_seed = int(rng.integers(0, 2**31))
        try:
            boot_fit = fit_fn(perturbed, seed=boot_seed, **fit_kwargs)
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
        - ``modularity_ci``: DataFrame with layer, estimate, lower, upper
        - ``community_count_ci``: DataFrame with layer, estimate, lower, upper
        - ``mean_node_stability``: DataFrame with layer, mean_stability
        - ``node_stability``: list of per-layer stability arrays
        - ``co_assignment``: list of per-layer co-assignment matrices
    """
    if boot_result.n_boot == 0:
        raise ValueError("No completed bootstrap replicates.")

    lower_q = alpha / 2
    upper_q = 1 - alpha / 2

    n_layers = len(boot_result.modularity_samples)
    point = boot_result.point_estimate

    # Modularity CIs
    mod_rows = []
    for i in range(n_layers):
        lc = point["layer_communities"][i]
        est = lc.modularity
        samples = boot_result.modularity_samples[i]
        valid = samples[~np.isnan(samples)]
        if len(valid) > 0:
            lo, hi = np.quantile(valid, [lower_q, upper_q])
        else:
            lo, hi = np.nan, np.nan
        mod_rows.append(
            {
                "layer": i + 1,
                "estimate": est if est is not None else np.nan,
                "lower": lo,
                "upper": hi,
            }
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
        "modularity_ci": pd.DataFrame(mod_rows),
        "community_count_ci": pd.DataFrame(count_rows),
        "mean_node_stability": pd.DataFrame(stab_rows),
        "node_stability": boot_result.node_stability,
        "co_assignment": boot_result.co_assignment,
    }
