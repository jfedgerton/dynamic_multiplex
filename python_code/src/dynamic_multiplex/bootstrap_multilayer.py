"""Bootstrap confidence intervals for multilayer community detection.

Provides parametric network-bootstrap uncertainty quantification by
redrawing the full edge set B times, re-running community detection, and
computing co-assignment probabilities, node-pair co-assignment intervals,
community-count reproducibility, and node-level stability measures.

.. warning::
   The percentile community-count confidence interval was replaced in
   version 1.1.0 with a descriptive reproducibility summary: a large
   simulation study found its coverage collapsed under community-size skew
   (~0.62) with no observable diagnostic to flag the failure. A modularity
   CI was removed in the same release because its empirical coverage is
   never close to nominal at any network size (community detection maximizes
   modularity, so the bootstrap interval concentrates around an upwardly
   biased value). The validated interval in this module is the node-pair
   co-assignment interval (``co_assignment_ci``). Raw ``modularity_samples``
   remain available for descriptive use.
"""

from __future__ import annotations

from dataclasses import dataclass

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

    # Agreement of every replicate with the point estimate (partition level:
    # NMI and ARI per layer; node level: Jaccard of each node's community).
    point_meta = [np.asarray(m) for m in point_estimate["meta_communities"]]
    nmi_rows: list[list[float]] = []
    ari_rows: list[list[float]] = []
    node_jaccard_accum = [np.zeros(n_nodes) for _ in range(n_layers)]

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


def community_est(
    boot_result: BootstrapResult,
) -> dict:
    """Summarize bootstrap community-count reproducibility (meta-communities).

    Reports, for each layer, the meta-community count from the observed
    network together with its *bootstrap reproducibility*: the proportion of
    bootstrap replicates in which the fitted number of communities equals
    the observed-network count. Also returns mean node stability, per-node
    stability, and the co-assignment matrices.

    .. note::
       **Why this is not a confidence interval.** Earlier versions returned
       a percentile ``community_count_ci``. It was replaced in 1.1.0 with a
       reproducibility summary because the interval's coverage is not robust
       to model misspecification. In a large simulation study the nominal
       95% community-count interval covered the truth at or above nominal on
       well-specified planted-partition networks (~0.99 for n >= 100 nodes),
       but coverage collapsed to ~0.62 when community sizes were strongly
       skewed, and no observable diagnostic reliably separated the
       trustworthy cases from the rest. Rather than ship an interval that
       silently undercovers, the function now reports how often the
       community count reproduces under resampling. This is a descriptive
       stability measure, not a calibrated interval: it makes no claim about
       the probability that any range contains the true community count. For
       a validated interval, use :func:`co_assignment_ci`, whose node-pair
       coverage held across the same misspecification stress tests. The raw
       per-replicate community counts are intentionally not exposed anywhere
       in the package output; only this reproducibility summary is returned.

    Parameters
    ----------
    boot_result : BootstrapResult
        Output from ``bootstrap_multilayer``.

    Returns
    -------
    dict
        Dictionary with keys:
        - ``community_count``: DataFrame with layer, estimate (observed-network
          community count), reproducibility (share of bootstrap replicates
          whose community count equals estimate, in [0, 1])
        - ``report``: list of one plain-language sentence per layer
        - ``mean_node_stability``: DataFrame with layer, mean_stability
        - ``node_stability``: list of per-layer stability arrays
        - ``co_assignment``: list of per-layer co-assignment matrices

    See Also
    --------
    co_assignment_ci : Wilson intervals for node-pair co-assignment.
    """
    if boot_result.n_boot == 0:
        raise ValueError("No completed bootstrap replicates.")

    n_layers = len(boot_result.modularity_samples)
    point = boot_result.point_estimate

    # Per-layer community count and bootstrap reproducibility. Reproducibility
    # is precomputed in bootstrap_multilayer from the raw per-replicate counts,
    # which are deliberately not exposed.
    count_rows = []
    for i in range(n_layers):
        est = len(np.unique(point["meta_communities"][i]))
        reproducibility = float(boot_result.community_count_reproducibility[i])
        count_rows.append(
            {"layer": i + 1, "estimate": est, "reproducibility": reproducibility}
        )

    # One plain-language sentence per layer.
    report = [
        f"Layer {r['layer']}: community count (K = {r['estimate']}) "
        f"reproduced in {round(100 * r['reproducibility'])}% of bootstrap resamples."
        for r in count_rows
    ]

    # Mean node stability per layer.
    stab_rows = []
    for i in range(n_layers):
        stab_rows.append(
            {
                "layer": i + 1,
                "mean_stability": float(np.mean(boot_result.node_stability[i])),
            }
        )

    return {
        "community_count": pd.DataFrame(count_rows),
        "report": report,
        "mean_node_stability": pd.DataFrame(stab_rows),
        "node_stability": boot_result.node_stability,
        "co_assignment": boot_result.co_assignment,
    }


def co_assignment_ci(
    boot_result: BootstrapResult,
    alpha: float = 0.05,
    method: str = "wilson",
    calibration_table=None,
) -> list[dict]:
    """Descriptive Monte Carlo interval for node-pair co-assignment.

    For every pair of nodes in every layer, returns the bootstrap
    co-assignment share (fraction of replicates in which the pair landed in
    the same meta-community) with a Wilson score interval that treats the
    replicates as binomial draws.

    .. warning::
       Read this as a diagnostic, not a calibrated confidence interval. In the
       package's simulation study its coverage of the fresh-data co-assignment
       propensity was near nominal only for pairs whose share is close to 0 or
       1 and fell to 0.03-0.15 for ambiguous pairs; no conditioning or
       alternative bootstrap repaired this, so ``method="calibrated"`` was
       removed in 1.3.0. The validated reliability product is
       :func:`partition_stability`.

    Parameters
    ----------
    boot_result : BootstrapResult
    alpha : float
        Significance level (0.05 gives 95 percent Wilson intervals).
    method : str
        ``"wilson"``. ``"calibrated"`` raises an informative error.
    calibration_table
        Ignored; retained for backward compatibility.

    Returns
    -------
    list[dict]
        One dict per layer with ``estimate``, ``lower``, ``upper`` (n x n
        arrays; diagonals are 1) and ``method``.
    """
    if method == "calibrated":
        raise ValueError(
            'method="calibrated" was removed in dynamic_multiplex 1.3.0: pair-level '
            "intervals could not be calibrated for ambiguous pairs. Use "
            "partition_stability() for the calibrated reliability score and the "
            'decided / undetermined pair flags, or method="wilson".'
        )
    if method != "wilson":
        raise ValueError('method must be "wilson"')
    if boot_result.n_boot == 0:
        raise ValueError("No completed bootstrap replicates.")

    from statistics import NormalDist  # stdlib; avoids a scipy dependency

    b = boot_result.n_boot
    z = NormalDist().inv_cdf(1 - alpha / 2)
    z2 = z**2
    layer_cis: list[dict] = []
    for phat in boot_result.co_assignment:
        denom = 1 + z2 / b
        center = (phat + z2 / (2 * b)) / denom
        half = z * np.sqrt(phat * (1 - phat) / b + z2 / (4 * b**2)) / denom
        lower = np.clip(center - half, 0.0, 1.0)
        upper = np.clip(center + half, 0.0, 1.0)
        np.fill_diagonal(lower, 1.0)
        np.fill_diagonal(upper, 1.0)
        layer_cis.append({"estimate": phat, "lower": lower, "upper": upper, "method": "wilson"})
    return layer_cis


# ---------------------------------------------------------------------------
# partition stability: score, calibrated accuracy floor, node floors, pairs
# ---------------------------------------------------------------------------

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


def _load_stability_table(calibration_table=None) -> pd.DataFrame:
    """Bundled default: ``dynamic_multiplex/data/stability_calibration_table.csv``,
    written by ``replication/post/12_stability.R``."""
    if calibration_table is None:
        from importlib.resources import files

        res = files("dynamic_multiplex").joinpath("data", "stability_calibration_table.csv")
        if not res.is_file():
            raise FileNotFoundError(
                "The bundled stability calibration table is missing; pass calibration_table."
            )
        tab = pd.read_csv(res)
    elif isinstance(calibration_table, (str, bytes)):
        tab = pd.read_csv(calibration_table)
    else:
        tab = pd.DataFrame(calibration_table)
    need = ["level", "stab_lo", "stab_hi", "n_calib", "acc_median", "acc_q05"]
    missing = [c for c in need if c not in tab.columns]
    if missing:
        raise ValueError(f"calibration table is missing columns: {missing}")
    return tab


def _floor_lookup(s: float, rows: pd.DataFrame) -> dict:
    if rows.empty:
        return {"floor": float("nan"), "median": float("nan"), "bin": None}
    rows = rows.sort_values("stab_lo").reset_index(drop=True)
    edges = np.append(rows["stab_lo"].to_numpy(), 1.0)
    idx = int(np.clip(np.searchsorted(edges, s, side="right") - 1, 0, len(rows) - 1))
    while idx > 0 and (pd.isna(rows.loc[idx, "acc_q05"]) or rows.loc[idx, "n_calib"] == 0):
        idx -= 1
    return {
        "floor": float(rows.loc[idx, "acc_q05"]),
        "median": float(rows.loc[idx, "acc_median"]),
        "bin": f"[{rows.loc[idx, 'stab_lo']:.1f}, {rows.loc[idx, 'stab_hi']:.1f})",
    }


def partition_stability(
    boot_result: BootstrapResult,
    metric: str = "nmi",
    decided: tuple[float, float] = (0.1, 0.9),
    calibration_table=None,
) -> dict:
    """Stability score and calibrated accuracy floor for the tracked partition.

    Summarises a ``bootstrap_multilayer`` result into one reliability report:

    * ``stability``: mean agreement (NMI by default, or ARI) between each
      bootstrap replicate's meta-partition and the point-estimate partition,
      averaged over layers.
    * ``floor``: calibrated 5th-percentile accuracy for the stability bin.
      In the package's simulation study, fits were binned by stability and
      the 5th percentile of accuracy against the planted partition recorded
      per bin on a calibration half of the configurations; on the held-out
      half, accuracy exceeded the floor in 95 percent of fits in every bin.
      Read as "with this stability, partition accuracy was at least ``floor``
      95 percent of the time in calibration".
    * ``node``: per-layer DataFrames of node-level Jaccard stability with a
      floor reported only above 0.9 (uninformative below).
    * ``pairs``: per-layer integer arrays, 1 = decidedly together
      (co-assignment share >= ``decided[1]``), -1 = decidedly apart
      (<= ``decided[0]``), 0 = undetermined. No interval is attached to pairs.

    The calibration is within the simulated planted-partition family; on
    networks that model does not describe, a stable but wrong partition could
    receive a floor it does not deserve.

    Parameters
    ----------
    boot_result : BootstrapResult
        From ``bootstrap_multilayer`` (1.3.0 or later).
    metric : str
        ``"nmi"`` (default) or ``"ari"``.
    decided : tuple[float, float]
        Thresholds for decidedly apart / together.
    calibration_table
        Optional path or DataFrame overriding the bundled table.

    Returns
    -------
    dict
        Keys ``stability``, ``stability_by_layer``, ``stability_mc_se``,
        ``metric``, ``floor``, ``floor_median``, ``bin``, ``node``, ``pairs``,
        ``pair_summary``, ``report``.
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

    tab = _load_stability_table(calibration_table)
    level = "partition_nmi" if metric == "nmi" else "partition_ari"
    fl = _floor_lookup(s, tab[tab["level"] == level])
    node_tab = tab[tab["level"] == "node_jaccard"]
    node = []
    for t in range(n_layers):
        st = np.asarray(boot_result.node_jaccard_stability[t])
        fr = np.array([_floor_lookup(v, node_tab)["floor"] if v >= 0.9 else np.nan for v in st])
        node.append(pd.DataFrame({"node": np.arange(1, len(st) + 1), "stability": st, "floor": fr}))

    pairs = []
    for P in boot_result.co_assignment:
        M = np.zeros(P.shape, dtype=int)
        M[P >= decided[1]] = 1
        M[P <= decided[0]] = -1
        np.fill_diagonal(M, 1)
        pairs.append(M)
    iu = np.triu_indices(pairs[0].shape[0], k=1)
    allp = np.concatenate([M[iu] for M in pairs])
    pair_summary = {
        "together": float((allp == 1).mean()),
        "apart": float((allp == -1).mean()),
        "undetermined": float((allp == 0).mean()),
    }
    report = (
        f"Partition stability {s:.2f} ({metric.upper()} over {boot_result.n_boot} bootstrap "
        f"replicates and {n_layers} layers; MC s.e. {mc_se:.3f}). In calibration, fits with "
        f"stability in {fl['bin']} had accuracy of at least {fl['floor']:.2f} in 95% of cases "
        f"(median {fl['median']:.2f}). {100 * (1 - pair_summary['undetermined']):.0f}% of node "
        f"pairs are decided ({100 * pair_summary['together']:.0f}% together, "
        f"{100 * pair_summary['apart']:.0f}% apart); {100 * pair_summary['undetermined']:.0f}% "
        "are undetermined."
    )
    return {
        "stability": s, "stability_by_layer": by_layer, "stability_mc_se": mc_se, "metric": metric,
        "floor": fl["floor"], "floor_median": fl["median"], "bin": fl["bin"],
        "node": node, "pairs": pairs, "pair_summary": pair_summary, "report": report,
    }
