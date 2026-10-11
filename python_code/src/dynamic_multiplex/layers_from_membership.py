"""Two-mode (membership) input: per-slice weighted layers from group memberships.

Membership panels (states in IGOs, legislators in caucuses) are usually
projected to a node-by-node network whose tie weight counts shared groups.
``layers_from_membership`` builds that projection, or one of two size-aware
alternatives, slice by slice, and returns layers the ``fit_multilayer_*``
functions accept. No weighting is safe by default; see the function notes.
"""
from __future__ import annotations

import warnings

import networkx as nx
import numpy as np
import pandas as pd

from .multilayer_utils import fit_layer_communities

_WEIGHTINGS = ("count", "size", "kumar")


def _as_wide(entry, t: int) -> pd.DataFrame:
    """Validated wide 0/1 membership table (rows nodes, columns groups) for slice t."""
    if not isinstance(entry, pd.DataFrame):
        raise ValueError(
            f"In slice {t}, the membership is not a pandas DataFrame. Pass a wide 0/1 table "
            "with node labels as the index, or a long table with columns 'node' and 'group' "
            "(arrays carry no row labels)."
        )
    if len(entry) == 0:
        raise ValueError(f"In slice {t}, the membership table is empty.")
    cols = list(entry.columns)
    if "node" in cols or "group" in cols:
        if sorted(cols) != ["group", "node"]:
            raise ValueError(f"In slice {t}, a long table must have exactly the columns 'node' and 'group'.")
        if entry.isna().to_numpy().any():
            raise ValueError(f"In slice {t}, the long table has missing values.")
        dup = entry.duplicated(["node", "group"])
        if dup.any():
            raise ValueError(f"In slice {t}, the long table has {int(dup.sum())} duplicate (node, group) rows.")
        e = entry.reset_index(drop=True)              # the long table's own index is irrelevant
        return pd.crosstab(e["node"], e["group"])
    if entry.index.hasnans:
        raise ValueError(f"In slice {t}, some row labels (node labels) are missing.")
    if entry.index.duplicated().any():
        raise ValueError(f"In slice {t}, there are duplicate row labels: {list(entry.index[entry.index.duplicated()][:10])}.")
    if entry.columns.duplicated().any():
        raise ValueError(f"In slice {t}, there are duplicate column labels: {list(entry.columns[entry.columns.duplicated()][:10])}.")
    if entry.isna().to_numpy().any():
        raise ValueError(f"In slice {t}, the membership table has missing values.")
    vals = entry.to_numpy()
    if not np.isin(vals, [0, 1]).all():
        raise ValueError(f"In slice {t}, membership entries must be 0/1.")
    return entry


def _prepare_slice(entry, t: int, min_group_size: int):
    """Validated slice: node labels, kept groups, 0/1 matrix B, group sizes, groups dropped.

    Groups below ``min_group_size`` are dropped with a warning; a node left
    with no group is an error.
    """
    wide = _as_wide(entry, t)
    labels = list(wide.index)
    B_all = wide.to_numpy(dtype=float)
    size_all = B_all.sum(axis=0)
    keep = size_all >= min_group_size
    n_drop = int((~keep).sum())
    if n_drop:
        warnings.warn(f"In slice {t}, dropped {n_drop} group(s) with fewer than {min_group_size} members.",
                      UserWarning, stacklevel=3)
    B = B_all[:, keep]
    size = size_all[keep]
    groups = [g for g, k in zip(wide.columns, keep) if k]
    orphan = B.sum(axis=1) == 0
    if orphan.any():
        bad = [labels[i] for i in np.flatnonzero(orphan)[:10]]
        raise ValueError(
            f"In slice {t}, {int(orphan.sum())} node(s) belong to no group of at least {min_group_size} "
            f"members, for example {bad}. Drop absent nodes from the slice before building layers."
        )
    return labels, groups, B, size, n_drop


def _kumar_targets(B: np.ndarray, labels: np.ndarray) -> np.ndarray:
    """Kumar reweighting targets r_j for a fixed node partition.

    r_j = [sum over the C communities of 1 / (1 + n_jc)] * (s_j + C) / M, with
    n_jc the members of group j in community c, s_j the size of group j and M
    the number of groups (hypernetx 2.4.3 ``hypergraph_modularity.kumar``).
    """
    B = np.asarray(B, dtype=float)
    _, inv = np.unique(np.asarray(labels), return_inverse=True)
    C = int(inv.max()) + 1
    H = np.zeros((B.shape[0], C))
    H[np.arange(B.shape[0]), inv] = 1.0
    n = B.T @ H                                          # M x C member counts
    return (1.0 / (1.0 + n)).sum(axis=1) * (B.sum(axis=0) + C) / B.shape[1]


def _projection(B: np.ndarray, group_weights: np.ndarray, size: np.ndarray, weighting: str) -> np.ndarray:
    if weighting == "count":
        W = B @ B.T
    else:
        W = (B * (group_weights / (size - 1.0))) @ B.T
    np.fill_diagonal(W, 0.0)
    return W


def _upper_edges(W: np.ndarray):
    iu, ju = np.nonzero(np.triu(W, k=1))
    return iu, ju, W[iu, ju]


def _cluster(W: np.ndarray, algorithm: str, resolution_parameter: float, seed) -> np.ndarray:
    """Same stage-1 routine (and settings) the fitters use, on one weighted slice."""
    g = nx.Graph()
    g.add_nodes_from(range(W.shape[0]))
    iu, ju, w = _upper_edges(W)
    g.add_weighted_edges_from(zip(iu.tolist(), ju.tolist(), w.tolist()))
    g.graph["dm_raw_ids"] = True                        # membership keyed by node index
    fit = fit_layer_communities([g], algorithm=algorithm, resolution_parameter=resolution_parameter,
                                directed=False, seed=seed)[0]
    return np.array([fit.membership[i] for i in range(W.shape[0])])


def layers_from_membership(
    memberships: list,
    weighting: str,
    algorithm: str | None = None,
    resolution_parameter: float = 1.0,
    seed: int | None = 123,
    max_iter: int = 50,
    tol: float = 0.01,
    min_group_size: int = 2,
    rescale: str = "none",
) -> list[nx.Graph]:
    """Build one weighted layer per time slice from group memberships.

    Let B be the N x M 0/1 membership matrix of a slice and s_j the size of
    group j.

    * ``"count"``: W = B Bᵀ (number of shared groups), diagonal 0.
    * ``"size"``: W_ik = sum_j B_ij B_kj / (s_j - 1), so each group adds total
      weight proportional to its size rather than its size squared.
    * ``"kumar"``: Kumar et al.'s hypergraph reweighting (as in hypernetx
      2.4.3). Group weights w_j start at 1; the slice is projected as in
      ``"size"`` with weights w_j and clustered with ``algorithm``; each group
      gets the target r_j = [sum_c 1 / (1 + n_jc)] (s_j + C) / M and
      w_j <- (w_j + r_j) / 2, until max_j |w_j - r_j| / 2 <= ``tol``. The
      projection with the final w_j is returned (the reference returns only
      the partition). The clustering calls the same per-layer routine as the
      fitters, so ``algorithm`` and ``resolution_parameter`` must be the ones
      you pass to the fitter.

    Every slice is built on its own, with the same ``seed``: the weights of
    slice t never depend on any other slice, so no later information enters
    an earlier layer.

    No weighting is safe in general. In small Python simulations (hypotheses,
    not established results): the count projection failed badly when a few
    very large groups were present, where size-weighted and Kumar both
    recovered the planted blocs; Kumar beat size-weighting only under noise
    groups of the same size as real ones; under small groups bridging blocs,
    count was best and both size-aware weightings failed; and stage-1 Louvain
    reproduced Kumar's own partition less reliably than Leiden. Nothing has
    been checked on real data beyond a pilot.

    Parameters
    ----------
    memberships : list
        One entry per slice: a wide 0/1 ``pandas.DataFrame`` (index = node
        labels, columns = groups) or a long ``DataFrame`` with exactly the
        columns ``"node"`` and ``"group"``.
    weighting : str
        ``"count"``, ``"size"`` or ``"kumar"``. Required; there is no default.
    algorithm : str or None
        ``"leiden"`` or ``"louvain"``; required for ``"kumar"`` and not
        accepted otherwise.
    resolution_parameter : float
        Resolution for Kumar's clustering (> 0).
    seed : int or None
        Seed for Kumar's clustering, reused for every slice. Default 123.
    max_iter : int
        Maximum number of Kumar reweighting passes (default 50). Failing to
        converge raises an error.
    tol : float
        Convergence tolerance on the largest half-change in a group weight.
    min_group_size : int
        Groups with fewer members (default 2) are dropped, with a warning.
    rescale : str
        ``"none"`` (default) or ``"mean"`` (each slice divided by its mean
        nonzero tie). Matters only where tie scale matters, such as the
        multislice coupling ``omega``; modularity-based two-stage fits are
        scale-invariant.

    Returns
    -------
    list[networkx.Graph]
        One weighted undirected graph per slice with the node labels of that
        slice. Pass ``allow_unequal_nodes=True`` to the fitter when node sets
        differ across slices. ``graph.graph["dm_membership"]`` holds
        ``weighting``, ``group_weights`` (final w_j by group), ``iterations``,
        ``final_diff``, ``groups_dropped``, ``rescale`` and ``partition``
        (Kumar's final partition, node label -> community; None otherwise).

    Notes
    -----
    Projected ties are not independent across pairs (one group creates many
    ties), so the package's bootstrap stability score has not been checked
    for layers built this way.
    """
    if not isinstance(memberships, list):
        raise ValueError("`memberships` must be a list with one entry per time slice.")
    if len(memberships) == 0:
        raise ValueError("`memberships` must contain at least one slice.")
    if weighting not in _WEIGHTINGS:
        raise ValueError('`weighting` must be "count", "size" or "kumar"; there is no default.')
    if weighting == "kumar":
        if algorithm not in ("leiden", "louvain"):
            raise ValueError('weighting="kumar" needs `algorithm` ("leiden" or "louvain"), the one you pass to the fitter.')
    elif algorithm is not None:
        raise ValueError('`algorithm` is used only by weighting="kumar"; leave it as None.')
    if isinstance(resolution_parameter, bool) or not isinstance(resolution_parameter, (int, float)) \
            or not resolution_parameter > 0:
        raise ValueError("`resolution_parameter` must be a number > 0.")
    if isinstance(max_iter, bool) or not isinstance(max_iter, (int, np.integer)) or max_iter < 1:
        raise ValueError("`max_iter` must be an integer >= 1.")
    if isinstance(tol, bool) or not isinstance(tol, (int, float)) or not tol > 0:
        raise ValueError("`tol` must be a number > 0.")
    if isinstance(min_group_size, bool) or not isinstance(min_group_size, (int, np.integer)) or min_group_size < 2:
        raise ValueError("`min_group_size` must be an integer >= 2.")
    if rescale not in ("none", "mean"):
        raise ValueError('`rescale` must be "none" or "mean".')

    out = []
    for t, entry in enumerate(memberships, start=1):
        labels, groups, B, size, n_drop = _prepare_slice(entry, t, min_group_size)

        w = np.ones(B.shape[1])
        iterations, final_diff, partition = 0, None, None
        if weighting == "kumar":
            comm = _cluster(_projection(B, w, size, "size"), algorithm, resolution_parameter, seed)
            while True:
                r = _kumar_targets(B, comm)
                final_diff = float(np.max(0.5 * np.abs(w - r)))
                w = 0.5 * w + 0.5 * r
                iterations += 1
                comm = _cluster(_projection(B, w, size, "size"), algorithm, resolution_parameter, seed)
                if final_diff <= tol:
                    break
                if iterations >= max_iter:
                    raise RuntimeError(
                        f"In slice {t}, Kumar reweighting did not converge within max_iter = {max_iter} passes "
                        f"(largest half-change {final_diff:.3g} > tol = {tol})."
                    )
            partition = {labels[i]: int(comm[i]) for i in range(len(labels))}

        W = _projection(B, w, size, "count" if weighting == "count" else "size")
        iu, ju, vals = _upper_edges(W)
        if not len(vals) or not vals.sum() > 0:
            raise ValueError(f"In slice {t}, the projection has zero total weight.")
        if rescale == "mean":
            vals = vals / vals.mean()

        g = nx.Graph()
        g.add_nodes_from(labels)
        g.add_weighted_edges_from((labels[i], labels[j], float(v)) for i, j, v in zip(iu, ju, vals))
        g.graph["dm_membership"] = {
            "weighting": weighting,
            "group_weights": pd.Series(w, index=groups, name="weight"),
            "iterations": iterations,
            "final_diff": final_diff,
            "groups_dropped": n_drop,
            "rescale": rescale,
            "partition": partition,
        }
        out.append(g)
    return out
