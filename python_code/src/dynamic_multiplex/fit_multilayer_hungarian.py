from __future__ import annotations

import numpy as np

from .multilayer_utils import (
    _sort_key,
    fit_layer_communities,
    make_layer_links,
    prepare_multilayer_graphs,
)


def _assign_max_overlap(overlap: np.ndarray) -> np.ndarray:
    """Return, for each row, the column it is matched to so that total overlap
    is maximised (Hungarian / optimal assignment). Falls back to a greedy
    matcher if SciPy is unavailable."""
    d = max(overlap.shape)
    square = np.zeros((d, d), dtype=int)
    square[: overlap.shape[0], : overlap.shape[1]] = overlap
    try:
        from scipy.optimize import linear_sum_assignment

        rows, cols = linear_sum_assignment(square.max() - square)
        col_for_row = np.empty(d, dtype=int)
        col_for_row[rows] = cols
        return col_for_row
    except Exception:  # pragma: no cover - greedy fallback
        used: set[int] = set()
        col_for_row = np.full(d, -1, dtype=int)
        for i in np.argsort(-square.max(axis=1)):
            for c in np.argsort(-square[i]):
                c = int(c)
                if c not in used:
                    col_for_row[i] = c
                    used.add(c)
                    break
        return col_for_row


def _match_communities_hungarian(memberships: list[dict]) -> list[np.ndarray]:
    """Align per-layer memberships by Hungarian matching of consecutive layers
    on the community overlap (contingency) matrix. Births get fresh labels.

    Each membership is a dict ``{node_id: community}``; overlap is counted on
    the node ids present in both layers, so layers with different node sets
    are matched actor by actor. Returns one array per layer in sorted node-id
    order (the order used by ``meta_communities`` elsewhere in the package).
    """
    orders = [sorted(m.keys(), key=_sort_key) for m in memberships]
    out: list[dict] = [dict(memberships[0])]
    first = np.array([out[0][n] for n in orders[0]], dtype=int)
    next_free = int(first.max()) + 1 if first.size else 1

    for t in range(1, len(memberships)):
        prev = out[t - 1]
        cur = memberships[t]
        cl = np.unique(np.array(list(cur.values()), dtype=int))
        pl = np.unique(np.array(list(prev.values()), dtype=int))

        shared = [n for n in orders[t] if n in prev]
        overlap = np.zeros((len(cl), len(pl)), dtype=int)
        if shared:
            cs = np.array([cur[n] for n in shared], dtype=int)
            ps = np.array([prev[n] for n in shared], dtype=int)
            for i, ci in enumerate(cl):
                for j, pj in enumerate(pl):
                    overlap[i, j] = int(np.sum((cs == ci) & (ps == pj)))

        assignment = _assign_max_overlap(overlap)

        relabel: dict[int, int] = {}
        for i, ci in enumerate(cl):
            col = int(assignment[i])
            if col < len(pl) and overlap[i, col] > 0:
                relabel[int(ci)] = int(pl[col])
        for ci in cl:
            if int(ci) not in relabel:
                relabel[int(ci)] = next_free
                next_free += 1
        next_free = max(next_free, max(relabel.values()) + 1)

        out.append({n: relabel[int(v)] for n, v in cur.items()})

    return [np.array([out[t][n] for n in orders[t]], dtype=int) for t in range(len(out))]


def fit_multilayer_hungarian(
    layers,
    algorithm: str = "leiden",
    resolution_parameter: float = 1.0,
    directed: bool = False,
    objective: str | None = None,
    seed: int | None = 123,
    allow_unequal_nodes: bool = False,
):
    """Fit multilayer communities with Hungarian snapshot matching.

    Detects communities independently in each layer, then tracks them across
    time by matching community labels between consecutive layers with the
    Hungarian (optimal assignment) algorithm on the community overlap matrix.
    Unlike the coupling-based fits (``fit_multilayer_jaccard``,
    ``fit_multilayer_overlap``, ``fit_multilayer_identity_ties``), this is a
    two-stage snapshot-and-match tracker: communities are found per layer and
    aligned post hoc, not jointly optimised on a coupled supra-graph.

    Parameters
    ----------
    allow_unequal_nodes : bool
        If True, layers may hold different node sets; they must then be
        networkx graphs with node labels. Community overlap is counted on the
        nodes present in both layers of each consecutive pair. Default False.

    Returns
    -------
    dict
        Keys ``layer_communities`` (per-layer detection), ``meta_communities``
        (one array per layer with labels aligned across consecutive layers by
        Hungarian matching; see ``extract_meta_membership``), ``meta_ids``,
        ``interlayer_ties`` (``None`` -- no supra-graph is built), ``method``
        (``"hungarian"``), and ``layer_links``.
    """
    graph_layers = prepare_multilayer_graphs(
        layers, directed=directed, require_same_nodes=not allow_unequal_nodes
    )
    links = make_layer_links(len(graph_layers), None)

    fit = fit_layer_communities(
        graph_layers,
        algorithm=algorithm,
        resolution_parameter=resolution_parameter,
        directed=directed,
        objective=objective,
        seed=seed,
    )

    meta = _match_communities_hungarian([dict(f.membership) for f in fit])
    meta_ids = sorted({int(v) for layer in meta for v in layer})

    return {
        "algorithm": algorithm,
        "layer_communities": fit,
        "meta_communities": meta,
        "meta_ids": meta_ids,
        "layer_links": links,
        "interlayer_ties": None,
        "method": "hungarian",
        "directed": directed,
        "node_labels": graph_layers[0].graph.get("node_labels"),
    }
