from __future__ import annotations

import time

import numpy as np
import pandas as pd

from .fit import (
    fit_multilayer_identity_ties,
    fit_multilayer_jaccard,
    fit_multilayer_overlap,
)
from .multilayer_utils import make_layer_links, weighted_jaccard
from .simulation import simulate_evolving_multilayer


def _community_sets(membership_vec: np.ndarray) -> dict[int, set[int]]:
    out: dict[int, set[int]] = {}
    for idx, cid in enumerate(membership_vec, start=1):
        out.setdefault(int(cid), set()).add(idx)
    return out


def _truth_interlayer_ties(membership: np.ndarray, links: pd.DataFrame, min_similarity: float = 0.05) -> pd.DataFrame:
    rows = []
    for _, link in links.iterrows():
        f, t, w = int(link["from"]), int(link["to"]), float(link["weight"])
        c_from = _community_sets(membership[:, f - 1])
        c_to = _community_sets(membership[:, t - 1])
        for cf, sf in c_from.items():
            for ct, st in c_to.items():
                sim = weighted_jaccard(sf, st)
                ws = sim * w
                if ws >= min_similarity:
                    rows.append({"from_layer": f, "to_layer": t, "from_community": cf, "to_community": ct, "weighted_similarity": ws})
    return pd.DataFrame(rows)


def _score(pred: pd.DataFrame, truth: pd.DataFrame, threshold: float = 0.05) -> dict:
    def keys(df: pd.DataFrame) -> set[str]:
        if df.empty:
            return set()
        return set((df["from_layer"].astype(str) + "::" + df["to_layer"].astype(str) + "::" + df["from_community"].astype(str) + "::" + df["to_community"].astype(str)).tolist())

    p = pred[pred["weighted_similarity"] >= threshold] if not pred.empty else pred
    t = truth[truth["weighted_similarity"] >= threshold] if not truth.empty else truth

    pk = keys(p)
    tk = keys(t)

    tp = len(pk.intersection(tk))
    fp = len(pk - tk)
    fn = len(tk - pk)

    precision = 0.0 if tp + fp == 0 else tp / (tp + fp)
    recall = 0.0 if tp + fn == 0 else tp / (tp + fn)
    f1 = 0.0 if precision + recall == 0 else 2 * precision * recall / (precision + recall)

    return {
        "tp": tp,
        "fp": fp,
        "fn": fn,
        "precision": precision,
        "recall": recall,
        "f1": f1,
        "predicted_ties": len(pk),
        "truth_ties": len(tk),
    }


def _adjusted_rand_index(x: list[int], y: list[int]) -> float:
    x = np.asarray(x)
    y = np.asarray(y)
    n = len(x)
    if n <= 1:
        return float("nan")

    _, xi = np.unique(x, return_inverse=True)
    _, yi = np.unique(y, return_inverse=True)
    cont = np.zeros((xi.max() + 1, yi.max() + 1), dtype=int)
    for i in range(n):
        cont[xi[i], yi[i]] += 1

    def comb2(v: np.ndarray) -> float:
        return float(np.sum(v * (v - 1) / 2))

    a = comb2(cont)
    b = comb2(cont.sum(axis=1))
    c = comb2(cont.sum(axis=0))
    d = n * (n - 1) / 2
    expected = (b * c) / d if d else 0
    max_idx = (b + c) / 2
    denom = max_idx - expected
    if denom == 0:
        return 0.0
    return float((a - expected) / denom)


def _mean_layer_ari(fit: dict, truth_membership: np.ndarray) -> float:
    vals = []
    for t in range(truth_membership.shape[1]):
        vals.append(_adjusted_rand_index(fit["layer_communities"][t]["membership"], truth_membership[:, t]))
    return float(np.nanmean(vals))


def _mean_layer_modularity(fit: dict) -> float:
    vals = [m["modularity"] for m in fit["layer_communities"] if m["modularity"] is not None]
    return float(np.nan if len(vals) == 0 else np.mean(vals))


def benchmark_interlayer_tie_recovery(
    n_reps: int = 50,
    simulation_args: dict | None = None,
    algorithm: str = "louvain",
    min_similarity: float = 0.05,
    threshold: float = 0.05,
    seed: int | None = None,
):
    rng = np.random.default_rng(seed)
    sim_defaults = dict(n_nodes=150, n_layers=6, n_communities=4, persistence=0.85, p_in=0.20, p_out=0.03, directed=False)
    sim_par = {**sim_defaults, **(simulation_args or {})}

    one_step = make_layer_links(sim_par["n_layers"])
    all_pairs = pd.DataFrame([(f, t, 1.0) for f in range(1, sim_par["n_layers"]) for t in range(f + 1, sim_par["n_layers"] + 1)], columns=["from", "to", "weight"])

    rows = []
    for rep in range(1, n_reps + 1):
        sim = simulate_evolving_multilayer(**sim_par, seed=int(rng.integers(0, 1_000_000_000)))
        truth = _truth_interlayer_ties(sim["membership"], one_step, min_similarity=min_similarity)

        t0 = time.perf_counter()
        fit_j = fit_multilayer_jaccard(
            sim["layers"],
            algorithm=algorithm,
            layer_links=one_step,
            min_similarity=min_similarity,
            directed=sim_par["directed"],
        )
        tj = time.perf_counter() - t0

        t0 = time.perf_counter()
        fit_o = fit_multilayer_overlap(
            sim["layers"],
            algorithm=algorithm,
            layer_links=one_step,
            min_similarity=min_similarity,
            directed=sim_par["directed"],
        )
        to = time.perf_counter() - t0

        t0 = time.perf_counter()
        fit_i = fit_multilayer_identity_ties(
            sim["layers"],
            algorithm=algorithm,
            layer_links=one_step,
            directed=sim_par["directed"],
        )
        ti = time.perf_counter() - t0

        t0 = time.perf_counter()
        fit_a = fit_multilayer_jaccard(
            sim["layers"],
            algorithm=algorithm,
            layer_links=all_pairs,
            min_similarity=min_similarity,
            directed=sim_par["directed"],
        )
        ta = time.perf_counter() - t0

        identity = fit_i["interlayer_ties"].copy()
        identity["from_community"] = [fit_i["layer_communities"][fl - 1]["membership"][n - 1] + 1 for fl, n in zip(identity["from_layer"], identity["node"])]
        identity["to_community"] = [fit_i["layer_communities"][tl - 1]["membership"][n - 1] + 1 for tl, n in zip(identity["to_layer"], identity["node"])]
        identity["weighted_similarity"] = identity["layer_weight"]
        identity = identity[["from_layer", "to_layer", "from_community", "to_community", "weighted_similarity"]].drop_duplicates()

        preds = {
            "bayesnet_jaccard": fit_j["interlayer_ties"][ ["from_layer", "to_layer", "from_community", "to_community", "weighted_similarity"] ],
            "bayesnet_overlap": fit_o["interlayer_ties"][ ["from_layer", "to_layer", "from_community", "to_community", "weighted_similarity"] ],
            "bayesnet_identity": identity,
            "multinet_all_to_all": fit_a["interlayer_ties"][ ["from_layer", "to_layer", "from_community", "to_community", "weighted_similarity"] ],
        }

        meta = {
            "bayesnet_jaccard": (tj, fit_j),
            "bayesnet_overlap": (to, fit_o),
            "bayesnet_identity": (ti, fit_i),
            "multinet_all_to_all": (ta, fit_a),
        }

        for model, df in preds.items():
            sc = _score(df, truth, threshold=threshold)
            elapsed, fit_obj = meta[model]
            sc.update(
                model=model,
                replicate=rep,
                elapsed_sec=elapsed,
                f1_per_second=np.nan if elapsed <= 0 else sc["f1"] / elapsed,
                mean_layer_ari=_mean_layer_ari(fit_obj, sim["membership"]),
                mean_layer_modularity=_mean_layer_modularity(fit_obj),
            )
            rows.append(sc)

    raw = pd.DataFrame(rows)
    summary = raw.groupby("model", as_index=False).mean(numeric_only=True)
    return {"raw_results": raw, "summary": summary}
