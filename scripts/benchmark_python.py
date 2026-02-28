"""
Benchmark dynamic_multiplex against competing temporal community detection approaches.

Compares accuracy (NMI, ARI) and runtime on synthetic planted-partition
networks with evolving community structure.

Competitors:
  1. dynamic_multiplex (adjacent-layer coupling) -- our method
  2. Independent Louvain (no temporal coupling)
  3. Full-coupling Louvain (all layer pairs connected)
  4. Aggregated Louvain (collapse layers into one network)

Usage:
    pip install -e ./python_code[louvain,dev]
    pip install scikit-learn
    python scripts/benchmark_python.py
"""

from __future__ import annotations

import itertools
import time
from collections import defaultdict
from math import log

import numpy as np
import pandas as pd

try:
    from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
except ImportError:
    # Fallback NMI/ARI implementations
    def _entropy(labels):
        n = len(labels)
        if n == 0:
            return 0.0
        counts = defaultdict(int)
        for lb in labels:
            counts[lb] += 1
        return -sum((c / n) * log(c / n) for c in counts.values())

    def normalized_mutual_info_score(a, b):
        n = len(a)
        if n == 0:
            return 0.0
        joint = defaultdict(int)
        ca, cb = defaultdict(int), defaultdict(int)
        for x, y in zip(a, b):
            joint[(x, y)] += 1
            ca[x] += 1
            cb[y] += 1
        mi = 0.0
        for (x, y), nxy in joint.items():
            mi += (nxy / n) * log((nxy * n) / (ca[x] * cb[y]))
        ha, hb = _entropy(a), _entropy(b)
        denom = (ha + hb) / 2
        return mi / denom if denom > 0 else 0.0

    def adjusted_rand_score(a, b):
        n = len(a)
        cont = defaultdict(int)
        ra, rb = defaultdict(int), defaultdict(int)
        for x, y in zip(a, b):
            cont[(x, y)] += 1
            ra[x] += 1
            rb[y] += 1
        comb2 = lambda k: k * (k - 1) / 2
        index = sum(comb2(v) for v in cont.values())
        sum_a = sum(comb2(v) for v in ra.values())
        sum_b = sum(comb2(v) for v in rb.values())
        expected = sum_a * sum_b / comb2(n) if comb2(n) > 0 else 0
        max_idx = (sum_a + sum_b) / 2
        denom = max_idx - expected
        return (index - expected) / denom if denom > 0 else 0.0


import community as community_louvain
import networkx as nx

from dynamic_multiplex import (
    fit_multilayer_jaccard,
    fit_multilayer_overlap,
    fit_multilayer_weighted_jaccard,
    fit_multilayer_weighted_overlap,
    fit_multilayer_identity_ties,
)


# ---------------------------------------------------------------------------
# Simulation
# ---------------------------------------------------------------------------

def simulate_evolving_networks(
    n_nodes: int,
    n_layers: int,
    n_communities: int,
    p_in: float = 0.3,
    p_out: float = 0.05,
    p_switch: float = 0.05,
    seed: int | None = None,
):
    """Generate temporal networks with slowly evolving planted partition."""
    rng = np.random.default_rng(seed)

    memberships = [rng.integers(1, n_communities + 1, size=n_nodes)]
    for _ in range(1, n_layers):
        prev = memberships[-1].copy()
        switch_mask = rng.random(n_nodes) < p_switch
        prev[switch_mask] = rng.integers(1, n_communities + 1, size=switch_mask.sum())
        memberships.append(prev)

    layers = []
    for t in range(n_layers):
        mat = np.zeros((n_nodes, n_nodes))
        for i in range(n_nodes - 1):
            for j in range(i + 1, n_nodes):
                prob = p_in if memberships[t][i] == memberships[t][j] else p_out
                tie = rng.binomial(1, prob)
                mat[i, j] = tie
                mat[j, i] = tie
        layers.append(mat)

    return layers, memberships


# ---------------------------------------------------------------------------
# Competing methods
# ---------------------------------------------------------------------------

def _memberships_from_fit(fit_result, n_layers):
    """Extract membership dicts from a dynamic_multiplex fit."""
    return [fit_result["layer_communities"][i].membership for i in range(n_layers)]


def run_dynamic_multiplex(layers, fit_type="jaccard"):
    """Our method: adjacent-layer coupling."""
    fit_fn = {
        "jaccard": fit_multilayer_jaccard,
        "overlap": fit_multilayer_overlap,
        "weighted_jaccard": fit_multilayer_weighted_jaccard,
        "weighted_overlap": fit_multilayer_weighted_overlap,
        "identity": fit_multilayer_identity_ties,
    }[fit_type]
    fit = fit_fn(layers, algorithm="louvain")
    return _memberships_from_fit(fit, len(layers))


def run_independent_louvain(layers):
    """Baseline: Louvain on each layer independently (no temporal coupling)."""
    results = []
    for mat in layers:
        g = nx.from_numpy_array(mat)
        partition = community_louvain.best_partition(g)
        results.append({k + 1: v + 1 for k, v in partition.items()})
    return results


def run_full_coupling(layers):
    """All-to-all layer coupling via Jaccard."""
    n = len(layers)
    links = [{"from": i + 1, "to": j + 1, "weight": 1.0}
             for i in range(n) for j in range(i + 1, n)]
    fit = fit_multilayer_jaccard(layers, algorithm="louvain", layer_links=links)
    return _memberships_from_fit(fit, n)


def run_aggregated_louvain(layers):
    """Collapse all layers into a single network."""
    agg = sum(layers)
    g = nx.from_numpy_array(agg)
    partition = community_louvain.best_partition(g)
    membership = {k + 1: v + 1 for k, v in partition.items()}
    return [membership] * len(layers)


# ---------------------------------------------------------------------------
# Evaluation
# ---------------------------------------------------------------------------

def evaluate(detected_memberships, true_memberships, n_nodes):
    """Compute mean NMI and ARI across layers."""
    nmis, aris = [], []
    for det, true in zip(detected_memberships, true_memberships):
        det_vec = [det.get(i + 1, 0) for i in range(n_nodes)]
        true_vec = list(true)
        nmis.append(normalized_mutual_info_score(true_vec, det_vec))
        aris.append(adjusted_rand_score(true_vec, det_vec))
    return float(np.mean(nmis)), float(np.mean(aris))


# ---------------------------------------------------------------------------
# Main benchmark
# ---------------------------------------------------------------------------

def run_benchmark():
    configs = list(itertools.product(
        [50, 100],           # n_nodes
        [5, 10],             # n_layers
        [3, 5],              # n_communities
        [0.0, 0.05, 0.15],   # p_switch
    ))

    methods = {
        "DynMux_Jaccard": lambda L: run_dynamic_multiplex(L, "jaccard"),
        "DynMux_Overlap": lambda L: run_dynamic_multiplex(L, "overlap"),
        "DynMux_WtJaccard": lambda L: run_dynamic_multiplex(L, "weighted_jaccard"),
        "DynMux_WtOverlap": lambda L: run_dynamic_multiplex(L, "weighted_overlap"),
        "DynMux_Identity": lambda L: run_dynamic_multiplex(L, "identity"),
        "Independent": run_independent_louvain,
        "FullCoupling": run_full_coupling,
        "Aggregated": run_aggregated_louvain,
    }

    n_reps = 5
    results = []

    total = len(configs) * n_reps
    print(f"Running {total} simulation scenarios x {len(methods)} methods ...")

    for idx, (n_nodes, n_layers, n_comms, p_switch) in enumerate(configs):
        for rep in range(n_reps):
            seed = idx * 1000 + rep
            layers, true_mem = simulate_evolving_networks(
                n_nodes=n_nodes, n_layers=n_layers, n_communities=n_comms,
                p_in=0.3, p_out=0.05, p_switch=p_switch, seed=seed,
            )

            for method_name, method_fn in methods.items():
                t0 = time.perf_counter()
                try:
                    detected = method_fn(layers)
                except Exception as e:
                    print(f"  SKIP {method_name}: {e}")
                    continue
                elapsed = time.perf_counter() - t0
                nmi, ari = evaluate(detected, true_mem, n_nodes)

                results.append({
                    "n_nodes": n_nodes,
                    "n_layers": n_layers,
                    "n_communities": n_comms,
                    "p_switch": p_switch,
                    "rep": rep,
                    "method": method_name,
                    "nmi": round(nmi, 4),
                    "ari": round(ari, 4),
                    "runtime_s": round(elapsed, 4),
                })

            done = idx * n_reps + rep + 1
            if done % 5 == 0 or done == total:
                print(f"  [{done}/{total}]")

    df = pd.DataFrame(results)
    outpath = "scripts/benchmark_python_results.csv"
    df.to_csv(outpath, index=False)
    print(f"\nResults saved to {outpath}")

    # Summary table
    summary = (
        df.groupby(["method", "p_switch"])
        .agg(mean_nmi=("nmi", "mean"), mean_ari=("ari", "mean"),
             mean_runtime=("runtime_s", "mean"))
        .round(3)
        .reset_index()
    )
    print("\n=== Summary (mean across configs and reps) ===")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    run_benchmark()
