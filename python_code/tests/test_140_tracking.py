"""1.4.0: chain repair, filtered tracking, Louvain resolution. Seed 123.

Written before the implementation. Fixtures
- twostage_regression_py.json: default output of the four two-stage fitters
  recorded from 1.3.2 (commit b55c66c).
- chain_panel_edges.csv: 60 nodes, 3 blocks of 20, 10 layers, block 1 splits
  at layer 6 (shared with the R tests).
"""
import importlib.util
import itertools
import json
import warnings

import networkx as nx
import numpy as np
import pandas as pd
import pytest

from _helpers_140 import FIX, chain_layers, sbm_layers
from dynamic_multiplex import (
    bootstrap_multilayer,
    fit_multilayer_jaccard,
    fit_multilayer_overlap,
    fit_multilayer_weighted_jaccard,
    fit_multilayer_weighted_overlap,
)
from dynamic_multiplex.multilayer_utils import fit_layer_communities, prepare_multilayer_graphs, repair_meta_chains

pytestmark = pytest.mark.skipif(
    importlib.util.find_spec("leidenalg") is None or importlib.util.find_spec("community") is None,
    reason="needs leidenalg and python-louvain",
)

FITTERS = {
    "jaccard": fit_multilayer_jaccard,
    "overlap": fit_multilayer_overlap,
    "weighted_jaccard": fit_multilayer_weighted_jaccard,
    "weighted_overlap": fit_multilayer_weighted_overlap,
}
ALGS = ["leiden", "louvain"]


def ids(meta):
    return np.unique(np.concatenate([np.asarray(m) for m in meta]))


def same_partition(a, b):
    """Equal as partitions of node-layers (labels may differ)."""
    x = np.concatenate([np.asarray(m) for m in a]); y = np.concatenate([np.asarray(m) for m in b])
    return np.array_equal(x[:, None] == x[None, :], y[:, None] == y[None, :])


# 1. defaults reproduce 1.3.2 exactly ------------------------------------------
def test_default_output_matches_132():
    rec = json.load(open(FIX / "twostage_regression_py.json"))
    assert rec["recorded_from"] == "b55c66c"
    panels = {"sbm": sbm_layers(), "chain": chain_layers()}
    for key, expected in rec["cases"].items():
        panel, fn, alg = key.split("|")
        for kwargs in ({}, {"repair_chains": False, "tracking": "smoothed"}):
            got = FITTERS[fn](panels[panel], algorithm=alg, seed=123, **kwargs)["meta_communities"]
            assert [m.tolist() for m in got] == expected, (key, kwargs)


# 2. the repair rule on a hand-built case (same case and answers in the R tests)
HAND_META = {"L1C1": 1, "L2C1": 1, "L3C1": 2, "L4C1": 2,
             "L1C2": 3, "L2C2": 3, "L3C2": 4, "L4C2": 4,
             "L1C3": 6, "L3C3": 5}
HAND_TIES = pd.DataFrame(
    [(1, 2, 1, 1, .8), (1, 2, 2, 2, .9), (1, 2, 3, 1, .9),
     (2, 3, 1, 1, .7), (2, 3, 2, 2, .15), (2, 3, 2, 3, .1),
     (3, 4, 1, 1, .9), (3, 4, 2, 2, .6), (3, 4, 3, 2, .3)],
    columns=["from_layer", "to_layer", "from_community", "to_community", "similarity"])
HAND_LINKS = pd.DataFrame({"from": [1, 2, 3], "to": [2, 3, 4], "weight": 1.0})


def test_repair_rule_hand_case_floor_02():
    meta, log = repair_meta_chains(HAND_META, HAND_TIES, HAND_LINKS, min_similarity=0.2)
    assert meta == {"L1C1": 1, "L2C1": 1, "L3C1": 1, "L4C1": 1,
                    "L1C2": 2, "L2C2": 2, "L1C3": 3, "L3C2": 4, "L4C2": 4, "L3C3": 5}
    assert log.to_dict("records") == [{"from_layer": 2, "from_community": 1, "to_layer": 3,
                                       "to_community": 1, "similarity": 0.7}]


def test_repair_rule_hand_case_floor_01():
    meta, log = repair_meta_chains(HAND_META, HAND_TIES, HAND_LINKS, min_similarity=0.1)
    assert meta == {"L1C1": 1, "L2C1": 1, "L3C1": 1, "L4C1": 1,
                    "L1C2": 2, "L2C2": 2, "L3C2": 2, "L4C2": 2, "L1C3": 3, "L3C3": 4}
    assert len(log) == 2


# 3. identical layers keep one id per block ------------------------------------
@pytest.mark.parametrize("fn,alg,tracking", list(itertools.product(FITTERS, ALGS, ["smoothed", "filtered"])))
def test_identical_layers_keep_three_ids_with_repair(fn, alg, tracking):
    base = chain_layers()[0]
    assert len(set(fit_layer_communities(prepare_multilayer_graphs([base, base]), algorithm=alg, seed=123)[0]
                   .membership.values())) == 3
    for T in (4, 12, 20):
        out = FITTERS[fn]([base.copy() for _ in range(T)], algorithm=alg, seed=123,
                          repair_chains=True, tracking=tracking)
        M = out["meta_communities"]
        assert len(ids(M)) == 3, (fn, alg, tracking, T)
        assert all(np.array_equal(M[0], m) for m in M)


# 4. filtered output never depends on later layers ------------------------------
@pytest.mark.parametrize("fn,alg,repair", list(itertools.product(FITTERS, ALGS, [False, True])))
def test_filtered_append_invariance(fn, alg, repair):
    L = chain_layers()
    full = FITTERS[fn](L, algorithm=alg, seed=123, tracking="filtered", repair_chains=repair)["meta_communities"]
    for t in (2, 5, 7):
        trunc = FITTERS[fn](L[:t], algorithm=alg, seed=123, tracking="filtered",
                            repair_chains=repair)["meta_communities"]
        for s in range(t):
            assert np.array_equal(full[s], trunc[s]), (fn, alg, repair, t, s)


def test_filtered_ignores_future_link_weights():
    L = chain_layers()
    base = [{"from": i, "to": i + 1, "weight": 1.0} for i in range(1, 10)]
    heavy = [dict(d) for d in base]; heavy[7]["weight"] = 5.0          # link 8 -> 9
    for repair in (False, True):
        a = fit_multilayer_jaccard(L, layer_links=base, seed=123, tracking="filtered", repair_chains=repair)
        b = fit_multilayer_jaccard(L, layer_links=heavy, seed=123, tracking="filtered", repair_chains=repair)
        for s in range(8):
            assert np.array_equal(a["meta_communities"][s], b["meta_communities"][s]), (repair, s)


# 5. repair only joins groups that never coexist ------------------------------
@pytest.mark.parametrize("fn,alg", list(itertools.product(FITTERS, ALGS)))
def test_repair_never_merges_coexisting_groups(fn, alg):
    for L in (chain_layers(), sbm_layers()):
        pre = FITTERS[fn](L, algorithm=alg, seed=123)
        post = FITTERS[fn](L, algorithm=alg, seed=123, repair_chains=True)
        layers_of = {}
        for t, (a, b) in enumerate(zip(pre["meta_communities"], post["meta_communities"]), start=1):
            for x, y in zip(a, b):
                layers_of.setdefault(int(y), {}).setdefault(int(x), set()).add(t)
        for y, groups in layers_of.items():
            spans = list(groups.values())
            for i, j in itertools.combinations(range(len(spans)), 2):
                assert not (spans[i] & spans[j]), (fn, alg, y)


# 6. similarity floor -----------------------------------------------------------
def test_repair_respects_floor_and_reports_joins():
    L = chain_layers()
    out = fit_multilayer_jaccard(L, seed=123, repair_chains=True)
    log = out["chain_repairs"]
    assert set(log.columns) == {"from_layer", "from_community", "to_layer", "to_community", "similarity"}
    assert len(log) > 0 and (log["similarity"] >= 0.2).all()
    strict = fit_multilayer_jaccard(L, seed=123, repair_chains=True, repair_min_similarity=0.6)["chain_repairs"]
    assert (strict["similarity"] >= 0.6).all() and len(strict) <= len(log)
    off = fit_multilayer_jaccard(L, seed=123, repair_chains=True, repair_min_similarity=1.0)["chain_repairs"]
    assert (off["similarity"] == 1.0).all()


def test_repair_reduces_fragmentation_on_chain_panel():
    L = chain_layers()
    before = len(ids(fit_multilayer_jaccard(L, seed=123)["meta_communities"]))
    after = len(ids(fit_multilayer_jaccard(L, seed=123, repair_chains=True)["meta_communities"]))
    assert after < before


# 7. fit object fields and argument checks -------------------------------------
def test_output_fields_and_argument_checks():
    L = sbm_layers()
    out = fit_multilayer_jaccard(L, seed=123)
    assert out["tracking"] == "smoothed" and out["repair_chains"] is False
    assert out["repair_min_similarity"] == 0.2 and len(out["chain_repairs"]) == 0
    with pytest.raises(ValueError, match="tracking"):
        fit_multilayer_jaccard(L, tracking="causal")
    for bad in (-0.1, 1.5):
        with pytest.raises(ValueError, match="repair_min_similarity"):
            fit_multilayer_jaccard(L, repair_chains=True, repair_min_similarity=bad)


# 8. unequal node sets and gapped layer links -----------------------------------
def test_unequal_nodes_and_gapped_links():
    L = chain_layers()
    L[3] = L[3].copy(); L[3].remove_node(60)
    links = [{"from": 1, "to": 3}, {"from": 2, "to": 4}, {"from": 3, "to": 5}, {"from": 4, "to": 6},
             {"from": 5, "to": 7}, {"from": 6, "to": 8}, {"from": 7, "to": 9}, {"from": 8, "to": 10}]
    for tracking in ("smoothed", "filtered"):
        out = fit_multilayer_jaccard(L, layer_links=links, seed=123, allow_unequal_nodes=True,
                                     repair_chains=True, tracking=tracking)
        assert [len(m) for m in out["meta_communities"]] == [60, 60, 60, 59, 60, 60, 60, 60, 60, 60]
    full = fit_multilayer_jaccard(L, layer_links=links, seed=123, allow_unequal_nodes=True,
                                  repair_chains=True, tracking="filtered")["meta_communities"]
    trunc_links = [d for d in links if d["to"] <= 6]
    trunc = fit_multilayer_jaccard(L[:6], layer_links=trunc_links, seed=123, allow_unequal_nodes=True,
                                   repair_chains=True, tracking="filtered")["meta_communities"]
    for s in range(6):
        assert np.array_equal(full[s], trunc[s]), s


# 9. Louvain resolution moves in the same direction as igraph's gamma (no code change)
def test_louvain_resolution_direction_matches_gamma():
    # python-louvain uses `resolution` as gamma in its move step, like igraph and R.
    # (Below 1 its behaviour is erratic; not tested or changed here.)
    g = prepare_multilayer_graphs(chain_layers()[:3])
    counts = [sum(len(set(f.membership.values())) for f in fit_layer_communities(g, algorithm="louvain",
              resolution_parameter=r, seed=123)) for r in (1.0, 2.0, 4.0)]
    assert counts[0] <= counts[1] <= counts[2] and counts[0] < counts[2], counts


# 10. bootstrap passes the new settings through to the refits ------------------
def test_bootstrap_passes_repair_and_tracking_through():
    L = [nx.to_numpy_array(g, nodelist=range(1, 61)) for g in chain_layers()[:3]]
    for kwargs in ({"repair_chains": True, "repair_min_similarity": 0.3}, {"tracking": "filtered"}):
        res = bootstrap_multilayer(L, n_boot=2, seed=123, **kwargs)
        pe = res.point_estimate
        assert pe["repair_chains"] == kwargs.get("repair_chains", False)
        assert pe["repair_min_similarity"] == kwargs.get("repair_min_similarity", 0.2)
        assert pe["tracking"] == kwargs.get("tracking", "smoothed")
    with pytest.raises(ValueError, match="repair_chains|tracking"):
        bootstrap_multilayer(L, fit_type="identity", n_boot=2, seed=123, repair_chains=True)


# 11. other fitters unchanged ----------------------------------------------------
def test_identity_and_hungarian_match_132():
    from dynamic_multiplex import fit_multilayer_hungarian, fit_multilayer_identity_ties
    rec = json.load(open(FIX / "other_regression_py.json"))
    assert rec["recorded_from"] == "b55c66c"
    L = chain_layers()
    for alg in ALGS:
        got = fit_multilayer_hungarian(L, algorithm=alg, seed=123)["meta_communities"]
        assert [m.tolist() for m in got] == rec["cases"][f"chain|hungarian|{alg}"], alg
    got = fit_multilayer_identity_ties(L, omega=1.0, seed=123)["meta_communities"]
    assert [np.asarray(m).tolist() for m in got] == rec["cases"]["chain|identity|leiden"]


# 12. Codex review fixes ------------------------------------------------------------
def test_filtered_requires_a_seed():
    # Unseeded stage 1 on future layers advances the RNG that filtered stage 2
    # draws from, so append-invariance needs a seed.
    from dynamic_multiplex import bootstrap_multilayer, fit_multilayer_from_partitions
    L = chain_layers()
    with pytest.raises(ValueError, match="seed"):
        fit_multilayer_jaccard(L, seed=None, tracking="filtered")
    parts = [dict(f.membership) for f in fit_multilayer_jaccard(L, seed=123)["layer_communities"]]
    with pytest.raises(ValueError, match="seed"):
        fit_multilayer_from_partitions(parts, seed=None, tracking="filtered")
    with pytest.raises(ValueError, match="seed"):
        bootstrap_multilayer(L, n_boot=2, seed=None, tracking="filtered")
    fit_multilayer_jaccard(L, seed=None)                     # smoothed unseeded still allowed


def test_bootstrap_forwards_and_validates_repair_args():
    from dynamic_multiplex import bootstrap_multilayer
    L = sbm_layers()
    res = bootstrap_multilayer(L, n_boot=2, seed=123, repair_min_similarity=0.8)
    assert res.point_estimate["repair_min_similarity"] == 0.8
    assert res.point_estimate["repair_chains"] is False
    with pytest.raises(ValueError, match="repair_min_similarity"):
        bootstrap_multilayer(L, n_boot=2, seed=123, repair_min_similarity=1.5)
    with pytest.raises(ValueError, match="repair_chains"):
        bootstrap_multilayer(L, n_boot=2, seed=123, repair_chains=1)
    with pytest.raises(ValueError, match="tracking"):
        bootstrap_multilayer(L, n_boot=2, seed=123, tracking="smooth")
