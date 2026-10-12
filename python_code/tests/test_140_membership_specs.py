"""1.4.0 (Python only): bipartite and all-or-nothing membership specifications, and
tracking of per-slice partitions over time. The R package deliberately has no
counterpart (see NEWS). Written before the implementation. Seed 123.
"""
import importlib.util
import itertools

import networkx as nx
import numpy as np
import pandas as pd
import pytest
from _helpers_140 import nmi

from _helpers_140 import chain_layers
from test_140_twomode import hand_wide, planted_long, truth
from dynamic_multiplex import (
    fit_multilayer_from_partitions,
    fit_multilayer_hungarian,
    fit_multilayer_jaccard,
    fit_multilayer_overlap,
    partition_membership_bipartite,
    partition_membership_strict,
)
from dynamic_multiplex.membership_partitions import _barber_modularity

pytestmark = pytest.mark.skipif(
    importlib.util.find_spec("leidenalg") is None or importlib.util.find_spec("community") is None,
    reason="needs leidenalg and python-louvain",
)
needs_hnx = pytest.mark.skipif(importlib.util.find_spec("hypernetx") is None, reason="hypernetx not installed")


def labels_of(part, nodes):
    return [part.nodes[v] for v in nodes]


# Barber bipartite modularity ------------------------------------------------------
def test_barber_modularity_hand_case():
    B = hand_wide().to_numpy(dtype=float)          # 5 nodes x 3 groups, m = 9 memberships
    node_c = np.array([0, 0, 0, 1, 1]); group_c = np.array([0, 1, 1])
    k = B.sum(1); d = B.sum(0); m = B.sum()
    exp = sum((B[i, j] - k[i] * d[j] / m) for i in range(5) for j in range(3) if node_c[i] == group_c[j]) / m
    assert _barber_modularity(B, node_c, group_c) == pytest.approx(exp, abs=1e-12)


def test_bipartite_recovers_planted_blocs_and_reports_barber_q():
    long = planted_long(123)
    part = partition_membership_bipartite([long])[0]
    assert set(part.nodes) == set(range(40))
    assert nmi(truth(), labels_of(part, range(40))) == pytest.approx(1.0)
    assert set(part.groups) == set(long["group"])
    # the reported quality is Barber's modularity of the returned partition
    wide = pd.crosstab(long["node"], long["group"])
    B = wide.to_numpy(dtype=float)
    q = _barber_modularity(B, np.array([part.nodes[v] for v in wide.index]),
                           np.array([part.groups[g] for g in wide.columns]))
    assert part.quality == pytest.approx(q, abs=1e-10)


def test_bipartite_beats_the_planted_assignment():
    long = planted_long(123, n_noise=15)
    part = partition_membership_bipartite([long])[0]
    wide = pd.crosstab(long["node"], long["group"])
    B = wide.to_numpy(dtype=float)
    node_c = truth()[wide.index.to_numpy()]
    group_c = np.array([np.bincount(node_c[B[:, j] > 0]).argmax() for j in range(B.shape[1])])
    assert part.quality >= _barber_modularity(B, node_c, group_c) - 1e-12


def test_bipartite_is_leiden_only_and_validates_input():
    with pytest.raises(ValueError, match="Leiden"):
        partition_membership_bipartite([planted_long(123)], algorithm="louvain")
    bad = hand_wide(); bad.iloc[0, 0] = 2
    with pytest.raises(ValueError, match="0/1"):
        partition_membership_bipartite([bad])
    with pytest.raises(ValueError, match="resolution_parameter"):
        partition_membership_bipartite([hand_wide()], resolution_parameter=0)


# all-or-nothing (strict) hypergraph modularity ----------------------------------------
@needs_hnx
def test_strict_recovers_planted_blocs_and_improves_on_its_start():
    part = partition_membership_strict([planted_long(123)])[0]
    assert nmi(truth(), labels_of(part, range(40))) == pytest.approx(1.0)
    assert part.quality >= part.initial_quality - 1e-12


@needs_hnx
def test_strict_matches_the_reference_pipeline():
    import random
    import hypernetx as hnx
    from hypernetx.algorithms.clustering import hypergraph_modularity as hmod
    long = planted_long(123, n_noise=15)
    part = partition_membership_strict([long])[0]
    random.seed(123); np.random.seed(123)
    H = hnx.Hypergraph(long.groupby("group")["node"].apply(list).to_dict())
    A = hmod.last_step(H, hmod.kumar(H), wdc=hmod.strict)
    ref = {v: k for k, block in enumerate(A) for v in block}
    nodes = sorted(ref)
    assert nmi([ref[v] for v in nodes], labels_of(part, nodes)) == pytest.approx(1.0)
    assert part.quality == pytest.approx(hmod.modularity(H, A, wdc=hmod.strict), abs=1e-10)


def test_strict_without_hypernetx_fails_loudly(monkeypatch):
    import builtins
    real = builtins.__import__

    def fake(name, *a, **k):
        if name.startswith("hypernetx"):
            raise ImportError("no hypernetx")
        return real(name, *a, **k)
    monkeypatch.setattr(builtins, "__import__", fake)
    with pytest.raises(ImportError, match="hypergraph"):
        partition_membership_strict([planted_long(123)])


# shared properties --------------------------------------------------------------------
SPECS = [(partition_membership_bipartite, {})] + (
    [(partition_membership_strict, {})] if importlib.util.find_spec("hypernetx") else [])


@pytest.mark.parametrize("fn,kw", SPECS)
def test_slice_independence_and_seed(fn, kw):
    S = [planted_long(123 + t, n_noise=10, prefix=f"t{t}_") for t in range(3)]
    full = fn(S, **kw)
    for t in range(3):
        alone = fn([S[t]], **kw)[0]
        assert alone.nodes == full[t].nodes and alone.quality == full[t].quality
    assert fn(S, **kw)[1].nodes == full[1].nodes


# tracking supplied partitions ------------------------------------------------------------
@pytest.mark.parametrize("alg", ["leiden", "louvain"])
@pytest.mark.parametrize("kw", [{}, {"repair_chains": True}, {"tracking": "filtered", "repair_chains": True}])
def test_tracking_stage1_partitions_reproduces_the_fitter(alg, kw):
    L = chain_layers()
    fit = fit_multilayer_jaccard(L, algorithm=alg, seed=123, **kw)
    parts = [dict(f.membership) for f in fit["layer_communities"]]
    out = fit_multilayer_from_partitions(parts, method="jaccard", algorithm=alg, seed=123, **kw)
    assert [m.tolist() for m in out["meta_communities"]] == [m.tolist() for m in fit["meta_communities"]]
    ov = fit_multilayer_overlap(L, algorithm=alg, seed=123)
    out = fit_multilayer_from_partitions([dict(f.membership) for f in ov["layer_communities"]],
                                         method="overlap", algorithm=alg, seed=123)
    assert [m.tolist() for m in out["meta_communities"]] == [m.tolist() for m in ov["meta_communities"]]


def test_tracking_hungarian_reproduces_the_fitter():
    L = chain_layers()
    fit = fit_multilayer_hungarian(L, seed=123)
    out = fit_multilayer_from_partitions([dict(f.membership) for f in fit["layer_communities"]], method="hungarian")
    assert [m.tolist() for m in out["meta_communities"]] == [m.tolist() for m in fit["meta_communities"]]


def test_tracking_membership_partitions_with_turnover():
    S = [planted_long(123, prefix="a"), planted_long(124, prefix="b"), planted_long(125, prefix="c")]
    S[1] = S[1].assign(node=S[1]["node"].replace({39: 99}))
    parts = partition_membership_bipartite(S)
    out = fit_multilayer_from_partitions(parts, method="jaccard", repair_chains=True)
    assert [len(m) for m in out["meta_communities"]] == [40, 40, 40]
    assert out["nodes"][1][-1] == 99 and 39 not in out["nodes"][1]
    joint_truth = np.concatenate([truth(), truth(), truth()])
    assert nmi(joint_truth, np.concatenate(out["meta_communities"])) == pytest.approx(1.0)


def test_tracking_validates_input():
    with pytest.raises(ValueError, match="method"):
        fit_multilayer_from_partitions([{1: 1}, {1: 1}], method="weighted_jaccard")
    with pytest.raises(ValueError, match="at least two"):
        fit_multilayer_from_partitions([{1: 1}])
    with pytest.raises(ValueError, match="partition"):
        fit_multilayer_from_partitions([{1: 1}, [1, 2]])


@needs_hnx
def test_strict_restores_rng_even_when_the_import_draws(monkeypatch):
    # Some environments draw from Python's `random` while importing hypernetx's
    # dependencies (Codex review). Simulate that with an import hook.
    import builtins
    import random
    real = builtins.__import__

    def drawing_import(name, *a, **k):
        if name.startswith("hypernetx"):
            random.random(); np.random.random()
        return real(name, *a, **k)
    monkeypatch.setattr(builtins, "__import__", drawing_import)
    random.seed(987); np.random.seed(987)
    py, npst = random.getstate(), np.random.get_state()
    partition_membership_strict([hand_wide()])
    after = np.random.get_state()
    assert random.getstate() == py
    assert npst[0] == after[0] and np.array_equal(npst[1], after[1]) and npst[2:] == after[2:]
