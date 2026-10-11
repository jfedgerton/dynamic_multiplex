"""1.4.0: layers_from_membership(), the two-mode (membership) input constructor.

Written before the implementation. Weightings exactly as in section 3 of the
project brief (claude/twomode_input_brief.md). Seed 123.

hypernetx 2.4.3 requires pandas < 3, so the hypernetx reference tests run only
where it is installed (a separate venv here); they are skipped otherwise.
"""
import importlib.util
import itertools

import networkx as nx
import numpy as np
import pandas as pd
import pytest
from _helpers_140 import nmi

from dynamic_multiplex import fit_multilayer_jaccard, layers_from_membership
from dynamic_multiplex.layers_from_membership import _kumar_targets
from dynamic_multiplex.multilayer_utils import fit_layer_communities

pytestmark = pytest.mark.skipif(
    importlib.util.find_spec("leidenalg") is None or importlib.util.find_spec("community") is None,
    reason="needs leidenalg and python-louvain",
)
needs_hnx = pytest.mark.skipif(importlib.util.find_spec("hypernetx") is None, reason="hypernetx not installed")
ALGS = ["leiden", "louvain"]


# fixtures ------------------------------------------------------------------------
def hand_wide():
    # g1 = {a, b, c}, g2 = {c, d}, g3 = {b, c, d, e}
    return pd.DataFrame({"g1": [1, 1, 1, 0, 0], "g2": [0, 0, 1, 1, 0], "g3": [0, 1, 1, 1, 1]},
                        index=list("abcde"))


def planted_long(seed, n_blocs=4, bloc=10, n_home=30, n_noise=0, prefix="g"):
    """Every node in at least one home group; home groups of 3-6 from one bloc."""
    rng = np.random.default_rng(seed)
    nodes = np.arange(n_blocs * bloc)
    rows = []
    for k in range(n_home):
        b = k % n_blocs
        members = rng.choice(nodes[b * bloc:(b + 1) * bloc], size=rng.integers(3, 7), replace=False)
        rows += [(int(v), f"{prefix}{k}") for v in members]
    covered = {v for v, _ in rows}
    for v in nodes:                                   # make sure nobody is left out
        if v not in covered:
            b = v // bloc
            mate = int(rng.choice([u for u in nodes[b * bloc:(b + 1) * bloc] if u != v]))
            rows += [(int(v), f"{prefix}fix{v}"), (mate, f"{prefix}fix{v}")]
    for k in range(n_noise):
        members = rng.choice(nodes, size=rng.integers(3, 7), replace=False)
        rows += [(int(v), f"{prefix}n{k}") for v in members]
    return pd.DataFrame(rows, columns=["node", "group"]).drop_duplicates()


def truth(n_blocs=4, bloc=10):
    return np.repeat(np.arange(n_blocs), bloc)


def edges(g):
    return sorted((min(u, v), max(u, v), d["weight"]) for u, v, d in g.edges(data=True))


# 1. hand-computed weights ---------------------------------------------------------
def test_hand_count():
    g = layers_from_membership([hand_wide()], weighting="count")[0]
    assert edges(g) == [("a", "b", 1.0), ("a", "c", 1.0), ("b", "c", 2.0), ("b", "d", 1.0), ("b", "e", 1.0),
                        ("c", "d", 2.0), ("c", "e", 1.0), ("d", "e", 1.0)]


def test_hand_size():
    g = layers_from_membership([hand_wide()], weighting="size")[0]
    exp = {("a", "b"): .5, ("a", "c"): .5, ("b", "c"): .5 + 1 / 3, ("b", "d"): 1 / 3, ("b", "e"): 1 / 3,
           ("c", "d"): 1 + 1 / 3, ("c", "e"): 1 / 3, ("d", "e"): 1 / 3}
    got = {(u, v): w for u, v, w in edges(g)}
    assert set(got) == set(exp)
    for k in exp:
        assert got[k] == pytest.approx(exp[k], abs=1e-12), k


def test_hand_kumar_targets():
    B = hand_wide().to_numpy()
    labels = np.array([0, 0, 0, 1, 1])               # {a, b, c} and {d, e}
    r = _kumar_targets(B, labels)
    assert r == pytest.approx([1.25 * 5 / 3, 4 / 3, 4 / 3], abs=1e-12)


# 2. count equals the dyadic number of shared groups --------------------------------
def test_count_equals_dyadic_shared_groups():
    long = planted_long(123, n_noise=10)
    pairs = long.merge(long, on="group")
    pairs = pairs[pairs["node_x"] < pairs["node_y"]].groupby(["node_x", "node_y"]).size()
    g = layers_from_membership([long], weighting="count")[0]
    got = {(min(u, v), max(u, v)): w for u, v, w in g.edges(data="weight")}
    assert got == {k: float(v) for k, v in pairs.items()}


# 3. Kumar against the reference implementation -------------------------------------
@needs_hnx
def test_kumar_targets_match_hypernetx_formula():
    import hypernetx as hnx
    long = planted_long(123, n_noise=8)
    H = hnx.Hypergraph(long.groupby("group")["node"].apply(list).to_dict())
    wide = pd.crosstab(long["node"], long["group"]).clip(upper=1)
    labels = np.arange(len(wide)) % 3                 # an arbitrary fixed partition
    CH = [set(wide.index[labels == c]) for c in range(3)]
    ref = [sum(1 / (1 + H.size(e, c)) for c in CH) * (H.size(e) + len(CH)) / len(H.edges)
           for e in wide.columns]
    assert _kumar_targets(wide.to_numpy(), labels) == pytest.approx(ref, abs=1e-12)


@needs_hnx
@pytest.mark.parametrize("alg", ALGS)
def test_kumar_partition_matches_hypernetx_on_easy_slice(alg):
    import hypernetx as hnx
    from hypernetx.algorithms.clustering import hypergraph_modularity as hmod
    long = planted_long(123)
    H = hnx.Hypergraph(long.groupby("group")["node"].apply(list).to_dict())
    ref = {v: k for k, part in enumerate(hmod.kumar(H)) for v in part}
    g = layers_from_membership([long], weighting="kumar", algorithm=alg)[0]
    part = g.graph["dm_membership"]["partition"]
    nodes = sorted(ref)
    assert nmi([ref[v] for v in nodes], [part[v] for v in nodes]) == pytest.approx(1.0)


# 4. stage 1 on the returned weights reproduces Kumar's own final partition ----------
@pytest.mark.parametrize("alg", ALGS)
def test_stage1_on_kumar_weights_reproduces_kumar_partition(alg):
    g = layers_from_membership([planted_long(123)], weighting="kumar", algorithm=alg)[0]
    part = g.graph["dm_membership"]["partition"]
    gs = nx.relabel_nodes(g, str)                    # string labels: membership keys are the labels
    fit = fit_layer_communities([gs], algorithm=alg, seed=123)[0].membership
    nodes = sorted(part)
    assert nmi([part[v] for v in nodes], [fit[str(v)] for v in nodes]) == pytest.approx(1.0)
    assert nmi(truth(), [part[v] for v in range(40)]) == pytest.approx(1.0)


# 5. slice t depends only on slice t ----------------------------------------------
@pytest.mark.parametrize("weighting,alg", [("count", None), ("size", None), ("kumar", "leiden"), ("kumar", "louvain")])
def test_slice_independence(weighting, alg):
    S = [planted_long(123 + t, n_noise=15, prefix=f"t{t}_") for t in range(4)]
    kw = {"weighting": weighting, "algorithm": alg}
    full = layers_from_membership(S, **kw)
    for t in range(4):
        alone = layers_from_membership([S[t]], **kw)[0]
        prefix = layers_from_membership(S[:t + 1], **kw)[t]
        for g in (alone, prefix):
            assert edges(g) == edges(full[t]), (weighting, alg, t)
            assert g.graph["dm_membership"]["group_weights"].equals(full[t].graph["dm_membership"]["group_weights"])


# 7. unequal node sets, group birth and death, accepted by the fitters ------------
def test_unequal_nodes_and_group_turnover_fit():
    S = [planted_long(123), planted_long(124), planted_long(125)]
    S[1] = S[1].assign(node=S[1]["node"].replace({39: 99}))   # node 39 leaves, node 99 enters
    S[2] = pd.concat([S[2], pd.DataFrame({"node": [40, 41, 42], "group": ["new"] * 3})])  # three new nodes
    G = layers_from_membership(S, weighting="size")
    assert [g.number_of_nodes() for g in G] == [40, 40, 43]
    assert 39 not in G[1] and 99 in G[1]
    out = fit_multilayer_jaccard(G, allow_unequal_nodes=True, seed=123)
    assert [len(m) for m in out["meta_communities"]] == [g.number_of_nodes() for g in G]


# 8. validation ---------------------------------------------------------------------
def test_validation_errors():
    W = hand_wide()
    bad = W.copy(); bad.iloc[0, 0] = 2
    with pytest.raises(ValueError, match="0/1"):
        layers_from_membership([bad], weighting="count")
    bad = W.astype(float); bad.iloc[0, 0] = np.nan
    with pytest.raises(ValueError, match="missing"):
        layers_from_membership([bad], weighting="count")
    bad = W.copy(); bad.index = list("abcda")
    with pytest.raises(ValueError, match="duplicate"):
        layers_from_membership([bad], weighting="count")
    bad = W.copy(); bad.columns = ["g1", "g1", "g3"]
    with pytest.raises(ValueError, match="duplicate"):
        layers_from_membership([bad], weighting="count")
    long = planted_long(123)
    with pytest.raises(ValueError, match="duplicate"):
        layers_from_membership([pd.concat([long, long.iloc[:1]])], weighting="count")
    with pytest.raises(ValueError, match="node.*group|columns"):
        layers_from_membership([long.rename(columns={"group": "org"})[["node", "org"]].assign(x=1)], weighting="count")
    with pytest.raises(ValueError, match="row labels|DataFrame"):
        layers_from_membership([W.to_numpy()], weighting="count")
    with pytest.raises(ValueError, match="list"):
        layers_from_membership(W, weighting="count")
    with pytest.raises(ValueError, match="at least one"):
        layers_from_membership([], weighting="count")
    with pytest.raises(ValueError, match="empty"):
        layers_from_membership([W.iloc[:0]], weighting="count")
    with pytest.raises(TypeError):
        layers_from_membership([W])                    # no default weighting
    with pytest.raises(ValueError, match="weighting"):
        layers_from_membership([W], weighting="hyper")
    with pytest.raises(ValueError, match="algorithm"):
        layers_from_membership([W], weighting="kumar")
    with pytest.raises(ValueError, match="algorithm"):
        layers_from_membership([W], weighting="count", algorithm="leiden")
    for kw, pat in [({"resolution_parameter": 0}, "resolution_parameter"), ({"max_iter": 0}, "max_iter"),
                    ({"tol": 0}, "tol"), ({"min_group_size": 1}, "min_group_size"), ({"rescale": "max"}, "rescale")]:
        with pytest.raises(ValueError, match=pat):
            layers_from_membership([W], weighting="kumar", algorithm="leiden", **kw)


def test_dropped_groups_and_orphaned_nodes():
    W = hand_wide(); W["solo"] = [1, 0, 0, 0, 0]
    with pytest.warns(UserWarning, match="dropped 1 group"):
        g = layers_from_membership([W], weighting="count")[0]
    assert g.graph["dm_membership"]["groups_dropped"] == 1
    orphan = hand_wide(); orphan.loc["f"] = 0; orphan["solo"] = [0, 0, 0, 0, 0, 1]
    with pytest.warns(UserWarning, match="dropped 1 group"):
        with pytest.raises(ValueError, match=r"slice 1.*'f'"):
            layers_from_membership([orphan], weighting="size")


def test_kumar_nonconvergence_fails_loudly():
    long = planted_long(123, n_noise=40)
    with pytest.raises(RuntimeError, match="did not converge"):
        layers_from_membership([long], weighting="kumar", algorithm="leiden", max_iter=1, tol=1e-15)


# 9. attributes, seed, rescale -------------------------------------------------------
@pytest.mark.parametrize("weighting,alg", [("count", None), ("size", None), ("kumar", "leiden"), ("kumar", "louvain")])
def test_attributes_and_seed(weighting, alg):
    long = planted_long(123, n_noise=15)
    a = layers_from_membership([long], weighting=weighting, algorithm=alg)[0]
    b = layers_from_membership([long], weighting=weighting, algorithm=alg)[0]
    assert edges(a) == edges(b)
    att = a.graph["dm_membership"]
    assert set(att) == {"weighting", "group_weights", "iterations", "final_diff", "groups_dropped",
                        "rescale", "partition"}
    assert att["weighting"] == weighting and att["rescale"] == "none"
    assert len(att["group_weights"]) == long["group"].nunique()
    if weighting == "kumar":
        assert 1 <= att["iterations"] <= 50 and att["final_diff"] <= 0.01 and att["partition"] is not None
    else:
        assert att["iterations"] == 0 and att["partition"] is None
        assert (att["group_weights"] == 1.0).all()


@pytest.mark.parametrize("weighting,alg", [("count", None), ("size", None), ("kumar", "leiden")])
def test_rescale_mean(weighting, alg):
    long = planted_long(123, n_noise=15)
    raw = layers_from_membership([long], weighting=weighting, algorithm=alg)[0]
    sc = layers_from_membership([long], weighting=weighting, algorithm=alg, rescale="mean")[0]
    w_raw = np.array([w for *_, w in edges(raw)]); w_sc = np.array([w for *_, w in edges(sc)])
    assert w_sc.mean() == pytest.approx(1.0, abs=1e-12)
    assert w_sc == pytest.approx(w_raw / w_raw.mean(), abs=1e-12)
    assert sc.graph["dm_membership"]["rescale"] == "mean"
