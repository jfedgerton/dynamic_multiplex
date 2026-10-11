"""1.4.0: nodes may enter and leave, with or without a change in layer size, in every
function. Node identity is the label, never the position. Seed 123."""
import importlib.util

import networkx as nx
import numpy as np
import pytest

from dynamic_multiplex import (
    bootstrap_multilayer,
    fit_multilayer_hungarian,
    fit_multilayer_identity_ties,
    fit_multilayer_jaccard,
    fit_multilayer_overlap,
    fit_multilayer_weighted_jaccard,
    fit_multilayer_weighted_overlap,
    partition_stability,
)
from dynamic_multiplex.multilayer_utils import prepare_multilayer_graphs

pytestmark = pytest.mark.skipif(
    importlib.util.find_spec("leidenalg") is None or importlib.util.find_spec("community") is None,
    reason="needs leidenalg and python-louvain",
)
FITTERS = [fit_multilayer_jaccard, fit_multilayer_overlap, fit_multilayer_weighted_jaccard,
           fit_multilayer_weighted_overlap, fit_multilayer_identity_ties, fit_multilayer_hungarian]


def blocks(nodes, seed):
    rng = np.random.default_rng(seed)
    nodes = list(nodes)
    g = nx.Graph(); g.add_nodes_from(nodes)
    for i in range(len(nodes)):
        for j in range(i + 1, len(nodes)):
            p = 0.8 if (nodes[i] - 1) // 5 == (nodes[j] - 1) // 5 else 0.03
            if rng.random() < p:
                g.add_edge(nodes[i], nodes[j], weight=1.0)
    return g


SAME_SIZE = [blocks(range(1, 21), 1), blocks(range(6, 26), 2), blocks(range(11, 31), 3)]   # 20 nodes each, members change
DIFF_SIZE = [blocks(range(1, 21), 1), blocks(range(1, 26), 2), blocks(range(6, 24), 3)]


@pytest.mark.parametrize("layers", [SAME_SIZE, DIFF_SIZE], ids=["same-size", "diff-size"])
def test_all_fitters_require_and_accept_labels(layers):
    for f in FITTERS:
        with pytest.raises(ValueError, match="same node set"):
            f(layers, seed=123)
        out = f(layers, seed=123, allow_unequal_nodes=True)
        assert [len(m) for m in out["meta_communities"]] == [g.number_of_nodes() for g in layers]


def test_matrices_rejected_with_allow_unequal_nodes_even_at_equal_size():
    mats = [nx.to_numpy_array(g) for g in SAME_SIZE]
    with pytest.raises(ValueError, match="node labels"):
        prepare_multilayer_graphs(mats, require_same_nodes=False)
    with pytest.raises(ValueError, match="node labels"):
        fit_multilayer_jaccard(mats, allow_unequal_nodes=True, seed=123)


def test_identity_ties_link_a_node_only_to_itself():
    out = fit_multilayer_identity_ties(SAME_SIZE, seed=123, allow_unequal_nodes=True)
    ties = out["interlayer_ties"]
    expected = {(1, 2, v) for v in range(6, 21)} | {(2, 3, v) for v in range(11, 26)}
    assert {(int(a), int(b), int(n)) for a, b, n in zip(ties["from_layer"], ties["to_layer"], ties["node"])} == expected


@pytest.mark.parametrize("layers", [SAME_SIZE, DIFF_SIZE], ids=["same-size", "diff-size"])
@pytest.mark.parametrize("fit_type", ["jaccard", "identity"])
def test_bootstrap_with_entering_and_leaving_nodes(layers, fit_type):
    with pytest.raises(ValueError, match="different node sets"):
        bootstrap_multilayer(layers, fit_type=fit_type, n_boot=2, seed=123)
    res = bootstrap_multilayer(layers, fit_type=fit_type, n_boot=4, seed=123, allow_unequal_nodes=True)
    assert res.n_boot == 4
    assert [len(v) for v in res.node_stability] == [g.number_of_nodes() for g in layers]
    p = partition_stability(res)
    assert 0 <= p["stability"] <= 1 and len(p["stability_by_layer"]) == 3


def test_bootstrap_equal_nodes_unchanged_by_flag():
    mats = [nx.to_numpy_array(g, nodelist=range(1, 21)) for g in (blocks(range(1, 21), 1), blocks(range(1, 21), 2))]
    graphs = [nx.relabel_nodes(nx.from_numpy_array(m), {i: i + 1 for i in range(20)}) for m in mats]
    a = partition_stability(bootstrap_multilayer(mats, n_boot=4, seed=123))
    b = partition_stability(bootstrap_multilayer(graphs, n_boot=4, seed=123, allow_unequal_nodes=True))
    assert a["stability"] == b["stability"]
