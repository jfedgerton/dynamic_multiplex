"""Input policy, node alignment and R-parity checks (package version 1.3.2)."""
import networkx as nx
import numpy as np
import pytest

from dynamic_multiplex import (
    fit_multilayer_hungarian,
    fit_multilayer_identity_ties,
    fit_multilayer_jaccard,
    fit_multilayer_overlap,
    fit_multilayer_weighted_jaccard,
    fit_multilayer_weighted_overlap,
)
from dynamic_multiplex.multilayer_utils import (
    _collapse_directed,
    layer_node_strengths,
    prepare_multilayer_graphs,
)

COMMUNITY_FITTERS = [
    fit_multilayer_jaccard,
    fit_multilayer_overlap,
    fit_multilayer_weighted_jaccard,
    fit_multilayer_weighted_overlap,
]


def test_self_loops_dropped_from_matrix_and_graph_input():
    L = np.array([[7, 1], [1, 0]])
    g = nx.from_numpy_array(L)
    for layers in ([L, L], [g, g]):
        prepared = prepare_multilayer_graphs(layers)
        assert nx.number_of_selfloops(prepared[0]) == 0
        assert list(prepared[0].edges(data="weight")) == [(0, 1, 1)]
    # strengths no longer count the loop
    assert layer_node_strengths(prepare_multilayer_graphs([g, g]))[0] == {1: 1.0, 2: 1.0}
    # the caller's graph is not modified
    assert nx.number_of_selfloops(g) == 1


def test_zero_weight_edges_dropped():
    g = nx.Graph()
    g.add_edge(0, 1, weight=0.0)
    g.add_edge(1, 2, weight=1.0)
    prepared = prepare_multilayer_graphs([g, g])
    assert list(prepared[0].edges()) == [(1, 2)]


def test_zero_vertex_layer_rejected():
    with pytest.raises(ValueError, match="no vertices"):
        fit_multilayer_jaccard([np.zeros((0, 0)), np.zeros((0, 0))])


def test_directed_louvain_collapse_sums_reciprocal_weights():
    d = nx.DiGraph()
    d.add_edge(0, 1, weight=0.2)
    d.add_edge(1, 0, weight=0.8)
    u = _collapse_directed(d)
    assert u[0][1]["weight"] == pytest.approx(1.0)


def test_identity_ties_unequal_integer_labels_consistent():
    # node 1 is the only shared actor; it must keep one id in both layers
    out = fit_multilayer_identity_ties(
        [nx.empty_graph([0, 1]), nx.empty_graph([1, 2])], allow_unequal_nodes=True
    )
    assert out["interlayer_ties"]["node"].tolist() == [1]
    assert set(out["layer_communities"][0].membership) == {0, 1}
    assert set(out["layer_communities"][1].membership) == {1, 2}
    m1, m2 = out["meta_communities"]
    assert m1[1] == m2[0]  # node 1: position 1 in layer 1, position 0 in layer 2


def test_identity_ties_unequal_negative_labels():
    out = fit_multilayer_identity_ties(
        [nx.empty_graph([0, 1]), nx.empty_graph([-1, 0])], allow_unequal_nodes=True
    )
    assert out["interlayer_ties"]["node"].tolist() == [0]
    m1, m2 = out["meta_communities"]
    assert m1[0] == m2[1]  # node 0: position 0 in layer 1, position 1 in layer 2


@pytest.mark.parametrize("fitter", COMMUNITY_FITTERS)
def test_community_fitters_accept_unequal_labelled_graphs(fitter):
    g1 = nx.Graph([("a", "b"), ("b", "c"), ("a", "c")])
    g2 = nx.Graph([("b", "c"), ("c", "d"), ("b", "d")])
    out = fitter([g1, g2], allow_unequal_nodes=True)
    assert [len(m) for m in out["meta_communities"]] == [3, 3]
    with pytest.raises(ValueError, match="same node set"):
        fitter([g1, g2])


@pytest.mark.parametrize("fitter", COMMUNITY_FITTERS + [fit_multilayer_identity_ties, fit_multilayer_hungarian])
def test_unequal_matrices_rejected(fitter):
    with pytest.raises(ValueError, match="node labels"):
        fitter([np.zeros((2, 2)), np.zeros((3, 3))], allow_unequal_nodes=True)


def test_hungarian_matches_by_label_with_unequal_nodes():
    g1 = nx.disjoint_union(nx.complete_graph(4), nx.complete_graph(4))
    g1 = nx.relabel_nodes(g1, dict(enumerate("abcdefgh")))
    # layer 2 drops "a", adds "z" to the second clique
    g2 = g1.copy()
    g2.remove_node("a")
    g2.add_edges_from([("z", n) for n in "efgh"])
    out = fit_multilayer_hungarian([g1, g2], allow_unequal_nodes=True, seed=123)
    lab1 = dict(zip(sorted(g1.nodes()), out["meta_communities"][0]))
    lab2 = dict(zip(sorted(g2.nodes()), out["meta_communities"][1]))
    for n in "bcdefgh":
        assert lab1[n] == lab2[n]
    assert lab2["z"] == lab2["e"]


def test_default_algorithm_is_leiden():
    A = np.array([[0, 1], [1, 0]])
    assert fit_multilayer_jaccard([A, A])["algorithm"] == "leiden"
