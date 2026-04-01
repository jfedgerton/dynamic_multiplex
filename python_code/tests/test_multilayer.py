import numpy as np

from dynamic_multiplex import (
    fit_multilayer_identity_ties,
    fit_multilayer_jaccard,
    fit_multilayer_overlap,
    fit_multilayer_weighted_jaccard,
    fit_multilayer_weighted_overlap,
    simulate_and_fit_multilayer,
)


def test_simulate_and_fit_louvain_jaccard_smoke():
    out = simulate_and_fit_multilayer(
        n_nodes=24,
        n_layers=3,
        n_communities=3,
        fit_type="jaccard",
        algorithm="louvain",
        seed=42,
    )
    assert len(out["layers"]) == 3
    assert not out["fit"]["interlayer_ties"].empty


def test_identity_ties_row_count_matches_nodes_times_links():
    layers = [np.eye(6), np.eye(6), np.eye(6)]
    layer_links = [{"from": 1, "to": 2, "weight": 1.0}, {"from": 2, "to": 3, "weight": 0.7}]
    out = fit_multilayer_identity_ties(layers, algorithm="louvain", layer_links=layer_links)
    assert out["interlayer_ties"].shape[0] == 12


def test_identity_ties_variable_node_sets():
    import networkx as nx

    g1 = nx.Graph()
    g1.add_nodes_from(["A", "B", "C"])
    g1.add_edges_from([("A", "B"), ("A", "C"), ("B", "C")])

    g2 = nx.Graph()
    g2.add_nodes_from(["B", "C", "D"])
    g2.add_edges_from([("B", "C"), ("B", "D"), ("C", "D")])

    g3 = nx.Graph()
    g3.add_nodes_from(["A", "C", "D"])
    g3.add_edges_from([("A", "C"), ("A", "D"), ("C", "D")])

    out = fit_multilayer_identity_ties([g1, g2, g3], algorithm="louvain")
    ties = out["interlayer_ties"]

    # Layer 1->2: shared = B, C (2 ties)
    # Layer 2->3: shared = C, D (2 ties)
    # Total: 4 ties
    assert ties.shape[0] == 4

    ties_12 = ties[(ties["from_layer"] == 1) & (ties["to_layer"] == 2)]
    assert set(ties_12["node"]) == {"B", "C"}

    ties_23 = ties[(ties["from_layer"] == 2) & (ties["to_layer"] == 3)]
    assert set(ties_23["node"]) == {"C", "D"}


def test_custom_layer_links_are_respected():
    layers = [np.eye(5), np.eye(5), np.eye(5), np.eye(5)]
    layer_links = [{"from": 1, "to": 3, "weight": 0.5}]
    out = fit_multilayer_jaccard(layers, algorithm="louvain", layer_links=layer_links)
    links = out["layer_links"]
    assert links.shape[0] == 1
    assert int(links.iloc[0]["from"]) == 1
    assert int(links.iloc[0]["to"]) == 3


def test_jaccard_self_loops_undirected():
    layers = [np.eye(6), np.eye(6), np.eye(6)]
    out = fit_multilayer_jaccard(layers, algorithm="louvain")
    ties = out["interlayer_ties"]
    loops = ties[(ties["from_layer"] == ties["to_layer"]) & (ties["from_community"] == ties["to_community"])]

    assert not loops.empty
    # Undirected self-loops use similarity=2 (Louvain aggregation convention)
    assert (loops["similarity"] == 2.0).all()
    assert (loops["weighted_similarity"] == 2.0 * loops["layer_weight"]).all()


def test_overlap_self_loop_multiplier_scales_weights():
    layers = [np.eye(6), np.eye(6), np.eye(6)]
    out = fit_multilayer_overlap(layers, algorithm="louvain", self_loop_multiplier=2.0)
    ties = out["interlayer_ties"]
    loops = ties[(ties["from_layer"] == ties["to_layer"]) & (ties["from_community"] == ties["to_community"])]

    assert not loops.empty
    # weighted_similarity = self_sim(2) * layer_weight * multiplier(2)
    assert (loops["weighted_similarity"] == 4.0 * loops["layer_weight"]).all()


def test_weighted_jaccard_uses_node_strength_weights():
    layer1 = np.array([
        [0, 10, 1],
        [10, 0, 0],
        [1, 0, 0],
    ], dtype=float)
    layer2 = np.array([
        [0, 1, 10],
        [1, 0, 0],
        [10, 0, 0],
    ], dtype=float)

    out = fit_multilayer_weighted_jaccard([layer1, layer2], algorithm="louvain", add_self_loops=False)
    ties = out["interlayer_ties"]

    assert np.isclose(ties.iloc[0]["similarity"], 13 / 31)


def test_weighted_overlap_uses_node_strength_weights():
    layer1 = np.array([
        [0, 10, 1],
        [10, 0, 0],
        [1, 0, 0],
    ], dtype=float)
    layer2 = np.array([
        [0, 1, 10],
        [1, 0, 0],
        [10, 0, 0],
    ], dtype=float)

    out = fit_multilayer_weighted_overlap([layer1, layer2], algorithm="louvain", add_self_loops=False)
    ties = out["interlayer_ties"]

    assert np.isclose(ties.iloc[0]["similarity"], 13 / 22)


def test_negative_weights_rejected_louvain():
    import pytest

    layer = np.array([
        [0, 1, -1],
        [1, 0, 1],
        [-1, 1, 0],
    ], dtype=float)

    with pytest.raises(ValueError, match="negative edge weights"):
        fit_multilayer_jaccard([layer, layer], algorithm="louvain")


def test_negative_weights_rejected_leiden_undirected():
    import pytest

    layer = np.array([
        [0, 1, -1],
        [1, 0, 1],
        [-1, 1, 0],
    ], dtype=float)

    with pytest.raises(ValueError, match="negative edge weights"):
        fit_multilayer_jaccard([layer, layer], algorithm="leiden")
