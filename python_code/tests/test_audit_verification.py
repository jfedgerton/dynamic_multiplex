"""
Verification tests added during the code audit.

These tests check:
1. Manual Jaccard/overlap computation on a known toy graph
2. Temporal adjacency correctness (no off-by-one)
3. Directed graph modularity returns None
4. Community ID stability in interlayer ties
5. Supra-structure dimensions check
6. Cross-method consistency on identical input
"""

import numpy as np
import pytest

from dynamic_multiplex import (
    fit_multilayer_identity_ties,
    fit_multilayer_jaccard,
    fit_multilayer_overlap,
    fit_multilayer_weighted_jaccard,
    fit_multilayer_weighted_overlap,
)
from dynamic_multiplex.multilayer_utils import (
    LayerCommunityFit,
    community_overlap_edges,
    make_layer_links,
    weighted_jaccard,
    weighted_jaccard_similarity,
    weighted_overlap,
    weighted_overlap_similarity,
)


# -----------------------------------------------------------------------
# 1. Manual Jaccard / overlap on known sets
# -----------------------------------------------------------------------


def test_jaccard_manual_computation():
    """Verify Jaccard on hand-computed example: {1,2,3} vs {2,3,4}."""
    a = [1, 2, 3]
    b = [2, 3, 4]
    # intersection = {2,3} => 2, union = {1,2,3,4} => 4
    assert weighted_jaccard(a, b) == pytest.approx(2 / 4)


def test_overlap_manual_computation():
    """Verify overlap on hand-computed example: {1,2,3} vs {2,3,4}."""
    a = [1, 2, 3]
    b = [2, 3, 4]
    # intersection = {2,3} => 2, min_size = min(3,3) = 3
    assert weighted_overlap(a, b) == pytest.approx(2 / 3)


def test_jaccard_disjoint_sets():
    """Disjoint sets should have Jaccard = 0."""
    assert weighted_jaccard([1, 2], [3, 4]) == 0.0


def test_jaccard_identical_sets():
    """Identical sets should have Jaccard = 1."""
    assert weighted_jaccard([1, 2, 3], [1, 2, 3]) == 1.0


def test_overlap_subset():
    """If A is a subset of B, overlap = 1."""
    assert weighted_overlap([1, 2], [1, 2, 3]) == 1.0


def test_weighted_jaccard_manual():
    """Weighted Jaccard with known node strengths."""
    a = [1, 2, 3]
    b = [2, 3, 4]
    wa = {1: 5.0, 2: 3.0, 3: 1.0}
    wb = {2: 2.0, 3: 4.0, 4: 6.0}
    # intersection = {2, 3}
    #   min(3,2) + min(1,4) = 2 + 1 = 3
    # union = {1,2,3,4}
    #   max(5,0) + max(3,2) + max(1,4) + max(0,6) = 5 + 3 + 4 + 6 = 18
    assert weighted_jaccard_similarity(a, b, wa, wb) == pytest.approx(3 / 18)


def test_weighted_overlap_manual():
    """Weighted overlap with known node strengths."""
    a = [1, 2, 3]
    b = [2, 3, 4]
    wa = {1: 5.0, 2: 3.0, 3: 1.0}
    wb = {2: 2.0, 3: 4.0, 4: 6.0}
    # intersection weight = min(3,2) + min(1,4) = 3
    # a_weight = 5+3+1 = 9, b_weight = 2+4+6 = 12
    # min_weight = 9
    assert weighted_overlap_similarity(a, b, wa, wb) == pytest.approx(3 / 9)


# -----------------------------------------------------------------------
# 2. Temporal adjacency — default layer links correctness
# -----------------------------------------------------------------------


def test_default_layer_links_chain_topology():
    """Default links should be 1->2, 2->3, ..., (n-1)->n."""
    links = make_layer_links(5)
    assert list(links["from"]) == [1, 2, 3, 4]
    assert list(links["to"]) == [2, 3, 4, 5]
    assert all(links["weight"] == 1.0)


def test_layer_links_no_skip():
    """Default links should NOT skip any layer."""
    links = make_layer_links(4)
    # Verify every consecutive pair is present
    for i in range(1, 4):
        row = links[(links["from"] == i) & (links["to"] == i + 1)]
        assert len(row) == 1, f"Missing link {i} -> {i+1}"


# -----------------------------------------------------------------------
# 3. Directed graph modularity should be None
# -----------------------------------------------------------------------


def test_directed_modularity_is_none():
    """For directed graphs, modularity should be None."""
    mat = np.array([
        [0, 1, 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1],
        [1, 0, 0, 0],
    ], dtype=float)
    layers = [mat, mat]
    out = fit_multilayer_jaccard(layers, algorithm="louvain", directed=True)
    for lc in out["layer_communities"]:
        assert lc.modularity is None, "Directed Louvain should return modularity=None"


# -----------------------------------------------------------------------
# 4. Interlayer ties connect only specified layer pairs
# -----------------------------------------------------------------------


def test_interlayer_ties_only_specified_pairs():
    """Interlayer ties should only appear between layer pairs in layer_links."""
    # Create 4 layers but only link 1->3
    rng = np.random.default_rng(42)
    layers = [rng.binomial(1, 0.3, (10, 10)).astype(float) for _ in range(4)]
    # Make symmetric
    for i in range(4):
        layers[i] = np.triu(layers[i], 1)
        layers[i] = layers[i] + layers[i].T

    custom_links = [{"from": 1, "to": 3, "weight": 1.0}]
    out = fit_multilayer_jaccard(layers, algorithm="louvain", layer_links=custom_links,
                                 add_self_loops=False)
    ties = out["interlayer_ties"]
    if not ties.empty:
        # All from_layer should be 1 and to_layer should be 3
        assert set(ties["from_layer"].unique()) == {1}
        assert set(ties["to_layer"].unique()) == {3}


# -----------------------------------------------------------------------
# 5. Supra-structure dimension checks
# -----------------------------------------------------------------------


def test_identity_ties_dimensions():
    """Identity ties: n_nodes * n_links rows expected."""
    n_nodes = 8
    layers = [np.eye(n_nodes)] * 5
    out = fit_multilayer_identity_ties(layers, algorithm="louvain")
    # Default links: 4 links for 5 layers
    expected_rows = n_nodes * 4
    assert out["interlayer_ties"].shape[0] == expected_rows


def test_community_overlap_edges_all_pairs_enumerated():
    """community_overlap_edges should create rows for all community pairs."""
    # Two layers, each with 2 communities, 1 link => 2*2 = 4 potential edges
    fit = [
        LayerCommunityFit(
            membership={1: 1, 2: 1, 3: 2, 4: 2},
            modularity=0.5,
            communities={1: [1, 2], 2: [3, 4]},
        ),
        LayerCommunityFit(
            membership={1: 1, 2: 1, 3: 2, 4: 2},
            modularity=0.5,
            communities={1: [1, 2], 2: [3, 4]},
        ),
    ]
    import pandas as pd
    links = pd.DataFrame({"from": [1], "to": [2], "weight": [1.0]})
    edges = community_overlap_edges(fit, links, metric="jaccard", min_similarity=0.0)
    # 2 from_communities x 2 to_communities = 4
    assert len(edges) == 4


def test_self_loops_directed_vs_undirected():
    """Self-loop similarity should be 1.0 for directed, 2.0 for undirected."""
    mat = np.array([
        [0, 1, 0, 0],
        [1, 0, 0, 0],
        [0, 0, 0, 1],
        [0, 0, 1, 0],
    ], dtype=float)
    layers = [mat, mat]

    out_undir = fit_multilayer_jaccard(layers, algorithm="louvain", directed=False)
    ties_u = out_undir["interlayer_ties"]
    loops_u = ties_u[(ties_u["from_layer"] == ties_u["to_layer"]) &
                     (ties_u["from_community"] == ties_u["to_community"])]
    if not loops_u.empty:
        assert (loops_u["similarity"] == 2.0).all()

    # For directed, create a directed graph
    mat_dir = np.array([
        [0, 1, 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1],
        [0, 0, 0, 0],
    ], dtype=float)
    layers_dir = [mat_dir, mat_dir]
    out_dir = fit_multilayer_jaccard(layers_dir, algorithm="louvain", directed=True)
    ties_d = out_dir["interlayer_ties"]
    loops_d = ties_d[(ties_d["from_layer"] == ties_d["to_layer"]) &
                     (ties_d["from_community"] == ties_d["to_community"])]
    if not loops_d.empty:
        assert (loops_d["similarity"] == 1.0).all()


# -----------------------------------------------------------------------
# 6. Cross-method consistency
# -----------------------------------------------------------------------


def test_all_methods_produce_valid_output_structure():
    """All five fit methods should produce structurally valid outputs."""
    rng = np.random.default_rng(99)
    n = 20
    mat = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            if (i < 10 and j < 10) or (i >= 10 and j >= 10):
                mat[i, j] = mat[j, i] = rng.binomial(1, 0.7)
            else:
                mat[i, j] = mat[j, i] = rng.binomial(1, 0.05)
    layers = [mat, mat]

    fns = [
        fit_multilayer_jaccard,
        fit_multilayer_overlap,
        fit_multilayer_weighted_jaccard,
        fit_multilayer_weighted_overlap,
        fit_multilayer_identity_ties,
    ]
    for fn in fns:
        out = fn(layers, algorithm="louvain")
        assert "layer_communities" in out
        assert "layer_links" in out
        assert "interlayer_ties" in out
        assert len(out["layer_communities"]) == 2
        # Each community fit should have membership covering all nodes
        for lc in out["layer_communities"]:
            assert len(lc.membership) == n
            # All node IDs should be 1-indexed
            assert min(lc.membership.keys()) >= 1
            assert max(lc.membership.keys()) <= n


# -----------------------------------------------------------------------
# 7. Edge cases
# -----------------------------------------------------------------------


def test_empty_graph_layers():
    """Layers with no edges should not crash."""
    layers = [np.zeros((5, 5)), np.zeros((5, 5))]
    out = fit_multilayer_jaccard(layers, algorithm="louvain")
    assert len(out["layer_communities"]) == 2


def test_single_node_layers():
    """Single-node layers should not crash."""
    layers = [np.zeros((1, 1)), np.zeros((1, 1))]
    out = fit_multilayer_jaccard(layers, algorithm="louvain")
    assert len(out["layer_communities"]) == 2
