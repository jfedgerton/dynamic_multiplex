"""
Tests for node-universe alignment, arbitrary node labels, and seed reproducibility.

Covers:
1. Unequal node-universe rejection across all fit methods
2. String and non-contiguous integer node labels via NetworkX graphs
3. Seed reproducibility for Louvain community detection
"""

import networkx as nx
import numpy as np
import pytest

from dynamic_multiplex import (
    fit_multilayer_identity_ties,
    fit_multilayer_jaccard,
    fit_multilayer_overlap,
    fit_multilayer_weighted_jaccard,
    fit_multilayer_weighted_overlap,
)

ALL_FIT_FNS = [
    fit_multilayer_jaccard,
    fit_multilayer_overlap,
    fit_multilayer_weighted_jaccard,
    fit_multilayer_weighted_overlap,
    fit_multilayer_identity_ties,
]


# -----------------------------------------------------------------------
# 1. Unequal node-universe rejection
# -----------------------------------------------------------------------


class TestUnequalNodeUniverseRejection:
    def test_different_matrix_sizes_rejected(self):
        """Layers with different matrix sizes should raise ValueError."""
        layer_small = np.eye(4)
        layer_large = np.eye(6)
        for fn in ALL_FIT_FNS:
            with pytest.raises(ValueError, match="same node"):
                fn([layer_small, layer_large], algorithm="louvain")

    def test_different_networkx_node_sets_rejected(self):
        """NetworkX graphs with different node sets should raise ValueError."""
        g1 = nx.Graph()
        g1.add_edges_from([(0, 1), (1, 2)])
        g2 = nx.Graph()
        g2.add_edges_from([(0, 1), (1, 3)])  # node 3 not in g1, node 2 not in g2

        for fn in ALL_FIT_FNS:
            with pytest.raises(ValueError, match="same node"):
                fn([g1, g2], algorithm="louvain")

    def test_same_size_matrices_accepted(self):
        """Layers with same-sized matrices should work fine."""
        layers = [np.eye(5), np.eye(5)]
        out = fit_multilayer_jaccard(layers, algorithm="louvain")
        assert len(out["layer_communities"]) == 2

    def test_three_layers_one_different_rejected(self):
        """Three layers where one differs should raise ValueError."""
        layers = [np.eye(5), np.eye(5), np.eye(4)]
        with pytest.raises(ValueError, match="same node"):
            fit_multilayer_jaccard(layers, algorithm="louvain")


# -----------------------------------------------------------------------
# 2. String and non-contiguous node labels
# -----------------------------------------------------------------------


class TestArbitraryNodeLabels:
    def _make_string_labeled_graph(self, labels, edges, weights=None):
        g = nx.Graph()
        g.add_nodes_from(labels)
        for i, (u, v) in enumerate(edges):
            w = weights[i] if weights else 1.0
            g.add_edge(u, v, weight=w)
        return g

    def test_string_node_labels_jaccard(self):
        """Graphs with string node labels should produce valid 1-indexed output."""
        labels = ["alice", "bob", "carol", "dave"]
        edges = [("alice", "bob"), ("bob", "carol"), ("carol", "dave"), ("alice", "carol")]
        g1 = self._make_string_labeled_graph(labels, edges)
        g2 = self._make_string_labeled_graph(labels, edges)

        out = fit_multilayer_jaccard([g1, g2], algorithm="louvain")
        for lc in out["layer_communities"]:
            assert len(lc.membership) == 4
            assert min(lc.membership.keys()) >= 1
            assert max(lc.membership.keys()) <= 4

    def test_string_labels_all_methods(self):
        """All fit methods should handle string-labeled graphs."""
        labels = ["x", "y", "z"]
        edges = [("x", "y"), ("y", "z")]
        g1 = self._make_string_labeled_graph(labels, edges)
        g2 = self._make_string_labeled_graph(labels, edges)

        for fn in ALL_FIT_FNS:
            out = fn([g1, g2], algorithm="louvain")
            assert len(out["layer_communities"]) == 2
            for lc in out["layer_communities"]:
                assert len(lc.membership) == 3

    def test_non_contiguous_integer_labels(self):
        """Graphs with non-contiguous integer labels (e.g., 10, 20, 30)."""
        g1 = nx.Graph()
        g1.add_nodes_from([10, 20, 30, 40])
        g1.add_edges_from([(10, 20, {"weight": 1.0}), (30, 40, {"weight": 1.0})])
        g2 = nx.Graph()
        g2.add_nodes_from([10, 20, 30, 40])
        g2.add_edges_from([(10, 20, {"weight": 1.0}), (30, 40, {"weight": 1.0})])

        out = fit_multilayer_jaccard([g1, g2], algorithm="louvain")
        for lc in out["layer_communities"]:
            assert len(lc.membership) == 4
            assert min(lc.membership.keys()) >= 1
            assert max(lc.membership.keys()) <= 4

    def test_identity_ties_with_string_labels(self):
        """Identity ties should produce correct node range with string labels."""
        labels = ["a", "b", "c", "d", "e"]
        edges = [("a", "b"), ("c", "d"), ("d", "e")]
        g1 = self._make_string_labeled_graph(labels, edges)
        g2 = self._make_string_labeled_graph(labels, edges)

        out = fit_multilayer_identity_ties([g1, g2], algorithm="louvain")
        ties = out["interlayer_ties"]
        # 5 nodes * 1 link = 5 rows
        assert ties.shape[0] == 5
        assert set(ties["node"]) == {1, 2, 3, 4, 5}

    def test_weighted_methods_with_string_labels(self):
        """Weighted methods should compute correct node strengths after remapping."""
        labels = ["p", "q", "r"]
        edges = [("p", "q"), ("q", "r")]
        weights = [3.0, 7.0]
        g1 = self._make_string_labeled_graph(labels, edges, weights)
        g2 = self._make_string_labeled_graph(labels, edges, weights)

        out = fit_multilayer_weighted_jaccard([g1, g2], algorithm="louvain", add_self_loops=False)
        assert "interlayer_ties" in out
        assert "weighting" in out
        assert out["weighting"] == "node_strength"


# -----------------------------------------------------------------------
# 3. Seed reproducibility
# -----------------------------------------------------------------------


class TestSeedReproducibility:
    def _make_layers(self, seed=42):
        rng = np.random.default_rng(seed)
        n = 30
        memberships = rng.integers(1, 4, size=n)
        layers = []
        for _ in range(3):
            mat = np.zeros((n, n))
            for i in range(n - 1):
                for j in range(i + 1, n):
                    prob = 0.4 if memberships[i] == memberships[j] else 0.05
                    tie = rng.binomial(1, prob)
                    mat[i, j] = tie
                    mat[j, i] = tie
            layers.append(mat)
        return layers

    def test_seed_produces_identical_results(self):
        """Same seed should produce identical community assignments."""
        layers = self._make_layers()

        out1 = fit_multilayer_jaccard(layers, algorithm="louvain", seed=123)
        out2 = fit_multilayer_jaccard(layers, algorithm="louvain", seed=123)

        for lc1, lc2 in zip(out1["layer_communities"], out2["layer_communities"]):
            assert lc1.membership == lc2.membership
            assert lc1.communities == lc2.communities

    def test_different_seeds_may_differ(self):
        """Different seeds should generally produce different results (not guaranteed
        but very likely for stochastic Louvain on non-trivial graphs)."""
        layers = self._make_layers()

        out1 = fit_multilayer_jaccard(layers, algorithm="louvain", seed=1)
        out2 = fit_multilayer_jaccard(layers, algorithm="louvain", seed=999)

        # Just check they run without error; exact equality is not guaranteed
        assert len(out1["layer_communities"]) == 3
        assert len(out2["layer_communities"]) == 3

    def test_seed_reproducibility_all_methods(self):
        """All fit methods should support seed parameter."""
        layers = self._make_layers()

        for fn in ALL_FIT_FNS:
            out1 = fn(layers, algorithm="louvain", seed=42)
            out2 = fn(layers, algorithm="louvain", seed=42)
            for lc1, lc2 in zip(out1["layer_communities"], out2["layer_communities"]):
                assert lc1.membership == lc2.membership

    def test_no_seed_runs_without_error(self):
        """Passing seed=None (default) should work fine."""
        layers = self._make_layers()
        out = fit_multilayer_jaccard(layers, algorithm="louvain", seed=None)
        assert len(out["layer_communities"]) == 3

    def test_seed_in_simulate_and_fit(self):
        """simulate_and_fit_multilayer with seed should be reproducible."""
        from dynamic_multiplex import simulate_and_fit_multilayer

        out1 = simulate_and_fit_multilayer(n_nodes=20, n_layers=2, seed=77)
        out2 = simulate_and_fit_multilayer(n_nodes=20, n_layers=2, seed=77)

        for lc1, lc2 in zip(
            out1["fit"]["layer_communities"], out2["fit"]["layer_communities"]
        ):
            assert lc1.membership == lc2.membership
