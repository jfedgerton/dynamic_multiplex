import numpy as np
import pytest

from dynamic_multiplex import (
    partition_stability,
    bootstrap_multilayer,
    extract_meta_membership,
    fit_multilayer_identity_ties,
    fit_multilayer_jaccard,
)


def _make_planted_layers(n_nodes=30, n_layers=3, seed=42):
    """Create small planted-partition layers for testing."""
    rng = np.random.default_rng(seed)
    memberships = rng.integers(1, 4, size=n_nodes)

    layers = []
    for _ in range(n_layers):
        mat = np.zeros((n_nodes, n_nodes))
        for i in range(n_nodes - 1):
            for j in range(i + 1, n_nodes):
                prob = 0.4 if memberships[i] == memberships[j] else 0.05
                tie = rng.binomial(1, prob)
                mat[i, j] = tie
                mat[j, i] = tie
        layers.append(mat)
    return layers


class TestBootstrapMultilayer:
    def test_basic_smoke(self):
        layers = _make_planted_layers()
        result = bootstrap_multilayer(layers, fit_type="jaccard", n_boot=5, seed=1)

        assert result.n_boot == 5
        assert len(result.co_assignment) == 3
        assert len(result.node_stability) == 3
        assert len(result.modularity_samples) == 3
        assert len(result.community_count_reproducibility) == 3
        assert result.point_estimate is not None

    def test_co_assignment_shape_and_range(self):
        layers = _make_planted_layers(n_nodes=20, n_layers=2)
        result = bootstrap_multilayer(layers, fit_type="jaccard", n_boot=10, seed=2)

        for co in result.co_assignment:
            assert co.shape == (20, 20)
            assert np.all(co >= 0.0)
            assert np.all(co <= 1.0)
            # Diagonal should be 1.0
            np.testing.assert_array_equal(np.diag(co), np.ones(20))
            # Symmetric
            np.testing.assert_array_almost_equal(co, co.T)

    def test_node_stability_range(self):
        layers = _make_planted_layers(n_nodes=20, n_layers=2)
        result = bootstrap_multilayer(layers, fit_type="overlap", n_boot=10, seed=3)

        for stab in result.node_stability:
            assert stab.shape == (20,)
            assert np.all(stab >= 0.0)
            assert np.all(stab <= 1.0)

    def test_modularity_samples_length(self):
        layers = _make_planted_layers(n_nodes=20, n_layers=2)
        result = bootstrap_multilayer(layers, fit_type="jaccard", n_boot=8, seed=4)

        for mod_s in result.modularity_samples:
            assert len(mod_s) == result.n_boot

    def test_community_count_reproducibility_range(self):
        layers = _make_planted_layers(n_nodes=20, n_layers=2)
        result = bootstrap_multilayer(layers, fit_type="jaccard", n_boot=5, seed=5)

        r = result.community_count_reproducibility
        assert len(r) == 2
        assert all(0.0 <= x <= 1.0 for x in r)

    def test_all_fit_types(self):
        layers = _make_planted_layers(n_nodes=20, n_layers=2)
        for fit_type in ["jaccard", "overlap", "weighted_jaccard",
                         "weighted_overlap", "identity"]:
            result = bootstrap_multilayer(
                layers, fit_type=fit_type, n_boot=3, seed=6
            )
            assert result.n_boot == 3

    def test_invalid_fit_type(self):
        layers = _make_planted_layers(n_nodes=10, n_layers=2)
        with pytest.raises(ValueError, match="fit_type"):
            bootstrap_multilayer(layers, fit_type="invalid", n_boot=1)

    def test_custom_layer_links(self):
        layers = _make_planted_layers(n_nodes=15, n_layers=4)
        layer_links = [{"from": 1, "to": 3, "weight": 0.5}]
        result = bootstrap_multilayer(
            layers, fit_type="jaccard", n_boot=3, seed=7,
            layer_links=layer_links
        )
        assert result.n_boot == 3
        # Point estimate should respect custom links
        assert result.point_estimate["layer_links"].shape[0] == 1


class TestEdgeResampling:
    def test_edge_bootstrap_runs(self):
        layers = _make_planted_layers(n_nodes=20, n_layers=2)
        boot = bootstrap_multilayer(
            layers, fit_type="jaccard", n_boot=6, seed=41)
        assert boot.n_boot == 6
        assert (boot.co_assignment[0] >= 0).all()
        assert (boot.co_assignment[0] <= 1).all()

    def test_no_resample_argument(self):
        # the legacy weights scheme was removed entirely in 1.1.0
        with pytest.raises(TypeError):
            bootstrap_multilayer(
                _make_planted_layers(n_nodes=15, n_layers=2),
                fit_type="jaccard", n_boot=3, resample="weights")


class TestMetaCommunities:
    def test_fit_returns_meta_communities(self):
        layers = _make_planted_layers(n_nodes=20, n_layers=4)
        fit = fit_multilayer_jaccard(layers, algorithm="leiden")
        assert "meta_communities" in fit
        assert fit["meta_communities"] is not None
        assert len(fit["meta_communities"]) == 4
        assert all(len(v) == 20 for v in fit["meta_communities"])
        # Meta communities cannot exceed the total per-layer communities.
        total_layer_comms = sum(
            len(lc.communities) for lc in fit["layer_communities"]
        )
        all_meta = np.concatenate(fit["meta_communities"])
        assert len(np.unique(all_meta)) <= total_layer_comms

    def test_extract_meta_membership(self):
        layers = _make_planted_layers(n_nodes=20, n_layers=3)
        fit = fit_multilayer_jaccard(layers, algorithm="leiden")
        mm = extract_meta_membership(fit)
        assert mm is fit["meta_communities"]
        assert len(mm) == 3

    def test_extract_meta_membership_missing_raises(self):
        with pytest.raises(ValueError, match="meta_communities"):
            extract_meta_membership({"layer_communities": []})

    def test_custom_layer_links_flow_through(self):
        layers = _make_planted_layers(n_nodes=20, n_layers=4)
        links = [{"from": 1, "to": 3, "weight": 1.0},
                 {"from": 2, "to": 4, "weight": 1.0}]
        fit = fit_multilayer_jaccard(
            layers, algorithm="leiden", layer_links=links
        )
        assert len(fit["meta_communities"]) == 4
        assert all(np.all(v >= 1) for v in fit["meta_communities"])

    def test_bootstrap_and_stability_run_on_meta(self):
        # 1.4.0: interval outputs removed; the bootstrap feeds only the stability score s
        layers = _make_planted_layers(n_nodes=20, n_layers=3)
        boot = bootstrap_multilayer(
            layers, fit_type="jaccard", algorithm="leiden",
            n_boot=6, seed=123,
        )
        assert boot.n_boot == 6
        ps = partition_stability(boot)
        assert 0 <= ps["stability"] <= 1 and len(ps["stability_by_layer"]) == 3

    def test_identity_fit_falls_back_without_error(self):
        layers = _make_planted_layers(n_nodes=20, n_layers=3)
        fit = fit_multilayer_identity_ties(layers, algorithm="leiden")
        # Identity now uses the node-level Mucha multislice; meta_ids is None.
        assert fit["meta_ids"] is None
        assert len(fit["meta_communities"]) == 3
        assert all(len(v) == 20 for v in fit["meta_communities"])

    def test_identity_bootstrap_runs(self):
        layers = _make_planted_layers(n_nodes=20, n_layers=2)
        boot = bootstrap_multilayer(
            layers, fit_type="identity", algorithm="leiden",
            n_boot=3, seed=9,
        )
        assert boot.n_boot == 3


class TestMultisliceIdentity:
    def test_identity_meta_shapes(self):
        layers = _make_planted_layers(n_nodes=20, n_layers=4)
        fit = fit_multilayer_identity_ties(layers, algorithm="leiden")
        assert fit["meta_ids"] is None
        assert len(fit["meta_communities"]) == 4
        assert all(len(v) == 20 for v in fit["meta_communities"])

    def test_multislice_not_degenerate(self):
        # A single supra-graph detection should recover more than one
        # meta-community overall (not collapse to one giant community).
        layers = _make_planted_layers(n_nodes=30, n_layers=4)
        fit = fit_multilayer_identity_ties(layers, algorithm="leiden")
        all_meta = np.concatenate(fit["meta_communities"])
        assert len(np.unique(all_meta)) > 1

    def test_custom_layer_links_change_partition(self):
        # All-to-all coupling ties nodes across every layer pair, so the meta
        # partition differs from the adjacent-only default coupling.
        layers = _make_planted_layers(n_nodes=25, n_layers=4)
        fit_default = fit_multilayer_identity_ties(layers, algorithm="leiden")

        n = 4
        all_to_all = [
            {"from": a, "to": b, "weight": 1.0}
            for a in range(1, n + 1)
            for b in range(a + 1, n + 1)
        ]
        fit_full = fit_multilayer_identity_ties(
            layers, algorithm="leiden", layer_links=all_to_all
        )

        default_meta = [v.tolist() for v in fit_default["meta_communities"]]
        full_meta = [v.tolist() for v in fit_full["meta_communities"]]
        assert default_meta != full_meta
