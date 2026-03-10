import numpy as np
import pytest

from dynamic_multiplex import bootstrap_multilayer, community_ci


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
        assert len(result.community_count_samples) == 3
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

    def test_community_count_samples_positive(self):
        layers = _make_planted_layers(n_nodes=20, n_layers=2)
        result = bootstrap_multilayer(layers, fit_type="jaccard", n_boot=5, seed=5)

        for count_s in result.community_count_samples:
            assert len(count_s) == result.n_boot
            assert np.all(count_s >= 1)

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


class TestCommunityCi:
    def test_basic_output_structure(self):
        layers = _make_planted_layers(n_nodes=20, n_layers=3)
        boot = bootstrap_multilayer(layers, fit_type="jaccard", n_boot=10, seed=10)
        ci = community_ci(boot)

        assert "modularity_ci" in ci
        assert "community_count_ci" in ci
        assert "mean_node_stability" in ci
        assert "node_stability" in ci
        assert "co_assignment" in ci

        # Check DataFrame shapes
        assert ci["modularity_ci"].shape[0] == 3
        assert ci["community_count_ci"].shape[0] == 3
        assert ci["mean_node_stability"].shape[0] == 3

    def test_ci_columns(self):
        layers = _make_planted_layers(n_nodes=20, n_layers=2)
        boot = bootstrap_multilayer(layers, fit_type="jaccard", n_boot=10, seed=11)
        ci = community_ci(boot)

        for df_name in ["modularity_ci", "community_count_ci"]:
            df = ci[df_name]
            assert set(df.columns) == {"layer", "estimate", "lower", "upper"}

    def test_lower_le_upper(self):
        layers = _make_planted_layers(n_nodes=20, n_layers=2)
        boot = bootstrap_multilayer(layers, fit_type="jaccard", n_boot=20, seed=12)
        ci = community_ci(boot)

        mod = ci["modularity_ci"]
        valid = mod.dropna(subset=["lower", "upper"])
        assert (valid["lower"] <= valid["upper"]).all()

        count = ci["community_count_ci"]
        assert (count["lower"] <= count["upper"]).all()

    def test_custom_alpha(self):
        layers = _make_planted_layers(n_nodes=20, n_layers=2)
        boot = bootstrap_multilayer(layers, fit_type="jaccard", n_boot=20, seed=13)

        ci_90 = community_ci(boot, alpha=0.10)
        ci_50 = community_ci(boot, alpha=0.50)

        # Wider alpha should give narrower intervals
        for i in range(2):
            width_90 = (ci_90["modularity_ci"].iloc[i]["upper"] -
                        ci_90["modularity_ci"].iloc[i]["lower"])
            width_50 = (ci_50["modularity_ci"].iloc[i]["upper"] -
                        ci_50["modularity_ci"].iloc[i]["lower"])
            if not np.isnan(width_90) and not np.isnan(width_50):
                assert width_50 <= width_90 + 1e-10

    def test_mean_stability_in_range(self):
        layers = _make_planted_layers(n_nodes=20, n_layers=2)
        boot = bootstrap_multilayer(layers, fit_type="jaccard", n_boot=10, seed=14)
        ci = community_ci(boot)

        stab = ci["mean_node_stability"]
        assert (stab["mean_stability"] >= 0).all()
        assert (stab["mean_stability"] <= 1).all()

    def test_no_completed_replicates_raises(self):
        layers = _make_planted_layers(n_nodes=20, n_layers=2)
        boot = bootstrap_multilayer(layers, fit_type="jaccard", n_boot=10, seed=15)
        # Manually set n_boot to 0 to simulate failure
        boot.n_boot = 0
        with pytest.raises(ValueError, match="No completed"):
            community_ci(boot)
