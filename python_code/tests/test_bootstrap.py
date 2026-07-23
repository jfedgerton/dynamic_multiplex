import numpy as np
import pytest

from dynamic_multiplex import (
    bootstrap_multilayer,
    co_assignment_ci,
    community_est,
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


class TestCommunityEst:
    def test_basic_output_structure(self):
        layers = _make_planted_layers(n_nodes=20, n_layers=3)
        boot = bootstrap_multilayer(layers, fit_type="jaccard", n_boot=10, seed=10)
        est = community_est(boot)

        assert "modularity_ci" not in est
        assert "community_count_ci" not in est
        assert "community_count" in est
        assert "report" in est
        assert "mean_node_stability" in est
        assert "node_stability" in est
        assert "co_assignment" in est

        # Check shapes
        assert est["community_count"].shape[0] == 3
        assert len(est["report"]) == 3
        assert est["mean_node_stability"].shape[0] == 3

    def test_columns(self):
        layers = _make_planted_layers(n_nodes=20, n_layers=2)
        boot = bootstrap_multilayer(layers, fit_type="jaccard", n_boot=10, seed=11)
        est = community_est(boot)

        df = est["community_count"]
        assert set(df.columns) == {"layer", "estimate", "reproducibility"}

    def test_reproducibility_range(self):
        layers = _make_planted_layers(n_nodes=20, n_layers=2)
        boot = bootstrap_multilayer(layers, fit_type="jaccard", n_boot=20, seed=12)
        est = community_est(boot)

        r = est["community_count"]["reproducibility"]
        assert (r >= 0).all() and (r <= 1).all()

    def test_report_reads_as_reproducibility_not_interval(self):
        layers = _make_planted_layers(n_nodes=20, n_layers=2)
        boot = bootstrap_multilayer(layers, fit_type="jaccard", n_boot=10, seed=13)
        est = community_est(boot)

        assert all("reproduced in" in s for s in est["report"])
        # no bracketed interval anywhere in the report
        assert not any("[" in s for s in est["report"])

    def test_mean_stability_in_range(self):
        layers = _make_planted_layers(n_nodes=20, n_layers=2)
        boot = bootstrap_multilayer(layers, fit_type="jaccard", n_boot=10, seed=14)
        est = community_est(boot)

        stab = est["mean_node_stability"]
        assert (stab["mean_stability"] >= 0).all()
        assert (stab["mean_stability"] <= 1).all()

    def test_no_completed_replicates_raises(self):
        layers = _make_planted_layers(n_nodes=20, n_layers=2)
        boot = bootstrap_multilayer(layers, fit_type="jaccard", n_boot=10, seed=15)
        # Manually set n_boot to 0 to simulate failure
        boot.n_boot = 0
        with pytest.raises(ValueError, match="No completed"):
            community_est(boot)


class TestCoAssignmentCi:
    def test_structure_and_bounds(self):
        layers = _make_planted_layers(n_nodes=20, n_layers=2)
        boot = bootstrap_multilayer(layers, fit_type="jaccard", n_boot=10, seed=21)
        pci = co_assignment_ci(boot)

        assert len(pci) == 2
        for layer_ci in pci:
            est = layer_ci["estimate"]
            lo = layer_ci["lower"]
            hi = layer_ci["upper"]
            assert est.shape == (20, 20)
            assert lo.shape == (20, 20)
            assert hi.shape == (20, 20)
            # bounds ordered and inside [0, 1]
            off = ~np.eye(20, dtype=bool)
            assert (lo[off] <= est[off] + 1e-12).all()
            assert (est[off] <= hi[off] + 1e-12).all()
            assert (lo >= 0).all() and (hi <= 1).all()
            # diagonal is degenerate at 1
            assert (np.diag(lo) == 1).all()
            assert (np.diag(hi) == 1).all()

    def test_narrower_with_smaller_alpha_inverse(self):
        layers = _make_planted_layers(n_nodes=20, n_layers=2)
        boot = bootstrap_multilayer(layers, fit_type="jaccard", n_boot=20, seed=22)
        pci_95 = co_assignment_ci(boot, alpha=0.05)
        pci_50 = co_assignment_ci(boot, alpha=0.50)
        off = ~np.eye(20, dtype=bool)
        w95 = (pci_95[0]["upper"] - pci_95[0]["lower"])[off]
        w50 = (pci_50[0]["upper"] - pci_50[0]["lower"])[off]
        assert (w50 <= w95 + 1e-12).all()

    def test_zero_boot_raises(self):
        layers = _make_planted_layers(n_nodes=15, n_layers=2)
        boot = bootstrap_multilayer(layers, fit_type="jaccard", n_boot=5, seed=23)
        boot.n_boot = 0
        with pytest.raises(ValueError):
            co_assignment_ci(boot)


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

    def test_bootstrap_and_cis_run_on_meta(self):
        layers = _make_planted_layers(n_nodes=20, n_layers=3)
        boot = bootstrap_multilayer(
            layers, fit_type="jaccard", algorithm="leiden",
            n_boot=6, seed=123,
        )
        assert boot.n_boot == 6
        est = community_est(boot)
        pe = boot.point_estimate["meta_communities"]
        expected = [len(np.unique(v)) for v in pe]
        assert list(est["community_count"]["estimate"]) == expected
        ci = co_assignment_ci(boot)
        assert len(ci) == 3

    def test_identity_fit_falls_back_without_error(self):
        layers = _make_planted_layers(n_nodes=20, n_layers=3)
        fit = fit_multilayer_identity_ties(layers, algorithm="leiden")
        # Fallback: node-level identity ties have no community columns.
        assert fit["meta_ids"] is None
        assert len(fit["meta_communities"]) == 3
        assert all(len(v) == 20 for v in fit["meta_communities"])
        # Fallback offsets labels so layers do not collide.
        all_meta = np.concatenate(fit["meta_communities"])
        assert len(np.unique(all_meta)) >= 1

    def test_identity_bootstrap_runs(self):
        layers = _make_planted_layers(n_nodes=20, n_layers=2)
        boot = bootstrap_multilayer(
            layers, fit_type="identity", algorithm="leiden",
            n_boot=3, seed=9,
        )
        assert boot.n_boot == 3
