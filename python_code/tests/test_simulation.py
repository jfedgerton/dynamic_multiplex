from bayesnet_py import simulate_and_fit_multilayer


def test_simulation_shape():
    out = simulate_and_fit_multilayer(n_nodes=30, n_layers=3, n_communities=3, seed=123)
    assert len(out["layers"]) == 3
    assert len(out["true_membership"]) == 30
    assert "fit" in out


def test_directed_option():
    out = simulate_and_fit_multilayer(
        n_nodes=20, n_layers=3, fit_type="identity", algorithm="leiden", directed=True, seed=1
    )
    assert out["directed"] is True
    layer = out["layers"][0]
    assert (layer != layer.T).any()
