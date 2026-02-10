import numpy as np

from dynamic_multiplex import (
    fit_multilayer_identity_ties,
    fit_multilayer_jaccard,
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


def test_custom_layer_links_are_respected():
    layers = [np.eye(5), np.eye(5), np.eye(5), np.eye(5)]
    layer_links = [{"from": 1, "to": 3, "weight": 0.5}]
    out = fit_multilayer_jaccard(layers, algorithm="louvain", layer_links=layer_links)
    links = out["layer_links"]
    assert links.shape[0] == 1
    assert int(links.iloc[0]["from"]) == 1
    assert int(links.iloc[0]["to"]) == 3
