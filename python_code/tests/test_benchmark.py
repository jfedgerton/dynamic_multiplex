from bayesnet_py import benchmark_interlayer_tie_recovery


def test_benchmark_metrics_exist():
    res = benchmark_interlayer_tie_recovery(
        n_reps=2,
        simulation_args={"n_nodes": 40, "n_layers": 4, "n_communities": 3},
        algorithm="louvain",
        seed=42,
    )
    assert "raw_results" in res
    assert "summary" in res

    cols = set(res["summary"].columns)
    required = {
        "precision",
        "recall",
        "f1",
        "elapsed_sec",
        "f1_per_second",
        "mean_layer_ari",
        "mean_layer_modularity",
        "predicted_ties",
        "truth_ties",
    }
    assert required.issubset(cols)
