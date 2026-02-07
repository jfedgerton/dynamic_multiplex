from bayesnet_py import benchmark_interlayer_tie_recovery


def main() -> None:
    bench = benchmark_interlayer_tie_recovery(
        n_reps=100,
        simulation_args={
            "n_nodes": 200,
            "n_layers": 8,
            "n_communities": 5,
            "persistence": 0.90,
            "p_in": 0.20,
            "p_out": 0.02,
            "directed": True,
        },
        algorithm="leiden",
        min_similarity=0.05,
        threshold=0.05,
        seed=2026,
    )

    print(bench["summary"].sort_values("f1", ascending=False))


if __name__ == "__main__":
    main()
