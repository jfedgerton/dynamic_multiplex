# Manuscript-ready simulation benchmark script
#
# Compares bayesnet interlayer tie recovery against common all-to-all coupling
# baselines in evolving multiplex systems and reports both quality + speed.

library(bayesnet)

set.seed(2026)

bench <- benchmark_interlayer_tie_recovery(
  n_reps = 100,
  simulation_args = list(
    n_nodes = 200,
    n_layers = 8,
    n_communities = 5,
    persistence = 0.90,
    p_in = 0.20,
    p_out = 0.02,
    directed = TRUE
  ),
  algorithm = "leiden",
  min_similarity = 0.05,
  threshold = 0.05
)

summary_tbl <- bench$summary[order(-bench$summary$f1), ]
print(summary_tbl)

cat("\nRecommended manuscript table columns:\n")
print(summary_tbl[, c(
  "model", "precision", "recall", "f1",
  "elapsed_sec", "f1_per_second", "mean_layer_ari", "mean_layer_modularity"
)])

# Optional manuscript figures with ggplot2 if installed
if (requireNamespace("ggplot2", quietly = TRUE)) {
  df <- bench$raw_results

  p_f1 <- ggplot2::ggplot(df, ggplot2::aes(x = model, y = f1, fill = model)) +
    ggplot2::geom_boxplot(alpha = 0.8, outlier.size = 0.5) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::labs(
      title = "Interlayer tie recovery in evolving multiplex networks",
      x = "Model",
      y = "F1 score"
    )

  p_speed <- ggplot2::ggplot(df[df$elapsed_sec > 0, ], ggplot2::aes(x = model, y = elapsed_sec, fill = model)) +
    ggplot2::geom_boxplot(alpha = 0.8, outlier.size = 0.5) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::labs(
      title = "Convergence speed proxy (runtime)",
      x = "Model",
      y = "Elapsed seconds"
    )

  p_ari <- ggplot2::ggplot(df, ggplot2::aes(x = model, y = mean_layer_ari, fill = model)) +
    ggplot2::geom_boxplot(alpha = 0.8, outlier.size = 0.5) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::labs(
      title = "Within-layer community recovery (ARI)",
      x = "Model",
      y = "Mean layer ARI"
    )

  print(p_f1)
  print(p_speed)
  print(p_ari)
}
