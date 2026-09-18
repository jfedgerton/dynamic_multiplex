#!/usr/bin/env Rscript
# 22_coverage_heatmap.R -- Ungated empirical coverage by network configuration.
# Heatmap tiles: within-community tie probability (p_in, x) by
# between-community tie probability (p_out, y); facet_grid rows = K
# (number of communities), columns = n (network size). Fill diverges around
# the nominal 0.95 (blue above, white at nominal, red below).
# Output: manuscript/figures/fig_coverage_heatmap.pdf/.png
set.seed(123)
suppressMessages({ library(ggplot2) })
# --- paths -------------------------------------------------------------
# DM_ROOT is the project root holding manuscript/output, manuscript/tables,
# and manuscript/figures. Set it before running, e.g.
#   export DM_ROOT=/path/to/dynamic_multiplex
# If unset, the current working directory is used, so the script also runs
# correctly when invoked from the project root.
ROOT <- Sys.getenv("DM_ROOT", unset = getwd())
OUT  <- file.path(ROOT, "manuscript", "output")
TAB  <- file.path(ROOT, "manuscript", "tables")
FIG  <- file.path(ROOT, "manuscript", "figures")
if (!dir.exists(OUT)) {
  stop("Simulation output not found at: ", OUT,
       "\n  Set DM_ROOT to the project root, or run from that directory.",
       call. = FALSE)
}
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
fs <- list.files(file.path(OUT, "coverage3_grid"), "^cov_task.*csv$", full.names = TRUE)
d <- do.call(rbind, lapply(fs, read.csv))
cat("rows:", nrow(d), "\n")
agg <- aggregate(cov_P_mean ~ p_in + p_out + n + K, data = d, FUN = mean)
cnt <- aggregate(cbind(nsims = cov_P_mean) ~ p_in + p_out + n + K, data = d, FUN = length)
agg <- merge(agg, cnt)
cat("cells:", nrow(agg), " n values:", paste(sort(unique(agg$n)), collapse = ","),
    " K values:", paste(sort(unique(agg$K)), collapse = ","), "\n")
agg$n_lab <- factor(paste0("n = ", agg$n), levels = paste0("n = ", sort(unique(agg$n))))
agg$K_lab <- factor(paste0("K = ", agg$K), levels = paste0("K = ", sort(unique(agg$K))))
p <- ggplot(agg, aes(factor(p_in), factor(p_out), fill = cov_P_mean)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", cov_P_mean)), size = 1.9) +
  facet_grid(K_lab ~ n_lab) +
  scale_fill_gradient2(low = "#c0392b", mid = "white", high = "#2166ac",
                       midpoint = 0.95, limits = c(min(agg$cov_P_mean), 1),
                       name = "Empirical\ncoverage") +
  labs(x = "Within-community tie probability", y = "Between-community tie probability") +
  theme_bw(base_size = 9) +
  theme(legend.position = "right", panel.grid = element_blank())
ggsave(file.path(FIG, "fig_coverage_heatmap.pdf"), p, width = 6.5, height = 5.2)
ggsave(file.path(FIG, "fig_coverage_heatmap.png"), p, width = 6.5, height = 5.2, dpi = 300)
cat("done 22\n")
