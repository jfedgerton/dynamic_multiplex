#!/usr/bin/env Rscript
# 22_coverage_heatmap.R -- Ungated empirical coverage by network configuration.
# Heatmap tiles: within-community tie probability (p_in, x) by
# between-community tie probability (p_out, y); facet_grid rows = K
# (number of communities), columns = n (network size). Fill diverges around
# the nominal 0.95 (blue above, white at nominal, red below).
# Output: manuscript/figures/fig_coverage_heatmap.pdf/.png
set.seed(123)
suppressMessages({ library(ggplot2) })
OUT <- "/storage/group/LiberalArts/default/jfe4_collab/dynamic_multiplex/manuscript/output"
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
ggsave("manuscript/figures/fig_coverage_heatmap.pdf", p, width = 6.5, height = 5.2)
ggsave("manuscript/figures/fig_coverage_heatmap.png", p, width = 6.5, height = 5.2, dpi = 300)
cat("done 22\n")
