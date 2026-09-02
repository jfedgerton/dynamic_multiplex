#!/usr/bin/env Rscript
# 23_coverage_alt_figs.R -- Candidate Section 5.1 figures to replace the
# sparse p_in x p_out heatmap (density is a paired 3-level factor, not a grid):
#   (1) fig_coverage_by_config: coverage by n, separation regime, faceted by K
#   (2) fig_calibration: reliability diagram from the 20-bin calibration files
#   (3) fig_width_ecdf: ECDF of interval widths by n with the 0.05 gate line
# Data: coverage3_grid cov_task*.csv / calib_task*.csv (jfe4_collab storage).
set.seed(123)
suppressMessages(library(ggplot2))
OUT <- "/storage/group/LiberalArts/default/jfe4_collab/dynamic_multiplex/manuscript/output/coverage3_grid"
fs <- list.files(OUT, "^cov_task.*csv$", full.names = TRUE)
d <- do.call(rbind, lapply(fs, read.csv))
cat("cov rows:", nrow(d), "\n")
d$sep <- factor(d$density, levels = c("weak", "default", "strong"),
                labels = c("Weak (0.20 / 0.10)", "Moderate (0.30 / 0.05)", "Strong (0.50 / 0.02)"))
agg <- aggregate(cov_P_mean ~ n + K + sep, data = d, FUN = mean)
agg$K_lab <- factor(paste0("K = ", agg$K), levels = paste0("K = ", sort(unique(agg$K))))
LEG <- expression("Community separation (" * p['in'] * " / " * p[out] * ")")
p1 <- ggplot(agg, aes(factor(n), cov_P_mean, group = sep, color = sep,
                      linetype = sep, shape = sep)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "grey40") +
  geom_line(linewidth = 0.6) + geom_point(size = 2.2) +
  facet_wrap(~K_lab) +
  scale_color_brewer(palette = "Dark2", name = LEG) +
  scale_linetype_manual(values = c("solid", "longdash", "dotdash"), name = LEG) +
  scale_shape_manual(values = c(16, 17, 15), name = LEG) +
  scale_y_continuous(limits = c(0.4, 1.0)) +
  labs(x = "Number of nodes", y = "Empirical coverage") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")
ggsave("manuscript/figures/fig_coverage_by_config.pdf", p1, width = 7.5, height = 3.6)
ggsave("manuscript/figures/fig_coverage_by_config.png", p1, width = 7.5, height = 3.6, dpi = 200)
cs <- list.files(OUT, "^calib_task.*csv$", full.names = TRUE)
cal <- do.call(rbind, lapply(cs, read.csv))
cat("calib rows:", nrow(cal), "\n")
ca <- aggregate(cbind(n_pairs, n_true) ~ bin + bin_lo + bin_hi, data = cal, FUN = sum)
ca$mid <- (ca$bin_lo + ca$bin_hi) / 2
ca$emp <- ca$n_true / ca$n_pairs
p2 <- ggplot(ca, aes(mid, emp)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  geom_line(color = "#2166ac", linewidth = 0.6) +
  geom_point(aes(size = n_pairs), color = "#2166ac", alpha = 0.8) +
  scale_size_continuous(trans = "log10", range = c(1, 4)) +
  labs(x = "Bootstrap co-membership probability (bin midpoint)",
       y = "Empirical co-membership frequency", size = "Node pairs") +
  theme_bw(base_size = 11)
ggsave("manuscript/figures/fig_calibration.pdf", p2, width = 5.2, height = 4.4)
ggsave("manuscript/figures/fig_calibration.png", p2, width = 5.2, height = 4.4, dpi = 200)
qq <- do.call(rbind, lapply(split(d$width_P_mean, d$n), function(w) {
  data.frame(p = seq(0, 1, length.out = 401),
             w = quantile(w, seq(0, 1, length.out = 401), na.rm = TRUE))
}))
qq$n <- rep(sort(unique(d$n)), each = 401)
p3 <- ggplot(qq, aes(w, p, color = factor(n))) +
  geom_step(linewidth = 0.6) +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "grey40") +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Mean width of 95% co-membership intervals (per simulation)",
       y = "Cumulative share of simulations", color = "Nodes") +
  theme_bw(base_size = 11)
ggsave("manuscript/figures/fig_width_ecdf.pdf", p3, width = 6.2, height = 3.8)
ggsave("manuscript/figures/fig_width_ecdf.png", p3, width = 6.2, height = 3.8, dpi = 200)
cat("done 23\n")
