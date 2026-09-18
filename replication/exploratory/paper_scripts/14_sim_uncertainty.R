#!/usr/bin/env Rscript
# ============================================================================
# 14_sim_uncertainty.R -- Uncertainty on the regime comparison (nmi_joint).
# Reads replication/extended/output/dynamic/dyn_cfg*.csv (72 configs x 30 reps,
# 16 methods each; one regime per file). Produces:
#   manuscript/tables/tab_regime_ci.csv       per-method mean nmi_joint + 95% CI, by regime
#   manuscript/tables/tab_regime_paired.csv   paired Wilcoxon DynMux vs baselines
#   manuscript/tables/tab_regime_paired.tex   booktabs version
#   manuscript/figures/fig_regime_ci.pdf      point + 95% CI, faceted by regime
# ============================================================================
set.seed(123)
suppressMessages({ library(ggplot2) })
D   <- "replication/extended/output/dynamic"
OUT <- "manuscript/tables"; FIG <- "manuscript/figures"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

files <- list.files(D, "dyn_cfg.*csv$", full.names = TRUE)
d <- do.call(rbind, lapply(files, function(f) {
  x <- read.csv(f, stringsAsFactors = FALSE); x$cfg <- basename(f); x
}))
cat("rows:", nrow(d), " configs:", length(files), "\n")
cat("methods:\n"); print(sort(unique(d$method)))
cat("regimes:", paste(unique(d$regime), collapse = ", "), "\n")
d$unit <- paste(d$cfg, d$rep, sep = "_")   # one simulated network instance

rn <- c(birthdeath = "Births and Deaths", churnswitch = "Churn and Switching",
        regimeshift = "Regime Shift", seasonality = "Seasonality")

## (A) per-method mean nmi_joint + 95% normal-approx CI, by regime -------------
ci <- do.call(rbind, lapply(split(d, list(d$regime, d$method), drop = TRUE), function(g) {
  x <- g$nmi_joint[is.finite(g$nmi_joint)]; n <- length(x); se <- if (n > 1) sd(x)/sqrt(n) else NA
  data.frame(regime = g$regime[1], method = g$method[1], n = n,
             mean = mean(x), lo = mean(x) - 1.96*se, hi = mean(x) + 1.96*se,
             stringsAsFactors = FALSE)
}))
write.csv(ci, file.path(OUT, "tab_regime_ci.csv"), row.names = FALSE)

## (B) paired Wilcoxon: main DynMux Leiden specs vs baselines -----------------
um <- unique(d$method)
refs <- intersect(c("DynMux Overlap", "DynMux Jaccard"), um)
if (!length(refs)) refs <- grep("^DynMux (Overlap|Jaccard)$", um, value = TRUE)
baselines <- intersect(c("Pooled Leiden", "Cross-sectional + Hungarian",
                         "multinet GLouvain", "Dynamic SBM"), um)
cat("refs:", paste(refs, collapse = " | "), "\n")
cat("baselines:", paste(baselines, collapse = " | "), "\n")

paired <- list()
for (rg in unique(d$regime)) {
  sub <- d[d$regime == rg, ]
  w <- reshape(sub[, c("unit","method","nmi_joint")], idvar = "unit",
               timevar = "method", direction = "wide")
  names(w) <- sub("nmi_joint.", "", names(w), fixed = TRUE)
  for (rf in refs) for (bl in baselines) {
    if (!(rf %in% names(w) && bl %in% names(w))) next
    a <- w[[rf]]; b <- w[[bl]]; ok <- is.finite(a) & is.finite(b)
    if (sum(ok) < 3) next
    t <- suppressWarnings(wilcox.test(a[ok], b[ok], paired = TRUE, exact = FALSE))
    paired[[length(paired)+1]] <- data.frame(regime = rg, dynmux = rf, baseline = bl,
      n_pairs = sum(ok), diff = mean(a[ok]-b[ok]), p = t$p.value, stringsAsFactors = FALSE)
  }
}
paired <- do.call(rbind, paired)
paired$p_bh <- ave(paired$p, paired$regime, FUN = function(p) p.adjust(p, "BH"))
paired <- paired[order(paired$regime, -paired$diff), ]
write.csv(paired, file.path(OUT, "tab_regime_paired.csv"), row.names = FALSE)

## (C) LaTeX booktabs for the paired table -----------------------------------
esc <- function(s) gsub("([&_])", "\\\\\\1", s)
con <- c("\\begin{table}[htbp]", "\\centering",
  "\\caption{Paired Wilcoxon signed-rank tests of DynMux (Leiden) against baselines on joint NMI, matched by simulated network within each regime. $p$-values Benjamini--Hochberg adjusted within regime.}",
  "\\label{tab:regime_paired}", "\\begin{tabular}{lllrrr}", "\\toprule",
  "Regime & DynMux & Baseline & $n$ & $\\Delta$NMI & $p_{\\text{BH}}$ \\\\", "\\midrule")
for (i in seq_len(nrow(paired))) { r <- paired[i, ]
  con <- c(con, sprintf("%s & %s & %s & %d & %.3f & %s \\\\",
    esc(unname(rn[r$regime])), esc(sub("DynMux ", "", r$dynmux)), esc(r$baseline),
    r$n_pairs, r$diff, formatC(r$p_bh, format = "e", digits = 1))) }
con <- c(con, "\\bottomrule", "\\end{tabular}", "\\end{table}")
writeLines(con, file.path(OUT, "tab_regime_paired.tex"))

## (D) error-bar figure: mean nmi_joint +/- 95% CI, faceted by regime --------
ci$regime_lab <- factor(unname(rn[ci$regime]), levels = unname(rn))
ci$is_dyn <- grepl("DynMux", ci$method)
lev <- c("DynMux Jaccard","DynMux Overlap","DynMux weighted Jaccard","DynMux weighted Overlap",
         "DynMux multislice (adjacent)","DynMux multislice (custom)",
         "DynMux Jaccard (Louvain)","DynMux Overlap (Louvain)","DynMux weighted Jaccard (Louvain)",
         "DynMux weighted Overlap (Louvain)","DynMux multislice (adjacent, Louvain)",
         "Cross-sectional + Hungarian","Pooled Leiden","Pooled Louvain","multinet GLouvain","Dynamic SBM")
u <- unique(ci$method); lev <- c(intersect(lev, u), setdiff(u, lev))
ci$method <- factor(ci$method, levels = rev(lev))
p <- ggplot(ci, aes(method, mean, colour = is_dyn)) +
  geom_point(size = 1.6) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.3) +
  facet_wrap(~ regime_lab, ncol = 2) + coord_flip() +
  scale_colour_manual(values = c(`TRUE` = "#1B9E77", `FALSE` = "grey55"), guide = "none") +
  labs(x = NULL, y = "Mean joint NMI (95% CI across 30 reps)") +
  theme_bw(base_size = 11)
ggsave(file.path(FIG, "fig_regime_ci.pdf"), p, width = 11, height = 9)
cat("DONE sim uncertainty\n")
