#!/usr/bin/env Rscript
# =============================================================================
# replication/post/18_mechanism_paired.R
# Appendix material from the mechanism tests (sim/03): the pairwise companion
# to post/11 for the second set of simulations.
#
#   fig_mech_paired_nmi.pdf     paired Delta joint NMI, 3 conditions x 2 levels
#   fig_mech_paired_kmae.pdf    paired Delta K MAE,     3 conditions x 2 levels
#   tab_app_mech_paired_ci.tex  paired 95% CIs, both metrics (booktabs tabular)
#
#   Delta Joint NMI = DynMux - baseline   (positive favors DynMux)
#   Delta K MAE     = baseline - DynMux   (positive favors DynMux; lower is better)
#
# Baselines: multislice (adjacent identity links) and cross-sectional + Hungarian,
# the two methods run in sim/03. Pairing unit = one simulated network (config x
# rep); every method saw the same network, so the difference is within-network.
# The "level" split mirrors the low/high intensity split of post/11 and is the
# condition-specific quantity that sets how hard the condition is:
#   node turnover   presence toggle 0.5 (low) / 0.7 (high)
#   core stability  periphery rotation probability 0.3 (low) / 0.5 (high)
#   long series     layers T = 30 (low) / 60 (high)
#
# Usage:  DM_ROOT=/path/to/dynamic_multiplex Rscript replication/post/18_mechanism_paired.R
# =============================================================================
set.seed(123)
suppressMessages({ library(ggplot2) })

ROOT <- Sys.getenv("DM_ROOT", unset = getwd())
MECH <- file.path(ROOT, "output", "mechanism")
TAB  <- file.path(ROOT, "manuscript", "tables")
FIG  <- file.path(ROOT, "manuscript", "figures")
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

write_tex <- function(header, body, align, path) {
  stopifnot(length(body) > 0, nchar(align) > 0)
  writeLines(c(sprintf("\\begin{tabular}{%s}", align),
               "\\toprule", header, "\\midrule",
               body,
               "\\bottomrule", "\\end{tabular}"), path)
  cat("wrote", basename(path), "(", length(body), "rows )\n")
}

# --- read ----------------------------------------------------------------
fs <- list.files(MECH, "^mech_cfg.*csv$", full.names = TRUE)
if (!length(fs)) stop("No mechanism output in ", MECH, " -- run sim/03 first.", call. = FALSE)
d <- do.call(rbind, lapply(fs, read.csv, stringsAsFactors = FALSE))
stopifnot(all(c("regime", "n", "r", "frac", "core", "p_rot", "K", "T", "rep", "method",
                "nmi_joint", "k_mae") %in% names(d)))
cat("files:", length(fs), " rows:", nrow(d), "\n")
stopifnot(length(fs) == 50L, nrow(d) == 50L * 10L * 3L)

ref       <- "DynMux Jaccard"
baselines <- c("Multislice adjacent", "Cross-sectional + Hungarian")
stopifnot(all(c(ref, baselines) %in% unique(d$method)))
d$unit <- paste(d$regime, d$n, d$r, d$frac, d$core, d$p_rot, d$K, d$T, d$rep, sep = "|")

cn <- c(turnover = "Node Turnover", coreperiphery = "Core Stability", longT = "Long Series")
cn_tex <- c(turnover = "Node turnover", coreperiphery = "Core stability", longT = "Long series")
blab <- c("Multislice adjacent" = "Multislice\nadjacent",
          "Cross-sectional + Hungarian" = "Hungarian\nmatching")
bl_tex <- c("Multislice adjacent" = "Multislice (adjacent links)",
            "Cross-sectional + Hungarian" = "Hungarian matching")
stopifnot(all(d$regime %in% names(cn)))

# condition-specific low / high level (see header)
d$level <- NA_character_
d$level[d$regime == "turnover"]      <- ifelse(d$frac[d$regime == "turnover"] >= 0.7, "high", "low")
d$level[d$regime == "coreperiphery"] <- ifelse(d$p_rot[d$regime == "coreperiphery"] >= 0.5, "high", "low")
d$level[d$regime == "longT"]         <- ifelse(d$T[d$regime == "longT"] >= 60, "high", "low")
stopifnot(!anyNA(d$level))
stopifnot(all(sort(unique(d$frac[d$regime == "turnover"])) == c(0.5, 0.7)),
          all(sort(unique(d$p_rot[d$regime == "coreperiphery"])) == c(0.3, 0.5)),
          all(sort(unique(d$T[d$regime == "longT"])) == c(30, 60)))

# --- paired differences within condition x level -------------------------
rows <- list()
for (rg in names(cn)) for (lv in c("low", "high")) for (mt in c("nmi_joint", "k_mae")) {
  sub <- d[d$regime == rg & d$level == lv, c("unit", "method", mt)]
  stopifnot(nrow(sub) > 0)
  w <- reshape(sub, idvar = "unit", timevar = "method", direction = "wide")
  names(w) <- sub(paste0(mt, "."), "", names(w), fixed = TRUE)
  for (bl in baselines) {
    stopifnot(all(c(ref, bl) %in% names(w)))
    dd <- if (mt == "nmi_joint") w[[ref]] - w[[bl]] else w[[bl]] - w[[ref]]
    dd <- dd[is.finite(dd)]
    stopifnot(length(dd) > 1)
    n <- length(dd); se <- sd(dd) / sqrt(n)
    rows[[length(rows) + 1]] <- data.frame(
      regime = rg, level = lv, metric = mt, baseline = bl,
      n = n, diff = mean(dd), lo = mean(dd) - 1.96 * se, hi = mean(dd) + 1.96 * se,
      win = mean(dd > 0), stringsAsFactors = FALSE)
  }
}
r <- do.call(rbind, rows)
stopifnot(nrow(r) == 3 * 2 * 2 * length(baselines))   # condition x level x metric x baseline
# n per cell: turnover 9 configs x 10 reps = 90 per level; core 8 x 10 = 80; longT 8 x 10 = 80
stopifnot(all(r$n[r$regime == "turnover"] == 90), all(r$n[r$regime != "turnover"] == 80))
cat("paired comparisons:", nrow(r), " n per comparison:", paste(unique(r$n), collapse = ","), "\n")
cat("CIs excluding zero:", sum(r$lo > 0 | r$hi < 0), "of", nrow(r), "\n")

# --- figures: one per metric ------------------------------------------------
r$cond_lab  <- factor(unname(cn[r$regime]), levels = unname(cn))
r$level_lab <- factor(ifelse(r$level == "low", "Low intensity", "High intensity"),
                      levels = c("Low intensity", "High intensity"))
r$baseline_lab <- factor(unname(blab[r$baseline]), levels = rev(unname(blab[baselines])))

mk <- function(mt, xlab) {
  s <- r[r$metric == mt, ]
  ggplot(s, aes(diff, baseline_lab)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
    geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.25, colour = "#1b9e77") +
    geom_point(size = 1.8, colour = "#1b9e77") +
    facet_grid(cond_lab ~ level_lab) +
    labs(x = xlab, y = NULL) +
    theme_bw(base_size = 9) +
    theme(legend.position  = "none",
          panel.grid.minor = element_blank(),
          strip.text       = element_text(face = "bold", size = 8.5),
          axis.text.y      = element_text(size = 8))
}
p1 <- mk("nmi_joint", expression(paste(Delta, " Joint NMI (positive favors DynMux), 95% CI")))
p2 <- mk("k_mae",     expression(paste(Delta, " ", italic(K), " MAE (positive favors DynMux), 95% CI")))
ggsave(file.path(FIG, "fig_mech_paired_nmi.pdf"),  p1, width = 6.5, height = 5.0)
ggsave(file.path(FIG, "fig_mech_paired_nmi.png"),  p1, width = 6.5, height = 5.0, dpi = 300)
ggsave(file.path(FIG, "fig_mech_paired_kmae.pdf"), p2, width = 6.5, height = 5.0)
ggsave(file.path(FIG, "fig_mech_paired_kmae.png"), p2, width = 6.5, height = 5.0, dpi = 300)
cat("wrote fig_mech_paired_nmi / _kmae (.pdf/.png)\n")

# --- CI table: both metrics side by side -----------------------------------
fmt3 <- function(x) sub("-", "$-$", sprintf("%.3f", round(x, 3) + 0), fixed = TRUE)  # +0 kills "-0.000"
r$cell <- sprintf("%s [%s, %s]", fmt3(r$diff), fmt3(r$lo), fmt3(r$hi))
w <- reshape(r[, c("regime", "level", "baseline", "metric", "cell")],
             idvar = c("regime", "level", "baseline"), timevar = "metric", direction = "wide")
names(w) <- sub("cell.", "", names(w), fixed = TRUE)
stopifnot(all(c("nmi_joint", "k_mae") %in% names(w)), nrow(w) == 3 * 2 * length(baselines))
w <- w[order(factor(w$regime,   levels = names(cn)),
             factor(w$level,    levels = c("low", "high")),
             factor(w$baseline, levels = baselines)), ]
grp  <- sprintf("%s (%s)", unname(cn_tex[w$regime]), ifelse(w$level == "low", "Low", "High"))
show <- c(TRUE, grp[-1] != grp[-length(grp)])
write_tex(header = "Condition (intensity) & Baseline & Joint NMI & $K$ MAE \\\\",
          body   = sprintf("%s & %s & %s & %s \\\\", ifelse(show, grp, ""), unname(bl_tex[w$baseline]),
                           w$nmi_joint, w$k_mae),
          align  = "llcc",
          path   = file.path(TAB, "tab_app_mech_paired_ci.tex"))
print(r[, c("regime", "level", "metric", "baseline", "n", "diff", "lo", "hi", "win")], digits = 3, row.names = FALSE)
cat("done 18_mechanism_paired\n")
