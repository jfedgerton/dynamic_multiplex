#!/usr/bin/env Rscript
# =============================================================================
# replication/post/15_omega_selection.R  --  appendix tables from sim/05 and
# sim/06, plus the empirical stability table from empirical/10
#   manuscript/tables/tab_app_omega_sweep.tex      joint NMI by regime for
#          multislice at omega in {0.25, 0.5, 1, 2, 4}, adjacent and same links,
#          with DynMux Jaccard on the same series
#   manuscript/figures/fig_omega_sweep.pdf
#   manuscript/tables/tab_app_selection_ci.tex     paired differences (rule vs always-DynMux / always-multislice / oracle) with 95% CIs
#   manuscript/tables/tab_app_selection_rule.tex   by regime: share of series
#          in which the method with the higher stability is the more accurate,
#          accuracy of the rule vs always-DynMux vs always-multislice
#   manuscript/tables/tab_app_empirical_stability.tex  s, floor, decided share
#          for the four networks x two methods
# Usage: DM_ROOT=... Rscript replication/post/15_omega_selection.R
# =============================================================================
ROOT <- Sys.getenv("DM_ROOT", unset = getwd())
suppressPackageStartupMessages(library(ggplot2))
TAB <- file.path(ROOT, "manuscript", "tables");  dir.create(TAB, showWarnings = FALSE, recursive = TRUE)
FIG <- file.path(ROOT, "manuscript", "figures"); dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
write_tex <- function(header, body, align, path) writeLines(c(sprintf("\\begin{tabular}{%s}", align), "\\toprule", header, "\\midrule", body, "\\bottomrule", "\\end{tabular}"), path)
rn <- c(birthdeath = "Births \\& deaths", churnswitch = "Gradual switching", regimeshift = "Abrupt rewiring", seasonality = "Recurring structure")
regs <- names(rn)
rd <- function(dir, pat) { f <- list.files(file.path(ROOT, "output", dir), pat, full.names = TRUE); if (!length(f)) return(NULL); do.call(rbind, lapply(f, read.csv, stringsAsFactors = FALSE)) }

# ---- omega sweep -----------------------------------------------------------------
om <- rd("omega", "^omega_cfg.*csv$")
if (!is.null(om)) {
  cat("omega rows:", nrow(om), " configs:", length(unique(paste(om$regime, om$n, om$r, om$intensity))), "\n")
  a <- aggregate(nmi_joint ~ variant + regime, om, mean)
  vars <- c("DynMux Jaccard", sprintf("Multislice adjacent omega=%g", c(0.25, 0.5, 1, 2, 4)), sprintf("Multislice custom links omega=%g", c(0.25, 0.5, 1, 2, 4)))
  lab <- c("DynMux (Jaccard)", sprintf("Multislice, adjacent, $\\omega = %g$", c(0.25, 0.5, 1, 2, 4)), sprintf("Multislice, same links, $\\omega = %g$", c(0.25, 0.5, 1, 2, 4)))
  body <- vapply(seq_along(vars), function(i) { x <- a[a$variant == vars[i], ]
    sprintf("%s & %s \\\\", lab[i], paste(sprintf("%.3f", x$nmi_joint[match(regs, x$regime)]), collapse = " & ")) }, character(1))
  write_tex(paste("Method &", paste(rn, collapse = " & "), "\\\\"), body, "lcccc", file.path(TAB, "tab_app_omega_sweep.tex"))
  os <- om[grepl("^Multislice", om$variant), ]
  os$omega <- as.numeric(sub(".*omega=", "", os$variant)); os$links <- ifelse(grepl("adjacent", os$variant), "adjacent links", "same links as DynMux")
  oa <- aggregate(nmi_joint ~ omega + links + regime, os, mean); oa$regime_lab <- factor(rn[oa$regime], levels = rn)
  ref <- aggregate(nmi_joint ~ regime, om[om$variant == "DynMux Jaccard", ], mean); ref$regime_lab <- factor(rn[ref$regime], levels = rn)
  g <- ggplot(oa, aes(omega, nmi_joint, colour = links)) + geom_line() + geom_point() +
    geom_hline(data = ref, aes(yintercept = nmi_joint), linetype = 2, colour = "grey40") +
    scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4)) + facet_wrap(~ regime_lab, ncol = 2) +
    labs(x = expression(omega), y = "Joint NMI", colour = "Multislice links") + theme_bw(base_size = 10) + theme(legend.position = "bottom")
  ggsave(file.path(FIG, "fig_omega_sweep.pdf"), g, width = 6.5, height = 5.5)
  cat("wrote tab_app_omega_sweep.tex, fig_omega_sweep.pdf\n")
} else cat("(no omega output)\n")

# ---- selection rule ---------------------------------------------------------------
se <- rd("selection", "^sel_cfg.*csv$")
if (!is.null(se)) {
  se$key <- paste(se$regime, se$n, se$r, se$intensity, se$rep)
  D <- se[se$method == "DynMux Jaccard", ]; M <- se[se$method == "Multislice (custom links)", ]
  D <- D[order(D$key), ]; M <- M[M$key %in% D$key, ]; M <- M[order(M$key), ]; D <- D[D$key %in% M$key, ]
  stopifnot(identical(D$key, M$key))
  pick_D <- D$stability >= M$stability
  rule_acc <- ifelse(pick_D, D$acc_joint, M$acc_joint); best <- pmax(D$acc_joint, M$acc_joint)
  correct <- (pick_D & D$acc_joint >= M$acc_joint) | (!pick_D & M$acc_joint >= D$acc_joint)
  rows <- lapply(c(regs, "all"), function(rg) { i <- if (rg == "all") rep(TRUE, nrow(D)) else D$regime == rg
    data.frame(regime = rg, n = sum(i), correct = mean(correct[i]), pickD = mean(pick_D[i]), rule = mean(rule_acc[i]),
               dyn = mean(D$acc_joint[i]), ms = mean(M$acc_joint[i]), oracle = mean(best[i])) })
  sr <- do.call(rbind, rows)
  write_tex("Regime & Series & Rule picks the more accurate method & Rule picks DynMux & Accuracy: rule & always DynMux & always multislice & oracle \\\\",
            sprintf("%s & %d & %.3f & %.2f & %.3f & %.3f & %.3f & %.3f \\\\", c(rn, all = "All regimes")[sr$regime], sr$n, sr$correct, sr$pickD, sr$rule, sr$dyn, sr$ms, sr$oracle),
            "lccccccc", file.path(TAB, "tab_app_selection_rule.tex"))
  cat("\n--- selection rule ---\n"); print(sr, row.names = FALSE, digits = 3)
  # paired differences: rule minus always-DynMux, rule minus always-multislice,
  # oracle minus rule; mean over series with 1.96 * sd / sqrt(n) intervals
  pd <- function(x) c(est = mean(x), lo = mean(x) - 1.96 * sd(x) / sqrt(length(x)), hi = mean(x) + 1.96 * sd(x) / sqrt(length(x)))
  ci <- do.call(rbind, lapply(c(regs, "all"), function(rg) { i <- if (rg == "all") rep(TRUE, nrow(D)) else D$regime == rg
    rbind(data.frame(regime = rg, comp = "Rule $-$ always multislice", t(pd(rule_acc[i] - M$acc_joint[i]))),
          data.frame(regime = rg, comp = "Rule $-$ always DynMux",     t(pd(rule_acc[i] - D$acc_joint[i]))),
          data.frame(regime = rg, comp = "Oracle $-$ rule",            t(pd(best[i] - rule_acc[i])))) }))
  stopifnot(nrow(ci) == 15L, sum(D$regime %in% regs) == nrow(D))
  fmt <- function(v) gsub("-", "$-$", sprintf("%.3f", v))
  write_tex("Regime & Comparison & Mean difference [95\\% CI] \\\\",
            sprintf("%s & %s & %s [%s, %s] \\\\", ifelse(duplicated(ci$regime), "", c(rn, all = "All regimes")[ci$regime]), ci$comp, fmt(ci$est), fmt(ci$lo), fmt(ci$hi)),
            "llc", file.path(TAB, "tab_app_selection_ci.tex"))
  cat("\n--- selection rule, paired differences ---\n"); print(ci, row.names = FALSE, digits = 3)
} else cat("(no selection output)\n")

# ---- empirical stability ---------------------------------------------------------
es <- rd("empirical", "_stability\\.csv$")
if (!is.null(es)) {
  net_lab <- c(atop = "Alliances (ATOP)", dca = "Defense cooperation", igo = "IGO co-membership", trade = "Trade")
  es <- es[order(match(es$net, names(net_lab)), es$method), ]
  write_tex("Network & Method & Stability $s$ & Floor & Decided pairs & Nodes with $s_i \\geq 0.9$ & Mean $K$ \\\\",
            sprintf("%s & %s & %.3f & %.2f & %.0f\\%% & %.0f\\%% & %.1f \\\\", ifelse(duplicated(es$net), "", net_lab[es$net]), es$method, es$stability, es$floor,
                    100 * (1 - es$pairs_undetermined), 100 * es$node_share_ge_0.9, es$mean_K),
            "llccccc", file.path(TAB, "tab_app_empirical_stability.tex"))
  cat("\n--- empirical stability ---\n"); print(es[, c("net", "method", "stability", "floor", "pairs_undetermined", "node_share_ge_0.9")], row.names = FALSE, digits = 3)
} else cat("(no empirical stability output)\n")
cat("done 15_omega_selection\n")
