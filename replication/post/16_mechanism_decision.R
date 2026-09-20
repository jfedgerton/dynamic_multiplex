#!/usr/bin/env Rscript
# =============================================================================
# replication/post/16_mechanism_decision.R  --  mechanism tests (sim/03) and
# the main-text decision-tree table (Section "When to use which method")
#   manuscript/tables/tab_decision_tree.tex   MAIN TEXT. One row per simulated
#          scenario (Table 2 regimes, mechanism regimes, coupling regimes):
#          joint NMI of DynMux (Jaccard), multislice (adjacent links; same
#          links as DynMux where the generator supplies them) and
#          cross-sectional + Hungarian; best method and its margin over the
#          runner-up. The prose rule is written from this table.
#   manuscript/tables/tab_app_mechanism.tex   APPENDIX. joint NMI, layer NMI,
#          K MAE by mechanism regime x method; paired DynMux minus multislice
#          difference and the share of series DynMux wins
#   manuscript/tables/tab_app_mechanism_cells.tex   APPENDIX. the paired
#          difference by design cell (turnover x n, core x rotation, K x T)
# Inputs: output/mechanism/mech_cfg*.csv (sim/03), output/regime/dyn_cfg*.csv
#         (sim/01), output/coupling/coup_cfg*.csv (sim/04)
# Usage: DM_ROOT=... Rscript replication/post/16_mechanism_decision.R
# =============================================================================
ROOT <- Sys.getenv("DM_ROOT", unset = getwd())
TAB <- file.path(ROOT, "manuscript", "tables"); dir.create(TAB, showWarnings = FALSE, recursive = TRUE)
write_tex <- function(header, body, align, path) writeLines(c(sprintf("\\begin{tabular}{%s}", align), "\\toprule", header, "\\midrule", body, "\\bottomrule", "\\end{tabular}"), path)
rd <- function(dir, pat) { f <- list.files(file.path(ROOT, "output", dir), pat, full.names = TRUE); if (!length(f)) return(NULL); do.call(rbind, lapply(f, read.csv, stringsAsFactors = FALSE)) }
fmt <- function(x) ifelse(is.na(x), "--", sprintf("%.3f", x))

# ---- 1. mechanism tests (appendix) ---------------------------------------------------
m <- rd("mechanism", "^mech_cfg.*csv$")
stopifnot(!is.null(m), nrow(m) > 0)
cat("mechanism rows:", nrow(m), " configs:", length(unique(paste(m$regime, m$n, m$r, m$frac, m$core, m$p_rot, m$K, m$T))), "\n")
mreg  <- c(turnover = "Node turnover", coreperiphery = "Core/periphery rotation", longT = "Long series, small communities")
mmeth <- c("DynMux Jaccard", "Multislice adjacent", "Cross-sectional + Hungarian")
mlab  <- c("DynMux (Jaccard)", "Multislice, adjacent", "Cross-sectional + Hungarian")
stopifnot(all(mmeth %in% m$method))
m$key <- paste(m$regime, m$n, m$r, m$frac, m$core, m$p_rot, m$K, m$T, m$rep)
J <- m[m$method == "DynMux Jaccard", ]; S <- m[m$method == "Multislice adjacent", ]
S <- S[match(J$key, S$key), ]; stopifnot(identical(J$key, S$key))
J$diff <- J$nmi_joint - S$nmi_joint; J$win <- as.numeric(J$diff > 0)
agg <- aggregate(cbind(nmi_joint, nmi_layer, k_mae) ~ regime + method, m, mean)
pd  <- aggregate(cbind(diff, win) ~ regime, J, mean)
body <- unlist(lapply(names(mreg), function(rg) {
  rows <- vapply(seq_along(mmeth), function(i) { x <- agg[agg$regime == rg & agg$method == mmeth[i], ]
    sprintf("%s & %s & %s & %s & %s \\\\", if (i == 1) mreg[rg] else "", mlab[i], fmt(x$nmi_joint), fmt(x$nmi_layer), fmt(x$k_mae)) }, character(1))
  p <- pd[pd$regime == rg, ]
  c(rows, sprintf("& \\emph{DynMux $-$ multislice, paired} & %+.3f & \\multicolumn{2}{l}{DynMux more accurate in %.0f\\%% of series} \\\\ \\addlinespace", p$diff, 100 * p$win)) }))
write_tex("Regime & Method & Joint NMI & Layer NMI & $K$ MAE \\\\", body, "llccc", file.path(TAB, "tab_app_mechanism.tex"))

cell_rows <- c()
a <- aggregate(cbind(diff, win) ~ frac + n, J[J$regime == "turnover", ], mean)
cell_rows <- c(cell_rows, sprintf("Node turnover & turnover %.1f, $n = %d$ & %+.3f & %.0f\\%% \\\\", a$frac, a$n, a$diff, 100 * a$win))
a <- aggregate(cbind(diff, win) ~ core + p_rot, J[J$regime == "coreperiphery", ], mean)
cell_rows <- c(cell_rows, sprintf("Core/periphery & core %.1f, rotation %.1f & %+.3f & %.0f\\%% \\\\", a$core, a$p_rot, a$diff, 100 * a$win))
a <- aggregate(cbind(diff, win) ~ K + T, J[J$regime == "longT", ], mean)
cell_rows <- c(cell_rows, sprintf("Long series & $K = %d$, $T = %d$ & %+.3f & %.0f\\%% \\\\", a$K, a$T, a$diff, 100 * a$win))
write_tex("Regime & Cell & DynMux $-$ multislice & DynMux wins \\\\", cell_rows, "llcc", file.path(TAB, "tab_app_mechanism_cells.tex"))
cat("\n--- mechanism (joint NMI) ---\n"); print(reshape(agg[, c("regime", "method", "nmi_joint")], idvar = "regime", timevar = "method", direction = "wide"), row.names = FALSE, digits = 3)
cat("--- paired diff / win share ---\n"); print(pd, row.names = FALSE, digits = 3)

# ---- 2. decision-tree table (main text) ----------------------------------------------
rows <- list()
add <- function(group, scen, j, s, h, note = "") rows[[length(rows) + 1L]] <<- data.frame(group = group, scen = scen, j = j, s = s, h = h, note = note, stringsAsFactors = FALSE)

r1 <- rd("regime", "^dyn_cfg.*csv$")
if (!is.null(r1)) {
  a1 <- aggregate(nmi_joint ~ regime + method, r1, function(x) mean(x, na.rm = TRUE))
  g1 <- function(rg, me) { v <- a1$nmi_joint[a1$regime == rg & a1$method == me]; if (length(v)) v else NA_real_ }
  rn <- c(birthdeath = "Communities are born and die", regimeshift = "Abrupt rewiring at a change-point", churnswitch = "Gradual switching, nodes persist")
  for (rg in names(rn)) add("Table 2 regimes", rn[rg], g1(rg, "DynMux Jaccard"), g1(rg, "DynMux multislice (adjacent)"), g1(rg, "Cross-sectional + Hungarian"))
  add("Table 2 regimes", "Recurring structure, period known", g1("seasonality", "DynMux Jaccard"), g1("seasonality", "DynMux multislice (custom)"), g1("seasonality", "Cross-sectional + Hungarian"), "multislice with period links")
  add("Table 2 regimes", "Recurring structure, period unknown", g1("seasonality", "DynMux Jaccard"), g1("seasonality", "DynMux multislice (adjacent)"), g1("seasonality", "Cross-sectional + Hungarian"), "multislice with adjacent links")
} else cat("(no regime output: Table 2 rows omitted)\n")

gm <- function(rg, me) { v <- agg$nmi_joint[agg$regime == rg & agg$method == me]; if (length(v)) v else NA_real_ }
for (rg in names(mreg)) add("Mechanism tests", mreg[rg], gm(rg, "DynMux Jaccard"), gm(rg, "Multislice adjacent"), gm(rg, "Cross-sectional + Hungarian"))

r4 <- rd("coupling", "^coup_cfg.*csv$")
if (!is.null(r4)) {
  a4 <- aggregate(nmi_joint_break ~ regime + method, r4, function(x) mean(x, na.rm = TRUE))
  g4 <- function(rg, me) { v <- a4$nmi_joint_break[a4$regime == rg & a4$method == me]; if (length(v)) v else NA_real_ }
  cn <- c(balanced = "Stable communities, edge noise only", split = "A community splits", merge = "Two communities merge", sizeskew = "Unequal community sizes", shrink = "A community shrinks", nested = "Nested communities")
  for (rg in names(cn)) add("Lineage regimes", cn[rg], g4(rg, "DynMux Jaccard"), g4(rg, "DynMux multislice"), NA_real_)
} else cat("(no coupling output: lineage rows omitted)\n")

dt <- do.call(rbind, rows)
best_of <- function(j, s, h) { v <- c("DynMux" = j, "Multislice" = s, "Hungarian" = h); v <- v[!is.na(v)]
  if (!length(v)) return(c("--", "--")); v <- sort(v, decreasing = TRUE); c(names(v)[1], if (length(v) > 1) sprintf("%+.3f", v[1] - v[2]) else "") }
bo <- t(mapply(best_of, dt$j, dt$s, dt$h)); dt$best <- bo[, 1]; dt$margin <- bo[, 2]
body <- unlist(lapply(unique(dt$group), function(g) { d <- dt[dt$group == g, ]
  c(sprintf("\\multicolumn{6}{l}{\\emph{%s}} \\\\", g),
    sprintf("%s & %s & %s & %s & %s & %s \\\\", d$scen, fmt(d$j), fmt(d$s), fmt(d$h), d$best, d$margin), "\\addlinespace") }))
write_tex("Scenario & DynMux (Jaccard) & Multislice & Hungarian & Best & Margin \\\\", body, "lcccll", file.path(TAB, "tab_decision_tree.tex"))
cat("\n--- decision-tree table (joint NMI) ---\n"); print(dt[, c("group", "scen", "j", "s", "h", "best", "margin")], row.names = FALSE, digits = 3)
cat("wrote tab_decision_tree.tex, tab_app_mechanism.tex, tab_app_mechanism_cells.tex\ndone 16_mechanism_decision\n")
