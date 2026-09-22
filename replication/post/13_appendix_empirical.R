#!/usr/bin/env Rscript
# =============================================================================
# replication/post/13_appendix_empirical.R
# Appendix tables for the empirical application: set-level recovery of the
# Braumoeller-coded international orders, from empirical/09_score_orders.R.
#
# Input:  $DM_ROOT/output/empirical/order_recovery.csv
#         long: net, order, kind, year, method, J, prec, rec, nB
# Output (bare booktabs tabulars, no caption/label):
#   tables/tab_order_recovery_summary.tex   net x method: mean J, precision, recall
#   tables/tab_order_recovery_<net>.tex     one per network; rows = orders,
#                                           columns = methods, cell = prec / rec
#                                           (means over the order's active years)
# Requires \usepackage{booktabs}.
#
# Usage:  DM_ROOT=/path/to/dynamic_multiplex Rscript replication/post/13_appendix_empirical.R
# =============================================================================
set.seed(123)

ROOT <- Sys.getenv("DM_ROOT", unset = getwd())
EMP  <- file.path(ROOT, "output", "empirical")
TAB  <- file.path(ROOT, "manuscript", "tables")
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)

write_tex <- function(header, body, align, path) {
  stopifnot(length(body) > 0, nchar(align) > 0)
  writeLines(c(sprintf("\\begin{tabular}{%s}", align),
               "\\toprule", header, "\\midrule",
               body,
               "\\bottomrule", "\\end{tabular}"), path)
  cat("wrote", basename(path), "(", length(body), "rows )\n")
}

orf <- file.path(EMP, "order_recovery.csv")
if (!file.exists(orf)) stop("Missing ", orf, " -- run empirical/09 first.", call. = FALSE)
o <- read.csv(orf, stringsAsFactors = FALSE)
stopifnot(all(c("net", "order", "kind", "year", "method", "J", "prec", "rec", "nB", "K") %in% names(o)))
cat("rows:", nrow(o), "\n")

net_lab  <- c(atop = "Alliances (ATOP)", dca = "Defense cooperation (DCA)",
              igo = "IGO co-membership", trade = "Trade")
meth_ord <- c("Jaccard", "Overlap", "multislice", "multinet", "Hungarian", "Pooled")
meth_lab <- c(Jaccard = "DynMux Jaccard", Overlap = "DynMux Overlap",
              multislice = "Multislice adj.", multinet = "Multislice full",
              Hungarian = "Hungarian", Pooled = "Pooled Leiden")
order_ord <- c("Concert GPs", "Concert Europe", "Interim Europe", "Bismarck Europe",
               "Wilhelm Europe", "Indep. after WWI", "Mandates", "League",
               "PW Liberal", "PW Communist", "Warsaw Pact", "PW Other", "PostCW Liberal")
stopifnot(all(o$net %in% names(net_lab)), all(o$method %in% meth_ord),
          all(o$order %in% order_ord))

# --- summary: net x method ------------------------------------------------
s <- aggregate(cbind(J, prec, rec, K) ~ net + method, data = o, FUN = mean)
s <- s[order(factor(s$net, levels = names(net_lab)), factor(s$method, levels = meth_ord)), ]
print(s)
grp  <- unname(net_lab[s$net])
show <- c(TRUE, grp[-1] != grp[-length(grp)])
write_tex(
  header = "Network & Method & Jaccard $J$ & Precision & Recall & $K$ \\\\",
  body   = sprintf("%s & %s & %.3f & %.3f & %.3f & %.1f \\\\",
                   ifelse(show, grp, ""), unname(meth_lab[s$method]), s$J, s$prec, s$rec, s$K),
  align  = "llcccc",
  path   = file.path(TAB, "tab_order_recovery_summary.tex"))

# --- per network: orders x methods, cell = precision / recall -------------
om <- aggregate(cbind(prec, rec, nB) ~ net + order + method, data = o, FUN = mean)
yrs <- aggregate(cbind(years = year) ~ net + order, data = o, FUN = function(y) length(unique(y)))
for (nt in names(net_lab)) {
  x <- om[om$net == nt, ]
  if (!nrow(x)) { cat("skip:", nt, "(no rows)\n"); next }
  ords <- order_ord[order_ord %in% unique(x$order)]
  body <- vapply(ords, function(od) {
    cells <- vapply(meth_ord, function(m) {
      r <- x[x$order == od & x$method == m, ]
      if (!nrow(r)) "---" else sprintf("%.2f / %.2f", r$prec, r$rec)
    }, character(1))
    ny <- yrs$years[yrs$net == nt & yrs$order == od]
    sprintf("%s & %d & %s \\\\", od, ny, paste(cells, collapse = " & "))
  }, character(1))
  write_tex(
    header = sprintf("Order & Years & %s \\\\", paste(unname(meth_lab[meth_ord]), collapse = " & ")),
    body   = unname(body),
    align  = paste0("lr", strrep("c", length(meth_ord))),
    path   = file.path(TAB, sprintf("tab_order_recovery_%s.tex", nt)))
}

cat("done 13_appendix_empirical\n")
