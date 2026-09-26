# =============================================================================
# replication/compare_tables.R
# Compare regenerated LaTeX tables with the versions printed in the paper.
#
# Usage: Rscript replication/compare_tables.R <new_tables_dir> <reference_tables_dir> <report.md>
#   new_tables_dir        manuscript/tables written by the local rerun
#   reference_tables_dir  the tables/ folder from the Overleaf export (what the paper prints)
#
# For each reference table the script
#   1. replaces every number with '#' and compares the remaining text
#      (a difference here is STRUCTURAL: a changed row, header, or label);
#   2. compares the numbers cell by cell, in order, and reports the largest
#      absolute difference and every cell that differs by more than TOL.
# The runtime table records the machine's own timings and is expected to differ.
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3)
new_dir <- args[1]; ref_dir <- args[2]; out_md <- args[3]
stopifnot(dir.exists(new_dir), dir.exists(ref_dir))
TOL <- 0.005
EXPECTED_DIFFERENT <- c("tab_app_runtime.tex")

norm <- function(x) {
  x <- gsub("\\$-\\$", "-", x)          # LaTeX minus
  x <- gsub("\\$\\+\\$", "+", x)
  x <- gsub("\\\\cellcolor\\{[^}]*\\}", "", x)
  gsub("[ \t]+", " ", x)
}
num_re <- "-?[0-9]+(\\.[0-9]+)?"

ref_files <- sort(list.files(ref_dir, pattern = "\\.tex$"))
rows <- list(); detail <- character()
for (f in ref_files) {
  rf <- file.path(ref_dir, f); nf <- file.path(new_dir, f)
  if (!file.exists(nf)) { rows[[f]] <- data.frame(table = f, status = "NOT REGENERATED", n_cells = NA, n_diff = NA, max_abs_diff = NA); next }
  r <- norm(readLines(rf, warn = FALSE)); n <- norm(readLines(nf, warn = FALSE))
  r <- r[nzchar(trimws(r))]; n <- n[nzchar(trimws(n))]
  skel_r <- gsub(num_re, "#", r); skel_n <- gsub(num_re, "#", n)
  nr <- regmatches(r, gregexpr(num_re, r)); nn <- regmatches(n, gregexpr(num_re, n))
  vr <- as.numeric(unlist(nr)); vn <- as.numeric(unlist(nn))
  structural <- !identical(skel_r, skel_n) || length(vr) != length(vn)
  if (structural) {
    status <- "STRUCTURAL"
    dl <- setdiff(union(skel_r, skel_n), intersect(skel_r, skel_n))
    detail <- c(detail, sprintf("\n### %s: structural difference\n", f),
                paste0("    ", head(dl, 12)))
    rows[[f]] <- data.frame(table = f, status = status, n_cells = length(vr), n_diff = NA, max_abs_diff = NA); next
  }
  d <- abs(vr - vn); big <- which(d > TOL + 1e-9)
  status <- if (f %in% EXPECTED_DIFFERENT) "EXPECTED DIFFERENT (timings)" else if (length(big)) "NUMERIC DIFF" else if (any(d > 0)) "ROUNDING ONLY" else "IDENTICAL"
  rows[[f]] <- data.frame(table = f, status = status, n_cells = length(vr), n_diff = length(big), max_abs_diff = if (length(d)) max(d) else 0)
  if (length(big) && !(f %in% EXPECTED_DIFFERENT)) {
    # label each differing cell with its line (first 60 characters)
    line_of <- rep(seq_along(nr), lengths(nr))
    detail <- c(detail, sprintf("\n### %s: %d cells differ by more than %.3f\n", f, length(big), TOL),
                sprintf("    line %d | paper %s | rerun %s | %s", line_of[big], vr[big], vn[big],
                        substr(trimws(r[line_of[big]]), 1, 60)))
  }
}
res <- do.call(rbind, rows)
extra <- setdiff(list.files(new_dir, pattern = "\\.tex$"), ref_files)
con <- file(out_md, "w")
writeLines(c("# Table comparison: local rerun vs paper", "",
             sprintf("New tables: `%s`  ", normalizePath(new_dir)),
             sprintf("Reference: `%s`  ", normalizePath(ref_dir)),
             sprintf("Tolerance: %.3f (two-decimal rounding)", TOL), "",
             "| Table | Status | Numbers | Cells > tol | Max abs diff |", "|---|---|---|---|---|",
             sprintf("| %s | %s | %s | %s | %s |", res$table, res$status, res$n_cells,
                     ifelse(is.na(res$n_diff), "", res$n_diff),
                     ifelse(is.na(res$max_abs_diff), "", formatC(res$max_abs_diff, format = "g", digits = 3))),
             "", if (length(extra)) c("Generated but not in the reference folder (not used by the paper):", paste0("- ", extra)) else "",
             "", "## Details", detail), con)
close(con)
print(table(res$status))
cat("report:", out_md, "\n")
