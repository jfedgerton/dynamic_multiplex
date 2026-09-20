# =============================================================================
# replication/tests/crosscheck_multislice_r_vs_python.R
# Cross-check of the two multislice implementations on the same networks:
#   R      dynamicmultiplex::fit_multilayer_identity_ties (generalized Louvain
#          with the Mucha per-slice null model, written for this package)
#   Python dynamic_multiplex.fit_multilayer_identity_ties (leidenalg multiplex
#          optimiser, the reference implementation)
# Writes output/tests/crosscheck_multislice.csv with, per network, the NMI
# between the two meta-partitions and each one's NMI to the planted truth.
# Both should agree closely on clear structure; on weak structure both are
# local optima of the same objective and may differ. Requires python with
# dynamic_multiplex, leidenalg and python-igraph installed.
# Usage: DM_ROOT=. Rscript replication/tests/crosscheck_multislice_r_vs_python.R
# =============================================================================
ROOT <- Sys.getenv("DM_ROOT", unset = getwd())
suppressMessages(pkgload::load_all(file.path(ROOT, "r_code"), quiet = TRUE)); suppressMessages(library(igraph))
outdir <- file.path(ROOT, "output", "tests"); dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
tmp <- file.path(outdir, "crosscheck_tmp"); dir.create(tmp, showWarnings = FALSE)
set.seed(123)
bl <- function(m, pin, pout) { n <- length(m); Pr <- ifelse(outer(m, m, "=="), pin, pout); diag(Pr) <- 0
  A <- matrix(0, n, n); up <- upper.tri(Pr); A[up] <- rbinom(sum(up), 1, Pr[up]); A + t(A) }
rows <- list(); k <- 0
for (n in c(100, 200)) for (sep in c("strong", "default", "weak")) for (sw in c(0, 0.1)) for (rep in 1:2) {
  k <- k + 1; K <- 4; T_ <- 6
  pin <- c(strong = 0.5, default = 0.3, weak = 0.2)[[sep]]; pout <- c(strong = 0.02, default = 0.05, weak = 0.1)[[sep]]
  mem <- sample(1:K, n, TRUE); truth <- list(); L <- list()
  for (t in 1:T_) { if (t > 1) { s <- runif(n) < sw; mem[s] <- sample(1:K, sum(s), TRUE) }; truth[[t]] <- mem; L[[t]] <- bl(mem, pin, pout) }
  for (t in 1:T_) write.table(L[[t]], file.path(tmp, sprintf("net%03d_layer%d.csv", k, t)), sep = ",", row.names = FALSE, col.names = FALSE)
  fit <- fit_multilayer_identity_ties(L, algorithm = "leiden", omega = 1, seed = 123)
  mR <- extract_meta_membership(fit)
  writeLines(as.character(unlist(mR)), file.path(tmp, sprintf("net%03d_R.txt", k)))
  writeLines(as.character(unlist(truth)), file.path(tmp, sprintf("net%03d_truth.txt", k)))
  rows[[k]] <- data.frame(net = k, n = n, separation = sep, p_switch = sw, rep = rep, T_ = T_)
}
meta <- do.call(rbind, rows); write.csv(meta, file.path(tmp, "meta.csv"), row.names = FALSE)
py <- sprintf('
import numpy as np, pandas as pd, sys
from dynamic_multiplex import fit_multilayer_identity_ties, extract_meta_membership
from dynamic_multiplex.bootstrap_multilayer import _nmi
tmp = %s; meta = pd.read_csv(tmp + "/meta.csv"); out = []
for _, r in meta.iterrows():
    k = int(r.net); T = int(r.T_)
    L = [np.loadtxt(f"{tmp}/net{k:03d}_layer{t}.csv", delimiter=",") for t in range(1, T + 1)]
    fit = fit_multilayer_identity_ties(L, algorithm="leiden", omega=1.0, seed=123)
    mP = np.concatenate([np.asarray(m) for m in extract_meta_membership(fit)])
    mR = np.loadtxt(f"{tmp}/net{k:03d}_R.txt", dtype=int); tr = np.loadtxt(f"{tmp}/net{k:03d}_truth.txt", dtype=int)
    out.append(dict(net=k, n=int(r.n), separation=r.separation, p_switch=float(r.p_switch), rep=int(r.rep),
                    nmi_R_vs_python=_nmi(mR, mP), nmi_R_truth=_nmi(mR, tr), nmi_python_truth=_nmi(mP, tr),
                    K_R=len(set(mR)), K_python=len(set(mP))))
pd.DataFrame(out).to_csv(tmp + "/python_result.csv", index=False)
', shQuote(tmp))
writeLines(py, file.path(tmp, "run.py"))
st <- system2("python3", file.path(tmp, "run.py"), stdout = TRUE, stderr = TRUE)
res <- read.csv(file.path(tmp, "python_result.csv"))
write.csv(res, file.path(outdir, "crosscheck_multislice.csv"), row.names = FALSE)
options(width = 160); print(res, row.names = FALSE, digits = 3)
cat(sprintf("\nmean NMI(R, python) = %.3f; strong/default structure: %.3f; R vs truth %.3f, python vs truth %.3f\n",
            mean(res$nmi_R_vs_python), mean(res$nmi_R_vs_python[res$separation != "weak"]), mean(res$nmi_R_truth), mean(res$nmi_python_truth)))
stopifnot(mean(res$nmi_R_vs_python[res$separation != "weak"]) > 0.95)
cat("CROSSCHECK PASSED\n")
