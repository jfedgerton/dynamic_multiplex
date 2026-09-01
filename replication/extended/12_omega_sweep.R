set.seed(123)
suppressMessages({ library(igraph); pkgload::load_all("r_code", quiet = TRUE) })
DATA <- "replication/extended/output/empirical_data"
OUT  <- "replication/extended/output/empirical"
OMEGA <- c(0.1, 0.25, 0.5, 1, 2, 4, 8)
nmi <- function(a, b) { k <- intersect(names(a), names(b)); if (length(k) < 2) return(NA_real_)
  igraph::compare(as.integer(a[k]), as.integer(b[k]), method = "nmi") }
metrics <- function(part, G) { Tn <- length(part); modv <- rep(NA_real_, Tn)
  for (t in seq_len(Tn)) { g <- G[[t]]; mem <- part[[t]]; if (is.null(mem) || vcount(g) == 0) next
    nm <- intersect(V(g)$name, names(mem)); if (length(nm) < 2) next
    gg <- induced_subgraph(g, which(V(g)$name %in% nm)); mm <- as.integer(mem[V(gg)$name])
    modv[t] <- tryCatch(modularity(gg, as.integer(factor(mm))), error = function(e) NA_real_) }
  pv <- sapply(seq_len(Tn - 1), function(t) nmi(part[[t]], part[[t + 1]]))
  c(mod = mean(modv, na.rm = TRUE), persist = mean(pv, na.rm = TRUE)) }
sweep <- function(net) {
  S <- readRDS(file.path(DATA, sprintf("%s_series.rds", net)))
  G <- lapply(S$graph_layers, function(g) delete_vertices(g, which(degree(g) == 0)))
  cat(sprintf("\n== %s omega sweep ==\n", toupper(net)))
  do.call(rbind, lapply(OMEGA, function(om) {
    f <- fit_multilayer_identity_ties(G, algorithm = "leiden", resolution_parameter = 1, omega = om)
    m <- extract_meta_membership(f)
    part <- lapply(seq_along(G), function(t) setNames(as.integer(m[[t]]), V(G[[t]])$name))
    r <- metrics(part, G)
    cat(sprintf("  omega=%-5g mod=%.3f persist=%.3f\n", om, r["mod"], r["persist"]))
    data.frame(network = net, omega = om, mod = round(r["mod"],3), persist = round(r["persist"],3),
               totcomm = length(unique(unlist(part)))) })) }
res <- rbind(sweep("atop"), sweep("dca"))
write.csv(res, file.path(OUT, "omega_sweep.csv"), row.names = FALSE); cat("\nSAVED omega_sweep.csv\n")
