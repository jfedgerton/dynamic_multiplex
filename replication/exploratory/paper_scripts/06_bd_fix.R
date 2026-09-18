#!/usr/bin/env Rscript
# Corrected BIRTHDEATH-ONLY rerun. Coupling/tracking methods fed the per-layer
# ALIVE node set (sim$active): absent nodes dropped, alive-isolates kept.
# Only birthdeath configs, only affected methods (DynMux x{Leiden,Louvain},
# Hungarian, multinet). Pooled + DynSBM excluded (immune / already correct).
# Identical grid + seeds; 72-task array self-skips non-birthdeath.
suppressMessages({ library(igraph) })
suppressMessages(pkgload::load_all("r_code", quiet = TRUE))
source("replication/extended/paper_scripts/gen06.R")
HAVE_CLUE     <- requireNamespace("clue",     quietly = TRUE)
HAVE_MULTINET <- requireNamespace("multinet", quietly = TRUE)
metric_cols <- c("nmi_layer","nmi_joint","nmi_change","k_mae",
                 "comembership_acc","mean_n_comm","total_n_comm")
cfgs <- expand.grid(
  regime    = c("seasonality","churnswitch","birthdeath","regimeshift"),
  n         = c(50L,100L,200L),
  r         = c(1.5,3,6),
  intensity = c("low","high"),
  stringsAsFactors = FALSE)
stopifnot(nrow(cfgs) == 72L)
set.seed(123); ord <- sample(nrow(cfgs))
TASK <- suppressWarnings(as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = NA)))
if (is.na(TASK)) TASK <- suppressWarnings(as.integer(commandArgs(TRUE)[1]))
stopifnot(!is.na(TASK), TASK >= 1L, TASK <= nrow(cfgs))
cfg <- cfgs[ord[TASK], ]
outdir  <- "replication/extended/output/dynamic_bdfix"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
outfile <- file.path(outdir, sprintf("dyn_cfg%02d.csv", TASK))
if (!identical(cfg$regime, "birthdeath")) { cat("[skip non-bd] task",TASK,cfg$regime,"\n"); quit(save="no") }
if (nzchar(Sys.getenv("SLURM_ARRAY_TASK_ID")) && file.exists(outfile)) { cat("[skip done]",TASK,"\n"); quit(save="no") }
alive_idx <- function(sim, t)
  if (!is.null(sim$active)) which(as.logical(sim$active[[t]])) else seq_len(nrow(sim$layers[[t]]))
AL <- function(sim) lapply(seq_along(sim$layers), function(t){
  idx <- alive_idx(sim, t)
  A   <- sim$layers[[t]][idx, idx, drop = FALSE]
  dimnames(A) <- list(as.character(idx), as.character(idx)); A })
EXPn <- function(m, sim, n) lapply(seq_along(m), function(t){
  idx <- alive_idx(sim, t); v <- rep(NA_integer_, n); v[idx] <- as.integer(m[[t]]); v })
leiden_layer <- function(mat){
  g <- igraph::graph_from_adjacency_matrix(mat, mode="undirected", weighted=TRUE, diag=FALSE)
  if (igraph::ecount(g)==0L) return(rep(1L, nrow(mat)))
  as.integer(igraph::membership(igraph::cluster_leiden(g, objective_function="modularity",
                                                       weights=igraph::E(g)$weight))) }
.greedy_assign <- function(M){ nr<-nrow(M); a<-rep(NA_integer_,nr); u<-integer(0)
  for(i in order(-apply(M,1,max))){ c<-order(-M[i,]); c<-c[!(c%in%u)]; a[i]<-c[1]; u<-c(u,c[1]) }; a }
dm <- function(fitfun, alg, links=TRUE){
  function(sim){
    n <- nrow(sim$layers[[1]])
    f <- if (links) fitfun(AL(sim), algorithm=alg, layer_links=sim$links)
         else        fitfun(AL(sim), algorithm=alg)
    EXPn(extract_meta_membership(f), sim, n)
  }
}
method_multinet <- function(sim){
  if (!HAVE_MULTINET) return(NULL)
  n <- nrow(sim$layers[[1]]); T_ <- length(sim$layers); nn <- multinet::ml_empty()
  for (t in seq_len(T_)){
    idx <- alive_idx(sim,t); if (!length(idx)) next
    A <- sim$layers[[t]][idx,idx,drop=FALSE]
    g <- igraph::graph_from_adjacency_matrix(A, mode="undirected", diag=FALSE)
    igraph::V(g)$name <- as.character(idx)
    multinet::add_igraph_layer_ml(nn, g, paste0("layer",t))
  }
  cm <- multinet::glouvain_ml(nn, gamma=1, omega=1)
  lapply(seq_len(T_), function(t){
    s <- cm[cm$layer==paste0("layer",t), , drop=FALSE]
    v <- rep(NA_integer_, n)
    if (nrow(s)) v[as.integer(s$actor)] <- as.integer(as.factor(s$cid)); v })
}
method_hungarian <- function(sim){
  n <- nrow(sim$layers[[1]]); T_ <- length(sim$layers)
  mems <- lapply(sim$layers, leiden_layer)
  out  <- vector("list", T_)
  a1 <- alive_idx(sim,1); v1 <- rep(NA_integer_,n); v1[a1] <- as.integer(mems[[1]][a1])
  out[[1]] <- v1; next_free <- max(v1, na.rm=TRUE) + 1L
  for (t in 2:T_){
    prev <- out[[t-1L]]; cur <- as.integer(mems[[t]])
    ac <- alive_idx(sim,t); ap <- alive_idx(sim,t-1L)
    avc <- logical(n); avc[ac] <- TRUE; avp <- logical(n); avp[ap] <- TRUE; both <- avc & avp
    cur_labs  <- sort(unique(cur[avc]))
    prev_labs <- sort(unique(prev[avp & !is.na(prev)]))
    M <- matrix(0, length(cur_labs), length(prev_labs))
    for(i in seq_along(cur_labs)) for(j in seq_along(prev_labs))
      M[i,j] <- sum(cur==cur_labs[i] & prev==prev_labs[j] & both, na.rm=TRUE)
    d <- max(nrow(M),ncol(M)); Ms <- matrix(0,d,d); Ms[seq_len(nrow(M)),seq_len(ncol(M))] <- M
    asg <- if (HAVE_CLUE) as.integer(clue::solve_LSAP(max(Ms)-Ms)) else .greedy_assign(Ms)
    mp <- rep(NA_integer_, length(cur_labs))
    for(i in seq_along(cur_labs)){ co<-asg[i]; if(co<=length(prev_labs) && M[i,co]>0) mp[i]<-prev_labs[co] }
    for(i in seq_along(cur_labs)) if(is.na(mp[i])){ mp[i]<-next_free; next_free<-next_free+1L }
    next_free <- max(next_free, max(mp,na.rm=TRUE)+1L); names(mp) <- as.character(cur_labs)
    v <- rep(NA_integer_, n); v[avc] <- mp[as.character(cur[avc])]
    out[[t]] <- v
  }
  out
}
METHODS <- list(
  "DynMux multislice (adjacent)"         = dm(fit_multilayer_identity_ties, "leiden", links=FALSE),
  "DynMux multislice (custom)"           = dm(fit_multilayer_identity_ties, "leiden", links=TRUE),
  "DynMux Jaccard"                       = dm(fit_multilayer_jaccard,        "leiden"),
  "DynMux Overlap"                       = dm(fit_multilayer_overlap,        "leiden"),
  "DynMux weighted Jaccard"              = dm(fit_multilayer_weighted_jaccard,"leiden"),
  "DynMux weighted Overlap"              = dm(fit_multilayer_weighted_overlap,"leiden"),
  "Cross-sectional + Hungarian"          = method_hungarian,
  "multinet GLouvain"                    = method_multinet,
  "DynMux multislice (adjacent, Louvain)"= dm(fit_multilayer_identity_ties, "louvain", links=FALSE),
  "DynMux Jaccard (Louvain)"             = dm(fit_multilayer_jaccard,        "louvain"),
  "DynMux Overlap (Louvain)"             = dm(fit_multilayer_overlap,        "louvain"),
  "DynMux weighted Jaccard (Louvain)"    = dm(fit_multilayer_weighted_jaccard,"louvain"),
  "DynMux weighted Overlap (Louvain)"    = dm(fit_multilayer_weighted_overlap,"louvain")
)
run_rep <- function(rep){
  seed <- 9000L + TASK*1000L + rep
  sim  <- simulate_regime(cfg, seed)
  rows <- vector("list", length(METHODS)); ri <- 0L
  for (mname in names(METHODS)){
    t0  <- proc.time()["elapsed"]
    det <- tryCatch(METHODS[[mname]](sim), error=function(e){
      message(sprintf("  [rep %d] %s errored: %s", rep, mname, conditionMessage(e))); NULL })
    el  <- as.numeric(proc.time()["elapsed"] - t0)
    m   <- if (is.null(det)) setNames(rep(NA_real_, length(metric_cols)), metric_cols) else eval_method(det, sim)
    ri  <- ri + 1L
    rows[[ri]] <- data.frame(
      regime=cfg$regime, n=cfg$n, r=cfg$r, intensity=cfg$intensity, rep=rep, method=mname,
      nmi_layer=round(m[["nmi_layer"]],4), nmi_joint=round(m[["nmi_joint"]],4),
      nmi_change=round(m[["nmi_change"]],4), k_mae=round(m[["k_mae"]],4),
      comembership_acc=round(m[["comembership_acc"]],4), runtime_s=round(el,3),
      mean_n_comm=round(m[["mean_n_comm"]],4), total_n_comm=m[["total_n_comm"]],
      stringsAsFactors=FALSE)
  }
  do.call(rbind, rows)
}
t0 <- proc.time()["elapsed"]
all_rows <- vector("list", REPS_FAST)
for (rep in seq_len(REPS_FAST)){
  all_rows[[rep]] <- tryCatch(run_rep(rep), error=function(e){
    message(sprintf("rep %d failed: %s", rep, conditionMessage(e))); NULL })
  if (rep %% 10L == 0L) cat(sprintf("  ... rep %d/%d (%.1f min)\n", rep, REPS_FAST,
                                    (proc.time()["elapsed"]-t0)/60))
}
res <- do.call(rbind, Filter(Negate(is.null), all_rows))
write.csv(res, outfile, row.names = FALSE)
cat(sprintf("[bdfix] task %d done: %d rows -> %s (%.1f min)\n",
            TASK, nrow(res), outfile, (proc.time()["elapsed"]-t0)/60))
