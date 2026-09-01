# =============================================================================
# 08_alliance_dca_empirical.R
# Empirical case study: all DynMux specs + baselines on two real dynamic
# international networks (ATOP formal alliances 1816-2018; Kinne DCAD defense
# cooperation 1980-2010). DynSBM excluded (does not scale to the 203-layer
# open population). DynMux runs on native per-year graphs with isolate states
# dropped; fixed-node baselines run on the union node set with an active mask.
# Usage: Rscript 08_alliance_dca_empirical.R <atop|dca>
# Requires prebuilt <net>_series.rds and <net>_union.rds (see build scripts).
# =============================================================================
args<-commandArgs(trailingOnly=TRUE); net<-if(length(args)>=1) args[1] else "atop"
DATA<-Sys.getenv("EMP_DATA", "/tmp"); OUT<-Sys.getenv("EMP_OUT", "/tmp")
suppressMessages({library(igraph)}); suppressMessages(pkgload::load_all("r_code",quiet=TRUE))
HAVE_CLUE<-requireNamespace("clue",quietly=TRUE); HAVE_MULTINET<-requireNamespace("multinet",quietly=TRUE)
S<-readRDS(file.path(DATA,sprintf("%s_series.rds",net))); yrs<-S$years
U<-readRDS(file.path(DATA,sprintf("%s_union.rds",net)));GLc<-lapply(seq_along(S$graph_layers),function(k){g<-S$graph_layers[[k]];igraph::delete_vertices(g,which(!U$present[[k]]))})
U<-readRDS(file.path(DATA,sprintf("%s_union.rds",net))); simU<-list(layers=U$layers,active=U$active,links=NULL)
leiden_layer<-function(mat){g<-igraph::graph_from_adjacency_matrix(mat,mode="undirected",weighted=TRUE,diag=FALSE);if(igraph::ecount(g)==0L)return(rep(1L,nrow(mat)));as.integer(igraph::membership(igraph::cluster_leiden(g,objective_function="modularity",weights=igraph::E(g)$weight)))}
method_pooled<-function(sim){L<-sim$layers;agg<-Reduce("+",L);g<-igraph::graph_from_adjacency_matrix(agg,mode="undirected",weighted=TRUE,diag=FALSE);mem<-igraph::membership(igraph::cluster_leiden(g,objective_function="modularity",weights=igraph::E(g)$weight));replicate(length(L),as.integer(mem),simplify=FALSE)}
.greedy<-function(M){nr<-nrow(M);a<-rep(NA_integer_,nr);u<-integer(0);for(i in order(-apply(M,1,max))){c<-order(-M[i,]);c<-c[!(c%in%u)];a[i]<-c[1];u<-c(u,c[1])};a}
match_hung<-function(mems){o<-vector("list",length(mems));o[[1]]<-as.integer(mems[[1]]);nf<-max(o[[1]])+1L
 for(t in 2:length(mems)){pv<-o[[t-1L]];cu<-as.integer(mems[[t]]);cl<-sort(unique(cu));pl<-sort(unique(pv));M<-matrix(0,length(cl),length(pl))
  for(i in seq_along(cl))for(j in seq_along(pl))M[i,j]<-sum(cu==cl[i]&pv==pl[j]);d<-max(nrow(M),ncol(M));Ms<-matrix(0,d,d);Ms[seq_len(nrow(M)),seq_len(ncol(M))]<-M
  as<-if(HAVE_CLUE)as.integer(clue::solve_LSAP(max(Ms)-Ms))else .greedy(Ms);mp<-rep(NA_integer_,length(cl))
  for(i in seq_along(cl)){co<-as[i];if(co<=length(pl)&&M[i,co]>0)mp[i]<-pl[co]};for(i in seq_along(cl))if(is.na(mp[i])){mp[i]<-nf;nf<-nf+1L}
  nf<-max(nf,max(mp)+1L);names(mp)<-as.character(cl);o[[t]]<-unname(mp[as.character(cu)])};lapply(o,as.integer)}
method_multinet<-function(sim){if(!HAVE_MULTINET)return(NULL);L<-sim$layers;n<-nrow(L[[1]]);T_<-length(L);nn<-multinet::ml_empty();ac<-as.character(seq_len(n))
 for(t in seq_len(T_)){g<-igraph::graph_from_adjacency_matrix(L[[t]],mode="undirected",diag=FALSE);igraph::V(g)$name<-ac;multinet::add_igraph_layer_ml(nn,g,paste0("layer",t))}
 cm<-multinet::glouvain_ml(nn,gamma=1,omega=1);lapply(seq_len(T_),function(t){s<-cm[cm$layer==paste0("layer",t),,drop=FALSE];v<-rep(NA_integer_,n);v[as.integer(s$actor)]<-as.integer(as.factor(s$cid));if(anyNA(v))v[is.na(v)]<-max(0L,v,na.rm=TRUE)+seq_len(sum(is.na(v)));v})}
dm<-function(ff,rp)function(){f<-ff(GLc,algorithm="leiden",resolution_parameter=rp);m<-extract_meta_membership(f);lapply(seq_along(GLc),function(t)setNames(as.integer(m[[t]]),igraph::V(GLc[[t]])$name))}
comp<-function(fn)function(){d<-fn(simU);if(is.null(d))return(NULL);lapply(seq_along(d),function(t){a<-U$present[[t]];setNames(as.integer(d[[t]][a]),U$names[a])})}
METHODS<-list(
 "DynMux multislice r1"=dm(fit_multilayer_identity_ties,1),"DynMux multislice r2"=dm(fit_multilayer_identity_ties,2),"DynMux multislice r4"=dm(fit_multilayer_identity_ties,4),
 "DynMux Jaccard r1"=dm(fit_multilayer_jaccard,1),"DynMux Jaccard r2"=dm(fit_multilayer_jaccard,2),"DynMux Jaccard r4"=dm(fit_multilayer_jaccard,4),
 "DynMux Overlap r1"=dm(fit_multilayer_overlap,1),"DynMux Overlap r2"=dm(fit_multilayer_overlap,2),"DynMux Overlap r4"=dm(fit_multilayer_overlap,4),
 "Pooled Leiden"=comp(method_pooled),"Cross-sectional + Hungarian"=comp(function(s)match_hung(lapply(s$layers,leiden_layer))),
 "multinet GLouvain"=comp(method_multinet))
outf<-file.path(OUT,sprintf("%s_partitions.rds",net)); results<-list()
if(file.exists(outf)) results<-readRDS(outf)$partitions
for(mn in names(METHODS)){ if(!is.null(results[[mn]])){cat(sprintf("[%s] %-28s cached\n",net,mn));next}
 gc();t0<-proc.time()["elapsed"]
 det<-tryCatch({setTimeLimit(elapsed=1500,transient=TRUE);r<-METHODS[[mn]]();setTimeLimit();r},error=function(e){setTimeLimit();message(sprintf("  %s FAILED: %s",mn,conditionMessage(e)));NULL})
 results[[mn]]<-det;ok<-!is.null(det);nc<-if(ok)length(unique(unlist(det)))else NA
 saveRDS(list(net=net,years=yrs,partitions=results),outf)
 cat(sprintf("[%s] %-28s %s ncomm=%s (%.1f min)\n",net,mn,if(ok)"ok"else"NA",nc,(proc.time()["elapsed"]-t0)/60));flush.console()}
cat("ALL_DONE\n")
