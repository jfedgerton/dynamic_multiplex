# Analyze saved empirical partitions: modularity, temporal persistence, community counts,
# break diagnostic (min year-to-year persistence), and cross-spec agreement. Net in {atop, dca}.
suppressMessages({library(igraph)})
DATA <- Sys.getenv("EMP_DATA","replication/extended/output/empirical_data")
OUT  <- Sys.getenv("EMP_OUT","replication/extended/output/empirical")
nmi <- function(a,b){ k<-intersect(names(a),names(b)); if(length(k)<2) return(NA_real_)
  igraph::compare(as.integer(a[k]), as.integer(b[k]), method="nmi") }
analyze <- function(net){
  S <- readRDS(file.path(DATA,sprintf("%s_series.rds",net))); yrs <- S$years
  U <- readRDS(file.path(DATA,sprintf("%s_union.rds",net)))
  P <- readRDS(file.path(OUT,sprintf("%s_partitions.rds",net)))$partitions
  G <- lapply(seq_along(S$graph_layers),function(k) igraph::delete_vertices(S$graph_layers[[k]],which(!U$present[[k]])))
  rows <- list(); traj <- list()
  for(m in names(P)){ part<-P[[m]]; if(is.null(part)) next; Tn<-length(part)
    modv<-rep(NA_real_,Tn); ncv<-rep(NA_real_,Tn)
    for(t in seq_len(Tn)){ g<-G[[t]]; mem<-part[[t]]; if(is.null(mem)||igraph::vcount(g)==0) next
      nm<-intersect(igraph::V(g)$name,names(mem)); if(length(nm)<2) next
      gg<-igraph::induced_subgraph(g, which(igraph::V(g)$name %in% nm)); mm<-as.integer(mem[igraph::V(gg)$name])
      ncv[t]<-length(unique(mm)); modv[t]<-tryCatch(igraph::modularity(gg,as.integer(factor(mm))),error=function(e)NA_real_) }
    pv<-sapply(seq_len(Tn-1),function(t) nmi(part[[t]],part[[t+1]])); traj[[m]]<-pv
    rows[[m]]<-data.frame(method=m,mean_mod=round(mean(modv,na.rm=TRUE),3),
      mean_persist=round(mean(pv,na.rm=TRUE),3),mean_ncomm=round(mean(ncv,na.rm=TRUE),1),
      total_ncomm=length(unique(unlist(part)))) }
  df<-do.call(rbind,rows); rownames(df)<-NULL
  cat("\n==============",toupper(net),"==============\n"); print(df,row.names=FALSE)
  for(sp in c("DynMux Overlap r2","DynMux multislice r2","Cross-sectional + Hungarian","Pooled Leiden")){
    pv<-traj[[sp]]; if(is.null(pv)||all(is.na(pv))) next; i<-which.min(pv)
    cat(sprintf("  min-persistence(break) %-27s ~%d->%d nmi=%.2f\n",sp,yrs[i],yrs[i+1],pv[i])) }
  pr<-function(a,b){pa<-P[[a]];pb<-P[[b]];if(is.null(pa)||is.null(pb))return(NA);mean(mapply(nmi,pa,pb),na.rm=TRUE)}
  cat(sprintf("  agreement: Overlap~multislice=%.2f  Overlap~Pooled=%.2f  multislice~Pooled=%.2f\n",
    pr("DynMux Overlap r2","DynMux multislice r2"),pr("DynMux Overlap r2","Pooled Leiden"),pr("DynMux multislice r2","Pooled Leiden")))
  list(net=net,df=df,traj=traj,years=yrs) }
res<-list(atop=analyze("atop"), dca=analyze("dca"))
saveRDS(res, file.path(OUT,"empirical_summary.rds")); cat("\nSAVED empirical_summary.rds\n")
