suppressMessages(library(igraph)); set.seed(123)
DATA<-"replication/extended/output/empirical_data"; OUT<-"replication/extended/output/empirical"; net<-"atop"
S<-readRDS(file.path(DATA,sprintf("%s_series.rds",net))); U<-readRDS(file.path(DATA,sprintf("%s_union.rds",net)))
P<-readRDS(file.path(OUT,sprintf("%s_partitions.rds",net)))$partitions; yrs<-S$years
G<-lapply(seq_along(S$graph_layers),function(k) igraph::delete_vertices(S$graph_layers[[k]],which(!U$present[[k]])))
mod_m<-function(ml){sapply(seq_along(G),function(t){g<-G[[t]];mem<-ml[[t]];if(is.null(mem)||vcount(g)==0)return(NA_real_);nm<-intersect(V(g)$name,names(mem));if(length(nm)<2)return(NA_real_);gg<-induced_subgraph(g,which(V(g)$name%in%nm));mm<-as.integer(mem[V(gg)$name]);tryCatch(modularity(gg,as.integer(factor(mm))),error=function(e)NA_real_)})}
mod_opt<-sapply(seq_along(G),function(t){g<-G[[t]];if(vcount(g)<2||ecount(g)==0)return(NA_real_);cl<-tryCatch(cluster_louvain(g),error=function(e)NULL);if(is.null(cl))return(NA_real_);modularity(cl)})
M<-list(Pooled=mod_m(P[["Pooled Leiden"]]),Overlap=mod_m(P[["DynMux Overlap r2"]]),Hungar=mod_m(P[["Cross-sectional + Hungarian"]]),OPTIMAL=mod_opt)
eras<-list(c(1816,1900),c(1901,1945),c(1946,1990),c(1991,2018))
cat(sprintf("%-12s %7s %7s %7s %7s\n","era","Pooled","Overlap","Hungar","OPTIMAL"))
for(e in eras){i<-which(yrs>=e[1]&yrs<=e[2]);cat(sprintf("%-12s %7.3f %7.3f %7.3f %7.3f\n",paste0(e[1],"-",e[2]),mean(M$Pooled[i],na.rm=TRUE),mean(M$Overlap[i],na.rm=TRUE),mean(M$Hungar[i],na.rm=TRUE),mean(M$OPTIMAL[i],na.rm=TRUE)))}
cat(sprintf("%-12s %7.3f %7.3f %7.3f %7.3f\n","ALL",mean(M$Pooled,na.rm=TRUE),mean(M$Overlap,na.rm=TRUE),mean(M$Hungar,na.rm=TRUE),mean(M$OPTIMAL,na.rm=TRUE)))
pl<-P[["Pooled Leiden"]]; cat("\nPooled ncomm layer1:",length(unique(pl[[1]])),"  top sizes:",paste(head(as.integer(sort(table(pl[[1]]),decreasing=TRUE)),8),collapse=" "),"\n")
same<-all(sapply(2:length(pl),function(t){k<-intersect(names(pl[[t]]),names(pl[[1]]));identical(unname(pl[[t]][k]),unname(pl[[1]][k]))})); cat("Pooled frozen (identical labels across layers):",same,"\n")
cat("median Optimal-Pooled gap:",round(median(M$OPTIMAL-M$Pooled,na.rm=TRUE),3),"\n")
