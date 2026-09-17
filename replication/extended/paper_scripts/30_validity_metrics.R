# Per-year validity metrics vs Braumoeller/Goodhart order codings, all layers/methods.
# Outputs ARI plus the Rand decomposition (pair recall / pair specificity) used for Fig 1A.
suppressMessages(library(igraph)); set.seed(123)
setwd("/storage/work/jfe4/dynamic_multiplex")
src <- readLines("replication/extended/paper_scripts/29_braumoeller_fullspan.R")
eval(parse(text=paste(src[1:(grep("^run<-function",src)-1)],collapse="\n")))
ch2<-function(x) x*(x-1)/2
metrics<-function(tv,mv){ tb<-table(tv,mv); n<-sum(tb)
 same_true<-sum(ch2(rowSums(tb))); diff_true<-ch2(n)-same_true
 same_both<-sum(ch2(tb)); same_pred<-sum(ch2(colSums(tb)))
 recall <- if(same_true>0) same_both/same_true else NA          # TPR-like
 spec   <- if(diff_true>0) 1-(same_pred-same_both)/diff_true else NA  # TNR-like
 c(recall=recall, spec=spec)}
meth<-c("DynMux Jaccard r1","DynMux Overlap r1","DynMux multislice r1","multinet GLouvain","Cross-sectional + Hungarian","Pooled Leiden")
mlab<-c("Jaccard","Overlap","multislice","multinet","Hungarian","Pooled")
rows<-list()
for(net in c("atop","igo","trade")){
 S<-readRDS(sprintf("replication/extended/output/empirical_data/%s_series.rds",net))
 U<-readRDS(sprintf("replication/extended/output/empirical_data/%s_union.rds",net))
 P<-readRDS(sprintf("replication/extended/output/empirical/%s_partitions.rds",net))$partitions
 yrs<-S$years; mask<-if(!is.null(U$present))U$present else U$active
 G<-lapply(seq_along(S$graph_layers),function(k) delete_vertices(S$graph_layers[[k]],which(!mask[[k]])))
 dg<-lapply(G,function(g) setNames(degree(g),V(g)$name))
 for(y in yrs){t<-which(yrs==y); d<-dg[[t]]; A<-P[[meth[1]]][[t]]
  al<-names(d)[d>0]; al<-al[al%in%names(A)]; if(length(al)<8)next
  tt<-truth(y,al); if(is.null(tt))next
  for(j in seq_along(meth)){Ai<-P[[meth[j]]][[t]]; if(!all(tt$n%in%names(Ai)))next
   mv<-as.integer(unlist(Ai[tt$n])); m<-metrics(tt$t,mv)
   rows[[length(rows)+1]]<-data.frame(net=net,year=y,era=tt$era,method=mlab[j],
     ari=round(ari(tt$t,mv),4), recall=round(m[["recall"]],4), spec=round(m[["spec"]],4),
     ncomm=length(unique(mv)), nnode=length(tt$n))}}}
D<-do.call(rbind,rows)
saveRDS(D,"replication/extended/paper_scripts/validity_metrics.rds")
write.csv(D,"/tmp/vm.csv",row.names=FALSE)
cat("rows:",nrow(D)," nets:",length(unique(D$net))," methods:",length(unique(D$method)),"\n")
print(table(D$net,D$era))
