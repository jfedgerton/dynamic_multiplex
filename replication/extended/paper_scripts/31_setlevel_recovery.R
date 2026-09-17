# Set-level bloc recovery: does a method produce a community that IS the bloc?
# For true bloc B, over detected communities C:  J = max |B n C| / |B u C|
# with precision |BnC|/|C| and recall |BnC|/|B| at the argmax, to separate
# "bloc split apart" (low recall) from "bloc buried in a bigger community" (low precision).
# Complements ARI, which is pair-based and dominated by the ~154-state residual bloc.
suppressMessages(library(igraph)); set.seed(123)
setwd("/storage/work/jfe4/dynamic_multiplex")
src<-readLines("replication/extended/paper_scripts/29_braumoeller_fullspan.R")
eval(parse(text=paste(src[1:(grep("^run<-function",src)-1)],collapse="\n")))
meth<-c("DynMux Jaccard r1","DynMux Overlap r1","DynMux multislice r1","multinet GLouvain","Cross-sectional + Hungarian","Pooled Leiden")
mlab<-c("Jaccard","Overlap","multislice","multinet","Hungarian","Pooled")
best<-function(B, mem, nodes){
  if(!length(B)) return(c(NA,NA,NA))
  bs<-nodes %in% B; out<-c(0,0,0)
  for(cc in unique(mem)){ cs<-mem==cc; i<-sum(bs & cs)
    j<-i/sum(bs | cs); if(j>out[1]) out<-c(j, i/sum(cs), i/sum(bs)) }
  out}
rows<-list()
for(net in c("atop","igo","trade")){
 S<-readRDS(sprintf("replication/extended/output/empirical_data/%s_series.rds",net))
 U<-readRDS(sprintf("replication/extended/output/empirical_data/%s_union.rds",net))
 P<-readRDS(sprintf("replication/extended/output/empirical/%s_partitions.rds",net))$partitions
 yrs<-S$years; mask<-if(!is.null(U$present))U$present else U$active
 G<-lapply(seq_along(S$graph_layers),function(k) delete_vertices(S$graph_layers[[k]],which(!mask[[k]])))
 dg<-lapply(G,function(g) setNames(degree(g),V(g)$name))
 for(y in 1946:1991){t<-which(yrs==y); if(!length(t))next; d<-dg[[t]]; A0<-P[[meth[1]]][[t]]
  al<-names(d)[d>0]; al<-al[al%in%names(A0)]
  lib<-intersect(inb(LIB,y),al); com<-intersect(inb(COM,y),al)
  if(length(lib)<3||length(com)<3) next
  for(j in seq_along(meth)){A<-P[[meth[j]]][[t]]; if(!all(al%in%names(A)))next
   mem<-as.integer(unlist(A[al]))
   bl<-best(lib,mem,al); bc<-best(com,mem,al)
   rows[[length(rows)+1]]<-data.frame(net=net,year=y,method=mlab[j],
     lib_J=bl[1],lib_prec=bl[2],lib_rec=bl[3],
     com_J=bc[1],com_prec=bc[2],com_rec=bc[3])}}}
D<-do.call(rbind,rows)
saveRDS(D,"replication/extended/paper_scripts/setlevel_recovery.rds")
cat("rows:",nrow(D),"years per net:\n"); print(table(D$net)/6)
for(nt in c("atop","igo","trade")){s<-D[D$net==nt,]
 cat(sprintf("\n### %s  (best-match set recovery, Cold War)\n",toupper(nt)))
 cat(sprintf("%-11s %-24s %-24s\n","","LIBERAL bloc","COMMUNIST bloc"))
 cat(sprintf("%-11s %7s %7s %7s %7s %7s %7s\n","method","J","prec","recall","J","prec","recall"))
 a<-aggregate(cbind(lib_J,lib_prec,lib_rec,com_J,com_prec,com_rec)~method,s,mean)
 a<-a[order(-a$lib_J),]
 for(i in 1:nrow(a)) cat(sprintf("%-11s %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\n",
   a$method[i],a$lib_J[i],a$lib_prec[i],a$lib_rec[i],a$com_J[i],a$com_prec[i],a$com_rec[i]))}
