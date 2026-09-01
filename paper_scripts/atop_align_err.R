suppressMessages(library(igraph)); set.seed(123)
DATA<-"replication/extended/output/empirical_data"; OUT<-"replication/extended/output/empirical"; net<-"atop"
S<-readRDS(file.path(DATA,sprintf("%s_series.rds",net))); U<-readRDS(file.path(DATA,sprintf("%s_union.rds",net)))
P<-readRDS(file.path(OUT,sprintf("%s_partitions.rds",net)))$partitions; yrs<-S$years
G<-lapply(seq_along(S$graph_layers),function(k) igraph::delete_vertices(S$graph_layers[[k]],which(!U$present[[k]])))
Tn<-length(G)
NB<-lapply(G,function(g){nm<-V(g)$name; al<-igraph::as_adj_list(g,mode="all"); setNames(lapply(al,function(ids) nm[as.integer(ids)]),nm)})
jac<-function(x,y){u<-length(union(x,y)); if(u==0) return(1); length(intersect(x,y))/u}
errfun<-function(ml){
 spb<-msb<-nb<-integer(Tn)
 for(t in 1:(Tn-1)){
  a<-ml[[t]]; b<-ml[[t+1]]; if(is.null(a)||is.null(b)) next
  common<-intersect(names(a),names(b)); if(length(common)<3) next
  la<-a[common]; lb<-b[common]; Na<-NB[[t]]; Nb<-NB[[t+1]]
  gA<-split(common,as.vector(la)); gB<-split(common,as.vector(lb))
  for(i in common){
   nti<-intersect(Na[[i]],common); ntj<-intersect(Nb[[i]],common)
   if(length(nti)<1 || length(ntj)<1) next
   J<-jac(nti,ntj)
   Ca<-setdiff(gA[[as.character(la[i])]],i); Cb<-setdiff(gB[[as.character(lb[i])]],i); coJ<-jac(Ca,Cb)
   nb[t]<-nb[t]+1
   if(J>=0.9 && coJ<=0.5) spb[t]<-spb[t]+1
   if(J<=0.3 && coJ>=0.9) msb[t]<-msb[t]+1
  }
 }
 list(spb=spb,msb=msb,nb=nb)
}
meth<-c("DynMux Jaccard r1","DynMux Overlap r1","DynMux multislice r1","multinet GLouvain","Cross-sectional + Hungarian","Pooled Leiden")
R<-lapply(meth,function(m) errfun(P[[m]])); names(R)<-meth
cat(sprintf("%-30s %6s %6s %6s %6s %7s\n","method","spur","miss","total","n","rate%"))
for(m in meth){r<-R[[m]]; sp<-sum(r$spb); ms<-sum(r$msb); n<-sum(r$nb); cat(sprintf("%-30s %6d %6d %6d %6d %7.2f\n",m,sp,ms,sp+ms,n,100*(sp+ms)/n))}
eras<-list(c(1816,1900),c(1901,1945),c(1946,2018))
cat("\n--- total error rate (%) by era ---\n"); cat(sprintf("%-30s %10s %10s %10s\n","method","1816-1900","1901-1945","1946-2018"))
for(m in meth){r<-R[[m]]; vals<-sapply(eras,function(e){i<-which(yrs>=e[1]&yrs<=e[2]); n<-sum(r$nb[i]); if(n==0)return(NA); 100*(sum(r$spb[i])+sum(r$msb[i]))/n}); cat(sprintf("%-30s %10.2f %10.2f %10.2f\n",m,vals[1],vals[2],vals[3]))}
