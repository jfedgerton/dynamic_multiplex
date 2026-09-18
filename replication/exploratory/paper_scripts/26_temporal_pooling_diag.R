# Diagnostic: does interlayer coupling smooth across real historical breaks?
# (1) consecutive-layer ARI by method -> does the partition register 1918/1945/1989?
# (2) bloc membership: are POL/HUN/CZR with RUS in the Cold War and with USA after?
suppressMessages(library(igraph))
S<-readRDS("replication/extended/output/empirical_data/atop_series.rds")
U<-readRDS("replication/extended/output/empirical_data/atop_union.rds")
P<-readRDS("replication/extended/output/empirical/atop_partitions.rds")$partitions
yrs<-S$years; mask<-if(!is.null(U$present))U$present else U$active
G<-lapply(seq_along(S$graph_layers),function(k) delete_vertices(S$graph_layers[[k]],which(!mask[[k]])))
deg<-lapply(G,function(g) setNames(degree(g),V(g)$name))
meth<-c("DynMux Jaccard r1","DynMux Overlap r1","DynMux multislice r1","multinet GLouvain","Cross-sectional + Hungarian","Pooled Leiden")
mlab<-c("Jaccard","Overlap","multislice","multinet","Hungarian","Pooled")
# ---- (1) consecutive-layer ARI on nodes active in both years ----
ari<-matrix(NA,length(yrs)-1,length(meth),dimnames=list(yrs[-1],mlab))
for(j in seq_along(meth)){m<-meth[j]
for(t in 2:length(yrs)){A<-P[[m]][[t-1]];B<-P[[m]][[t]];d1<-deg[[t-1]];d2<-deg[[t]]
ns<-intersect(names(A)[names(A)%in%names(d1)],names(B)[names(B)%in%names(d2)])
ns<-ns[sapply(ns,function(n) d1[[n]]>0 && d2[[n]]>0)]
if(length(ns)<5)next
ari[t-1,j]<-igraph::compare(as.integer(unlist(A[ns])),as.integer(unlist(B[ns])),method="adjusted.rand")}}
cat("=== MEAN CONSECUTIVE-LAYER ARI (higher = more persistent) ===\n")
print(round(colMeans(ari,na.rm=TRUE),3))
cat("\n=== ARI AT HISTORICAL BREAKS (year = transition INTO that year) ===\n")
brk<-c(1918,1919,1945,1946,1989,1990,1991,1992)
print(round(ari[as.character(brk[brk%in%rownames(ari)]),],2))
cat("\n=== MEAN ARI WITHIN STABLE SPANS vs AT BREAKS ===\n")
stab<-rownames(ari)[as.integer(rownames(ari))%in%c(1955:1985,1995:2015)]
brkw<-rownames(ari)[as.integer(rownames(ari))%in%c(1917:1920,1944:1947,1989:1992)]
cat("stable:"); print(round(colMeans(ari[stab,],na.rm=TRUE),3))
cat("breaks:"); print(round(colMeans(ari[brkw,],na.rm=TRUE),3))
# ---- (2) bloc membership of former Warsaw Pact states ----
key<-c("290"="POL","310"="HUN","316"="CZR","315"="CZE","265"="GDR")
sameC<-function(m,t,a,b){A<-P[[m]][[t]];d<-deg[[t]]
if(!(a%in%names(A)&&b%in%names(A)&&a%in%names(d)&&b%in%names(d)))return(NA)
if(d[[a]]<=0||d[[b]]<=0)return(NA); A[[a]]==A[[b]]}
cat("\n=== BLOC MEMBERSHIP: U=with USA, R=with RUS, B=both, n=neither, .=inactive ===\n")
for(kk in seq_along(key)){cat("\n--",key[kk],"--\n     ")
ys<-c(1950,1960,1970,1980,1985,1988,1990,1993,1996,1999,2002,2008,2015)
cat(paste(sprintf("%5d",ys),collapse=""),"\n")
for(j in seq_along(meth)){cat(sprintf("%-11s",mlab[j]))
for(y in ys){t<-which(yrs==y); u<-sameC(meth[j],t,"2",names(key)[kk]); r<-sameC(meth[j],t,"365",names(key)[kk])
s<-if(is.na(u)||is.na(r))"." else if(u&&r)"B" else if(u)"U" else if(r)"R" else "n"
cat(sprintf("%5s",s))}; cat("\n")}}
