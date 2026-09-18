# External validity: score ATOP community partitions against Braumoeller/Goodhart
# Cold War bloc coding (Only the Dead replication, Replication.R lines 1637-1795,
# coding credited to Andy Goodhart). Liberal/Communist/Other is exhaustive and
# mutually exclusive 1945-1991 (verified: 0 state-years in >1 bloc).
# NOTE: original coding has ccode 100 (Colombia) labelled "Bulgaria"; Bulgaria (355)
# is therefore absent from the Communist bloc. FIXBUG toggles the correction.
suppressMessages(library(igraph)); set.seed(123)
FIXBUG <- TRUE
LIB <- c("2"="1945:1991","20"="1945:1991","200"="1945:1991","205"="1945:1991","210"="1945:1991","211"="1945:1991","212"="1945:1991","220"="1945:1991","225"="1945:1991","230"="1982:1991","235"="1945:1991","255"="1990:1991","260"="1955:1990","305"="1945:1991","325"="1945:1991","350"="1952:1991","380"="1945:1991","385"="1945:1991","390"="1945:1991","395"="1945:1991","640"="1952:1991","666"="1948:1991","713"="1947:1991","732"="1948:1991","740"="1952:1991","900"="1945:1991","920"="1945:1991")
COM <- c("100"="1955:1991","265"="1955:1990","290"="1955:1991","310"="1955:1991","315"="1955:1991","339"="1955:1991","359"="1945:1991","360"="1955:1991","365"="1955:1991","366"="1945:1991","367"="1945:1991","368"="1945:1991","369"="1945:1991","370"="1945:1991","371"="1945:1991","372"="1945:1991","373"="1945:1991","701"="1945:1991","702"="1945:1991","703"="1945:1991","704"="1945:1991","705"="1945:1991","710"="1949:1991","712"="1945:1991","731"="1948:1991","812"="1949:1991","816"="1945:1991")
if(FIXBUG){ COM <- COM[names(COM)!="100"]; COM["355"] <- "1955:1991" }
inb <- function(tbl,y){ ks<-names(tbl); ks[sapply(tbl,function(r){p<-as.integer(strsplit(r,":")[[1]]); y>=p[1] && y<=p[2]})] }
S<-readRDS("replication/extended/output/empirical_data/atop_series.rds")
U<-readRDS("replication/extended/output/empirical_data/atop_union.rds")
P<-readRDS("replication/extended/output/empirical/atop_partitions.rds")$partitions
yrs<-S$years; mask<-if(!is.null(U$present))U$present else U$active
G<-lapply(seq_along(S$graph_layers),function(k) delete_vertices(S$graph_layers[[k]],which(!mask[[k]])))
deg<-lapply(G,function(g) setNames(degree(g),V(g)$name))
meth<-c("DynMux Jaccard r1","DynMux Overlap r1","DynMux multislice r1","multinet GLouvain","Cross-sectional + Hungarian","Pooled Leiden")
mlab<-c("Jaccard","Overlap","multislice","multinet","Hungarian","Pooled")
A3<-A2<-matrix(NA,0,length(meth)); yrsused<-c(); nl<-nc<-no<-c()
for(y in 1946:1991){t<-which(yrs==y); if(!length(t))next; d<-deg[[t]]
 A<-P[[meth[1]]][[t]]
 alive<-names(d)[d>0]; alive<-alive[alive%in%names(A)]
 lib<-intersect(inb(LIB,y),alive); com<-intersect(inb(COM,y),alive)
 oth<-setdiff(alive,c(lib,com))
 if(length(lib)<3||length(com)<3) next
 n3<-c(lib,com,oth); t3<-c(rep(1,length(lib)),rep(2,length(com)),rep(3,length(oth)))
 n2<-c(lib,com);     t2<-c(rep(1,length(lib)),rep(2,length(com)))
 r3<-sapply(meth,function(m){Ai<-P[[m]][[t]]; if(!all(n3%in%names(Ai)))return(NA)
   igraph::compare(t3,as.integer(unlist(Ai[n3])),method="adjusted.rand")})
 r2<-sapply(meth,function(m){Ai<-P[[m]][[t]]; if(!all(n2%in%names(Ai)))return(NA)
   igraph::compare(t2,as.integer(unlist(Ai[n2])),method="adjusted.rand")})
 A3<-rbind(A3,r3); A2<-rbind(A2,r2); yrsused<-c(yrsused,y)
 nl<-c(nl,length(lib)); nc<-c(nc,length(com)); no<-c(no,length(oth))}
rownames(A3)<-rownames(A2)<-yrsused; colnames(A3)<-colnames(A2)<-mlab
cat("years scored:",length(yrsused),"(",min(yrsused),"-",max(yrsused),")\n")
cat("mean bloc sizes in ATOP-active nodes: Liberal",round(mean(nl),1)," Communist",round(mean(nc),1)," Other",round(mean(no),1),"\n\n")
cat("=== MEAN ARI vs Braumoeller/Goodhart blocs ===\n")
cat(sprintf("%-11s %14s %14s\n","method","3-bloc(L/C/O)","2-bloc(L vs C)"))
for(j in seq_along(mlab)) cat(sprintf("%-11s %14.3f %14.3f\n",mlab[j],mean(A3[,j],na.rm=TRUE),mean(A2[,j],na.rm=TRUE)))
cat("\n=== bloc purity: share of bloc in its single largest community ===\n")
pur<-function(m,ns,t){Ai<-P[[m]][[t]]; v<-as.integer(unlist(Ai[ns])); max(table(v))/length(v)}
cat(sprintf("%-11s %8s %8s\n","method","Liberal","Commun."))
for(j in seq_along(meth)){pl<-pc<-c()
 for(i in seq_along(yrsused)){y<-yrsused[i]; t<-which(yrs==y); d<-deg[[t]]; A<-P[[meth[j]]][[t]]
  alive<-names(d)[d>0]; alive<-alive[alive%in%names(A)]
  lib<-intersect(inb(LIB,y),alive); com<-intersect(inb(COM,y),alive)
  if(length(lib)>2) pl<-c(pl,pur(meth[j],lib,t)); if(length(com)>2) pc<-c(pc,pur(meth[j],com,t))}
 cat(sprintf("%-11s %8.3f %8.3f\n",mlab[j],mean(pl),mean(pc)))}
saveRDS(list(A3=A3,A2=A2,years=yrsused),"replication/extended/paper_scripts/braumoeller_validity.rds")
