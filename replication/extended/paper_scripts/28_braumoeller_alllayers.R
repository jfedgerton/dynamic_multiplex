# Braumoeller/Goodhart Cold War bloc validation across all four empirical networks.
suppressMessages(library(igraph)); set.seed(123)
setwd("/storage/work/jfe4/dynamic_multiplex")
LIB <- c("2"="1945:1991","20"="1945:1991","200"="1945:1991","205"="1945:1991","210"="1945:1991","211"="1945:1991","212"="1945:1991","220"="1945:1991","225"="1945:1991","230"="1982:1991","235"="1945:1991","255"="1990:1991","260"="1955:1990","305"="1945:1991","325"="1945:1991","350"="1952:1991","380"="1945:1991","385"="1945:1991","390"="1945:1991","395"="1945:1991","640"="1952:1991","666"="1948:1991","713"="1947:1991","732"="1948:1991","740"="1952:1991","900"="1945:1991","920"="1945:1991")
COM <- c("265"="1955:1990","290"="1955:1991","310"="1955:1991","315"="1955:1991","339"="1955:1991","355"="1955:1991","359"="1945:1991","360"="1955:1991","365"="1955:1991","366"="1945:1991","367"="1945:1991","368"="1945:1991","369"="1945:1991","370"="1945:1991","371"="1945:1991","372"="1945:1991","373"="1945:1991","701"="1945:1991","702"="1945:1991","703"="1945:1991","704"="1945:1991","705"="1945:1991","710"="1949:1991","712"="1945:1991","731"="1948:1991","812"="1949:1991","816"="1945:1991")
inb <- function(tbl,y){ ks<-names(tbl); ks[sapply(tbl,function(r){p<-as.integer(strsplit(r,":")[[1]]); y>=p[1] && y<=p[2]})] }
meth<-c("DynMux Jaccard r1","DynMux Overlap r1","DynMux multislice r1","multinet GLouvain","Cross-sectional + Hungarian","Pooled Leiden")
mlab<-c("Jaccard","Overlap","multislice","multinet","Hungarian","Pooled")
run<-function(net){
 S<-readRDS(sprintf("replication/extended/output/empirical_data/%s_series.rds",net))
 U<-readRDS(sprintf("replication/extended/output/empirical_data/%s_union.rds",net))
 P<-readRDS(sprintf("replication/extended/output/empirical/%s_partitions.rds",net))$partitions
 yrs<-S$years; mask<-if(!is.null(U$present))U$present else U$active
 G<-lapply(seq_along(S$graph_layers),function(k) delete_vertices(S$graph_layers[[k]],which(!mask[[k]])))
 deg<-lapply(G,function(g) setNames(degree(g),V(g)$name))
 A3<-A2<-matrix(NA,0,length(meth)); ys<-c()
 for(y in 1946:1991){t<-which(yrs==y); if(!length(t))next; d<-deg[[t]]; A<-P[[meth[1]]][[t]]
  alive<-names(d)[d>0]; alive<-alive[alive%in%names(A)]
  lib<-intersect(inb(LIB,y),alive); com<-intersect(inb(COM,y),alive); oth<-setdiff(alive,c(lib,com))
  if(length(lib)<3||length(com)<3) next
  n3<-c(lib,com,oth); t3<-c(rep(1,length(lib)),rep(2,length(com)),rep(3,length(oth)))
  n2<-c(lib,com);     t2<-c(rep(1,length(lib)),rep(2,length(com)))
  f<-function(ns,tr) sapply(meth,function(m){Ai<-P[[m]][[t]]; if(!all(ns%in%names(Ai)))return(NA)
    igraph::compare(tr,as.integer(unlist(Ai[ns])),method="adjusted.rand")})
  A3<-rbind(A3,f(n3,t3)); A2<-rbind(A2,f(n2,t2)); ys<-c(ys,y)}
 if(!length(ys)){cat("\n###",net,": no scorable years\n"); return(invisible(NULL))}
 cat(sprintf("\n### %s  (%d years, %d-%d)\n", toupper(net), length(ys), min(ys), max(ys)))
 cat(sprintf("%-11s %14s %14s\n","method","3-bloc","2-bloc(L vs C)"))
 for(j in seq_along(mlab)) cat(sprintf("%-11s %14.3f %14.3f\n",mlab[j],mean(A3[,j],na.rm=TRUE),mean(A2[,j],na.rm=TRUE)))}
for(n in c("atop","igo","trade","dca")) run(n)
