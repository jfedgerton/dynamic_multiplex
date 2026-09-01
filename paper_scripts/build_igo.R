suppressMessages({library(peacesciencer);library(igraph)})
OUT<-Sys.getenv("EMP_DATA","replication/extended/output/empirical_data"); THRESH<-as.integer(Sys.getenv("IGO_THRESH","1"))
dy<-create_dyadyears(subset_years=1816:2014, directed=FALSE)
dy<-add_igos(dy)
dy$w<-ifelse(is.na(dy$dyadigos),0,dy$dyadigos)
E<-dy[dy$w>=THRESH, c("ccode1","ccode2","year","w")]
YEARS<-sort(unique(E$year)); cat("IGO years with data:",length(YEARS),"range",paste(range(YEARS),collapse="-"),"\n")
ucodes<-sort(unique(c(E$ccode1,E$ccode2))); N<-length(ucodes); idx<-setNames(seq_len(N),as.character(ucodes))
graph_layers<-vector("list",length(YEARS)); layers<-vector("list",length(YEARS)); active<-vector("list",length(YEARS))
for(k in seq_along(YEARS)){
 e<-E[E$year==YEARS[k],]
 A<-matrix(0,N,N,dimnames=list(as.character(ucodes),as.character(ucodes)))
 if(nrow(e)>0){ii<-idx[as.character(e$ccode1)];jj<-idx[as.character(e$ccode2)];A[cbind(ii,jj)]<-e$w;A[cbind(jj,ii)]<-e$w}
 layers[[k]]<-A; active[[k]]<-rowSums(A>0)>0
 g<-graph_from_adjacency_matrix(A,mode="undirected",weighted=TRUE,diag=FALSE); V(g)$name<-as.character(ucodes); graph_layers[[k]]<-g
}
saveRDS(list(years=YEARS,graph_layers=graph_layers),file.path(OUT,"igo_series.rds"))
saveRDS(list(layers=layers,active=active,names=as.character(ucodes)),file.path(OUT,"igo_union.rds"))
cat("BUILT igo | years",length(YEARS),"| union",N,"| mean active",round(mean(sapply(active,sum)),1),"| mean edges",round(mean(sapply(layers,function(m)sum(m>0)/2)),0),"| mean wt",round(mean(unlist(lapply(layers,function(m)m[m>0]))),1),"\n")
