suppressMessages({library(peacesciencer);library(igraph)})
OUT<-Sys.getenv("EMP_DATA","replication/extended/output/empirical_data")
dy<-create_dyadyears(subset_years=1870:2014, directed=FALSE)
dy<-add_cow_trade(dy)
f1<-ifelse(is.na(dy$flow1),0,dy$flow1); f2<-ifelse(is.na(dy$flow2),0,dy$flow2); dy$tot<-f1+f2
YEARS0<-sort(unique(dy$year)); keep<-rep(FALSE,nrow(dy))
for(y in YEARS0){i<-which(dy$year==y & dy$tot>0); if(length(i)<4)next; thr<-as.numeric(quantile(dy$tot[i],0.75)); keep[i[dy$tot[i]>=thr]]<-TRUE}
E<-dy[keep,c("ccode1","ccode2","year")]; E$w<-log(dy$tot[keep]+1)
YEARS<-sort(unique(E$year)); ucodes<-sort(unique(c(E$ccode1,E$ccode2))); N<-length(ucodes); idx<-setNames(seq_len(N),as.character(ucodes))
graph_layers<-vector("list",length(YEARS)); layers<-vector("list",length(YEARS)); active<-vector("list",length(YEARS))
for(k in seq_along(YEARS)){e<-E[E$year==YEARS[k],]
 A<-matrix(0,N,N,dimnames=list(as.character(ucodes),as.character(ucodes)))
 if(nrow(e)>0){ii<-idx[as.character(e$ccode1)];jj<-idx[as.character(e$ccode2)];A[cbind(ii,jj)]<-e$w;A[cbind(jj,ii)]<-e$w}
 layers[[k]]<-A; active[[k]]<-rowSums(A>0)>0
 g<-graph_from_adjacency_matrix(A,mode="undirected",weighted=TRUE,diag=FALSE); V(g)$name<-as.character(ucodes); graph_layers[[k]]<-g}
saveRDS(list(years=YEARS,graph_layers=graph_layers),file.path(OUT,"trade_series.rds"))
saveRDS(list(layers=layers,active=active,names=as.character(ucodes)),file.path(OUT,"trade_union.rds"))
cat("BUILT trade | years",length(YEARS),"| range",paste(range(YEARS),collapse="-"),"| union",N,"| mean active",round(mean(sapply(active,sum)),1),"| mean edges",round(mean(sapply(layers,function(m)sum(m>0)/2)),0),"\n")
