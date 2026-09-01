suppressMessages(library(igraph))
B<-"/storage/work/jfe4/dynamic_multiplex"
S<-readRDS(file.path(B,"replication/extended/output/empirical_data/dca_series.rds"))
P<-readRDS(file.path(B,"replication/extended/output/empirical/dca_partitions.rds"))$partitions
yrs<-S$years; cat("DCA years:",min(yrs),"-",max(yrs),"(n=",length(yrs),")\n")
mk<-"multinet GLouvain"
deg<-function(t){g<-S$graph_layers[[t]]; if(is.null(g))return(setNames(integer(0),character(0))); setNames(igraph::degree(g),V(g)$name)}
cat("\nyr nAct nComm maxSz  USA POL HUN GER RUS JPN ; RUSdeg\n")
for(y in c(1990,1995,2000,2005,2010)){t<-which(yrs==y); if(length(t)!=1)next
 A<-P[[mk]][[t]]; d<-deg(t); sizes<-sort(table(A),decreasing=TRUE)
 g<-function(cc){cc<-as.character(cc); if(!is.null(A)&&cc%in%names(A))A[[cc]] else NA}
 rdeg<-if("365"%in%names(d))d[["365"]] else 0
 cat(sprintf("%d %4d %4d %5d   %s %s %s %s %s %s ; RUSdeg=%d\n",y,sum(d>0),length(sizes),as.integer(sizes[1]),g(2),g(290),g(310),g(255),g(365),g(740),rdeg))}
nc<-c(); ls<-c()
for(t in seq_along(yrs)){A<-P[[mk]][[t]]; if(is.null(A)||length(A)==0)next; tb<-table(A); nc<-c(nc,length(tb)); ls<-c(ls,max(tb)/sum(tb))}
cat("\nmultinet mean nComm:",round(mean(nc),1),"| mean largest-comm share:",round(mean(ls),2),"\n")
