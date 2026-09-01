suppressMessages(library(igraph))
B<-"/storage/work/jfe4/dynamic_multiplex"
S<-readRDS(file.path(B,"replication/extended/output/empirical_data/dca_series.rds"))
P<-readRDS(file.path(B,"replication/extended/output/empirical/dca_partitions.rds"))$partitions
yrs<-S$years
for(mk in c("Pooled Leiden","multinet GLouvain")){
 cat("\n==== ",mk," ====\n"); cat("yr nComm maxSz  USA POL HUN GER RUS JPN\n")
 for(y in c(1995,2005,2010)){t<-which(yrs==y); A<-P[[mk]][[t]]; sizes<-sort(table(A),decreasing=TRUE)
  g<-function(cc){cc<-as.character(cc); if(!is.null(A)&&cc%in%names(A))A[[cc]] else NA}
  cat(sprintf("%d %4d %5d   %s %s %s %s %s %s\n",y,length(sizes),as.integer(sizes[1]),g(2),g(290),g(310),g(255),g(365),g(740)))}
 # is USA in same comm as POL/HUN/GER in aggregate? show comm sizes for 2005
 t<-which(yrs==2005); A<-P[[mk]][[t]]; cat("2005 comm sizes:",paste(sort(as.integer(table(A)),decreasing=TRUE),collapse=","),"\n")}
