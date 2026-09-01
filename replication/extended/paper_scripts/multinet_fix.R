# multinet_fix.R <atop|dca> : recompute ONLY multinet GLouvain, active-per-layer.
# Nodes drop in/out (actors added only in years they are active); present-but-
# isolate states get singletons post hoc; scored on the present mask.
suppressMessages(library(igraph))
args<-commandArgs(trailingOnly=TRUE); net<-if(length(args)) args[1] else "atop"
DATA<-Sys.getenv("EMP_DATA","/tmp"); OUT<-Sys.getenv("EMP_OUT","/tmp")
stopifnot(requireNamespace("multinet",quietly=TRUE))
U<-readRDS(file.path(DATA,sprintf("%s_union.rds",net)))
L<-U$layers; ac<-U$present; n<-nrow(L[[1]]); T_<-length(L); nm<-as.character(seq_len(n))
nn<-multinet::ml_empty()
for(t in seq_len(T_)){
  keep<-which(ac[[t]]); if(length(keep)==0) next
  A<-L[[t]][keep,keep,drop=FALSE]
  g<-igraph::graph_from_adjacency_matrix(A,mode="undirected",diag=FALSE)
  igraph::V(g)$name<-nm[keep]
  multinet::add_igraph_layer_ml(nn,g,paste0("layer",t))
}
cat("actor-layer entries:",sum(sapply(ac,sum))," (vs padded",n*T_,")\n")
cm<-multinet::glouvain_ml(nn,gamma=1,omega=1)
d<-lapply(seq_len(T_),function(t){
  s<-cm[cm$layer==paste0("layer",t),,drop=FALSE]
  v<-rep(NA_integer_,n)
  if(nrow(s)>0) v[as.integer(s$actor)]<-as.integer(as.factor(s$cid))
  if(anyNA(v)) v[is.na(v)]<-max(0L,v,na.rm=TRUE)+seq_len(sum(is.na(v)))
  v
})
part<-lapply(seq_len(T_),function(t){a<-U$present[[t]];setNames(as.integer(d[[t]][a]),U$names[a])})
outf<-file.path(OUT,sprintf("%s_partitions.rds",net))
R<-readRDS(outf); R$partitions[["multinet GLouvain"]]<-part; saveRDS(R,outf)
cat("multinet refit active-per-layer | net=",net,"| layers=",T_,"| ncomm=",length(unique(unlist(part))),"\n")
