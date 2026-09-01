coupling <- c("DynMux multislice (adjacent)","DynMux multislice (custom)","DynMux Jaccard","DynMux Overlap",
 "DynMux weighted Jaccard","DynMux weighted Overlap","Cross-sectional + Hungarian",
 "multinet GLouvain","DynMux multislice (adjacent, Louvain)","DynMux Jaccard (Louvain)",
 "DynMux Overlap (Louvain)","DynMux weighted Jaccard (Louvain)","DynMux weighted Overlap (Louvain)")
readdir<-"replication/extended/output/dynamic_prebdfix_backup"
writedir<-"replication/extended/output/dynamic"
fixdir<-"replication/extended/output/dynamic_bdfix"
ff<-list.files(fixdir,"csv$")
for(fn in ff){
  o<-read.csv(file.path(readdir,fn),check.names=FALSE)
  x<-read.csv(file.path(fixdir,fn),check.names=FALSE)
  stopifnot(all(o$regime=="birthdeath"))
  kept<-o[!(o$method %in% coupling),]
  m<-rbind(kept, x[,names(o)])
  write.csv(m, file.path(writedir,fn), row.names=FALSE)
}
cat("SPLICE DONE:",length(ff),"files\n")
o1<-read.csv(file.path(writedir,ff[1]),check.names=FALSE)
cat("example",ff[1],"rows:",nrow(o1)," methods:",length(unique(o1$method))," kept-nonfix:",paste(setdiff(unique(o1$method),coupling),collapse=";"),"\n")
