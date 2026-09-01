DATA<-"replication/extended/output/empirical_data"; OUT<-"replication/extended/output/empirical"; net<-"dca"
S<-readRDS(file.path(DATA,sprintf("%s_series.rds",net))); U<-readRDS(file.path(DATA,sprintf("%s_union.rds",net))); P<-readRDS(file.path(OUT,sprintf("%s_partitions.rds",net)))$partitions; yrs<-S$years
mask<-if(!is.null(U$present))U$present else U$active
G<-lapply(seq_along(S$graph_layers),function(k) igraph::delete_vertices(S$graph_layers[[k]],which(!mask[[k]])))
deg<-lapply(G,function(g) setNames(igraph::degree(g),igraph::V(g)$name))
cat("DCA years",min(yrs),max(yrs),"layers",length(yrs),"methods",length(P),"\n")
t0<-which(yrs==2000); if(length(t0)){cat("sample nodes 2000:",paste(head(sort(names(P[[1]][[t0]])),15),collapse=" "),"\n")}
eras<-list(ColdWar=c(1980,1989),LIO=c(1990,2010))
D<-c("ColdWar|2|740|+","ColdWar|2|260|+","ColdWar|365|290|+","ColdWar|365|265|+","ColdWar|2|365|-","ColdWar|710|365|-","ColdWar|345|365|-","ColdWar|265|260|-","LIO|2|290|+","LIO|2|310|+","LIO|2|255|+","LIO|290|365|-","LIO|310|365|-","LIO|316|365|-")
meth<-intersect(c("DynMux Jaccard r1","DynMux Overlap r1","DynMux multislice r1","multinet GLouvain","Cross-sectional + Hungarian","Pooled Leiden"),names(P))
mlab<-c("DynMux Jaccard r1"="Jaccard","DynMux Overlap r1"="Overlap","DynMux multislice r1"="multislice","multinet GLouvain"="multinet","Cross-sectional + Hungarian"="Hungarian","Pooled Leiden"="Pooled")
frac<-function(m,a,b,y0,y1){co<-c();for(y in y0:y1){t<-which(yrs==y);if(!length(t))next;A<-P[[m]][[t]];dg<-deg[[t]];if(!(a%in%names(A)&&b%in%names(A)&&a%in%names(dg)&&b%in%names(dg)))next;if(dg[[a]]<=0||dg[[b]]<=0)next;co<-c(co,as.integer(A[[a]]==A[[b]]))};if(length(co)<1)return(NA);mean(co)}
rows<-list()
for(d in D){p<-strsplit(d,"\\|")[[1]];er<-p[1];a<-p[2];b<-p[3];ty<-p[4];w<-eras[[er]];for(m in meth){f<-frac(m,a,b,w[1],w[2]);rows[[length(rows)+1]]<-data.frame(era=er,type=ty,method=unname(mlab[m]),frac=f,stringsAsFactors=FALSE)}}
FR<-do.call(rbind,rows); d2<-FR[!is.na(FR$frac),]
cat("\n== DCA continuous balanced accuracy ==\n"); cat(sprintf("%-11s %5s %8s %8s %8s\n","method","n","align","separate","balAcc"))
for(m in unname(mlab[meth])){s<-d2[d2$method==m,];al<-mean(s$frac[s$type=="+"]);se<-mean(1-s$frac[s$type=="-"]);cat(sprintf("%-11s %5d %8.2f %8.2f %8.2f\n",m,nrow(s),al,se,(al+se)/2))}
nmm<-c("2"="USA","740"="JPN","260"="FRG","365"="RUS","290"="POL","265"="GDR","710"="CHN","345"="YUG","255"="GER","310"="HUN","316"="CZR")
cat("\n== co-cluster fraction per dyad (NA = states not both in DCA) ==\n")
for(d in D){p<-strsplit(d,"\\|")[[1]];er<-p[1];a<-p[2];b<-p[3];ty<-p[4];w<-eras[[er]];fj<-frac("DynMux Jaccard r1",a,b,w[1],w[2]);fm<-frac("DynMux multislice r1",a,b,w[1],w[2]);fp<-frac("Pooled Leiden",a,b,w[1],w[2]);cat(sprintf("%-8s %s-%s(%s) J=%s ms=%s P=%s\n",er,ifelse(a%in%names(nmm),nmm[a],a),ifelse(b%in%names(nmm),nmm[b],b),ty,ifelse(is.na(fj),"NA",sprintf("%.2f",fj)),ifelse(is.na(fm),"NA",sprintf("%.2f",fm)),ifelse(is.na(fp),"NA",sprintf("%.2f",fp))))}
