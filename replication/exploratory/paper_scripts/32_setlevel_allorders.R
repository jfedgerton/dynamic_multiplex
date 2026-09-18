# Set-level order recovery, generalized to every coded order and every layer.
# For coded order B (states in year y) and detected communities C:
#   J = max_C |B n C|/|B u C| ; precision |BnC|/|C| ; recall |BnC|/|B|
# Needs only the target SET, so it works for partial/overlapping orders where ARI cannot.
suppressMessages(library(igraph)); set.seed(123)
setwd("/storage/work/jfe4/dynamic_multiplex")
src<-readLines("replication/extended/paper_scripts/29_braumoeller_fullspan.R")
eval(parse(text=paste(src[1:(grep("^run<-function",src)-1)],collapse="\n")))
CGP<-c("200"="1816:1852","220"="1816:1852","255"="1816:1852","300"="1816:1852","365"="1816:1852")
WP <-c("100"="1955:1991","265"="1955:1990","290"="1955:1991","310"="1955:1991","315"="1955:1991","339"="1955:1968","360"="1955:1991","365"="1955:1991")
MAN<-c("645"="1922:1939","652"="1922:1939","660"="1922:1939","663"="1922:1939","666"="1922:1939")
IAW<-c("205"="1921:1939","290"="1918:1939","305"="1918:1938","310"="1918:1939","315"="1918:1939","339"="1912:1938","345"="1919:1939","355"="1908:1938","366"="1918:1939","367"="1918:1939","368"="1918:1939","369"="1918:1920","371"="1918:1920","372"="1918:1921","373"="1918:1920","375"="1917:1939","651"="1922:1939","670"="1932:1938","678"="1918:1938","700"="1919:1938")
ORD<-list(
 list(id="Concert GPs",     y=c(1816,1852), kind="membership", f=function(y,al) intersect(inb(CGP,y),al)),
 list(id="Concert Europe",  y=c(1816,1852), kind="geographic", f=function(y,al) al[as.integer(al)>=200 & as.integer(al)<400]),
 list(id="Interim Europe",  y=c(1855,1870), kind="geographic", f=function(y,al) al[as.integer(al)>=200 & as.integer(al)<400]),
 list(id="Bismarck Europe", y=c(1871,1890), kind="geographic", f=function(y,al) al[as.integer(al)>=200 & as.integer(al)<400]),
 list(id="Wilhelm Europe",  y=c(1891,1914), kind="geographic", f=function(y,al) al[as.integer(al)>=200 & as.integer(al)<400]),
 list(id="Indep. after WWI",y=c(1908,1939), kind="membership", f=function(y,al) intersect(inb(IAW,y),al)),
 list(id="Mandates",        y=c(1922,1939), kind="membership", f=function(y,al) intersect(inb(MAN,y),al)),
 list(id="League",          y=c(1920,1941), kind="membership", f=function(y,al) intersect(inb(LEAGUE,y),al)),
 list(id="PW Liberal",      y=c(1945,1991), kind="membership", f=function(y,al) intersect(inb(LIB,y),al)),
 list(id="PW Communist",    y=c(1955,1991), kind="membership", f=function(y,al) intersect(inb(COM,y),al)),
 list(id="Warsaw Pact",     y=c(1955,1991), kind="membership", f=function(y,al) intersect(inb(WP,y),al)),
 list(id="PW Other",        y=c(1945,1991), kind="residual",   f=function(y,al) setdiff(al,c(inb(LIB,y),inb(COM,y)))),
 list(id="PostCW Liberal",  y=c(1992,2018), kind="membership", f=function(y,al) intersect(inb(PCW,y),al)))
meth<-c("DynMux Jaccard r1","DynMux Overlap r1","DynMux multislice r1","multinet GLouvain","Cross-sectional + Hungarian","Pooled Leiden")
mlab<-c("Jaccard","Overlap","multislice","multinet","Hungarian","Pooled")
best<-function(B,mem,nodes){ bs<-nodes %in% B; if(sum(bs)<3) return(NULL); out<-c(0,0,0)
 for(cc in unique(mem)){cs<-mem==cc; i<-sum(bs&cs); j<-i/sum(bs|cs); if(j>out[1]) out<-c(j,i/sum(cs),i/sum(bs))}
 c(out,sum(bs))}
rows<-list()
for(net in c("atop","igo","trade","dca")){
 S<-readRDS(sprintf("replication/extended/output/empirical_data/%s_series.rds",net))
 U<-readRDS(sprintf("replication/extended/output/empirical_data/%s_union.rds",net))
 P<-readRDS(sprintf("replication/extended/output/empirical/%s_partitions.rds",net))$partitions
 yrs<-S$years; mask<-if(!is.null(U$present))U$present else U$active
 G<-lapply(seq_along(S$graph_layers),function(k) delete_vertices(S$graph_layers[[k]],which(!mask[[k]])))
 dg<-lapply(G,function(g) setNames(degree(g),V(g)$name))
 for(o in ORD) for(y in o$y[1]:o$y[2]){t<-which(yrs==y); if(!length(t))next
  d<-dg[[t]]; A0<-P[[meth[1]]][[t]]; al<-names(d)[d>0]; al<-al[al%in%names(A0)]
  if(length(al)<8) next
  B<-o$f(y,al); if(length(B)<3) next
  for(j in seq_along(meth)){A<-P[[meth[j]]][[t]]; if(!all(al%in%names(A)))next
   r<-best(B,as.integer(unlist(A[al])),al); if(is.null(r))next
   rows[[length(rows)+1]]<-data.frame(net=net,order=o$id,kind=o$kind,year=y,method=mlab[j],
     J=r[1],prec=r[2],rec=r[3],nB=r[4],stringsAsFactors=FALSE)}}}
D<-do.call(rbind,rows)
saveRDS(D,"replication/extended/paper_scripts/setlevel_allorders.rds")
a<-aggregate(cbind(J,prec,rec,nB)~net+order+kind+method,D,mean)
n<-aggregate(year~net+order+method,D,length); a<-merge(a,n,by=c("net","order","method"))
write.csv(a,"/tmp/sl_all.csv",row.names=FALSE)
cat("obs:",nrow(D)," cells:",nrow(a),"\n"); print(table(D$net,D$order))
