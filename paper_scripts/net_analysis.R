suppressMessages(library(igraph)); set.seed(123)
args<-commandArgs(trailingOnly=TRUE); net<-if(length(args)>=1)args[1] else "atop"
DATA<-"replication/extended/output/empirical_data"; OUT<-"replication/extended/output/empirical"
S<-readRDS(file.path(DATA,sprintf("%s_series.rds",net))); U<-readRDS(file.path(DATA,sprintf("%s_union.rds",net)))
P<-readRDS(file.path(OUT,sprintf("%s_partitions.rds",net)))$partitions; yrs<-S$years
mask<-if(!is.null(U$present))U$present else U$active
G<-lapply(seq_along(S$graph_layers),function(k) igraph::delete_vertices(S$graph_layers[[k]],which(!mask[[k]])))
deg<-lapply(G,function(g) setNames(igraph::degree(g),V(g)$name)); Tn<-length(G)
NB<-lapply(G,function(g){nm<-V(g)$name; al<-igraph::as_adj_list(g,mode="all"); setNames(lapply(al,function(ids) nm[as.integer(ids)]),nm)})
jc<-function(x,y){u<-length(union(x,y)); if(u==0)return(1); length(intersect(x,y))/u}
errfun<-function(ml){sp<-ms<-n<-0; for(t in 1:(Tn-1)){a<-ml[[t]];b<-ml[[t+1]]; if(is.null(a)||is.null(b))next; cm<-intersect(names(a),names(b)); if(length(cm)<3)next; la<-a[cm];lb<-b[cm]; gA<-split(cm,as.vector(la));gB<-split(cm,as.vector(lb)); for(i in cm){nti<-intersect(NB[[t]][[i]],cm);ntj<-intersect(NB[[t+1]][[i]],cm); if(length(nti)<1||length(ntj)<1)next; J<-jc(nti,ntj); Ca<-setdiff(gA[[as.character(la[i])]],i);Cb<-setdiff(gB[[as.character(lb[i])]],i); coJ<-jc(Ca,Cb); n<-n+1; if(J>=0.9&&coJ<=0.5)sp<-sp+1; if(J<=0.3&&coJ>=0.9)ms<-ms+1}}; c(sp=sp,ms=ms,n=n)}
meth<-intersect(c("DynMux Jaccard r1","DynMux Overlap r1","DynMux multislice r1","multinet GLouvain","Cross-sectional + Hungarian","Pooled Leiden"),names(P))
cat("=====",net,"TEMPORAL ALIGNMENT ERROR =====\n"); cat(sprintf("%-14s %6s %6s %6s %7s\n","method","spur","miss","total","rate%"))
for(m in meth){r<-errfun(P[[m]]); cat(sprintf("%-14s %6d %6d %6d %7.2f\n",sub("DynMux ","",m),r["sp"],r["ms"],r["sp"]+r["ms"],100*(r["sp"]+r["ms"])/max(r["n"],1)))}
West<-c("2","20","200","220","260","325","211","210","212","385","390","395","235","350","640","230","740"); East<-c("365","290","265","315","310","360","355","712","40")
ac<-function(t,c,a,dg)(c%in%names(a))&&(c%in%names(dg))&&dg[[c]]>0
cls<-function(m){inc<-mis<-cor<-0; for(y in 1955:1989){t<-which(yrs==y); if(!length(t))next; a<-P[[m]][[t]];dg<-deg[[t]]; Aa<-West[sapply(West,function(c)ac(t,c,a,dg))];Ba<-East[sapply(East,function(c)ac(t,c,a,dg))]; if(length(Aa)<2||length(Ba)<2)next; cA<-names(sort(table(a[Aa]),dec=TRUE))[1];cB<-names(sort(table(a[Ba]),dec=TRUE))[1];mg<-cA==cB; for(s in c(Aa,Ba)){own<-if(s%in%Aa)cA else cB;riv<-if(s%in%Aa)cB else cA;cs<-as.character(a[[s]]); if(mg||cs==riv)inc<-inc+1 else if(cs==own)cor<-cor+1 else mis<-mis+1}}; c(inc=inc,mis=mis,tot=inc+mis+cor)}
cat("\n=====",net,"COLD WAR West-vs-East institutional recovery (1955-1989) =====\n"); cat(sprintf("%-14s %10s %8s %6s\n","method","incorrect","missing","n"))
for(m in meth){r<-cls(m); if(r["tot"]>0)cat(sprintf("%-14s %6d(%4.0f%%) %8d %6d\n",sub("DynMux ","",m),r["inc"],100*r["inc"]/r["tot"],r["mis"],r["tot"]))}
gp<-c("2"="USA","200"="UK","220"="FRA","255"="GER","300"="AUH","325"="ITA","365"="RUS","710"="CHN","740"="JPN","260"="WGE","640"="TUR","290"="POL","345"="YUG")
JAC<-if("DynMux Jaccard r1"%in%names(P))"DynMux Jaccard r1" else meth[1]
cat("\n=====",net,"great-power groupings (",JAC,") =====\n")
for(y in c(1875,1910,1940,1960,1985,2000,2010)){t<-which(yrs==y); if(!length(t))next; a<-P[[JAC]][[t]];dg<-deg[[t]]; pres<-names(gp)[sapply(names(gp),function(c)ac(t,c,a,dg))]; if(length(pres)<2){cat(y,": <2 active\n");next}; grps<-split(pres,sapply(pres,function(c)a[[c]])); grps<-grps[order(-sapply(grps,length))]; cat(sprintf("%d: %s\n",y,paste(sapply(grps,function(g)paste(gp[g],collapse="+")),collapse="  |  ")))}
