args<-commandArgs(TRUE); NET<-if(length(args)>=1)args[1] else "atop"
suppressMessages(library(igraph))
DATA<-"replication/extended/output/empirical_data"; OUT<-"replication/extended/output/empirical"
S<-readRDS(file.path(DATA,paste0(NET,"_series.rds")))
P<-readRDS(file.path(OUT,paste0(NET,"_partitions.rds")))$partitions
yrs<-S$years; ymin<-min(yrs); ymax<-max(yrs)
deg<-lapply(seq_along(yrs),function(t){g<-S$graph_layers[[t]]; if(is.null(g))return(setNames(integer(0),character(0))); setNames(igraph::degree(g),V(g)$name)})
nm<-c("2"="USA","20"="CAN","200"="UK","220"="FRA","255"="GER","260"="FRG","265"="GDR","300"="AUH","325"="ITA","365"="RUS","290"="POL","316"="CZR","310"="HUN","345"="YUG","375"="FIN","710"="CHN","740"="JPN","731"="PRK")
eras<-list(Concert=c(1816,1853),Bismarck=c(1854,1890),WWI=c(1891,1918),Interwar=c(1919,1938),WWII=c(1939,1945),ColdWar=c(1946,1989),LIO=c(1990,2018))
ab<-c(Jaccard="DynMux Jaccard r1",Overlap="DynMux Overlap r1",multislice="DynMux multislice r1",multinet="multinet GLouvain",Hungarian="Cross-sectional + Hungarian",Pooled="Pooled Leiden")
D<-c("Concert|365|300|+","Concert|300|255|+","Concert|365|255|+","Concert|220|365|-","Concert|200|365|-","Bismarck|255|300|+","Bismarck|255|325|+","Bismarck|220|255|-","Bismarck|300|365|-","WWI|255|300|+","WWI|220|365|+","WWI|200|220|+","WWI|200|740|+","WWI|255|220|-","WWI|255|200|-","WWI|255|365|-","WWI|325|200|-","Interwar|220|290|+","Interwar|255|325|+","Interwar|220|365|+","Interwar|220|255|-","Interwar|255|365|-","WWII|2|365|+","WWII|2|200|+","WWII|255|740|+","WWII|255|325|+","WWII|255|365|-","WWII|740|2|-","WWII|325|200|-","WWII|375|365|-","ColdWar|2|740|+","ColdWar|2|260|+","ColdWar|365|290|+","ColdWar|365|265|+","ColdWar|710|731|+","ColdWar|365|731|+","ColdWar|710|365|T","ColdWar|2|365|-","ColdWar|345|365|-","ColdWar|265|260|-","ColdWar|2|731|-","ColdWar|740|731|-","LIO|2|290|+","LIO|2|310|+","LIO|2|255|+","LIO|2|740|+","LIO|710|731|+","LIO|290|365|-","LIO|310|365|-","LIO|316|365|-","LIO|2|731|-","LIO|740|731|-")
sm<-function(m,a,b,t){A<-P[[ab[m]]][[t]]; na<-as.character(a); nb<-as.character(b); if(is.null(A)||!(na%in%names(A))||!(nb%in%names(A)))return(NA); A[[na]]==A[[nb]]}
act<-function(a,b,t){d<-deg[[t]]; na<-as.character(a); nb<-as.character(b); length(d)>0&&na%in%names(d)&&nb%in%names(d)&&d[[na]]>0&&d[[nb]]>0}
expsign<-function(a,b,ty,y){ if(a==710&&b==365) return(if(y<=1960)"+" else "-"); ty }
con<-file(paste0("paper_scripts/grid_",NET,".txt"),"w")
for(er in names(eras)){ e0<-max(eras[[er]][1],ymin); e1<-min(eras[[er]][2],ymax); if(e0>e1) next
 ds<-D[startsWith(D,paste0(er,"|"))]; hdr<-FALSE
 for(d in ds){p<-strsplit(d,"\\|")[[1]]; a<-as.integer(p[2]); b<-as.integer(p[3]); ty<-p[4]
  tyc<-if(ty=="T")"+/-" else ty; lab<-paste0(nm[p[2]],"-",nm[p[3]]," (",tyc,")")
  anyact<-FALSE; strs<-list()
  for(m in names(ab)){s<-character(0)
   for(y in e0:e1){t<-which(yrs==y)
    if(length(t)!=1||!act(a,b,t)){s<-c(s,"-");next}
    anyact<-TRUE; co<-sm(m,a,b,t); es<-expsign(a,b,ty,y)
    ok<-if(is.na(co))NA else ((es=="+"&&co)||(es=="-"&&!co))
    s<-c(s,if(is.na(ok))"-" else if(ok)"o" else "x")}
   strs[[m]]<-paste(s,collapse="")}
  if(!anyact)next
  if(!hdr){cat(sprintf("#ERA %s %d %d\n",er,e0,e1),file=con);hdr<-TRUE}
  for(m in names(ab)) cat(sprintf("%s %s %s\n",lab,m,strs[[m]]),file=con)}}
close(con); cat("wrote grid_",NET,".txt\n",sep="")
