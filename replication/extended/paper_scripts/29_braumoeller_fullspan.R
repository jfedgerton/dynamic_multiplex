# Full-span external validation, 1816-2018, against Braumoeller/Goodhart order codings
# (Only the Dead replication archive; coding credited to Andy Goodhart).
# Each era scored against whatever partition the coding supports.
# NOTE: 19th-c orders are coded as "all European states (ccode 200-399)", so those
# eras test a Europe/non-Europe split, not a fine bloc structure. Interpret accordingly.
# Colombia/Bulgaria coding bug (ccode 100 labelled Bulgaria) corrected here.
suppressMessages(library(igraph)); set.seed(123)
setwd("/storage/work/jfe4/dynamic_multiplex")
LIB<-c("2"="1945:1991","20"="1945:1991","200"="1945:1991","205"="1945:1991","210"="1945:1991","211"="1945:1991","212"="1945:1991","220"="1945:1991","225"="1945:1991","230"="1982:1991","235"="1945:1991","255"="1990:1991","260"="1955:1990","305"="1945:1991","325"="1945:1991","350"="1952:1991","380"="1945:1991","385"="1945:1991","390"="1945:1991","395"="1945:1991","640"="1952:1991","666"="1948:1991","713"="1947:1991","732"="1948:1991","740"="1952:1991","900"="1945:1991","920"="1945:1991")
COM<-c("265"="1955:1990","290"="1955:1991","310"="1955:1991","315"="1955:1991","339"="1955:1991","355"="1955:1991","359"="1945:1991","360"="1955:1991","365"="1955:1991","366"="1945:1991","367"="1945:1991","368"="1945:1991","369"="1945:1991","370"="1945:1991","371"="1945:1991","372"="1945:1991","373"="1945:1991","701"="1945:1991","702"="1945:1991","703"="1945:1991","704"="1945:1991","705"="1945:1991","710"="1949:1991","712"="1945:1991","731"="1948:1991","812"="1949:1991","816"="1945:1991")
LATINAM<-c("70"="1821:1856","90"="1821:1856","91"="1821:1856","92"="1821:1856","93"="1821:1856","94"="1821:1856","95"="1821:1856","100"="1820:1856","101"="1821:1856","130"="1822:1856","135"="1821:1856","140"="1822:1856","145"="1825:1856","150"="1816:1856","155"="1816:1856","160"="1816:1856","165"="1828:1856")
LEAGUE<-c("20"="1920:1939","40"="1920:1939","41"="1920:1939","42"="1924:1939","70"="1931:1939","90"="1920:1936","91"="1920:1936","92"="1920:1937","93"="1920:1936","94"="1920:1925","95"="1920:1939","100"="1920:1939","101"="1920:1938","130"="1934:1939","135"="1920:1939","140"="1920:1926","145"="1920:1939","150"="1920:1935","155"="1920:1938","160"="1920:1939","165"="1920:1939","200"="1920:1939","205"="1921:1939","210"="1920:1939","211"="1920:1939","212"="1920:1939","220"="1920:1941","225"="1920:1939","230"="1920:1939","235"="1920:1939","255"="1926:1933","290"="1920:1939","305"="1920:1938","310"="1922:1939","315"="1920:1939","325"="1920:1937","339"="1920:1939","345"="1920:1939","350"="1920:1939","355"="1920:1939","360"="1920:1940","365"="1934:1939","366"="1921:1939","367"="1921:1939","368"="1921:1939","375"="1920:1939","380"="1920:1939","385"="1920:1939","390"="1920:1939","450"="1920:1939","530"="1923:1936","560"="1920:1939","630"="1920:1939","640"="1932:1939","645"="1932:1939","651"="1937:1939","700"="1934:1939","710"="1920:1939","740"="1920:1933","750"="1920:1939","800"="1920:1939","900"="1920:1939","920"="1920:1939")
PCW<-c("2"="1992:2018","20"="1992:2018","200"="1992:2018","205"="1992:2018","210"="1992:2018","211"="1992:2018","212"="1992:2018","220"="1992:2018","225"="1992:2018","230"="1992:2018","235"="1992:2018","255"="1992:2018","260"="1992:2018","290"="2000:2018","305"="1992:2018","310"="2000:2018","316"="2000:2018","317"="2005:2018","325"="1992:2018","338"="2005:2018","339"="2010:2018","344"="2010:2018","349"="2005:2018","350"="1992:2018","352"="2005:2018","355"="2005:2018","360"="2005:2018","366"="2005:2018","367"="2005:2018","368"="2005:2018","375"="1995:2018","380"="1992:2018","385"="1992:2018","390"="1992:2018","395"="1992:2018","640"="1992:2018","666"="1992:2018","713"="1992:2018","732"="1992:2018","740"="1992:2018","900"="1992:2018","920"="1992:2018")
ari<-function(a,b){tb<-table(a,b); n<-sum(tb); if(n<2) return(NA)
 ch2<-function(x) x*(x-1)/2
 idx<-sum(ch2(tb)); ea<-sum(ch2(rowSums(tb))); eb<-sum(ch2(colSums(tb)))
 exp<-ea*eb/ch2(n); mx<-(ea+eb)/2
 if(mx==exp) return(0); (idx-exp)/(mx-exp)}
inb<-function(tb,y){k<-names(tb); k[sapply(tb,function(r){p<-as.integer(strsplit(r,":")[[1]]); y>=p[1]&&y<=p[2]})]}
eur<-function(al) al[as.integer(al)>=200 & as.integer(al)<400]
truth<-function(y,al){
 if(y>=1816&&y<=1852){e<-eur(al); l<-intersect(inb(LATINAM,y),al); o<-setdiff(al,c(e,l)); if(length(e)<3||length(l)<3)return(NULL); list(era="Concert 1816-52",n=c(e,l,o),t=c(rep(1,length(e)),rep(2,length(l)),rep(3,length(o))))}
 else if(y>=1855&&y<=1870){e<-eur(al); o<-setdiff(al,e); if(length(e)<3||length(o)<3)return(NULL); list(era="Interim 1855-70",n=c(e,o),t=c(rep(1,length(e)),rep(2,length(o))))}
 else if(y>=1871&&y<=1890){e<-eur(al); o<-setdiff(al,e); if(length(e)<3||length(o)<3)return(NULL); list(era="Bismarck 1871-90",n=c(e,o),t=c(rep(1,length(e)),rep(2,length(o))))}
 else if(y>=1891&&y<=1914){e<-eur(al); o<-setdiff(al,e); if(length(e)<3||length(o)<3)return(NULL); list(era="Wilhelm 1891-1914",n=c(e,o),t=c(rep(1,length(e)),rep(2,length(o))))}
 else if(y>=1920&&y<=1941){g<-intersect(inb(LEAGUE,y),al); o<-setdiff(al,g); if(length(g)<3||length(o)<3)return(NULL); list(era="League 1920-41",n=c(g,o),t=c(rep(1,length(g)),rep(2,length(o))))}
 else if(y>=1945&&y<=1991){l<-intersect(inb(LIB,y),al); c2<-intersect(inb(COM,y),al); o<-setdiff(al,c(l,c2)); if(length(l)<3||length(c2)<3)return(NULL); list(era="ColdWar 1945-91",n=c(l,c2,o),t=c(rep(1,length(l)),rep(2,length(c2)),rep(3,length(o))))}
 else if(y>=1992){p<-intersect(inb(PCW,y),al); o<-setdiff(al,p); if(length(p)<3||length(o)<3)return(NULL); list(era="PostCW 1992-2018",n=c(p,o),t=c(rep(1,length(p)),rep(2,length(o))))}
 else NULL}
run<-function(net){
 S<-readRDS(sprintf("replication/extended/output/empirical_data/%s_series.rds",net))
 U<-readRDS(sprintf("replication/extended/output/empirical_data/%s_union.rds",net))
 P<-readRDS(sprintf("replication/extended/output/empirical/%s_partitions.rds",net))$partitions
 mm<-names(P); yrs<-S$years; mask<-if(!is.null(U$present))U$present else U$active
 G<-lapply(seq_along(S$graph_layers),function(k) delete_vertices(S$graph_layers[[k]],which(!mask[[k]])))
 dg<-lapply(G,function(g) setNames(degree(g),V(g)$name))
 rows<-list()
 for(y in yrs){t<-which(yrs==y); d<-dg[[t]]; A<-P[[mm[1]]][[t]]
  al<-names(d)[d>0]; al<-al[al%in%names(A)]; if(length(al)<8)next
  tr<-truth(y,al); if(is.null(tr))next
  r<-sapply(mm,function(m){Ai<-P[[m]][[t]]; if(!all(tr$n%in%names(Ai)))return(NA)
    ari(tr$t, unlist(Ai[tr$n]))})
  rows[[length(rows)+1]]<-data.frame(net=net,year=y,era=tr$era,t(r),check.names=FALSE)}
 do.call(rbind,rows)}
all<-do.call(rbind,lapply(c("atop","igo","trade"),run))
saveRDS(all,"replication/extended/paper_scripts/braumoeller_fullspan.rds")
mm<-setdiff(names(all),c("net","year","era"))
for(nt in c("atop","igo","trade")){s<-all[all$net==nt,]
 cat("\n########",toupper(nt),"- years scored:",nrow(s),"\n")
 agg<-aggregate(s[,mm],list(era=s$era),function(x) round(mean(x,na.rm=TRUE),3))
 print(agg[,c("era",mm[c(4,7,1,10,11,12)])])}
cat("\n######## OVERALL MEAN ACROSS ALL ERAS ########\n")
for(nt in c("atop","igo","trade")){s<-all[all$net==nt,]
 cat(sprintf("%-6s ",nt)); print(round(colMeans(s[,mm],na.rm=TRUE),3))}
