suppressMessages(library(igraph)); set.seed(123)
DATA<-"replication/extended/output/empirical_data"; OUT<-"replication/extended/output/empirical"; net<-"atop"
S<-readRDS(file.path(DATA,sprintf("%s_series.rds",net))); U<-readRDS(file.path(DATA,sprintf("%s_union.rds",net)))
P<-readRDS(file.path(OUT,sprintf("%s_partitions.rds",net)))$partitions; yrs<-S$years; JAC<-"DynMux Jaccard r1"
G<-lapply(seq_along(S$graph_layers),function(k) igraph::delete_vertices(S$graph_layers[[k]],which(!U$present[[k]])))
deg<-lapply(G,function(g){setNames(igraph::degree(g),V(g)$name)})
act<-function(t,c){a<-P[[JAC]][[t]];dg<-deg[[t]];(c%in%names(a))&&(c%in%names(dg))&&dg[[c]]>0}
comm<-function(t,c){P[[JAC]][[t]][[c]]}
edge<-function(t,c1,c2){g<-G[[t]];if(!(c1%in%V(g)$name)||!(c2%in%V(g)$name))return(FALSE);are_adjacent(g,c1,c2)}
rng<-function(v){if(!length(v))return("(never)");b<-c();s<-v[1];p<-v[1];for(i in 2:length(v)){if(i<=length(v)&&v[i]==p+1){p<-v[i]}else{b<-c(b,if(s==p)s else paste0(s,"-",p));s<-v[i];p<-v[i]}};b<-c(b,if(s==p)s else paste0(s,"-",p));paste(b,collapse=", ")}
coclust<-function(codes){yy<-c();for(t in seq_along(yrs)){if(!all(sapply(codes,function(c)act(t,c))))next;labs<-sapply(codes,function(c)comm(t,c));if(length(unique(labs))==1)yy<-c(yy,yrs[t])};yy}
cases<-list(
 c(lab="Holy Alliance: Russia+Austria+Prussia", codes="365,300,255"),
 c(lab="Dreikaiserbund: Germany+AusHun+Russia", codes="255,300,365"),
 c(lab="Franco-Russian: France+Russia",         codes="220,365"),
 c(lab="Triple Alliance: Germany+AusHun+Italy", codes="255,300,325"),
 c(lab="NAZI-SOVIET: Germany+Russia",           codes="255,365"),
 c(lab="Rome-Berlin: Germany+Italy",            codes="255,325"),
 c(lab="Axis core: Germany+Italy+Japan",        codes="255,325,740"),
 c(lab="Axis+east: Germany+Hungary+Romania+Bulgaria+Finland", codes="255,310,360,355,375"),
 c(lab="Sino-Soviet: China+USSR",               codes="710,365"),
 c(lab="NATO core: USA+UK+WGermany",            codes="2,200,260"),
 c(lab="Warsaw Pact: USSR+Poland+EGermany",     codes="365,290,265"),
 c(lab="US-Japan: USA+Japan",                   codes="2,740"),
 c(lab="NATO enlargement: USA+Poland",          codes="2,290"),
 c(lab="NATO enlargement: USA+Hungary",         codes="2,310"))
cat("=== Jaccard r1 co-membership years (all listed states in ONE community) ===\n")
for(x in cases){cd<-strsplit(x[["codes"]],",")[[1]];yy<-coclust(cd);cat(sprintf("%-52s : %s\n",x[["lab"]],rng(yy)))}
gp<-c("2"="USA","200"="UK","220"="FRA","255"="GER","300"="AUH","325"="ITA","365"="RUS","710"="CHN","740"="JPN","260"="WGE","640"="TUR","290"="POL","345"="YUG")
snap<-c(1825,1850,1875,1890,1910,1914,1925,1938,1940,1943,1950,1960,1975,1985,1995,2005,2015)
cat("\n=== Great-power groupings by year (Jaccard r1; states sharing a community are grouped) ===\n")
for(y in snap){t<-which(yrs==y);if(!length(t))next;pres<-names(gp)[sapply(names(gp),function(c)act(t,c))];if(length(pres)<2){cat(y,": (fewer than 2 active)\n");next};labs<-sapply(pres,function(c)comm(t,c));grps<-split(pres,labs);grps<-grps[order(-sapply(grps,length))];out<-sapply(grps,function(g)paste(gp[g],collapse="+"));cat(sprintf("%d: %s\n",y,paste(out,collapse="  |  ")))}
