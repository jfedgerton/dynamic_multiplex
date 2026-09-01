suppressMessages(library(igraph)); set.seed(123)
DATA<-"replication/extended/output/empirical_data"; OUT<-"replication/extended/output/empirical"
S<-readRDS(file.path(DATA,"atop_series.rds")); U<-readRDS(file.path(DATA,"atop_union.rds"))
P<-readRDS(file.path(OUT,"atop_partitions.rds"))$partitions; yrs<-S$years
mask<-if(!is.null(U$present))U$present else U$active
G<-lapply(seq_along(S$graph_layers),function(k) delete_vertices(S$graph_layers[[k]],which(!mask[[k]])))
ab<-c(Jac="DynMux Jaccard r1",Ovl="DynMux Overlap r1",msl="DynMux multislice r1",mnt="multinet GLouvain",Hun="Cross-sectional + Hungarian",Pool="Pooled Leiden")
pairs<-list(c(2,290,"USA-POL"),c(2,310,"USA-HUN"),c(2,255,"USA-GER"))
edge<-function(a,b,t){g<-G[[t]];na<-as.character(a);nb<-as.character(b);if(!(na%in%V(g)$name)||!(nb%in%V(g)$name))return(NA);are_adjacent(g,na,nb)}
sm<-function(m,a,b,t){A<-P[[m]][[t]];na<-as.character(a);nb<-as.character(b);if(is.null(A)||!(na%in%names(A))||!(nb%in%names(A)))return(NA);A[[na]]==A[[nb]]}
f<-function(x)ifelse(is.na(x),".",ifelse(x,"Y","n"))
for(pr in pairs){a<-as.integer(pr[1]);b<-as.integer(pr[2]);cat("====",pr[3],"====\n");cat(sprintf("%-5s %-4s %-4s %-4s %-4s %-4s %-4s %-4s\n","yr","edge","Jac","Ovl","msl","mnt","Hun","Pool"))
for(t in which(yrs>=1990&yrs<=2018)){e<-edge(a,b,t);v<-sapply(ab,function(m)sm(m,a,b,t));cat(sprintf("%-5d %-4s %-4s %-4s %-4s %-4s %-4s %-4s\n",yrs[t],f(e),f(v[1]),f(v[2]),f(v[3]),f(v[4]),f(v[5]),f(v[6])))};cat("\n")}
t10<-which(yrs==2010);A<-P[["DynMux Jaccard r1"]][[t10]];cu<-A[["2"]];co<-names(A)[A==cu]
cat("2010 Jaccard USA community",cu,"size",length(co),"members:",paste(co,collapse=" "),"\n")
nato<-c(2,20,200,255,260,220,325,210,211,235,230,385,395,205,375,290,310,316,315,290,349,360,355)
cat("NATO present:",paste(sort(intersect(as.integer(co),nato)),collapse=" "),"\n")
cat("NATO missing:",paste(sort(setdiff(nato,as.integer(co))),collapse=" "),"\n")
