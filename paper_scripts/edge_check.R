suppressMessages(library(igraph))
DATA<-"replication/extended/output/empirical_data"
S<-readRDS(file.path(DATA,"atop_series.rds")); yrs<-S$years
G<-S$graph_layers
chk<-function(a,b,lab){yy<-c();for(t in seq_along(yrs)){g<-G[[t]];na<-as.character(a);nb<-as.character(b)
  if(!is.null(g)&&na%in%V(g)$name&&nb%in%V(g)$name&&are_adjacent(g,na,nb))yy<-c(yy,yrs[t])}
  if(length(yy)==0){cat(sprintf("%-10s: NO EDGE ever\n",lab))}else{cat(sprintf("%-10s: %d-%d (n=%d)\n",lab,min(yy),max(yy),length(yy)))}}
chk(710,731,"CHN-PRK"); chk(365,731,"RUS-PRK"); chk(2,740,"USA-JPN"); chk(740,731,"JPN-PRK"); chk(2,731,"USA-PRK"); chk(710,365,"CHN-RUS")
