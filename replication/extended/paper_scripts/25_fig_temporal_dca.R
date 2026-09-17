# fig_temporal_dca: DCA (Kinne DCAD) analog of 24_fig_temporal_free.R.
# Periods x methods; unified Correct/Misaligned/Inactive; blue/red (colorblind-safe).
# NOTE: pre-1992 DCAD is too sparse for community detection (4 edges in 1980, 47 in
# 1991, 635 in 2010), so the series is restricted to 1992-2010 in two periods.
suppressMessages({library(igraph); library(ggplot2); library(patchwork)}); set.seed(123)
setwd("/storage/work/jfe4/dynamic_multiplex")
S<-readRDS("replication/extended/output/empirical_data/dca_series.rds")
P<-readRDS("replication/extended/output/empirical/dca_partitions.rds")$partitions
yrs<-S$years; deg<-lapply(S$graph_layers,function(g) setNames(degree(g),V(g)$name))
nm<-c("2"="USA","200"="UK","220"="FRA","255"="GER","290"="POL","310"="HUN","365"="RUS","640"="TUR","666"="ISR","705"="KAZ","710"="CHN","750"="IND","770"="PAK")
eras<-list(`Post-Cold War`=c(1992,2001),`Post-9/11`=c(2002,2010))
D<-c("Post-Cold War|2|290|+","Post-Cold War|2|310|+","Post-Cold War|2|200|+","Post-Cold War|220|255|+","Post-Cold War|365|705|+",
"Post-Cold War|2|365|-","Post-Cold War|2|710|-","Post-Cold War|365|640|-","Post-Cold War|750|770|-",
"Post-9/11|2|290|+","Post-9/11|2|750|+","Post-9/11|2|666|+","Post-9/11|220|255|+","Post-9/11|365|710|+",
"Post-9/11|2|365|-","Post-9/11|2|710|-","Post-9/11|365|640|-","Post-9/11|750|770|-")
meth<-c("DynMux Jaccard r1","DynMux Overlap r1","DynMux multislice r1","multinet GLouvain","Cross-sectional + Hungarian","Pooled Leiden")
mlab<-c("DynMux Jaccard r1"="Jaccard","DynMux Overlap r1"="Overlap","DynMux multislice r1"="multislice","multinet GLouvain"="multinet","Cross-sectional + Hungarian"="Hungarian","Pooled Leiden"="Pooled")
rows<-list()
for(d in D){p<-strsplit(d,"\\|")[[1]];er<-p[1];a<-p[2];b<-p[3];ty<-p[4];lab<-paste0(nm[a],"-",nm[b]," (",ty,")");w<-eras[[er]]
for(y in w[1]:w[2]){t<-which(yrs==y);if(!length(t))next; for(m in meth){A<-P[[m]][[t]];dg<-deg[[t]]
act<-(a%in%names(A))&&(b%in%names(A))&&(a%in%names(dg))&&(b%in%names(dg))&&dg[[a]]>0&&dg[[b]]>0
oc<-if(!act)"na" else {co<-A[[a]]==A[[b]]; if((ty=="+"&&co)||(ty=="-"&&!co))"ok" else "bad"}
rows[[length(rows)+1]]<-data.frame(era=er,dyad=lab,method=unname(mlab[m]),year=y,oc=oc,stringsAsFactors=FALSE)}}}
DF<-do.call(rbind,rows); DF$era<-factor(DF$era,levels=names(eras)); DF$method<-factor(DF$method,levels=unname(mlab[meth]))
cat("ACTIVE SHARE by dyad:\n"); print(round(with(DF,tapply(oc!="na",list(era,dyad),mean)),2))
a<-DF[DF$oc!="na",]; a$ty<-ifelse(grepl("\\(\\+\\)",a$dyad),"pos","neg")
cat("\nCORRECT RATE overall:\n"); print(round(with(a,tapply(oc=="ok",method,mean)),3))
cat("\nCORRECT RATE by sign:\n"); print(round(with(a,tapply(oc=="ok",list(ty,method),mean)),3))
cat("\nCORRECT RATE by period:\n"); print(round(with(a,tapply(oc=="ok",list(era,method),mean)),2))
saveRDS(DF,"replication/extended/paper_scripts/temporal_DF_dca.rds")
DF$status<-factor(ifelse(DF$oc=="na","Inactive",ifelse(DF$oc=="ok","Correct","Misaligned")),levels=c("Correct","Misaligned","Inactive"))
pal<-c("Correct"="#0571b0","Misaligned"="#ca0020","Inactive"="#e6e6e6")
ers<-levels(DF$era); plots<-list(); hts<-numeric(0)
for(i in seq_along(ers)){
  d<-DF[DF$era==ers[i],]; d$dyad<-factor(d$dyad,levels=rev(unique(d$dyad)))
  p<-ggplot(d,aes(year,dyad,fill=status))+geom_tile(width=1)+facet_grid(era~method)+
    scale_fill_manual(values=pal,name=NULL,drop=FALSE)+
    scale_x_continuous(breaks=function(l){b<-scales::breaks_pretty(4)(l); b[b>l[1]+1.2 & b<l[2]-1.2]},expand=expansion(add=0.5))+
    labs(x=NULL,y=NULL)+theme_minimal(base_size=9)+
    theme(panel.grid=element_blank(),axis.text.x=element_text(size=6),panel.spacing.x=unit(6,"pt"))
  if(i>1) p<-p+theme(strip.text.x=element_blank())
  if(i==length(ers)) p<-p+labs(x="Year")
  plots[[i]]<-p; hts[i]<-length(unique(d$dyad))+ifelse(i==1,1.6,0.6)}
pc<-wrap_plots(plots,ncol=1,heights=hts)+plot_layout(guides="collect")+
  plot_annotation(title="Temporal alignment vs misalignment by method and period (DCA)",
    subtitle="each cell = one year; blue = correct, red = misaligned, gray = inactive; x-axis range varies by period")
pc<-pc & theme(legend.position="bottom")
ggsave("manuscript/figures/fig_temporal_dca.png",pc,width=13,height=7,dpi=150)
ggsave("manuscript/figures/fig_temporal_dca.pdf",pc,width=13,height=7)
cat("\nsaved fig_temporal_dca\n")
