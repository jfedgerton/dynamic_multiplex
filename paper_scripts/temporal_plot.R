source("paper_scripts/leakage_plot.R"); suppressMessages(library(ggplot2))
rows<-list()
for(d in D){p<-strsplit(d,"\\|")[[1]];er<-p[1];a<-p[2];b<-p[3];ty<-p[4];lab<-paste0(nm[a],"-",nm[b]," (",ty,")");w<-eras[[er]]
 for(y in w[1]:w[2]){t<-which(yrs==y);if(!length(t))next; for(m in meth){A<-P[[m]][[t]];dg<-deg[[t]]
  act<-(a%in%names(A))&&(b%in%names(A))&&(a%in%names(dg))&&(b%in%names(dg))&&dg[[a]]>0&&dg[[b]]>0
  oc<-if(!act)"na" else {co<-A[[a]]==A[[b]]; if((ty=="+"&&co)||(ty=="-"&&!co))"ok" else "bad"}
  rows[[length(rows)+1]]<-data.frame(era=er,dyad=lab,method=unname(mlab[m]),year=y,oc=oc,stringsAsFactors=FALSE)}}}
DF<-do.call(rbind,rows); DF$era<-factor(DF$era,levels=names(eras)); DF$method<-factor(DF$method,levels=unname(mlab[meth]))
DF$fill<-ifelse(DF$oc=="na","inactive",ifelse(DF$oc=="ok",paste(DF$era,"correct"),"misaligned"))
pal<-c("Concert correct"="#33a02c","Bismarck correct"="#6a3d9a","WWI correct"="#1b9e77","Interwar correct"="#7570b3","WWII correct"="#66a61e","ColdWar correct"="#2c7fb8","LIO correct"="#2c7fb8","misaligned"="#d7301f","inactive"="#e6e6e6")
saveRDS(DF,"paper_scripts/temporal_DF.rds")
for(er in names(eras)){d<-DF[DF$era==er,]; d$dyad<-factor(d$dyad,levels=rev(unique(d$dyad)))
 p<-ggplot(d,aes(year,dyad,fill=fill))+geom_tile()+facet_wrap(~method,nrow=1)+scale_fill_manual(values=pal,guide="none")+labs(title=paste0(er," (",eras[[er]][1],"-",eras[[er]][2],"): temporal alignment/misalignment"),subtitle="year-by-year; era color=correct, red=misaligned, gray=inactive. (+) should co-cluster, (-) should be separate.",x="Year",y=NULL)+theme_minimal(base_size=11)+theme(panel.grid=element_blank(),axis.text.x=element_text(size=7))
 ggsave(sprintf("manuscript/figures/fig_temporal_%s.png",er),p,width=13,height=1.4+0.4*length(unique(d$dyad)),dpi=150)}
DFc<-DF; DFc$dyad<-factor(DFc$dyad,levels=rev(unique(DFc$dyad)))
pc<-ggplot(DFc,aes(year,dyad,fill=fill))+geom_tile()+facet_grid(era~method,scales="free_y",space="free_y",switch="y")+scale_fill_manual(values=pal,name=NULL)+labs(title="Temporal alignment vs misalignment by method and era (ATOP)",subtitle="each cell=one year; era color=correct, red=misaligned, gray=inactive",x="Year",y=NULL)+theme_minimal(base_size=9)+theme(panel.grid=element_blank(),strip.text.y.left=element_text(angle=0),axis.text.x=element_text(size=6),legend.position="bottom")
ggsave("manuscript/figures/fig_temporal.png",pc,width=14,height=14,dpi=150); cat("saved temporal plots; rows",nrow(DF),"\n")
