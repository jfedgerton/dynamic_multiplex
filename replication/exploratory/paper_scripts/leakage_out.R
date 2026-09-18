source("replication/extended/paper_scripts/leakage_plot.R")
suppressMessages(library(ggplot2)); dir.create("manuscript/tables",showWarnings=FALSE)
for(er in levels(DF$era)){d<-DF[DF$era==er,]; d$dyad<-factor(as.character(d$dyad),levels=rev(unique(as.character(d$dyad))))
 p<-ggplot(d,aes(method,dyad,fill=fill))+geom_tile(color="white",linewidth=0.6)+scale_fill_manual(values=pal,guide="none")+labs(title=paste0(er," (",eras[[er]][1],"-",eras[[er]][2],"): alignment vs misalignment"),subtitle="(+) should co-cluster; (-) should be separate. Saturated=correct, contrast=misaligned, gray=inactive.",x=NULL,y=NULL)+theme_minimal(base_size=12)+theme(panel.grid=element_blank())
 ggsave(sprintf("manuscript/figures/fig_leakage_%s.png",er),p,width=7,height=1.2+0.4*length(unique(d$dyad)),dpi=150)}
cat("per-era plots done\n")
sym<-function(oc) if(oc=="ok")"$\\checkmark$" else if(oc=="bad")"$\\times$" else "--"
mc<-levels(DF$method); L<-c("% requires amssymb for \\checkmark","\\begin{table}[htbp]\\centering\\small","\\caption{Alliance co-membership ground truth by era: (+) pairs should be co-clustered, (--) pairs should be separate. Cell = whether the method is correct ($\\checkmark$) or misaligned ($\\times$); -- = inactive that era.}","\\label{tab:leakage}",paste0("\\begin{tabular}{ll",paste(rep("c",length(mc)),collapse=""),"}"),"\\toprule",paste0("Dyad & Truth & ",paste(mc,collapse=" & "),"\\\\"),"\\midrule")
for(er in levels(DF$era)){L<-c(L,paste0("\\addlinespace \\multicolumn{",2+length(mc),"}{l}{\\textbf{",er,"}}\\\\")); dz<-unique(as.character(DF$dyad[DF$era==er]))
 for(dd in dz){tr<-if(grepl("[(][+][)]",dd))"align" else "sep"; nmn<-sub(" [(][-+][)]","",dd); cells<-sapply(mc,function(x){oc<-DF$oc[DF$era==er & as.character(DF$dyad)==dd & DF$method==x][1]; sym(oc)}); L<-c(L,paste0(nmn," & ",tr," & ",paste(cells,collapse=" & "),"\\\\"))}}
L<-c(L,"\\bottomrule","\\end{tabular}","\\end{table}"); writeLines(L,"manuscript/tables/tab_leakage.tex"); cat("wrote tab_leakage.tex lines",length(L),"\n")
