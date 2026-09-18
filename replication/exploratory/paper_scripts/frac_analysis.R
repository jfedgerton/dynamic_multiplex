source("replication/extended/paper_scripts/leakage_plot.R")
frac<-function(m,a,b,y0,y1){co<-c();for(y in y0:y1){t<-which(yrs==y);if(!length(t))next;A<-P[[m]][[t]];dg<-deg[[t]];if(!(a%in%names(A)&&b%in%names(A)&&a%in%names(dg)&&b%in%names(dg)))next;if(dg[[a]]<=0||dg[[b]]<=0)next;co<-c(co,as.integer(A[[a]]==A[[b]]))};if(length(co)<1)return(c(f=NA_real_,n=0));c(f=mean(co),n=length(co))}
rows<-list()
for(d in D){p<-strsplit(d,"\\|")[[1]];er<-p[1];a<-p[2];b<-p[3];ty<-p[4];w<-eras[[er]];for(m in meth){r<-frac(m,a,b,w[1],w[2]);rows[[length(rows)+1]]<-data.frame(era=er,dyad=paste0(nm[a],"-",nm[b]),type=ty,method=unname(mlab[m]),frac=r[["f"]],stringsAsFactors=FALSE)}}
FR<-do.call(rbind,rows); FR$era<-factor(FR$era,levels=names(eras)); ml<-unname(mlab[meth]); d2<-FR[!is.na(FR$frac),]
cat("== continuous balanced accuracy (align=mean frac, sep=mean 1-frac) ==\n"); cat(sprintf("%-11s %8s %8s %8s\n","method","align","separate","balAcc"))
for(m in ml){s<-d2[d2$method==m,];al<-mean(s$frac[s$type=="+"]);se<-mean(1-s$frac[s$type=="-"]);cat(sprintf("%-11s %8.2f %8.2f %8.2f\n",m,al,se,(al+se)/2))}
cat("\n== per-era balAcc: Jaccard / Overlap / multislice ==\n")
for(er in names(eras)){v<-sapply(c("Jaccard","Overlap","multislice"),function(m){s<-d2[d2$era==er&d2$method==m,];(mean(s$frac[s$type=="+"])+mean(1-s$frac[s$type=="-"]))/2});cat(sprintf("%-10s %5.2f %5.2f %5.2f\n",er,v[1],v[2],v[3]))}
mc<-ml; sy<-function(f) if(is.na(f))"--" else sprintf("%.2f",f)
L<-c("\\begin{table}[htbp]\\centering\\small","\\caption{Co-clustering fraction: share of co-present years the two states share a community. (+) should be near 1, (--) near 0.}","\\label{tab:frac}",paste0("\\begin{tabular}{ll",paste(rep("c",length(mc)),collapse=""),"}"),"\\toprule",paste0("Dyad & Truth & ",paste(mc,collapse=" & "),"\\\\"),"\\midrule")
for(er in names(eras)){L<-c(L,paste0("\\addlinespace \\multicolumn{",2+length(mc),"}{l}{\\textbf{",er,"}}\\\\"));dz<-unique(FR$dyad[FR$era==er]);for(dd in dz){sub<-FR[FR$era==er&FR$dyad==dd,];ty<-sub$type[1];cells<-sapply(mc,function(x)sy(sub$frac[sub$method==x][1]));L<-c(L,paste0(dd," & ",ifelse(ty=="+","align","sep")," & ",paste(cells,collapse=" & "),"\\\\"))}}
L<-c(L,"\\bottomrule","\\end{tabular}","\\end{table}"); writeLines(L,"manuscript/tables/tab_frac.tex"); cat("\nwrote tab_frac.tex rows",nrow(FR),"\n")
