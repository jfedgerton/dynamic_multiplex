# =============================================================================
# 07_score_orders.R
# Set-level recovery of Braumoeller (2019) / Goodhart coded international
# orders by the detected communities, for every coded order, every layer
# (year) and every method, on all four empirical networks. This is the data
# behind the manuscript's Figure 3 (precision vs recall of order recovery by
# method and network) and the appendix precision/recall table.
#
# For coded order B (states coded into the order in year y, restricted to the
# states that are present with degree > 0 in that year's layer) and detected
# communities C:
#   J    = max_C |B n C| / |B u C|        (best-matching community)
#   prec = |B n C| / |C|                  at the argmax C
#   rec  = |B n C| / |B|                  at the argmax C
#   nB   = |B|
# Needs only the target SET, so it works for partial/overlapping orders where
# ARI cannot.
#
# Order codings: Braumoeller, Only the Dead replication archive (coding
# credited to Andy Goodhart). Colombia/Bulgaria coding bug (ccode 100 labelled
# Bulgaria) corrected in COM (355 in, 100 out). 19th-century orders are coded
# as "all European states (ccode 200-399)", so those eras test a
# Europe/non-Europe split, not a fine bloc structure.
#
# Consolidates: replication/extended/paper_scripts/29_braumoeller_fullspan.R
#               (codings + helpers) and 32_setlevel_allorders.R (scoring loop),
#               with no readLines()+eval().
#
# Usage:  Rscript 07_score_orders.R
# Env:    DM_ROOT  project root (default getwd())
# Input:  $DM_ROOT/output/empirical_data/<net>_series.rds, <net>_union.rds
#         $DM_ROOT/output/empirical/<net>_partitions.rds   (06_fit_networks.R)
# Output: $DM_ROOT/output/empirical/order_recovery.csv
#         $DM_ROOT/output/empirical/order_recovery.rds
#         long format: net, order, kind, year, method, J, prec, rec, nB
# =============================================================================
set.seed(123)
suppressMessages(library(igraph))

ROOT     <- Sys.getenv("DM_ROOT", unset = getwd())
EMP_DATA <- file.path(ROOT, "output", "empirical_data")
EMP_OUT  <- file.path(ROOT, "output", "empirical")
dir.create(EMP_OUT, recursive = TRUE, showWarnings = FALSE)
stopifnot(dir.exists(EMP_DATA), dir.exists(EMP_OUT))

NETS <- c("atop", "igo", "trade", "dca")

# ---------------------------------------------------------------------------
# Order codings (29_braumoeller_fullspan.R). Each entry: ccode = "start:end".
# ---------------------------------------------------------------------------
LIB<-c("2"="1945:1991","20"="1945:1991","200"="1945:1991","205"="1945:1991","210"="1945:1991","211"="1945:1991","212"="1945:1991","220"="1945:1991","225"="1945:1991","230"="1982:1991","235"="1945:1991","255"="1990:1991","260"="1955:1990","305"="1945:1991","325"="1945:1991","350"="1952:1991","380"="1945:1991","385"="1945:1991","390"="1945:1991","395"="1945:1991","640"="1952:1991","666"="1948:1991","713"="1947:1991","732"="1948:1991","740"="1952:1991","900"="1945:1991","920"="1945:1991")
COM<-c("265"="1955:1990","290"="1955:1991","310"="1955:1991","315"="1955:1991","339"="1955:1991","355"="1955:1991","359"="1945:1991","360"="1955:1991","365"="1955:1991","366"="1945:1991","367"="1945:1991","368"="1945:1991","369"="1945:1991","370"="1945:1991","371"="1945:1991","372"="1945:1991","373"="1945:1991","701"="1945:1991","702"="1945:1991","703"="1945:1991","704"="1945:1991","705"="1945:1991","710"="1949:1991","712"="1945:1991","731"="1948:1991","812"="1949:1991","816"="1945:1991")
LEAGUE<-c("20"="1920:1939","40"="1920:1939","41"="1920:1939","42"="1924:1939","70"="1931:1939","90"="1920:1936","91"="1920:1936","92"="1920:1937","93"="1920:1936","94"="1920:1925","95"="1920:1939","100"="1920:1939","101"="1920:1938","130"="1934:1939","135"="1920:1939","140"="1920:1926","145"="1920:1939","150"="1920:1935","155"="1920:1938","160"="1920:1939","165"="1920:1939","200"="1920:1939","205"="1921:1939","210"="1920:1939","211"="1920:1939","212"="1920:1939","220"="1920:1941","225"="1920:1939","230"="1920:1939","235"="1920:1939","255"="1926:1933","290"="1920:1939","305"="1920:1938","310"="1922:1939","315"="1920:1939","325"="1920:1937","339"="1920:1939","345"="1920:1939","350"="1920:1939","355"="1920:1939","360"="1920:1940","365"="1934:1939","366"="1921:1939","367"="1921:1939","368"="1921:1939","375"="1920:1939","380"="1920:1939","385"="1920:1939","390"="1920:1939","450"="1920:1939","530"="1923:1936","560"="1920:1939","630"="1920:1939","640"="1932:1939","645"="1932:1939","651"="1937:1939","700"="1934:1939","710"="1920:1939","740"="1920:1933","750"="1920:1939","800"="1920:1939","900"="1920:1939","920"="1920:1939")
PCW<-c("2"="1992:2018","20"="1992:2018","200"="1992:2018","205"="1992:2018","210"="1992:2018","211"="1992:2018","212"="1992:2018","220"="1992:2018","225"="1992:2018","230"="1992:2018","235"="1992:2018","255"="1992:2018","260"="1992:2018","290"="2000:2018","305"="1992:2018","310"="2000:2018","316"="2000:2018","317"="2005:2018","325"="1992:2018","338"="2005:2018","339"="2010:2018","344"="2010:2018","349"="2005:2018","350"="1992:2018","352"="2005:2018","355"="2005:2018","360"="2005:2018","366"="2005:2018","367"="2005:2018","368"="2005:2018","375"="1995:2018","380"="1992:2018","385"="1992:2018","390"="1992:2018","395"="1992:2018","640"="1992:2018","666"="1992:2018","713"="1992:2018","732"="1992:2018","740"="1992:2018","900"="1992:2018","920"="1992:2018")
# NOTE: 29 also codes LATINAM (Concert-era Latin American republics) and
# defines ari()/eur()/truth() for the era-level ARI validation. None of these
# feed the set-level scoring (32 uses only LIB, COM, LEAGUE, PCW and inb()),
# so they are not carried here.

# Additional codings (32_setlevel_allorders.R)
CGP<-c("200"="1816:1852","220"="1816:1852","255"="1816:1852","300"="1816:1852","365"="1816:1852")
WP <-c("100"="1955:1991","265"="1955:1990","290"="1955:1991","310"="1955:1991","315"="1955:1991","339"="1955:1968","360"="1955:1991","365"="1955:1991")
# NOTE: WP is reproduced exactly as coded in 32. Its entry "100" is Colombia's
# ccode; the Warsaw Pact member it stands in for is Bulgaria (355). This is the
# same Colombia/Bulgaria labelling bug that 29 corrects in COM; it is left
# uncorrected here because the sources score with it as-is.
MAN<-c("645"="1922:1939","652"="1922:1939","660"="1922:1939","663"="1922:1939","666"="1922:1939")
IAW<-c("205"="1921:1939","290"="1918:1939","305"="1918:1938","310"="1918:1939","315"="1918:1939","339"="1912:1938","345"="1919:1939","355"="1908:1938","366"="1918:1939","367"="1918:1939","368"="1918:1939","369"="1918:1920","371"="1918:1920","372"="1918:1921","373"="1918:1920","375"="1917:1939","651"="1922:1939","670"="1932:1938","678"="1918:1938","700"="1919:1938")

for (cod in list(LIB, COM, LEAGUE, PCW, CGP, WP, MAN, IAW))
  stopifnot(is.character(cod), !is.null(names(cod)), all(grepl("^[0-9]+:[0-9]+$", cod)))

# states in coding table `tb` whose span covers year y
inb<-function(tb,y){k<-names(tb); k[sapply(tb,function(r){p<-as.integer(strsplit(r,":")[[1]]); y>=p[1]&&y<=p[2]})]}

# The 13 coded orders: id, year span, kind, and the coded set B given (y, al)
# where al = states present with degree > 0 in that year's layer.
ORD<-list(
 list(id="Concert GPs",     y=c(1816,1852), kind="membership", f=function(y,al) intersect(inb(CGP,y),al)),
 list(id="Concert Europe",  y=c(1816,1852), kind="geographic", f=function(y,al) al[as.integer(al)>=200 & as.integer(al)<400]),
 list(id="Interim Europe",  y=c(1855,1870), kind="geographic", f=function(y,al) al[as.integer(al)>=200 & as.integer(al)<400]),
 list(id="Bismarck Europe", y=c(1871,1890), kind="geographic", f=function(y,al) al[as.integer(al)>=200 & as.integer(al)<400]),
 list(id="Wilhelm Europe",  y=c(1891,1914), kind="geographic", f=function(y,al) al[as.integer(al)>=200 & as.integer(al)<400]),
 list(id="Indep. after WWI",y=c(1908,1939), kind="membership", f=function(y,al) intersect(inb(IAW,y),al)),
 list(id="Mandates",        y=c(1922,1939), kind="membership", f=function(y,al) intersect(inb(MAN,y),al)),
 list(id="League",          y=c(1920,1941), kind="membership", f=function(y,al) intersect(inb(LEAGUE,y),al)),
 list(id="PW Liberal",      y=c(1945,1991), kind="membership", f=function(y,al) intersect(inb(LIB,y),al)),
 list(id="PW Communist",    y=c(1955,1991), kind="membership", f=function(y,al) intersect(inb(COM,y),al)),
 list(id="Warsaw Pact",     y=c(1955,1991), kind="membership", f=function(y,al) intersect(inb(WP,y),al)),
 list(id="PW Other",        y=c(1945,1991), kind="residual",   f=function(y,al) setdiff(al,c(inb(LIB,y),inb(COM,y)))),
 list(id="PostCW Liberal",  y=c(1992,2018), kind="membership", f=function(y,al) intersect(inb(PCW,y),al)))
stopifnot(length(ORD) == 13L)

# Methods scored (partition names as written by 06) and their figure labels
meth<-c("DynMux Jaccard r1","DynMux Overlap r1","DynMux multislice r1","multinet GLouvain","Cross-sectional + Hungarian","Pooled Leiden")
mlab<-c("Jaccard","Overlap","multislice","multinet","Hungarian","Pooled")
stopifnot(length(meth) == length(mlab))

# best-matching community for set B among membership `mem` over `nodes`:
# returns c(J, precision, recall, |B|), or NULL if fewer than 3 of B are in nodes
best<-function(B,mem,nodes){ bs<-nodes %in% B; if(sum(bs)<3) return(NULL); out<-c(0,0,0)
 for(cc in unique(mem)){cs<-mem==cc; i<-sum(bs&cs); j<-i/sum(bs|cs); if(j>out[1]) out<-c(j,i/sum(cs),i/sum(bs))}
 c(out,sum(bs))}

# ---------------------------------------------------------------------------
# Score every network x order x year x method
# ---------------------------------------------------------------------------
rows<-list()
for(net in NETS){
 series_f <- file.path(EMP_DATA, sprintf("%s_series.rds", net))
 union_f  <- file.path(EMP_DATA, sprintf("%s_union.rds",  net))
 part_f   <- file.path(EMP_OUT,  sprintf("%s_partitions.rds", net))
 for (f in c(series_f, union_f, part_f))
   if (!file.exists(f)) stop("Missing ", f, " -- run 05_build_networks.R and 06_fit_networks.R ", net, " first.", call. = FALSE)
 S<-readRDS(series_f); U<-readRDS(union_f); P<-readRDS(part_f)$partitions
 yrs<-S$years; mask<-if(!is.null(U$present))U$present else U$active
 stopifnot(length(S$graph_layers) == length(yrs), length(mask) == length(yrs), !is.null(P))
 stopifnot(meth[1] %in% names(P))                       # reference partition for the node set
 miss<-setdiff(meth,names(P))
 if(length(miss)){ if(identical(miss,"multinet GLouvain")) warning(sprintf("[%s] multinet GLouvain missing from partitions; scored as NA/skipped",net),call.=FALSE)
   else stop(sprintf("[%s] partitions missing methods: %s",net,paste(miss,collapse="; ")),call.=FALSE) }
 for(m in intersect(meth,names(P))) stopifnot(length(P[[m]]) == length(yrs))
 G<-lapply(seq_along(S$graph_layers),function(k) delete_vertices(S$graph_layers[[k]],which(!mask[[k]])))
 dg<-lapply(G,function(g) setNames(degree(g),V(g)$name))
 n0<-length(rows)
 for(o in ORD) for(y in o$y[1]:o$y[2]){t<-which(yrs==y); if(!length(t))next
  d<-dg[[t]]; A0<-P[[meth[1]]][[t]]; al<-names(d)[d>0]; al<-al[al%in%names(A0)]
  if(length(al)<8) next
  B<-o$f(y,al); if(length(B)<3) next
  for(j in seq_along(meth)){A<-P[[meth[j]]][[t]]; if(!all(al%in%names(A)))next
   r<-best(B,as.integer(unlist(A[al])),al); if(is.null(r))next
   rows[[length(rows)+1]]<-data.frame(net=net,order=o$id,kind=o$kind,year=y,method=mlab[j],
     J=r[1],prec=r[2],rec=r[3],nB=r[4],stringsAsFactors=FALSE)}}
 cat(sprintf("[%s] years=%d (%d-%d) rows=%d\n",net,length(yrs),min(yrs),max(yrs),length(rows)-n0))
 stopifnot(length(rows) > n0)
}
D<-do.call(rbind,rows); rownames(D)<-NULL
stopifnot(is.data.frame(D), nrow(D) > 0,
          identical(names(D), c("net","order","kind","year","method","J","prec","rec","nB")),
          all(D$J >= 0 & D$J <= 1), all(D$prec >= 0 & D$prec <= 1), all(D$rec >= 0 & D$rec <= 1),
          all(D$nB >= 3), !anyNA(D))

# ---------------------------------------------------------------------------
# Write
# ---------------------------------------------------------------------------
csv_f <- file.path(EMP_OUT, "order_recovery.csv")
rds_f <- file.path(EMP_OUT, "order_recovery.rds")
write.csv(D, csv_f, row.names = FALSE)
saveRDS(D, rds_f)
stopifnot(file.exists(csv_f), file.exists(rds_f))

a<-aggregate(cbind(J,prec,rec,nB)~net+order+kind+method,D,mean)
n<-aggregate(year~net+order+method,D,length); a<-merge(a,n,by=c("net","order","method"))
cat("obs:",nrow(D)," cells:",nrow(a),"\n"); print(table(D$net,D$order))
cat("\nmean recovery by net x method (across orders and years):\n")
print(aggregate(cbind(J,prec,rec)~net+method,D,function(x) round(mean(x),3)))
cat("SCORE_DONE ->", csv_f, "\n")
