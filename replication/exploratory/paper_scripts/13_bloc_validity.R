#!/usr/bin/env Rscript
# ============================================================================
# 13_bloc_validity.R  --  Historical-plausibility check, faceted by model.
# Focal years; states classified US-aligned / USSR-aligned / Other under each
# main-text model; black outline = historically-implausible co-membership
# (antagonist pair sharing a community that era). Set1 colours + hand-rolled
# hatch texture (no ggpattern -> no sf/units/proj dependency).
# ============================================================================
set.seed(123)
suppressMessages({ library(ggplot2) })
DATA <- "replication/extended/output/empirical_data"; OUT <- "replication/extended/output/empirical"
P   <- readRDS(file.path(OUT,  "atop_partitions.rds"))$partitions
yrs <- readRDS(file.path(DATA, "atop_series.rds"))$years

models <- c("DynMux Overlap r1","DynMux Jaccard r1","DynMux multislice r1",
            "Pooled Leiden","Cross-sectional + Hungarian","multinet GLouvain")
mlab <- c("Overlap (r1)","Jaccard (r1)","Multislice (r1)",
          "Pooled Leiden","Cross-sec + Hungarian","multinet GLouvain"); names(mlab) <- models
cc <- c(USA=2,UK=200,France=220,Germany=255,`W.Germany`=260,Italy=325,Canada=20,Japan=740,
        Poland=290,Hungary=310,Czechoslovakia=315,`E.Germany`=265,Romania=360,
        Bulgaria=355,USSR=365,China=710,`Austria-Hung`=300)
focal <- c(1914,1943,1960,1985)
eras <- list(list(yr=1914:1918,A=c(255,300,640),B=c(200,220,365,2)),
             list(yr=1939:1945,A=c(255,325,740),B=c(2,200,365,220)),
             list(yr=1947:1989,A=c(2,200,220,260,325,20,740),B=c(365,290,265,315,310,360,355)))
badcodes <- function(mem,yr){ out<-character(0); for(e in eras){ if(!(yr %in% e$yr)) next
  la<-mem[intersect(as.character(e$A),names(mem))]; lb<-mem[intersect(as.character(e$B),names(mem))]
  sh<-intersect(la,lb); out<-c(out,names(la)[la %in% sh],names(lb)[lb %in% sh]) }; unique(out) }

rows<-list()
for(m in models){ part<-P[[m]]; if(is.null(part)) next
  for(i in which(yrs %in% focal)){ mem<-part[[i]]; if(is.null(mem)) next; yr<-yrs[i]
    usa<-mem["2"]; usr<-mem["365"]; bc<-badcodes(mem,yr)
    for(nm in names(cc)){ code<-as.character(cc[nm]); if(!code %in% names(mem)) next; lab<-mem[code]
      side<-if(!is.na(usa)&&lab==usa) "US-aligned" else if(!is.na(usr)&&lab==usr) "USSR/Russia-aligned" else "Other"
      rows[[length(rows)+1]]<-data.frame(model=mlab[m],year=yr,country=nm,side=side,bad=code %in% bc) } } }
df<-do.call(rbind,rows)
df$model<-factor(df$model,levels=mlab)
ylev<-rev(names(cc))
df$xn<-match(as.character(df$year),as.character(focal))
df$yn<-match(df$country,ylev)
df$side<-factor(df$side,levels=c("US-aligned","USSR/Russia-aligned","Other"))

# hand-rolled hatch: 2 parallel "/" per non-US tile; add "\" for Other (crosshatch)
h<-0.45; nz<-df[df$side!="US-aligned",]; ot<-df[df$side=="Other",]
seg<-rbind(
  data.frame(model=nz$model, x=nz$xn-h, y=nz$yn,   xend=nz$xn,   yend=nz$yn+h),
  data.frame(model=nz$model, x=nz$xn,   y=nz$yn-h, xend=nz$xn+h, yend=nz$yn),
  data.frame(model=ot$model, x=ot$xn-h, y=ot$yn+h, xend=ot$xn+h, yend=ot$yn-h))
seg$model<-factor(seg$model,levels=mlab)

pal<-c(`US-aligned`="#377EB8",`USSR/Russia-aligned`="#E41A1C",`Other`="#BDBDBD")  # Set1 blue/red + grey
p<-ggplot(df,aes(xn,yn))+
  geom_tile(aes(fill=side),colour="grey85")+
  geom_segment(data=seg,aes(x=x,y=y,xend=xend,yend=yend),colour="white",linewidth=0.3,inherit.aes=FALSE)+
  geom_tile(data=subset(df,bad),fill=NA,colour="black",linewidth=1)+
  facet_wrap(~model,nrow=2)+
  scale_fill_manual(values=pal)+
  scale_x_continuous(breaks=seq_along(focal),labels=focal,expand=c(0,0))+
  scale_y_continuous(breaks=seq_along(ylev),labels=ylev,expand=c(0,0))+
  labs(x=NULL,y=NULL,fill=NULL,subtitle="Black outline = historically implausible co-membership")+
  theme_bw(base_size=10)+theme(legend.position="bottom",panel.grid=element_blank())
ggsave("manuscript/figures/fig_bloc_validity.pdf",p,width=12,height=8)
cat("wrote fig_bloc_validity.pdf\n")
