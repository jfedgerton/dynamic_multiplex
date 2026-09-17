# fig_temporal v2: per-era x scales (patchwork stack), unified Correct/Misaligned/Inactive legend,
# era strips on right (vertical, ggplot default). Reads cached temporal_DF.rds.
suppressMessages({library(ggplot2); library(patchwork)})
setwd("/storage/work/jfe4/dynamic_multiplex")
DF <- readRDS("replication/extended/paper_scripts/temporal_DF.rds")
DF$status <- factor(ifelse(DF$oc=="na","Inactive",ifelse(DF$oc=="ok","Correct","Misaligned")),
                    levels=c("Correct","Misaligned","Inactive"))
pal <- c("Correct"="#1b9e77","Misaligned"="#d7301f","Inactive"="#e6e6e6")
ers <- levels(DF$era)
plots <- list(); hts <- numeric(0)
for(i in seq_along(ers)){
  d <- DF[DF$era==ers[i],]
  d$dyad <- factor(d$dyad, levels=rev(unique(d$dyad)))
  p <- ggplot(d, aes(year, dyad, fill=status)) + geom_tile(width=1) +
    facet_grid(era~method) +
    scale_fill_manual(values=pal, name=NULL, drop=FALSE) +
    scale_x_continuous(expand=expansion(add=0.5)) +
    labs(x=NULL, y=NULL) +
    theme_minimal(base_size=9) +
    theme(panel.grid=element_blank(), axis.text.x=element_text(size=6))
  if(i>1) p <- p + theme(strip.text.x=element_blank())
  if(i==length(ers)) p <- p + labs(x="Year")
  plots[[i]] <- p
  hts[i] <- length(unique(d$dyad)) + ifelse(i==1,1.6,0.6)
}
pc <- wrap_plots(plots, ncol=1, heights=hts) + plot_layout(guides="collect") +
  plot_annotation(title="Temporal alignment vs misalignment by method and era (ATOP)",
    subtitle="each cell = one year; teal = correct, red = misaligned, gray = inactive; x-axis range varies by era")
pc <- pc & theme(legend.position="bottom")
ggsave("manuscript/figures/fig_temporal.png", pc, width=13, height=14, dpi=150)
ggsave("manuscript/figures/fig_temporal.pdf", pc, width=13, height=14)
cat("saved fig_temporal free-scale version\n")
