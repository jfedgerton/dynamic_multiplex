#!/usr/bin/env Rscript
# 15_sim_paired_fig.R -- Paired-difference figure matching the main table.
# facet_grid: rows = regimes, columns = metrics (Joint NMI, K MAE).
# Delta Joint NMI = DynMux minus baseline (positive favors DynMux).
# Delta K MAE    = baseline minus DynMux (positive favors DynMux; lower MAE is better).
# 95% CI of paired differences within simulated network.
# Produces: manuscript/figures/fig_regime_paired.pdf/.png
#           manuscript/tables/tab_regime_paired_fig.csv
set.seed(123)
suppressMessages({ library(ggplot2) })

D <- "replication/extended/output/dynamic"
d <- do.call(rbind, lapply(list.files(D, "dyn_cfg.*csv$", full.names=TRUE),
  function(f){ x <- read.csv(f); x$cfg <- basename(f); x }))
relab <- c("DynMux multislice (adjacent)" = "Multislice (adjacent)",
           "Cross-sectional + Hungarian"  = "Hungarian matching")
d$method <- ifelse(d$method %in% names(relab), relab[d$method], d$method)
d$unit <- paste(d$cfg, d$rep, sep = "_")

rn <- c(birthdeath="Births & Deaths", churnswitch="Churn & Switch",
        regimeshift="Regime Shift", seasonality="Seasonality")
refs <- c("DynMux Jaccard", "DynMux Overlap")
baselines <- c("Hungarian matching", "Multislice (adjacent)",
               "Dynamic SBM", "Pooled Leiden", "multinet GLouvain")
blab <- c("Hungarian matching"    = "Hungarian\nmatching",
          "Multislice (adjacent)" = "Multislice\nadjacent",
          "Dynamic SBM"           = "Dynamic\nSBM",
          "Pooled Leiden"         = "Pooled\nLeiden",
          "multinet GLouvain"     = "Multislice full\n(multinet)")

rows <- list()
for (rg in unique(d$regime)) {
  for (mt in c("nmi_joint", "k_mae")) {
    sub <- d[d$regime == rg, c("unit","method",mt)]
    w <- reshape(sub, idvar="unit", timevar="method", direction="wide")
    names(w) <- sub(paste0(mt,"."), "", names(w), fixed=TRUE)
    for (rf in refs) for (bl in baselines) {
      if (mt == "nmi_joint") dd <- w[[rf]] - w[[bl]] else dd <- w[[bl]] - w[[rf]]
      dd <- dd[is.finite(dd)]
      n <- length(dd); se <- sd(dd)/sqrt(n)
      rows[[length(rows)+1]] <- data.frame(regime=rg, metric=mt, ref=rf,
        baseline=bl, n=n, diff=mean(dd), lo=mean(dd)-1.96*se, hi=mean(dd)+1.96*se)
    }
  }
}
r <- do.call(rbind, rows)
write.csv(r, "manuscript/tables/tab_regime_paired_fig.csv", row.names=FALSE)

r$regime_lab <- factor(unname(rn[r$regime]), levels=unname(rn))
r$metric_lab <- factor(ifelse(r$metric=="nmi_joint", "Joint NMI", "K MAE"),
                       levels=c("Joint NMI","K MAE"))
r$baseline_lab <- factor(unname(blab[r$baseline]),
                         levels=rev(unname(blab[baselines])))
r$ref <- sub("DynMux ", "", r$ref)

p <- ggplot(r, aes(diff, baseline_lab, colour=ref)) +
  geom_vline(xintercept=0, linetype=2, colour="grey40") +
  geom_point(position=position_dodge(width=0.5), size=1.6) +
  geom_errorbarh(aes(xmin=lo, xmax=hi),
                 position=position_dodge(width=0.5), height=0.25) +
  facet_grid(regime_lab ~ metric_lab, scales="free_x") +
  scale_colour_manual(values=c(Jaccard="#1b9e77", Overlap="#d95f02")) +
  labs(x="Paired difference (positive favors DynMux), 95% CI",
       y=NULL, colour="DynMux coupling") +
  theme_bw(base_size=9) +
  theme(legend.position="bottom", panel.grid.minor=element_blank())
ggsave("manuscript/figures/fig_regime_paired.pdf", p, width=6.5, height=7.5)
ggsave("manuscript/figures/fig_regime_paired.png", p, width=6.5, height=7.5, dpi=300)
cat("done: fig_regime_paired (facet_grid regimes x metrics)\n")
