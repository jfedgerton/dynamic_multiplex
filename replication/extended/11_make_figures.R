set.seed(123)
suppressMessages({ library(ggplot2) })
FIGDIR <- "manuscript/figures"; dir.create(FIGDIR, recursive = TRUE, showWarnings = FALSE)
theme_set(theme_bw(base_size = 11))
ax_x <- theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, lineheight = 0.85))
wrap2 <- function(x) gsub("DynMux ", "DynMux\n", gsub("Cross-sectional \\+ ", "Cross-sec\n+ ", x))
rf <- list.files("replication/extended/output/dynamic", "dyn_cfg.*csv$", full.names = TRUE)
if (length(rf)) { d <- do.call(rbind, lapply(rf, read.csv)); agg <- aggregate(nmi_joint ~ method + regime, d, mean)
  rn <- c(birthdeath="Births and Deaths", churnswitch="Churn and Switching", regimeshift="Regime Shift", seasonality="Seasonality")
  agg$regime <- factor(rn[agg$regime], levels = rn); agg$dyn <- grepl("DynMux", agg$method)
  p1 <- ggplot(agg, aes(reorder(method, nmi_joint), nmi_joint, fill = dyn)) + geom_col() + facet_wrap(~ regime, ncol = 2) + coord_flip() +
    scale_fill_manual(values = c(`TRUE`="#1B9E77",`FALSE`="grey75"), guide="none") + labs(x = NULL, y = "Mean joint NMI")
  ggsave(file.path(FIGDIR, "fig_regime.pdf"), p1, width = 11, height = 9); cat("fig_regime\n") }
sf <- "replication/extended/output/empirical/empirical_summary.rds"
if (file.exists(sf)) { s <- readRDS(sf)
  keep <- c("DynMux Overlap r2","DynMux Jaccard r2","DynMux multislice r2","Cross-sectional + Hungarian","Pooled Leiden")
  build <- function(x, net) do.call(rbind, lapply(keep, function(m){ pv <- x$traj[[m]]; if (is.null(pv)) return(NULL); data.frame(network=net, method=m, year=x$years[seq_along(pv)+1], persistence=pv) }))
  pl <- rbind(build(s$atop,"ATOP (1816-2018)"), build(s$dca,"DCA (1980-2010)"))
  p2 <- ggplot(pl, aes(year, persistence, colour = method)) + geom_line(na.rm = TRUE) + facet_wrap(~ network, scales="free_x", ncol=1) +
    labs(x = NULL, y = "Year-to-year persistence (NMI)", colour = NULL) + scale_colour_brewer(palette="Dark2") + theme(legend.position="bottom")
  ggsave(file.path(FIGDIR, "fig_empirical_persistence.pdf"), p2, width = 9, height = 7); cat("fig_persistence\n")
  mk <- function(x, net){ z <- x$df[,c("method","mean_mod")]; z$network <- net; z }; mm <- rbind(mk(s$atop,"ATOP"), mk(s$dca,"DCA")); mm$dyn <- grepl("DynMux", mm$method); mm$mlab <- wrap2(mm$method)
  p3 <- ggplot(mm, aes(reorder(mlab, -mean_mod), mean_mod, fill = dyn)) + geom_col() + facet_wrap(~ network, ncol=1, scales="free_x") +
    scale_fill_manual(values = c(`TRUE`="#1B9E77",`FALSE`="grey75"), guide="none") + labs(x = NULL, y = "Mean per-year modularity") + ax_x
  ggsave(file.path(FIGDIR, "fig_empirical_quality.pdf"), p3, width = 10, height = 8); cat("fig_quality\n") }
of <- "replication/extended/output/empirical/omega_sweep.csv"
if (file.exists(of)) { om <- read.csv(of)
  p5 <- ggplot(om, aes(omega, mod, colour = toupper(network))) + geom_line() + geom_point() + scale_x_log10() + scale_colour_brewer(palette="Dark2") +
    geom_hline(yintercept = c(0.421,0.330), linetype="dashed", colour="grey50") + labs(x = "interlayer coupling omega (log)", y = "Mean per-year modularity", colour = "network", subtitle = "dashed = Overlap r1 reference")
  ggsave(file.path(FIGDIR, "fig_omega.pdf"), p5, width = 8, height = 5); cat("fig_omega\n") }
cat("DONE\n")
