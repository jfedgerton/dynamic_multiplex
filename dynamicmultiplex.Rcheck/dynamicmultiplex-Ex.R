pkgname <- "dynamicmultiplex"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "dynamicmultiplex-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('dynamicmultiplex')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("animate_multilayer_gif")
### * animate_multilayer_gif

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: animate_multilayer_gif
### Title: Animate multiplex networks as a GIF by layer
### Aliases: animate_multilayer_gif

### ** Examples

## Not run: 
##D if (requireNamespace("ggplot2", quietly = TRUE) &&
##D     requireNamespace("gganimate", quietly = TRUE)) {
##D   set.seed(123)
##D   layers <- lapply(1:3, function(i) {
##D     m <- matrix(rbinom(64, 1, 0.35), nrow = 8)
##D     m <- pmax(m, t(m))
##D     diag(m) <- 0
##D     m
##D   })
##D   fit <- fit_multilayer_jaccard(layers, algorithm = "louvain")
##D   gif_path <- animate_multilayer_gif(
##D     layers,
##D     fit = fit,
##D     output_file = tempfile(fileext = ".gif")
##D   )
##D }
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("animate_multilayer_gif", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("bootstrap_multilayer")
### * bootstrap_multilayer

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: bootstrap_multilayer
### Title: Bootstrap confidence intervals for multilayer community
###   detection
### Aliases: bootstrap_multilayer

### ** Examples

set.seed(123)
layers <- lapply(1:3, function(i) {
  m <- matrix(rbinom(64, 1, 0.35), nrow = 8)
  m <- pmax(m, t(m))
  diag(m) <- 0
  m
})
boot <- bootstrap_multilayer(
  layers,
  fit_type = "jaccard",
  algorithm = "louvain",
  n_boot = 5,
  seed = 123
)
names(boot)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("bootstrap_multilayer", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("community_ci")
### * community_ci

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: community_ci
### Title: Summarize bootstrap results into confidence intervals
### Aliases: community_ci

### ** Examples

set.seed(123)
layers <- lapply(1:3, function(i) {
  m <- matrix(rbinom(64, 1, 0.35), nrow = 8)
  m <- pmax(m, t(m))
  diag(m) <- 0
  m
})
boot <- bootstrap_multilayer(
  layers,
  fit_type = "jaccard",
  algorithm = "louvain",
  n_boot = 5,
  seed = 123
)
ci <- community_ci(boot, alpha = 0.05)
str(ci, max.level = 1)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("community_ci", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("fit_multilayer_identity_ties")
### * fit_multilayer_identity_ties

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fit_multilayer_identity_ties
### Title: Fit multilayer communities with identity interlayer ties
### Aliases: fit_multilayer_identity_ties

### ** Examples

set.seed(123)
layers <- lapply(1:3, function(i) {
  m <- matrix(rbinom(64, 1, 0.35), nrow = 8)
  m <- pmax(m, t(m))
  diag(m) <- 0
  m
})
fit <- fit_multilayer_identity_ties(layers, algorithm = "louvain")
names(fit)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fit_multilayer_identity_ties", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("fit_multilayer_jaccard")
### * fit_multilayer_jaccard

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fit_multilayer_jaccard
### Title: Fit multilayer communities and interlayer weighted Jaccard ties
### Aliases: fit_multilayer_jaccard

### ** Examples

set.seed(123)
layers <- lapply(1:3, function(i) {
  m <- matrix(rbinom(64, 1, 0.35), nrow = 8)
  m <- pmax(m, t(m))
  diag(m) <- 0
  m
})
fit <- fit_multilayer_jaccard(layers, algorithm = "louvain")
names(fit)
fit$layer_communities[[1]]$membership




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fit_multilayer_jaccard", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("fit_multilayer_overlap")
### * fit_multilayer_overlap

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fit_multilayer_overlap
### Title: Fit multilayer communities and interlayer weighted overlap ties
### Aliases: fit_multilayer_overlap

### ** Examples

set.seed(123)
layers <- lapply(1:3, function(i) {
  m <- matrix(rbinom(64, 1, 0.35), nrow = 8)
  m <- pmax(m, t(m))
  diag(m) <- 0
  m
})
fit <- fit_multilayer_overlap(layers, algorithm = "louvain")
names(fit)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fit_multilayer_overlap", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("fit_multilayer_weighted_jaccard")
### * fit_multilayer_weighted_jaccard

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fit_multilayer_weighted_jaccard
### Title: Fit multilayer communities and interlayer node-strength weighted
###   Jaccard ties
### Aliases: fit_multilayer_weighted_jaccard

### ** Examples

set.seed(123)
layers <- lapply(1:3, function(i) {
  w <- matrix(runif(64), 8, 8) * matrix(rbinom(64, 1, 0.35), 8, 8)
  w <- (w + t(w)) / 2
  diag(w) <- 0
  w
})
fit <- fit_multilayer_weighted_jaccard(layers, algorithm = "louvain")
names(fit)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fit_multilayer_weighted_jaccard", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("fit_multilayer_weighted_overlap")
### * fit_multilayer_weighted_overlap

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fit_multilayer_weighted_overlap
### Title: Fit multilayer communities and interlayer node-strength weighted
###   overlap ties
### Aliases: fit_multilayer_weighted_overlap

### ** Examples

set.seed(123)
layers <- lapply(1:3, function(i) {
  w <- matrix(runif(64), 8, 8) * matrix(rbinom(64, 1, 0.35), 8, 8)
  w <- (w + t(w)) / 2
  diag(w) <- 0
  w
})
fit <- fit_multilayer_weighted_overlap(layers, algorithm = "louvain")
names(fit)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fit_multilayer_weighted_overlap", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot_multilayer_alluvial")
### * plot_multilayer_alluvial

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot_multilayer_alluvial
### Title: Plot an alluvial view of community transitions over time
### Aliases: plot_multilayer_alluvial

### ** Examples

if (requireNamespace("ggplot2", quietly = TRUE) &&
    requireNamespace("ggalluvial", quietly = TRUE)) {
  set.seed(123)
  layers <- lapply(1:3, function(i) {
    m <- matrix(rbinom(64, 1, 0.35), nrow = 8)
    m <- pmax(m, t(m))
    diag(m) <- 0
    m
  })
  fit <- fit_multilayer_jaccard(layers, algorithm = "louvain")
  p <- plot_multilayer_alluvial(fit)
}




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot_multilayer_alluvial", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot_multilayer_series")
### * plot_multilayer_series

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot_multilayer_series
### Title: Plot multiplex networks as wrapped panels
### Aliases: plot_multilayer_series

### ** Examples

if (requireNamespace("ggplot2", quietly = TRUE)) {
  set.seed(123)
  layers <- lapply(1:3, function(i) {
    m <- matrix(rbinom(64, 1, 0.35), nrow = 8)
    m <- pmax(m, t(m))
    diag(m) <- 0
    m
  })
  fit <- fit_multilayer_jaccard(layers, algorithm = "louvain")
  p <- plot_multilayer_series(layers, fit = fit)
}




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot_multilayer_series", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("simulate_and_fit_multilayer")
### * simulate_and_fit_multilayer

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: simulate_and_fit_multilayer
### Title: Simulate multiplex layers and fit interlayer models
### Aliases: simulate_and_fit_multilayer

### ** Examples

sim <- simulate_and_fit_multilayer(
  n_nodes = 30,
  n_layers = 3,
  n_communities = 3,
  fit_type = "jaccard",
  algorithm = "louvain",
  seed = 123
)
names(sim)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("simulate_and_fit_multilayer", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
