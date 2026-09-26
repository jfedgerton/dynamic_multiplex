# =============================================================================
# 11_alignment_dyads.R
# Face-validity check for the empirical partitions: does the tracker put
# well-known ALIGNED dyads (USA-UK after 1941, USSR-Poland 1945-89, ...) in the
# same community, and well-known ANTI-ALIGNED dyads (USA-USSR in the Cold War,
# India-Pakistan, North-South Korea, ...) in different communities, in the
# years when that alignment held?
#
# For every network, method, hand-coded dyad and year inside the dyad's
# period, same = 1 if both states are present with degree > 0 in that year's
# layer and share a community, 0 if present and split, NA if either is absent.
# Aligned dyads should score near 1, anti-aligned dyads near 0. Nothing here
# is estimated: the dyad list and periods are hard-coded below (COW ccodes).
#
# Usage:  DM_ROOT=/path/to/dynamic_multiplex Rscript replication/empirical/11_alignment_dyads.R
# Input:  $DM_ROOT/output/empirical_data/<net>_series.rds, <net>_union.rds
#         $DM_ROOT/output/empirical/<net>_partitions.rds   (08_fit_networks.R)
# Output: $DM_ROOT/output/empirical/alignment_dyads.csv
#           long: net, dyad_id, dyad, type, a, b, year, method, same
#         $DM_ROOT/output/empirical/alignment_summary.csv
#           net, method, type, dyads, dyad_years, share_same
#         $DM_ROOT/manuscript/tables/tab_app_alignment.tex
#           share of period years in the same community, DynMux (Jaccard), by net
#         $DM_ROOT/manuscript/figures/fig_app_alignment.pdf (+ .png)
#           year x dyad tiles, DynMux (Jaccard), one panel per network
# =============================================================================
set.seed(123)
suppressMessages({ library(igraph); library(ggplot2) })

ROOT     <- Sys.getenv("DM_ROOT", unset = getwd())
EMP_DATA <- file.path(ROOT, "output", "empirical_data")
EMP_OUT  <- file.path(ROOT, "output", "empirical")
TAB      <- file.path(ROOT, "manuscript", "tables")
FIG      <- file.path(ROOT, "manuscript", "figures")
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
stopifnot(dir.exists(EMP_DATA), dir.exists(EMP_OUT))

NETS    <- c("atop", "igo", "trade", "dca")
NET_LAB <- c(atop = "Alliances (ATOP)", igo = "IGO membership", trade = "Trade", dca = "Defense cooperation (DCA)")
METH    <- c("DynMux Jaccard r1"           = "DynMux (Jaccard)",
             "DynMux multislice r1"        = "Multislice",
             "Cross-sectional + Hungarian" = "Hungarian",
             "Pooled Leiden"               = "Pooled Leiden")

# ---------------------------------------------------------------------------
# Hand-coded dyads. COW ccodes as character. Period = years the alignment or
# rivalry is uncontroversial. Germany: 255 before 1945 and from 1991; West
# Germany 260 and East Germany 265 for 1955-1990. Russia/USSR = 365 throughout.
# ---------------------------------------------------------------------------
D <- rbind(
  # ---- aligned: expect the same community --------------------------------
  data.frame(id = "A01", a = "2",   b = "200", dyad = "USA-UK",                  type = "aligned", y1 = 1941, y2 = 2014),
  data.frame(id = "A02", a = "2",   b = "20",  dyad = "USA-Canada",              type = "aligned", y1 = 1940, y2 = 2014),
  data.frame(id = "A03", a = "2",   b = "260", dyad = "USA-West Germany",        type = "aligned", y1 = 1955, y2 = 1990),
  data.frame(id = "A04", a = "2",   b = "255", dyad = "USA-Germany",             type = "aligned", y1 = 1991, y2 = 2014),
  data.frame(id = "A05", a = "2",   b = "740", dyad = "USA-Japan",               type = "aligned", y1 = 1952, y2 = 2014),
  data.frame(id = "A06", a = "2",   b = "732", dyad = "USA-South Korea",         type = "aligned", y1 = 1953, y2 = 2014),
  data.frame(id = "A07", a = "2",   b = "666", dyad = "USA-Israel",              type = "aligned", y1 = 1967, y2 = 2014),
  data.frame(id = "A08", a = "2",   b = "900", dyad = "USA-Australia",           type = "aligned", y1 = 1951, y2 = 2014),
  data.frame(id = "A09", a = "2",   b = "670", dyad = "USA-Saudi Arabia",        type = "aligned", y1 = 1945, y2 = 2014),
  data.frame(id = "A10", a = "200", b = "220", dyad = "UK-France",               type = "aligned", y1 = 1949, y2 = 2014),
  data.frame(id = "A11", a = "220", b = "260", dyad = "France-West Germany",     type = "aligned", y1 = 1963, y2 = 1990),
  data.frame(id = "A12", a = "200", b = "235", dyad = "UK-Portugal",             type = "aligned", y1 = 1816, y2 = 2014),
  data.frame(id = "A13", a = "365", b = "290", dyad = "USSR-Poland",             type = "aligned", y1 = 1945, y2 = 1989),
  data.frame(id = "A14", a = "365", b = "265", dyad = "USSR-East Germany",       type = "aligned", y1 = 1955, y2 = 1990),
  data.frame(id = "A15", a = "365", b = "315", dyad = "USSR-Czechoslovakia",     type = "aligned", y1 = 1948, y2 = 1989),
  data.frame(id = "A16", a = "365", b = "40",  dyad = "USSR-Cuba",               type = "aligned", y1 = 1962, y2 = 1991),
  data.frame(id = "A17", a = "365", b = "710", dyad = "USSR-China (treaty)",     type = "aligned", y1 = 1950, y2 = 1959),
  data.frame(id = "A18", a = "365", b = "750", dyad = "USSR-India",              type = "aligned", y1 = 1971, y2 = 1991),
  data.frame(id = "A19", a = "365", b = "710", dyad = "Russia-China (SCO)",      type = "aligned", y1 = 2001, y2 = 2014),
  data.frame(id = "A20", a = "710", b = "731", dyad = "China-North Korea",       type = "aligned", y1 = 1961, y2 = 2014),
  data.frame(id = "A21", a = "255", b = "300", dyad = "Germany-Austria-Hungary", type = "aligned", y1 = 1879, y2 = 1918),
  data.frame(id = "A22", a = "220", b = "365", dyad = "France-Russia",           type = "aligned", y1 = 1894, y2 = 1917),
  data.frame(id = "A23", a = "200", b = "740", dyad = "UK-Japan",                type = "aligned", y1 = 1902, y2 = 1921),
  data.frame(id = "A24", a = "255", b = "325", dyad = "Germany-Italy (Axis)",    type = "aligned", y1 = 1936, y2 = 1943),
  data.frame(id = "A25", a = "255", b = "740", dyad = "Germany-Japan (Axis)",    type = "aligned", y1 = 1940, y2 = 1945),
  # ---- anti-aligned: expect different communities ------------------------
  data.frame(id = "B01", a = "2",   b = "365", dyad = "USA-USSR (Cold War)",     type = "anti", y1 = 1947, y2 = 1989),
  data.frame(id = "B02", a = "2",   b = "710", dyad = "USA-China",               type = "anti", y1 = 1950, y2 = 1971),
  data.frame(id = "B03", a = "2",   b = "40",  dyad = "USA-Cuba",                type = "anti", y1 = 1961, y2 = 2014),
  data.frame(id = "B04", a = "2",   b = "731", dyad = "USA-North Korea",         type = "anti", y1 = 1950, y2 = 2014),
  data.frame(id = "B05", a = "2",   b = "630", dyad = "USA-Iran",                type = "anti", y1 = 1980, y2 = 2014),
  data.frame(id = "B06", a = "2",   b = "645", dyad = "USA-Iraq",                type = "anti", y1 = 1990, y2 = 2003),
  data.frame(id = "B07", a = "2",   b = "740", dyad = "USA-Japan (WWII)",        type = "anti", y1 = 1941, y2 = 1945),
  data.frame(id = "B08", a = "365", b = "710", dyad = "USSR-China (split)",      type = "anti", y1 = 1961, y2 = 1989),
  data.frame(id = "B09", a = "260", b = "265", dyad = "West-East Germany",       type = "anti", y1 = 1955, y2 = 1989),
  data.frame(id = "B10", a = "750", b = "770", dyad = "India-Pakistan",          type = "anti", y1 = 1947, y2 = 2014),
  data.frame(id = "B11", a = "666", b = "651", dyad = "Israel-Egypt",            type = "anti", y1 = 1948, y2 = 1978),
  data.frame(id = "B12", a = "666", b = "652", dyad = "Israel-Syria",            type = "anti", y1 = 1948, y2 = 2014),
  data.frame(id = "B13", a = "630", b = "645", dyad = "Iran-Iraq",               type = "anti", y1 = 1980, y2 = 1988),
  data.frame(id = "B14", a = "731", b = "732", dyad = "North-South Korea",       type = "anti", y1 = 1950, y2 = 2014),
  data.frame(id = "B15", a = "710", b = "713", dyad = "China-Taiwan",            type = "anti", y1 = 1949, y2 = 2014),
  data.frame(id = "B16", a = "710", b = "816", dyad = "China-Vietnam",           type = "anti", y1 = 1979, y2 = 1991),
  data.frame(id = "B17", a = "350", b = "640", dyad = "Greece-Turkey",           type = "anti", y1 = 1955, y2 = 2014),
  data.frame(id = "B18", a = "530", b = "520", dyad = "Ethiopia-Somalia",        type = "anti", y1 = 1977, y2 = 1991),
  data.frame(id = "B19", a = "220", b = "255", dyad = "France-Germany (pre-WWI)",type = "anti", y1 = 1871, y2 = 1914),
  data.frame(id = "B20", a = "220", b = "255", dyad = "France-Germany (interwar)",type = "anti", y1 = 1919, y2 = 1939),
  data.frame(id = "B21", a = "200", b = "255", dyad = "UK-Germany (pre-WWI)",    type = "anti", y1 = 1900, y2 = 1918),
  data.frame(id = "B22", a = "200", b = "255", dyad = "UK-Germany (Nazi era)",   type = "anti", y1 = 1933, y2 = 1945),
  data.frame(id = "B23", a = "200", b = "365", dyad = "UK-Russia (Great Game)",  type = "anti", y1 = 1830, y2 = 1906),
  data.frame(id = "B24", a = "365", b = "300", dyad = "Russia-Austria-Hungary",  type = "anti", y1 = 1908, y2 = 1918),
  data.frame(id = "B25", a = "740", b = "710", dyad = "Japan-China",             type = "anti", y1 = 1931, y2 = 1945),
  stringsAsFactors = FALSE)
stopifnot(!anyDuplicated(D$id), all(D$y1 <= D$y2), all(D$type %in% c("aligned", "anti")), all(D$a != D$b))
cat("dyads:", nrow(D), " aligned", sum(D$type == "aligned"), " anti", sum(D$type == "anti"), "\n")

# ---------------------------------------------------------------------------
# Co-membership by net x method x dyad x year
# ---------------------------------------------------------------------------
rows <- list()
for (net in NETS) {
  series_f <- file.path(EMP_DATA, sprintf("%s_series.rds", net))
  union_f  <- file.path(EMP_DATA, sprintf("%s_union.rds",  net))
  part_f   <- file.path(EMP_OUT,  sprintf("%s_partitions.rds", net))
  for (f in c(series_f, union_f, part_f)) if (!file.exists(f)) stop("Missing ", f, call. = FALSE)
  S <- readRDS(series_f); U <- readRDS(union_f); P <- readRDS(part_f)$partitions
  yrs <- S$years; mask <- if (!is.null(U$present)) U$present else U$active
  stopifnot(length(S$graph_layers) == length(yrs), length(mask) == length(yrs))
  miss <- setdiff(names(METH), names(P))
  if (length(miss)) stop(sprintf("[%s] partitions missing methods: %s", net, paste(miss, collapse = "; ")), call. = FALSE)
  for (m in names(METH)) stopifnot(length(P[[m]]) == length(yrs))
  # states present with degree > 0 in each year (same rule as 09_score_orders.R)
  G  <- lapply(seq_along(yrs), function(k) delete_vertices(S$graph_layers[[k]], which(!mask[[k]])))
  dg <- lapply(G, function(g) setNames(degree(g), V(g)$name))
  n0 <- length(rows)
  for (t in seq_along(yrs)) {
    y  <- yrs[t]; d <- dg[[t]]; al <- names(d)[d > 0]
    Dy <- D[D$y1 <= y & D$y2 >= y, ]
    if (!nrow(Dy)) next
    for (m in names(METH)) {
      A <- P[[m]][[t]]
      both <- Dy$a %in% al & Dy$b %in% al & Dy$a %in% names(A) & Dy$b %in% names(A)
      same <- rep(NA_integer_, nrow(Dy))
      same[both] <- as.integer(unlist(A[Dy$a[both]]) == unlist(A[Dy$b[both]]))
      rows[[length(rows) + 1]] <- data.frame(net = net, dyad_id = Dy$id, dyad = Dy$dyad, type = Dy$type,
                                             a = Dy$a, b = Dy$b, year = y, method = unname(METH[m]),
                                             same = same, stringsAsFactors = FALSE)
    }
  }
  cat(sprintf("[%s] years %d (%d-%d)  dyad-year rows %d\n", net, length(yrs), min(yrs), max(yrs),
              sum(sapply(rows[(n0 + 1):length(rows)], nrow))))
}
L <- do.call(rbind, rows)
stopifnot(nrow(L) > 0, all(L$same %in% c(0L, 1L, NA)),
          identical(names(L), c("net", "dyad_id", "dyad", "type", "a", "b", "year", "method", "same")))
write.csv(L, file.path(EMP_OUT, "alignment_dyads.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# Summary: share of observed dyad-years in the same community, by type
# ---------------------------------------------------------------------------
Lo <- L[!is.na(L$same), ]
S1 <- aggregate(same ~ net + method + type, Lo, mean)
S1$dyad_years <- aggregate(same ~ net + method + type, Lo, length)$same
S1$dyads      <- aggregate(dyad_id ~ net + method + type, Lo, function(x) length(unique(x)))$dyad_id
names(S1)[names(S1) == "same"] <- "share_same"
S1 <- S1[order(S1$net, S1$method, S1$type), c("net", "method", "type", "dyads", "dyad_years", "share_same")]
write.csv(S1, file.path(EMP_OUT, "alignment_summary.csv"), row.names = FALSE)
cat("\nshare of dyad-years in the same community (aligned should be high, anti low):\n")
W <- reshape(S1[, c("net", "method", "type", "share_same")], idvar = c("net", "method"), timevar = "type", direction = "wide")
W$gap <- W$share_same.aligned - W$share_same.anti
print(W[order(W$net, -W$gap), ], row.names = FALSE, digits = 3)

# ---------------------------------------------------------------------------
# Appendix table: per dyad, DynMux (Jaccard), share of period years same community
# ---------------------------------------------------------------------------
Lj <- Lo[Lo$method == "DynMux (Jaccard)", ]
Tj <- aggregate(same ~ dyad_id + dyad + type + net, Lj, mean)
Tw <- reshape(Tj, idvar = c("dyad_id", "dyad", "type"), timevar = "net", direction = "wide")
Tw <- merge(D[, c("id", "y1", "y2")], Tw, by.x = "id", by.y = "dyad_id")
Tw <- Tw[order(Tw$type, Tw$id), ]
cell <- function(x) ifelse(is.na(x), "--", sprintf("%.2f", x))
# Rows for a longtable; the \\begin{longtable}, caption and label live in appendix.tex
hdr <- "Dyad & Period & ATOP & IGO & Trade & DCA \\\\"
tex <- c("\\toprule", hdr, "\\midrule", "\\endfirsthead",
         "\\multicolumn{6}{l}{\\textit{Table~\\ref{tab:app_alignment} continued}} \\\\",
         "\\toprule", hdr, "\\midrule", "\\endhead",
         "\\midrule", "\\multicolumn{6}{r}{\\textit{continued on next page}} \\\\", "\\endfoot",
         "\\bottomrule", "\\endlastfoot",
         "\\multicolumn{6}{l}{\\emph{Aligned (expect near 1)}} \\\\")
for (ty in c("aligned", "anti")) {
  if (ty == "anti") tex <- c(tex, "\\addlinespace", "\\multicolumn{6}{l}{\\emph{Anti-aligned (expect near 0)}} \\\\")
  for (i in which(Tw$type == ty))
    tex <- c(tex, sprintf("%s & %d--%d & %s & %s & %s & %s \\\\", Tw$dyad[i], Tw$y1[i], Tw$y2[i],
                          cell(Tw[["same.atop"]][i]), cell(Tw[["same.igo"]][i]), cell(Tw[["same.trade"]][i]), cell(Tw[["same.dca"]][i])))
}
writeLines(tex, file.path(TAB, "tab_app_alignment.tex"))

# ---------------------------------------------------------------------------
# Appendix figure: year x dyad tiles, DynMux (Jaccard), one panel per network
# ---------------------------------------------------------------------------
Fj <- L[L$method == "DynMux (Jaccard)", ]
Fj$status  <- factor(ifelse(is.na(Fj$same), "Not both present", ifelse(Fj$same == 1L, "Same community", "Different communities")),
                     levels = c("Same community", "Different communities", "Not both present"))
Fj$net_lab <- factor(unname(NET_LAB[Fj$net]), levels = unname(NET_LAB))
ord <- D[order(D$type, D$y1), ]
Fj$dyad_f  <- factor(Fj$dyad, levels = rev(unique(ord$dyad)))
Fj$type_lab <- factor(ifelse(Fj$type == "aligned", "Aligned", "Anti-aligned"), levels = c("Aligned", "Anti-aligned"))
p <- ggplot(Fj, aes(year, dyad_f, fill = status)) +
  geom_tile(height = 0.85) +
  facet_grid(type_lab ~ net_lab, scales = "free", space = "free_y") +
  scale_fill_manual(values = c("Same community" = "#1b9e77", "Different communities" = "#d95f02", "Not both present" = "grey85"), name = NULL) +
  scale_x_continuous(breaks = seq(1850, 2000, by = 50)) +
  labs(x = "Year", y = NULL) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", panel.grid = element_blank(),
        axis.text.y = element_text(size = 7), strip.text = element_text(size = 9))
ggsave(file.path(FIG, "fig_app_alignment.pdf"), p, width = 10, height = 8.5)
ggsave(file.path(FIG, "fig_app_alignment.png"), p, width = 10, height = 8.5, dpi = 200)
cat("ALIGN_DONE ->", file.path(EMP_OUT, "alignment_dyads.csv"), "\n")
