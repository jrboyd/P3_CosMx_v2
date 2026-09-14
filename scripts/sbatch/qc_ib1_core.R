# Is IB1's tissue core salvageable? The cluster-i defect follows a radial gradient
# (12% at the core -> 61% at the periphery) while RB1 is flat. If the core reaches
# RB1-like marker support, IB1 can be analysed core-only instead of discarded.
#
# Reports, per radial quintile and per candidate filter, the marker-based support
# for the B compartment -- the metric that actually failed (IB1 NBC_MBC: 7.5% B,
# 1.3% strongB vs RB1 35.3% / 13.8%).
suppressPackageStartupMessages(library(Seurat))
md <- readRDS("output/04_merge/merged_01_merged.Rds")@meta.data
GOOD <- function(v) !v %in% c("all_markers_0", "below_background", "ambiguous")

for (s in c("IB1", "RB1")) {
  m <- md[md$merge_sample == s, ]
  cx <- median(m$x_slide_mm); cy <- median(m$y_slide_mm)
  m$r <- sqrt((m$x_slide_mm - cx)^2 + (m$y_slide_mm - cy)^2)
  m$q <- cut(m$r, breaks = quantile(m$r, 0:5/5), include.lowest = TRUE, labels = FALSE)
  nb <- m$celltype_coarse == "NBC_MBC"
  cat("\n##########", s, "-- B-compartment support by radial quintile ##########\n")
  out <- data.frame(
    quintile = 1:5,
    n_cells = as.vector(table(m$q)),
    pct_cluster_i = round(100 * tapply(m$inSitu_clusters == "i", m$q, mean), 1),
    pct_marker_informative = round(100 * tapply(GOOD(m$lineage), m$q, mean), 1),
    n_NBC_MBC = as.vector(table(m$q[nb])),
    NBC_pct_markerB = round(100 * tapply(m$lineage[nb] == "B", m$q[nb], mean), 1),
    NBC_pct_strongB = round(100 * tapply(m$strongB[nb], m$q[nb], mean), 1)
  )
  print(out, row.names = FALSE)

  cat("\n-- candidate filters for", s, "--\n")
  filters <- list(
    "as-is (all retained)" = rep(TRUE, nrow(m)),
    "drop cluster i" = m$inSitu_clusters != "i",
    "marker-informative only" = GOOD(m$lineage),
    "drop cluster i + marker-informative" = m$inSitu_clusters != "i" & GOOD(m$lineage),
    "inner 2 quintiles" = m$q <= 2,
    "inner 2 quintiles + drop cluster i" = m$q <= 2 & m$inSitu_clusters != "i"
  )
  res <- do.call(rbind, lapply(names(filters), function(fn) {
    k <- filters[[fn]]; nbk <- k & nb
    data.frame(filter = fn, n_kept = sum(k),
               pct_kept = round(100 * mean(k), 1),
               n_NBC = sum(nbk),
               NBC_pct_markerB = round(100 * mean(m$lineage[nbk] == "B"), 1),
               NBC_pct_strongB = round(100 * mean(m$strongB[nbk]), 1))
  }))
  print(res, row.names = FALSE)
}
cat("\nOK\n")
