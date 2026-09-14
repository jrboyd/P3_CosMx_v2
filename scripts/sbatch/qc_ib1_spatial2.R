# IB1 failure-mode analysis, corrected.
#
# v1 (qc_ib1_spatial.R) used x_FOV_px/y_FOV_px -- WITHIN-FOV pixel coordinates --
# so its "radial distance" was noise and its null result meaningless. This uses
# x_slide_mm/y_slide_mm (global) and adds a real distance-to-tissue-edge, plus the
# per-cell inSitu_clusters (a-j) that v1 discovered.
#
# Questions:
#  1. Is high cluster-i / dim membrane spatially PERIPHERAL (fixation gradient)
#     or patchy (imaging failure)? Tested against both radial distance and
#     distance to the tissue boundary (which handles non-convex tissue).
#  2. Characterize cluster i per cell: counts, counts/um2, area, channels.
#     low counts AND low counts/um2 -> RNA/probe failure (over-fixation)
#     normal counts spread over LARGER area -> segmentation dilution
#  3. Does cluster i drive IB1's unsupported NBC_MBC calls? Cross-tab against
#     celltype_coarse and the independent marker lineage columns.
#
# Usage: sbatch --mem=64G scripts/sbatch/run_rscript.sh scripts/sbatch/qc_ib1_spatial2.R

suppressPackageStartupMessages(library(Seurat))

seu <- readRDS("output/04_merge/merged_01_merged.Rds")
md <- seu@meta.data
rm(seu); invisible(gc())
md$cps <- md$nCount_RNA / pmax(md$Area.um2, 1)   # counts per um^2

for (cc in c("inSitu_clusters", "x_slide_mm", "y_slide_mm", "Area.um2", "nCount_RNA")) {
  stopifnot(cc %in% colnames(md))
}

qs <- function(v) round(quantile(v, c(.25, .5, .75), na.rm = TRUE), 3)

for (s in c("IB1", "RB1")) {
  m <- md[md$merge_sample == s, ]
  cat("\n##########################", s, "(", nrow(m), "retained cells )##########################\n")

  # ---- 1. per-CELL cluster i burden (not the FOV-level frac_i) --------------
  ct <- sort(table(m$inSitu_clusters), decreasing = TRUE)
  cat("\ninSitu_clusters per-cell composition (%):\n")
  print(round(100 * ct / sum(ct), 2))
  is_i <- m$inSitu_clusters == "i"
  cat("cluster i cells:", sum(is_i), sprintf("(%.1f%%)\n", 100 * mean(is_i)))

  # ---- 2. what IS cluster i? counts vs area vs channels --------------------
  keys <- c("nCount_RNA", "nFeature_RNA", "Area.um2", "cps", "complexity",
            "propNegative", "Mean.Membrane", "Mean.CD45", "Mean.DAPI",
            "inSitu_clusters_posterior_probability")
  cat("\nmedian per-cell QC, cluster i vs rest:\n")
  tab <- t(sapply(keys, function(k)
    c(cluster_i = round(median(m[[k]][is_i], na.rm = TRUE), 3),
      rest      = round(median(m[[k]][!is_i], na.rm = TRUE), 3),
      ratio     = round(median(m[[k]][is_i], na.rm = TRUE) /
                          median(m[[k]][!is_i], na.rm = TRUE), 3))))
  print(tab)

  # ---- 3. does cluster i drive the cell-type calls? ------------------------
  cat("\ncelltype_coarse composition WITHIN cluster i (%):\n")
  print(round(100 * sort(table(m$celltype_coarse[is_i]), decreasing = TRUE) / sum(is_i), 1))
  cat("\ncluster-i share OF each celltype_coarse (%):\n")
  print(round(100 * tapply(is_i, m$celltype_coarse, mean), 1))
  for (lc in c("lineage", "per_cell_marker_lineage", "strongB")) {
    if (lc %in% colnames(m)) {
      cat("\n", lc, " x cluster i (row % within cluster i):\n", sep = "")
      print(round(100 * sort(table(m[[lc]][is_i]), decreasing = TRUE) / sum(is_i), 1))
      cat(lc, "among NBC_MBC-called cells (%):\n")
      nb <- m$celltype_coarse == "NBC_MBC"
      print(round(100 * sort(table(m[[lc]][nb]), decreasing = TRUE) / sum(nb), 1))
    }
  }

  # ---- 4. spatial: radial AND distance-to-tissue-edge ----------------------
  x <- m$x_slide_mm; y <- m$y_slide_mm
  cat("\nslide coords: x", paste(round(range(x, na.rm = TRUE), 2), collapse = " to "),
      "mm | y", paste(round(range(y, na.rm = TRUE), 2), collapse = " to "), "mm\n")

  # occupancy grid at 0.1 mm; edge bins are occupied bins touching an empty bin,
  # so internal holes count as edge too (a fixation gradient is about surfaces).
  BIN <- 0.1
  ix <- floor((x - min(x, na.rm = TRUE)) / BIN) + 1
  iy <- floor((y - min(y, na.rm = TRUE)) / BIN) + 1
  occ <- matrix(FALSE, max(ix, na.rm = TRUE), max(iy, na.rm = TRUE))
  occ[cbind(ix, iy)] <- TRUE
  pad <- matrix(FALSE, nrow(occ) + 2, ncol(occ) + 2)
  pad[2:(nrow(occ) + 1), 2:(ncol(occ) + 1)] <- occ
  nb_empty <- matrix(FALSE, nrow(occ), ncol(occ))
  for (dx in -1:1) for (dy in -1:1) {
    if (dx == 0 && dy == 0) next
    nb_empty <- nb_empty | !pad[(2:(nrow(occ) + 1)) + dx, (2:(ncol(occ) + 1)) + dy]
  }
  edge <- which(occ & nb_empty, arr.ind = TRUE)
  cat("occupied bins:", sum(occ), "| edge bins:", nrow(edge), "\n")
  # distance from each cell's bin to the nearest edge bin (chunked, mm)
  cells_bin <- unique(cbind(ix, iy))
  d <- numeric(nrow(cells_bin))
  for (k in seq_len(nrow(cells_bin))) {
    d[k] <- BIN * sqrt(min((edge[, 1] - cells_bin[k, 1])^2 + (edge[, 2] - cells_bin[k, 2])^2))
  }
  key <- paste(cells_bin[, 1], cells_bin[, 2])
  m$d_edge <- d[match(paste(ix, iy), key)]
  cx <- median(x, na.rm = TRUE); cy <- median(y, na.rm = TRUE)
  m$r <- sqrt((x - cx)^2 + (y - cy)^2)

  for (metric in c("d_edge", "r")) {
    lbl <- if (metric == "d_edge") "distance to tissue edge (mm)" else "radial distance from centroid (mm)"
    cat("\n-- by ", lbl, " (quintiles) --\n", sep = "")
    q <- cut(m[[metric]], breaks = unique(quantile(m[[metric]], 0:5 / 5, na.rm = TRUE)),
             include.lowest = TRUE, labels = FALSE)
    out <- data.frame(
      bin = 1:max(q, na.rm = TRUE),
      median_dist = round(tapply(m[[metric]], q, median), 2),
      pct_cluster_i = round(100 * tapply(is_i, q, mean), 1),
      Membrane = round(tapply(m$Mean.Membrane, q, median), 0),
      CD45 = round(tapply(m$Mean.CD45, q, median), 0),
      DAPI = round(tapply(m$Mean.DAPI, q, median), 0),
      nCount = round(tapply(m$nCount_RNA, q, median), 0),
      cps = round(tapply(m$cps, q, median), 3)
    )
    print(out, row.names = FALSE)
    cat(sprintf("  spearman rho(%s, cluster_i) = %+.3f | vs Membrane %+.3f | vs nCount %+.3f\n",
                metric,
                cor(m[[metric]], as.numeric(is_i), method = "spearman", use = "complete.obs"),
                cor(m[[metric]], m$Mean.Membrane, method = "spearman", use = "complete.obs"),
                cor(m[[metric]], m$nCount_RNA, method = "spearman", use = "complete.obs")))
  }
  rm(m); invisible(gc())
}
cat("\nOK\n")
