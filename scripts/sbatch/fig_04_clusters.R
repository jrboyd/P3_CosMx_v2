# Step-4 joint clustering: sample distribution across clusters + annotation by the
# cell-typing strategies that actually discriminate (marker lineage, strongB,
# InSituType). celltype_coarse/fine from the reference transfer are near-constant
# across clusters; the transfer was disabled in step 3 and is no longer reported.
#
# Sample segregation is evaluated at BOTH resolutions -- low (res 0.4) is the
# primary check, since coarser clusters would mix samples if the high-resolution
# segregation were merely over-splitting. Every figure carries its resolution in
# its filename.
# Figures -> output/04_merge/figures, each with a companion CSV of its own data.
# Usage: sbatch --mem=128G scripts/sbatch/run_rscript.sh scripts/sbatch/fig_04_clusters.R

suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(data.table); setDTthreads(1)
})

OUT <- "output/04_merge"; FIG <- file.path(OUT, "figures")
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
MIN_N <- 100L   # a per-sample rate below this many cells is not plotted

## ---- design tokens (documented default palette, light surface) -------------
SURF <- "#fcfcfb"; INK <- "#0b0b0b"; INK2 <- "#52514e"; MUTED <- "#898781"
GRID <- "#e1e0d9"; AXIS <- "#c3c2b7"
S_IB1 <- "#2a78d6"; S_RB1 <- "#eb6834"          # categorical slots 1, 2
SEQ <- c("#cde2fb","#9ec5f4","#6da7ec","#3987e5","#256abf","#184f95","#0d366b")

theme_viz <- function(base = 11, xgrid = FALSE) {
  th <- theme_minimal(base_size = base) +
    theme(
      plot.background  = element_rect(fill = SURF, colour = NA),
      panel.background = element_rect(fill = SURF, colour = NA),
      panel.grid.major.y = element_line(colour = GRID, linewidth = 0.25),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line.x  = element_line(colour = AXIS, linewidth = 0.3),
      axis.ticks   = element_blank(),
      axis.text    = element_text(colour = MUTED, size = base - 2),
      axis.title   = element_text(colour = INK2, size = base - 1),
      plot.title   = element_text(colour = INK, face = "bold", size = base + 2),
      plot.subtitle= element_text(colour = INK2, size = base - 1),
      plot.caption = element_text(colour = MUTED, size = base - 3, hjust = 0),
      strip.text   = element_text(colour = INK, face = "bold", size = base - 1),
      legend.text  = element_text(colour = INK2, size = base - 2),
      legend.title = element_text(colour = INK2, size = base - 2),
      legend.position = "top", legend.justification = "left",
      legend.key.height = unit(9, "pt"))
  if (xgrid) th <- th + theme(panel.grid.major.x = element_line(colour = GRID, linewidth = 0.25))
  th
}
save_fig <- function(name, p, w, h) {
  f <- file.path(FIG, name); ggsave(f, p, width = w, height = h, dpi = 150, bg = SURF)
  cat("wrote", f, "\n")
}
save_csv <- function(name, d) { fwrite(d, file.path(FIG, name)); cat("wrote", file.path(FIG, name), "\n") }

## ---- data -----------------------------------------------------------------
seu <- readRDS(file.path(OUT, "merged_processed.Rds"))
md <- as.data.table(seu@meta.data, keep.rownames = "cell")
ump <- Embeddings(seu, "pearsonumap")
md[, `:=`(umap1 = ump[cell, 1], umap2 = ump[cell, 2])]
rm(seu, ump); invisible(gc())
md[, samp := factor(merge_sample, levels = c("IB1", "RB1"))]
md[, is_i := inSitu_clusters == "i"]
EXP_IB1 <- md[, mean(samp == "IB1")]
cat("cells:", nrow(md), "| expected IB1 share:", round(100*EXP_IB1, 2), "%\n")
cat("UMAP1 range:", paste(round(range(md$umap1), 2), collapse = " to "),
    "| UMAP2 range:", paste(round(range(md$umap2), 2), collapse = " to "), "\n")

RES <- c(high = "pearson_clusters_high", low = "pearson_clusters_low")

## ---- master summary -------------------------------------------------------
summ <- rbindlist(lapply(names(RES), function(rn) {
  d <- md[, .(
    n = .N, n_IB1 = sum(samp == "IB1"), n_RB1 = sum(samp == "RB1"),
    pct_IB1 = round(100*mean(samp == "IB1"), 2),
    pct_strongB = round(100*mean(strongB, na.rm = TRUE), 2),
    pct_strongB_IB1 = round(100*mean(strongB[samp == "IB1"], na.rm = TRUE), 2),
    pct_strongB_RB1 = round(100*mean(strongB[samp == "RB1"], na.rm = TRUE), 2),
    pct_markerB = round(100*mean(lineage == "B", na.rm = TRUE), 2),
    pct_marker_informative = round(100*mean(!lineage %in%
      c("all_markers_0","below_background","ambiguous"), na.rm = TRUE), 2),
    pct_cluster_i = round(100*mean(is_i), 2),
    pct_cluster_i_IB1 = round(100*mean(is_i[samp == "IB1"]), 2),
    pct_cluster_i_RB1 = round(100*mean(is_i[samp == "RB1"]), 2),
    med_nCount = median(nCount_RNA), med_nFeature = median(nFeature_RNA),
    med_Area_um2 = round(median(Area.um2), 1),
    med_PanCK = round(median(Mean.PanCK), 0), med_Membrane = round(median(Mean.Membrane), 0),
    med_CD45 = round(median(Mean.CD45), 0),
    lineage_dominant = names(sort(table(lineage), decreasing = TRUE))[1],
    lineage_dominant_pct = round(100*max(table(lineage))/.N, 1),
  ), by = .(cluster = get(RES[[rn]]))][, resolution := rn][]
  d[, log2_enrich_IB1 := round(log2((pct_IB1/100 + 1e-6)/EXP_IB1), 3)]
  d[, enriched_in := ifelse(log2_enrich_IB1 > 0, "IB1", "RB1")]
  setorder(d, -pct_IB1)[]
}), use.names = TRUE)
setcolorder(summ, c("resolution", "cluster", "n"))
save_csv("04_cluster_summary.csv", summ)
make_res_figs <- function(RN) {
RLAB <- c(low = "res 0.4", high = "res 0.8")[[RN]]
CLCOL <- RES[[RN]]
hi <- summ[resolution == RN]
ord <- hi[order(-pct_IB1), as.character(cluster)]

## ==== FIG 1: sample composition (IB1 at the base, so the reference line reads) ====
comp <- melt(hi[, .(cluster = as.character(cluster), n, IB1 = pct_IB1, RB1 = 100 - pct_IB1)],
             id.vars = c("cluster","n"), variable.name = "samp", value.name = "pct")
comp[, cluster := factor(cluster, levels = ord)]
comp[, samp := factor(samp, levels = c("RB1","IB1"))]   # IB1 last => drawn at base
save_csv(sprintf("04_cluster_sample_composition_%s.csv", RN), comp[order(cluster, samp)])
lab <- hi[log2_enrich_IB1 >= 1][, cluster := factor(as.character(cluster), levels = ord)]
p1 <- ggplot(comp, aes(cluster, pct, fill = samp)) +
  geom_col(width = 0.78, colour = SURF, linewidth = 0.35) +
  geom_hline(yintercept = 100*EXP_IB1, colour = INK, linetype = "22", linewidth = 0.4) +
  geom_text(data = lab, aes(cluster, 100, label = sprintf("%.0f%% IB1", pct_IB1)),
            inherit.aes = FALSE, vjust = -0.6, size = 3, colour = INK) +
  scale_fill_manual(values = c(IB1 = S_IB1, RB1 = S_RB1), name = NULL,
                    breaks = c("IB1","RB1")) +
  scale_y_continuous(expand = expansion(c(0, 0.10)), labels = function(x) paste0(x, "%")) +
  labs(title = sprintf("Joint clusters are strongly sample-segregated (%s)", RLAB),
       subtitle = sprintf(paste0("IB1 share of each joint cluster (%s, %d clusters). ",
                 "Dashed line = whole-merge IB1 share (%.1f%%).\nOnly IB1-enriched clusters ",
                 "(>=2x expected) are labelled; every other cluster sits below the line."),
                 RLAB, nrow(hi), 100*EXP_IB1),
       x = sprintf("joint cluster (%s)", CLCOL), y = "share of cluster",
       caption = "Data: 04_cluster_sample_composition.csv / 04_cluster_summary.csv") +
  theme_viz()
save_fig(sprintf("fig_04_cluster_sample_composition_%s.png", RN), p1, if (RN=="low") 7 else 11, 5.4)

## ==== FIG 2: signed enrichment ====
e <- hi[, .(cluster = factor(as.character(cluster),
            levels = hi[order(log2_enrich_IB1), as.character(cluster)]),
            log2_enrich_IB1, enriched_in, n)]
p2 <- ggplot(e, aes(log2_enrich_IB1, cluster, fill = enriched_in)) +
  geom_col(width = 0.72) +
  geom_vline(xintercept = 0, colour = AXIS, linewidth = 0.4) +
  scale_fill_manual(values = c(IB1 = S_IB1, RB1 = S_RB1), name = "enriched in") +
  scale_x_continuous(expand = expansion(0.05)) +
  labs(title = sprintf("Sample enrichment per joint cluster (%s)", RLAB),
       subtitle = "log2(observed IB1 share / expected IB1 share). 0 = matches the merge, +1 = 2x over-represented.",
       x = "log2 enrichment of IB1", y = "joint cluster",
       caption = "Data: 04_cluster_summary.csv") +
  theme_viz(xgrid = TRUE) + theme(panel.grid.major.y = element_blank())
save_fig(sprintf("fig_04_cluster_sample_enrichment_%s.png", RN), p2, 8, if (RN=="low") 3.6 else 7)

## ==== FIG 3/4: composition heatmaps (ramp spans the DATA range) ====
heat <- function(col, fname, ttl, sub, xlab_angle = 40) {
  d <- md[, .N, by = c(CLCOL, col)]
  setnames(d, c("cluster","grp","N"))
  d[, pct := 100*N/sum(N), by = cluster]
  ordc <- hi[order(-pct_strongB), as.character(cluster)]
  d[, cluster := factor(as.character(cluster), levels = rev(ordc))]
  d[, grp := factor(as.character(grp))]
  csv <- sub("^fig_04_", "04_", sub("\\.png$", ".csv", fname))
  save_csv(csv, d[order(cluster, -pct)])
  top <- ceiling(max(d$pct)/10)*10
  d[, lbl := ifelse(pct >= 15, sprintf("%.0f", pct), NA_character_)]
  p <- ggplot(d, aes(grp, cluster, fill = pct)) +
    geom_tile(colour = SURF, linewidth = 0.4) +
    geom_text(aes(label = lbl, colour = pct > 0.62*top), size = 2.6, na.rm = TRUE) +
    scale_fill_gradientn(colours = SEQ, limits = c(0, top), name = "% of cluster") +
    scale_colour_manual(values = c(`TRUE` = SURF, `FALSE` = INK), guide = "none") +
    labs(title = ttl, subtitle = sub, x = NULL,
         y = "joint cluster (ordered by % strongB)",
         caption = paste0("Labelled where >=15% of the cluster. Data: ", csv)) +
    theme_viz() +
    theme(panel.grid.major.y = element_blank(),
          axis.text.x = element_text(angle = xlab_angle, hjust = 1, colour = INK2))
  save_fig(fname, p, 9, if (RN=="low") 4.4 else 7.5)
}
heat("lineage", sprintf("fig_04_cluster_lineage_heatmap_%s.png", RN),
     sprintf("Marker lineage composition of each joint cluster (%s)", RLAB),
     "Per-cell marker lineage calls - independent of the reference transfer.")
heat("inSitu_clusters", sprintf("fig_04_cluster_insitu_heatmap_%s.png", RN),
     sprintf("InSituType cluster composition of each joint cluster (%s)", RLAB),
     "InSituType assignments a-j; cluster i is the low-information population.", 0)

## ==== FIG 5: per-sample marker profile (rates on >=100 cells only) ====
M_SB <- "% strongB (B-marker positive)"; M_CI <- "% InSituType cluster i (low information)"
prof <- rbindlist(list(
  hi[, .(cluster = as.character(cluster), measure = M_SB,
         IB1 = pct_strongB_IB1, RB1 = pct_strongB_RB1, n_IB1, n_RB1)],
  hi[, .(cluster = as.character(cluster), measure = M_CI,
         IB1 = pct_cluster_i_IB1, RB1 = pct_cluster_i_RB1, n_IB1, n_RB1)]))
prof <- melt(prof, id.vars = c("cluster","measure","n_IB1","n_RB1"),
             variable.name = "samp", value.name = "pct")
prof[, n_sample := ifelse(samp == "IB1", n_IB1, n_RB1)]
prof[, plotted := n_sample >= MIN_N]
prof[, `:=`(cluster = factor(cluster, levels = ord),
            samp = factor(samp, levels = c("IB1","RB1")),
            measure = factor(measure, levels = c(M_SB, M_CI)))]
save_csv(sprintf("04_cluster_marker_profile_%s.csv", RN), prof[order(measure, cluster, samp)])
cat("suppressed (per-sample n <", MIN_N, "):", prof[plotted == FALSE, .N], "of", nrow(prof), "bars\n")
base <- data.table(
  measure = factor(c(M_SB, M_SB, M_CI, M_CI), levels = c(M_SB, M_CI)),
  samp = factor(c("IB1","RB1","IB1","RB1"), levels = c("IB1","RB1")),
  y = c(100*md[samp == "IB1", mean(strongB, na.rm = TRUE)],
        100*md[samp == "RB1", mean(strongB, na.rm = TRUE)],
        100*md[samp == "IB1", mean(is_i)], 100*md[samp == "RB1", mean(is_i)]))
p5 <- ggplot(prof[plotted == TRUE], aes(cluster, pct, fill = samp)) +
  geom_col(position = position_dodge(width = 0.8, preserve = "single"), width = 0.74) +
  geom_hline(data = base, aes(yintercept = y, colour = samp),
             linetype = "22", linewidth = 0.4) +
  facet_wrap(~measure, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c(IB1 = S_IB1, RB1 = S_RB1), name = NULL) +
  scale_colour_manual(values = c(IB1 = S_IB1, RB1 = S_RB1), guide = "none") +
  scale_y_continuous(expand = expansion(c(0, 0.06)), labels = function(x) paste0(x, "%")) +
  labs(title = sprintf("Do joint clusters carry the same signal in both samples? (%s)", RLAB),
       subtitle = paste0("Rates computed within each sample. Dashed lines = that sample's ",
                         "whole-tissue baseline.\nBars omitted where the sample contributes ",
                         "fewer than ", MIN_N, " cells to the cluster."),
       x = sprintf("joint cluster (%s)", CLCOL), y = NULL,
       caption = "Data: 04_cluster_marker_profile.csv (column `plotted` flags omitted bars)") +
  theme_viz()
save_fig(sprintf("fig_04_cluster_marker_profile_%s.png", RN), p5, if (RN=="low") 7 else 11, 7)
} # end make_res_figs
for (RN in c("low", "high")) { cat("\n===== figures for resolution", RN, "=====\n"); make_res_figs(RN) }



## ==== HEADLINE: does coarsening the clustering mix the samples? ====
seg <- copy(summ)
seg[, resolution_lab := factor(ifelse(resolution == "low", "res 0.4 (7 clusters)",
                                      "res 0.8 (26 clusters)"),
        levels = c("res 0.4 (7 clusters)", "res 0.8 (26 clusters)"))]
seg[, clu := factor(paste(resolution, cluster), levels = seg[order(resolution, pct_IB1), paste(resolution, cluster)])]
TOL <- log2(1.5)
seg[, within := abs(log2_enrich_IB1) < TOL]
save_csv("04_sample_segregation_by_resolution.csv",
         seg[order(resolution, -pct_IB1),
             .(resolution, cluster, n, n_IB1, n_RB1, pct_IB1, log2_enrich_IB1,
               enriched_in, within_1.5x_of_expected = within, lineage_dominant,
               lineage_dominant_pct)])
sm <- seg[, .(clusters = .N, within = sum(within),
              pct_cells_within = round(100*sum(n[within])/sum(n), 1)), by = resolution_lab]
cat("
segregation summary:
"); print(sm)
pseg <- ggplot(seg, aes(pct_IB1, clu, colour = enriched_in)) +
  geom_segment(aes(x = 100*EXP_IB1, xend = pct_IB1, yend = clu), linewidth = 0.5) +
  geom_point(aes(size = n)) +
  geom_vline(xintercept = 100*EXP_IB1, colour = INK, linetype = "22", linewidth = 0.4) +
  facet_grid(resolution_lab ~ ., scales = "free_y", space = "free_y") +
  scale_colour_manual(values = c(IB1 = S_IB1, RB1 = S_RB1), name = "enriched in") +
  scale_size_area(max_size = 7, labels = function(x) format(x, big.mark = ","), name = "cells") +
  scale_x_continuous(limits = c(0, 100), labels = function(x) paste0(x, "%")) +
  scale_y_discrete(labels = function(v) sub("^(low|high) ", "", v)) +
  labs(title = "Sample segregation persists when the clustering is coarsened",
       subtitle = paste0("IB1 share per cluster at both resolutions; dashed line = whole-merge ",
         sprintf("expectation (%.1f%%).\nAt res 0.4, 0 of 7 clusters fall within 1.5x of it ",
                 100*EXP_IB1),
         "- so the segregation is not high-resolution over-splitting."),
       x = "IB1 share of cluster", y = "joint cluster",
       caption = "Data: 04_sample_segregation_by_resolution.csv") +
  theme_viz(xgrid = TRUE) + theme(panel.grid.major.y = element_blank(),
        strip.text.y = element_text(angle = 0))
save_fig("fig_04_sample_segregation_by_resolution.png", pseg, 8.6, 9.5)

## ==== FIG 6/7: binned UMAP (explicit tile geometry; blank-panel bug fixed) ====
NB <- 130
x0 <- min(md$umap1); y0 <- min(md$umap2)
wx <- diff(range(md$umap1))/NB; wy <- diff(range(md$umap2))/NB
md[, bx := pmin(NB, floor((umap1 - x0)/wx) + 1L)]
md[, byy := pmin(NB, floor((umap2 - y0)/wy) + 1L)]
bins <- md[, .(n = .N,
               pct_IB1 = 100*mean(samp == "IB1"),
               pct_strongB = 100*mean(strongB, na.rm = TRUE),
               pct_cluster_i = 100*mean(is_i)), by = .(bx, byy)][n >= 20]
bins[, `:=`(x = x0 + (bx - 0.5)*wx, y = y0 + (byy - 0.5)*wy)]
save_csv("04_umap_bins.csv", bins)
cat("bins >=20 cells:", nrow(bins), "covering", round(100*sum(bins$n)/nrow(md), 1), "% of cells\n")
cat("bin x range:", paste(round(range(bins$x), 2), collapse=" to "),
    "| y range:", paste(round(range(bins$y), 2), collapse=" to "),
    "| tile", round(wx, 3), "x", round(wy, 3), "\n")
stopifnot(nrow(bins) > 0, all(is.finite(bins$x)), all(is.finite(bins$y)))

umap_fig <- function(col, fname, ttl, sub, legname) {
  d <- copy(bins); d[, v := get(col)]
  p <- ggplot(d, aes(x, y, fill = v)) +
    geom_tile(width = wx, height = wy) +
    scale_fill_gradientn(colours = SEQ, limits = c(0, 100), name = legname) +
    coord_equal(ratio = wy/wx, expand = TRUE) +
    labs(title = ttl, subtitle = sub, x = "UMAP 1", y = "UMAP 2",
         caption = "Bins of >=20 cells on the joint pearson UMAP. Data: 04_umap_bins.csv") +
    theme_viz() +
    theme(panel.grid.major.y = element_blank(), axis.text = element_blank(),
          axis.line.x = element_blank())
  save_fig(fname, p, 7.6, 7.2)
}
umap_fig("pct_IB1", "fig_04_umap_sample_share.png",
         "Where each sample sits on the joint UMAP",
         sprintf("IB1 share per bin; whole-merge expectation is %.1f%%. Dark = IB1-dominated.", 100*EXP_IB1),
         "% IB1")
umap_fig("pct_strongB", "fig_04_umap_strongB.png",
         "B-marker-positive cells on the joint UMAP",
         "strongB rate per bin - marker evidence, independent of the reference transfer.",
         "% strongB")
cat("\nOK\n")
