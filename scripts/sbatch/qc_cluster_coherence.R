# Evaluate cell-type calls at CLUSTER level, not per cell.
#
# Per-cell marker absence is expected at this dropout rate, so the right standard
# is: does a cell called NBC_MBC cluster with cells that DO show B markers?
# Uses the unsupervised scPearsonPCA clusters (from process_sample) and the
# InSituType clusters -- neither derived from the reference transfer.
#
# Decisive question: for clusters that are ENRICHED in strongB cells, what does
# celltype_coarse call them? A strongB-enriched cluster labelled epithelial is a
# real error by the cluster-level standard; labelled NBC_MBC means the calls are
# fine and only the per-cell view looked bad.
suppressPackageStartupMessages(library(Seurat))
md <- readRDS("output/04_merge/merged_01_merged.Rds")@meta.data

for (s in c("IB1", "RB1")) {
  m <- md[md$merge_sample == s, ]
  base_sb <- mean(m$strongB, na.rm = TRUE)
  base_mb <- mean(m$lineage == "B", na.rm = TRUE)
  cat("\n################################", s, "################################\n")
  cat(sprintf("baseline: strongB %.1f%%  marker-B %.1f%%  | called NBC_MBC %.1f%%  epithelial %.1f%%\n",
              100*base_sb, 100*base_mb,
              100*mean(m$celltype_coarse=="NBC_MBC"), 100*mean(m$celltype_coarse=="epithelial")))

  for (cl in c("pearson_clusters_low", "pearson_clusters_high", "inSitu_clusters")) {
    cat("\n===== clustering:", cl, "=====\n")
    g <- as.character(m[[cl]])
    tab <- do.call(rbind, lapply(sort(unique(g)), function(k) {
      i <- g == k
      tc <- sort(table(m$celltype_coarse[i]), decreasing = TRUE)
      data.frame(
        cluster = k, n = sum(i),
        pct_strongB = round(100*mean(m$strongB[i], na.rm=TRUE), 1),
        enrich_sB = round(mean(m$strongB[i], na.rm=TRUE)/base_sb, 2),
        pct_markerB = round(100*mean(m$lineage[i]=="B", na.rm=TRUE), 1),
        PanCK = round(median(m$Mean.PanCK[i], na.rm=TRUE), 0),
        called = names(tc)[1],
        called_pct = round(100*tc[1]/sum(i), 0),
        pct_NBC = round(100*mean(m$celltype_coarse[i]=="NBC_MBC"), 0),
        pct_epi = round(100*mean(m$celltype_coarse[i]=="epithelial"), 0)
      )
    }))
    tab <- tab[order(-tab$pct_strongB), ]
    print(tab, row.names = FALSE)

    # the cluster-level verdict
    bcl <- tab$enrich_sB >= 2 & tab$n >= 500
    if (any(bcl)) {
      cat(sprintf("\n  strongB-ENRICHED clusters (>=2x baseline, n>=500): %d clusters, %d cells\n",
                  sum(bcl), sum(tab$n[bcl])))
      cat(sprintf("    of their cells: %.1f%% called NBC_MBC, %.1f%% called epithelial\n",
                  weighted.mean(tab$pct_NBC[bcl], tab$n[bcl]),
                  weighted.mean(tab$pct_epi[bcl], tab$n[bcl])))
      cat("    dominant label per enriched cluster: ",
          paste(tab$called[bcl], collapse=", "), "\n")
    } else cat("\n  no clusters reach 2x strongB enrichment\n")
  }
}
cat("\nOK\n")
