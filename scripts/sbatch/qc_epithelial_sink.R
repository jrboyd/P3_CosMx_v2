# Is `epithelial` acting as a SINK class in the reference transfer?
#
# epithelial is 0.47% of the tonsil atlas but 53.5% of RB1 / 37.6% of IB1 calls,
# and 78% of strongB cells land there. Test with readouts independent of BOTH the
# transfer and the (unknown-provenance) lineage/strongB columns:
#   * Mean.PanCK / Max.PanCK -- PROTEIN channel; real epithelium is PanCK+
#   * keratin/epithelial TRANSCRIPTS vs B/T transcripts per called label
# If epithelial-called cells are not PanCK-high and not keratin-high, the label
# is wrong regardless of how lineage was defined.
suppressPackageStartupMessages({library(Seurat); library(Matrix)})
md <- readRDS("output/04_merge/merged_01_merged.Rds")@meta.data

cat("=== reference composition (annotation_level_1, ALL cells) ===\n")
ref <- readRDS("output/00_tonsil_reference/tonsil_ref_processed.Rds")@meta.data
rt <- sort(table(ref$annotation_level_1), decreasing = TRUE)
print(data.frame(n = as.vector(rt), pct = round(100*as.vector(rt)/sum(rt), 2),
                 row.names = names(rt)))
rm(ref); invisible(gc())

EPI <- c("KRT5","KRT8","KRT18","KRT19","KRT14","KRT17","EPCAM","SFN","S100A7")
BT  <- c("MS4A1","CD79A","CD3E","TRAC","PTPRC","CD74")

for (s in c("IB1","RB1")) {
  m <- md[md$merge_sample == s, ]
  cat("\n##########", s, "##########\n")
  cat("celltype_coarse composition (%):\n")
  print(round(100*sort(table(m$celltype_coarse), decreasing=TRUE)/nrow(m), 2))
  cat("\nprotein channels + score by called label (medians):\n")
  labs <- names(sort(table(m$celltype_coarse), decreasing=TRUE))
  tab <- t(sapply(labs, function(L) {
    k <- m$celltype_coarse == L
    c(n = sum(k),
      PanCK = round(median(m$Mean.PanCK[k], na.rm=TRUE), 0),
      Membrane = round(median(m$Mean.Membrane[k], na.rm=TRUE), 0),
      CD45 = round(median(m$Mean.CD45[k], na.rm=TRUE), 0),
      score = round(median(m$celltype_coarse.score[k], na.rm=TRUE), 3),
      pct_strongB = round(100*mean(m$strongB[k], na.rm=TRUE), 1))
  }))
  print(tab)
}

# transcript-level: are epithelial-called cells keratin+ or B+?
for (s in c("IB1","RB1")) {
  cat("\n##########", s, "transcript detection by called label ##########\n")
  seu <- readRDS(file.path("output/03_filter_and_process", paste0(s, "_03_HCA_cell_typed.Rds")))
  dat <- LayerData(seu[["sketch"]], layer = "counts")
  lab <- as.character(seu@meta.data[colnames(dat), "celltype_coarse"])
  g <- intersect(c(EPI, BT), rownames(dat))
  for (L in c("epithelial", "NBC_MBC", "CD4_T")) {
    k <- lab == L
    if (sum(k) < 50) next
    pct <- round(100 * Matrix::rowSums(dat[g, k, drop=FALSE] > 0) / sum(k), 1)
    cat("\n", L, " (n=", sum(k), ") %% cells detecting:\n", sep="")
    print(pct[order(-pct)])
  }
  rm(seu, dat); invisible(gc())
}
cat("\nOK\n")
