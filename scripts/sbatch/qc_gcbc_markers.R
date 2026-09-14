# Where do GC-B-cell-like cells end up in celltype_coarse?
#
# The coarse transfer calls almost no GCBC (2 cells in RB1, 0 in IB1) even though
# these are tonsils. The question is not whether GCBC is used as a label, but
# whether cells that actually express GC markers land somewhere defensible
# (NBC_MBC) rather than somewhere nonsensical (epithelial, myeloid, FDC).
#
# Marker scoring uses the `sketch` assay (100k in-memory cells per sample) so no
# BPCells row-scan is needed; those are also the cells whose labels came straight
# from the anchor transfer, so sketch_celltype_* and celltype_* can be compared.
#
# Usage: sbatch --mem=96G scripts/sbatch/run_rscript.sh \
#          scripts/sbatch/qc_gcbc_markers.R IB1 RB1

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

prefixes <- commandArgs(trailingOnly = TRUE)
stopifnot(length(prefixes) > 0)
OUT_DIR <- "output/03_filter_and_process"

# canonical germinal-centre B markers (AID, BCL6 program, DZ/LZ, proliferation)
GC <- c("AICDA", "BCL6", "MEF2B", "RGS13", "GCSAM", "MME", "LMO2", "SERPINA9",
        "MYBL1", "CD38", "MKI67", "CXCR4", "CD83", "FAS", "EZH2", "POU2AF1")
# contrast sets: naive/memory B (where GC cells are apparently being absorbed)
NBC <- c("TCL1A", "IGHD", "FCER2", "SELL", "CR2", "CD27", "MS4A1", "CD19", "CD79A")
EPI <- c("KRT5", "KRT19", "EPCAM", "KRT8", "KRT18", "KRT14")

score_set <- function(dat, genes) {
  g <- intersect(genes, rownames(dat))
  m <- as.matrix(dat[g, , drop = FALSE])
  m <- t(scale(t(m)))                       # z-score each gene across cells
  m[!is.finite(m)] <- 0
  list(score = colMeans(m), used = g, missing = setdiff(genes, rownames(dat)))
}

for (prefix in prefixes) {
  f <- file.path(OUT_DIR, paste0(prefix, "_03_HCA_cell_typed.Rds"))
  cat("\n########## ", prefix, " ##########\n", sep = "")
  stopifnot(file.exists(f))
  seu <- readRDS(f)
  md <- seu@meta.data

  # ---- 1. metadata only: where do GC-ish FINE calls sit in COARSE? ----------
  gc_fine <- grep("GC|DZ|LZ|Cycling", unique(as.character(md$celltype_fine)),
                  value = TRUE)
  cat("\n-- fine labels treated as GC-like (n=", length(gc_fine), "):\n", sep = "")
  print(gc_fine)
  is_gc_fine <- as.character(md$celltype_fine) %in% gc_fine
  cat("\ncells with a GC-like fine call:", sum(is_gc_fine), "of", nrow(md), "\n")
  cat("their coarse assignment:\n")
  print(sort(table(as.character(md$celltype_coarse)[is_gc_fine]), decreasing = TRUE))

  # ---- 2. expression-based, independent of the transfer ---------------------
  dat <- LayerData(seu[["sketch"]], layer = "data")
  cells <- colnames(dat)
  sc_gc <- score_set(dat, GC); sc_nbc <- score_set(dat, NBC); sc_epi <- score_set(dat, EPI)
  cat("\n-- GC markers in panel (", length(sc_gc$used), "/", length(GC), "): ",
      paste(sc_gc$used, collapse = ", "), "\n", sep = "")
  if (length(sc_gc$missing)) cat("   absent from panel: ", paste(sc_gc$missing, collapse = ", "), "\n", sep = "")

  lab <- as.character(md[cells, "celltype_coarse"])
  cat("\nmean GC score by coarse label (sketch cells):\n")
  tab <- data.frame(
    n = as.vector(table(lab)[sort(unique(lab))]),
    GC = round(tapply(sc_gc$score, lab, mean)[sort(unique(lab))], 3),
    NBC = round(tapply(sc_nbc$score, lab, mean)[sort(unique(lab))], 3),
    EPI = round(tapply(sc_epi$score, lab, mean)[sort(unique(lab))], 3)
  )
  print(tab[order(-tab$GC), ])

  # ---- 3. the actual question: top GC-scoring cells -> which coarse label? ---
  for (q in c(0.99, 0.95)) {
    thr <- quantile(sc_gc$score, q)
    top <- sc_gc$score >= thr
    cat("\ntop ", round((1 - q) * 100), "% GC-score cells (n=", sum(top),
        "), coarse assignment:\n", sep = "")
    print(round(100 * sort(table(lab[top]), decreasing = TRUE) / sum(top), 1))
  }
  thr <- quantile(sc_gc$score, 0.99)
  top <- sc_gc$score >= thr
  cat("\ntop 1% GC-score cells, fine assignment (top 12):\n")
  print(head(sort(table(as.character(md[cells, "celltype_fine"])[top]), decreasing = TRUE), 12))

  rm(seu, dat); invisible(gc())
}
cat("\nOK\n")
