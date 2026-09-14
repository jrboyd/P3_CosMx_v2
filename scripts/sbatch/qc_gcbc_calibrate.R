# Calibrate the GC marker score before trusting qc_gcbc_markers.R.
#
# qc_gcbc_markers.R found the top GC-scoring query cells to be mostly epithelial,
# but 4 of 16 canonical GC genes are absent from the 6K panel and several of the
# remaining 12 are not GC-specific (MKI67/EZH2 = any cycling cell, RGS13 = mast,
# MME/CD10 = epithelium). Two checks:
#
#  1. REFERENCE calibration -- score the same 12 panel genes on the tonsil atlas,
#     where the true GCBC label is known. If atlas GCBC does not separate from
#     epithelial on this score, the score cannot answer the question and the
#     query result is an artifact of marker choice, not evidence of misassignment.
#
#  2. B-LINEAGE-restricted query check -- within B-lineage cells only, do the
#     GC-program-high cells sit in NBC_MBC? This is the question that actually
#     matters, with epithelium/mast/myeloid confounds removed.
#
# Usage: sbatch --mem=128G scripts/sbatch/run_rscript.sh \
#          scripts/sbatch/qc_gcbc_calibrate.R IB1 RB1

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

prefixes <- commandArgs(trailingOnly = TRUE)
OUT_DIR <- "output/03_filter_and_process"
REF_RDS <- "output/00_tonsil_reference/tonsil_ref_processed.Rds"
REF_ANNO_COL_COARSE <- "annotation_level_1"

GC <- c("AICDA", "BCL6", "MEF2B", "RGS13", "GCSAM", "MME", "LMO2", "SERPINA9",
        "MYBL1", "CD38", "MKI67", "CXCR4", "CD83", "FAS", "EZH2", "POU2AF1")
# GC genes with the fewest non-GC confounders, for a stricter score
GC_STRICT <- c("AICDA", "BCL6", "LMO2", "POU2AF1")
BLIN <- c("MS4A1", "CD79A", "CD79B", "CD19")

zscore_mean <- function(dat, genes) {
  g <- intersect(genes, rownames(dat))
  if (!length(g)) return(list(score = rep(NA_real_, ncol(dat)), used = g))
  m <- as.matrix(dat[g, , drop = FALSE])
  m <- t(scale(t(m)))
  m[!is.finite(m)] <- 0
  list(score = colMeans(m), used = g)
}

report <- function(dat, lab, tag) {
  s <- zscore_mean(dat, GC)
  ss <- zscore_mean(dat, GC_STRICT)
  keep <- sort(unique(lab))
  tab <- data.frame(
    n = as.vector(table(lab)[keep]),
    GC12 = round(tapply(s$score, lab, mean)[keep], 3),
    GC4 = round(tapply(ss$score, lab, mean)[keep], 3)
  )
  cat("\n[", tag, "] mean GC score by label (GC12 = 12 panel genes, GC4 = AICDA/BCL6/LMO2/POU2AF1):\n", sep = "")
  print(tab[order(-tab$GC4), ])
  for (nm in c("GC12", "GC4")) {
    sc <- if (nm == "GC12") s$score else ss$score
    top <- sc >= quantile(sc, 0.99)
    cat("\n[", tag, "] top 1% by ", nm, " -> label %:\n", sep = "")
    print(round(100 * sort(table(lab[top]), decreasing = TRUE) / sum(top), 1))
  }
  invisible(NULL)
}

# ------------------------------------------------- 1. reference calibration ----
cat("##################### REFERENCE (ground truth labels) #####################\n")
ref <- readRDS(REF_RDS)
rdat <- LayerData(ref[["sketch"]], layer = "data")
rlab <- as.character(ref@meta.data[colnames(rdat), REF_ANNO_COL_COARSE])
cat("reference sketch cells:", ncol(rdat), "| GC genes present:",
    length(intersect(GC, rownames(rdat))), "/", length(GC),
    "| strict:", paste(intersect(GC_STRICT, rownames(rdat)), collapse = ","), "\n")
report(rdat, rlab, "REF")
cat("\nVERDICT INPUT: if GCBC is not at/near the top above, the score is not GC-specific.\n")
rm(ref, rdat); invisible(gc())

# --------------------------------------- 2. B-lineage-restricted query check ----
for (prefix in prefixes) {
  f <- file.path(OUT_DIR, paste0(prefix, "_03_HCA_cell_typed.Rds"))
  cat("\n#####################", prefix, "#####################\n")
  seu <- readRDS(f)
  dat <- LayerData(seu[["sketch"]], layer = "data")
  md <- seu@meta.data
  lab <- as.character(md[colnames(dat), "celltype_coarse"])

  report(dat, lab, paste0(prefix, " all sketch cells"))

  # B lineage two ways: by transferred coarse label, and by marker positivity
  b_lab <- lab %in% c("NBC_MBC", "GCBC", "PC")
  bs <- zscore_mean(dat, BLIN)
  b_mrk <- bs$score >= quantile(bs$score, 0.75)
  cat("\nB-lineage by label:", sum(b_lab), "| by marker (top 25% B score):", sum(b_mrk),
      "| overlap:", sum(b_lab & b_mrk), "\n")

  if (sum(b_mrk) > 100) {
    report(dat[, b_mrk, drop = FALSE], lab[b_mrk],
           paste0(prefix, " B-marker-high cells only"))
  }
  rm(seu, dat); invisible(gc())
}
cat("\nOK\n")
