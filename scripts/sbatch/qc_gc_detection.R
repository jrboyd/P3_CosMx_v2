# Is the missing GC signal biology or detection sensitivity?
#
# qc_gcbc_calibrate.R showed the GC score cleanly identifies GCBC in the tonsil
# atlas (top 1% = 99.9% GCBC) but has almost no dynamic range in IB1/RB1. Two
# explanations: the tissue genuinely lacks a GC program, or CosMx does not detect
# these transcripts at scRNA-seq sensitivity. Compare per-gene DETECTION RATE
# (fraction of cells with count > 0) in the query against the reference's GCBC
# cells, with non-GC positive controls to show detection is otherwise sane.
#
# Usage: sbatch --mem=128G scripts/sbatch/run_rscript.sh \
#          scripts/sbatch/qc_gc_detection.R IB1 RB1

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

prefixes <- commandArgs(trailingOnly = TRUE)
OUT_DIR <- "output/03_filter_and_process"
REF_RDS <- "output/00_tonsil_reference/tonsil_ref_processed.Rds"

GC4 <- c("AICDA", "BCL6", "LMO2", "POU2AF1")
GC_REST <- c("RGS13", "MME", "CD38", "MKI67", "CXCR4", "CD83", "FAS", "EZH2")
CTRL <- c("MS4A1", "CD79A", "CD3E", "KRT5", "CD74", "PTPRC")
GENES <- c(GC4, GC_REST, CTRL)

det <- function(counts, genes) {
  g <- intersect(genes, rownames(counts))
  m <- counts[g, , drop = FALSE]
  data.frame(
    gene = g,
    pct_detected = round(100 * Matrix::rowSums(m > 0) / ncol(m), 2),
    mean_count = round(Matrix::rowMeans(m), 3),
    row.names = NULL
  )
}

cat("##################### REFERENCE #####################\n")
ref <- readRDS(REF_RDS)
rc <- LayerData(ref[["sketch"]], layer = "counts")
gcbc <- rownames(ref@meta.data)[ref@meta.data$annotation_level_1 == "GCBC"]
gcbc <- intersect(gcbc, colnames(rc))
cat("ref sketch cells:", ncol(rc), "| GCBC among them:", length(gcbc), "\n")
cat("median counts/cell (all ref sketch):", median(Matrix::colSums(rc)), "\n\n")
d_all <- det(rc, GENES); d_gc <- det(rc[, gcbc, drop = FALSE], GENES)
ref_tab <- data.frame(gene = d_all$gene,
                      ref_all_pct = d_all$pct_detected,
                      ref_GCBC_pct = d_gc$pct_detected,
                      ref_GCBC_mean = d_gc$mean_count)
print(ref_tab)
rm(ref, rc); invisible(gc())

for (prefix in prefixes) {
  cat("\n#####################", prefix, "#####################\n")
  seu <- readRDS(file.path(OUT_DIR, paste0(prefix, "_03_HCA_cell_typed.Rds")))
  qc <- LayerData(seu[["sketch"]], layer = "counts")
  cat("sketch cells:", ncol(qc), "| median counts/cell:", median(Matrix::colSums(qc)), "\n")
  bcell <- rownames(seu@meta.data)[seu@meta.data$celltype_coarse %in% c("NBC_MBC", "GCBC")]
  bcell <- intersect(bcell, colnames(qc))
  cat("B-compartment cells (NBC_MBC/GCBC):", length(bcell), "\n\n")
  d_all <- det(qc, GENES); d_b <- det(qc[, bcell, drop = FALSE], GENES)
  out <- merge(ref_tab,
               data.frame(gene = d_all$gene,
                          qry_all_pct = d_all$pct_detected,
                          qry_B_pct = d_b$pct_detected,
                          qry_B_mean = d_b$mean_count),
               by = "gene", all = TRUE)
  out$set <- ifelse(out$gene %in% GC4, "GC-strict",
                    ifelse(out$gene %in% GC_REST, "GC-other", "control"))
  print(out[order(out$set, -out$ref_GCBC_pct), ])
  rm(seu, qc); invisible(gc())
}
cat("\nREAD: if GC-strict genes are detected in a decent %% of ref GCBC cells but",
    "near 0%% of query cells while controls are comparable, the absent GC signal is",
    "a detection limit, not biology.\n")
cat("\nOK\n")
