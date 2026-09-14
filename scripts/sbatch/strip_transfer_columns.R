# Remove the disabled reference-transfer columns from the DOWNSTREAM objects.
#
# The transfer was disabled in step 3 (2026-08-06) because its labels anti-track
# marker evidence; scripts/04_script_merge_and_process.R now drops these columns
# on read, but objects saved before that still carry them and would mislead
# anyone who opened them.
#
# The step-3 objects (*_03_HCA_cell_typed.Rds) are deliberately NOT touched: they
# are the record, and keeping them means a re-run of step 4 can restore these
# columns without re-enabling the transfer.
#
# Usage: sbatch --mem=200G scripts/sbatch/run_rscript.sh scripts/sbatch/strip_transfer_columns.R

suppressPackageStartupMessages(library(Seurat))

DROP <- c("celltype_coarse", "celltype_fine", "celltype_coarse.score",
          "celltype_fine.score", "sketch_celltype_coarse", "sketch_celltype_fine",
          "sample_celltype_coarse", "sample_celltype_fine")
TARGETS <- c("output/04_merge/merged_processed.Rds",
             "output/04_merge/merged_01_merged.Rds",
             "output/05_integrate/merged_integrated.Rds",
             "output/05_integrate/merged_sketch_integrated.Rds",
             "output/05_integrate/merged_sketch_rpca.Rds",
             "output/05_integrate/merged_sketch_cca.Rds")

for (f in TARGETS) {
  if (!file.exists(f)) { cat("skip (absent):", f, "\n"); next }
  cat("\n==", f, "\n")
  o <- readRDS(f)
  present <- intersect(DROP, colnames(o@meta.data))
  if (!length(present)) { cat("  nothing to drop\n"); rm(o); invisible(gc()); next }
  cat("  dropping:", paste(present, collapse = ", "), "\n")
  for (cc in present) o[[cc]] <- NULL
  stopifnot(!any(DROP %in% colnames(o@meta.data)))
  tmp <- paste0(f, ".tmp")            # write-then-rename, never truncate in place
  saveRDS(o, tmp); stopifnot(file.size(tmp) > 0)
  file.rename(tmp, f)
  cat("  rewrote (", format(file.size(f), big.mark = ","), " bytes ), ",
      ncol(o@meta.data), " metadata columns remain\n", sep = "")
  rm(o); invisible(gc())
}
cat("\nstep-3 objects intentionally left intact:\n")
for (f in Sys.glob("output/03_filter_and_process/*_03_HCA_cell_typed.Rds")) cat("  ", f, "\n")
cat("\nOK\n")
