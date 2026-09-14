# Post-hoc repair of step-03 outputs produced before the label fix in
# scripts/03_script_filter_and_process.R.
#
# TransferData() inside ProjectData() rewrote '_' to '-' in the propagated
# celltype_coarse/celltype_fine columns ("Names of identity class contain
# underscores"), diverging them from the reference vocabulary and from the
# sketch_celltype_* columns. This maps them back, exactly as the patched script
# now does at the end of transfer_reference().
#
# Usage: sbatch --mem=64G scripts/sbatch/run_rscript.sh \
#          scripts/sbatch/remap_labels_03.R IB1 [RB1 ...]

suppressPackageStartupMessages(library(Seurat))

prefixes <- commandArgs(trailingOnly = TRUE)
stopifnot(length(prefixes) > 0)

OUT_DIR <- "output/03_filter_and_process"
REF_RDS <- "output/00_tonsil_reference/tonsil_ref_processed.Rds"
REF_ANNO_COL_COARSE <- "annotation_level_1"
REF_ANNO_COL_FINE <- "annotation_20230508"

restore_ref_labels <- function(x, ref_levels) {
  ref_levels <- unique(as.character(ref_levels))
  key <- setNames(ref_levels, gsub("_", "-", ref_levels, fixed = TRUE))
  key <- key[!names(key) %in% ref_levels]          # no collisions with real labels
  key <- key[!duplicated(names(key))]
  x <- as.character(x)
  hit <- x %in% names(key)
  x[hit] <- unname(key[x[hit]])
  x
}

# full reference vocabulary (not just the labels this sample happened to predict)
message("loading reference vocabulary ...")
ref_md <- readRDS(REF_RDS)@meta.data
coarse_levels <- unique(as.character(ref_md[[REF_ANNO_COL_COARSE]]))
fine_levels <- unique(as.character(ref_md[[REF_ANNO_COL_FINE]]))
rm(ref_md); invisible(gc())
message("reference labels: ", length(coarse_levels), " coarse, ", length(fine_levels), " fine")

for (prefix in prefixes) {
  f <- file.path(OUT_DIR, paste0(prefix, "_03_HCA_cell_typed.Rds"))
  cat("\n=====", prefix, "=====\n", f, "\n")
  stopifnot(file.exists(f))

  seu <- readRDS(f)
  before_c <- as.character(seu$celltype_coarse)
  before_f <- as.character(seu$celltype_fine)
  after_c <- restore_ref_labels(before_c, coarse_levels)
  after_f <- restore_ref_labels(before_f, fine_levels)

  cat("coarse cells relabelled:", sum(before_c != after_c), "of", length(after_c), "\n")
  cat("fine   cells relabelled:", sum(before_f != after_f), "of", length(after_f), "\n")
  chg_c <- unique(cbind(before_c, after_c)[before_c != after_c, , drop = FALSE])
  chg_f <- unique(cbind(before_f, after_f)[before_f != after_f, , drop = FALSE])
  if (nrow(chg_c)) { cat("coarse mapping:\n"); print(chg_c, quote = FALSE) }
  if (nrow(chg_f)) { cat("fine mapping:\n"); print(chg_f, quote = FALSE) }

  if (!nrow(chg_c) && !nrow(chg_f)) {
    cat("already consistent with the reference vocabulary -- not rewriting\n")
    next
  }

  seu$celltype_coarse <- after_c
  seu$celltype_fine <- after_f
  stopifnot(
    all(seu$celltype_coarse %in% coarse_levels),
    all(seu$celltype_fine %in% fine_levels),
    !any(is.na(seu$celltype_coarse)), !any(is.na(seu$celltype_fine))
  )
  # sketch_celltype_* come straight from TransferData and were never mangled
  sk <- !is.na(seu$sketch_celltype_coarse)
  cat("sketch-vs-projected agreement after remap: coarse ",
      round(mean(seu$celltype_coarse[sk] == seu$sketch_celltype_coarse[sk]) * 100, 1),
      "% | fine ",
      round(mean(seu$celltype_fine[sk] == seu$sketch_celltype_fine[sk]) * 100, 1), "%\n", sep = "")

  # write-then-rename so an interrupted job cannot truncate the original
  tmp <- paste0(f, ".tmp")
  saveRDS(seu, tmp)
  stopifnot(file.size(tmp) > 0)
  file.rename(tmp, f)
  cat("rewrote ", f, " (", format(file.size(f), big.mark = ","), " bytes)\n", sep = "")
  rm(seu); invisible(gc())
}
cat("\nOK\n")
