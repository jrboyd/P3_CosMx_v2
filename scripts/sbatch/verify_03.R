# QC check on a step-03 output object (scripts/03_script_filter_and_process.R).
# Usage: sbatch --mem=48G scripts/sbatch/run_rscript.sh scripts/sbatch/verify_03.R <prefix>
suppressPackageStartupMessages(library(Seurat))

prefix <- commandArgs(trailingOnly = TRUE)[1]
f <- file.path("output/03_filter_and_process", paste0(prefix, "_03_HCA_cell_typed.Rds"))
stopifnot(file.exists(f))
cat("== ", prefix, " ==\n", f, "  (", format(file.size(f), big.mark = ","), " bytes)\n", sep = "")

s <- readRDS(f)
cat("dims:", paste(dim(s), collapse = " x "), "| default assay:", DefaultAssay(s), "\n")
cat("assays:", paste(Assays(s), collapse = ","), "\n")
cat("sketch layers:", paste(Layers(s[["sketch"]]), collapse = ","),
    "| sketch cells:", length(Cells(s[["sketch"]])), "\n")
cat("reductions:", paste(Reductions(s), collapse = ","), "\n")
cat("pca.full:", paste(dim(Embeddings(s, "pca.full")), collapse = " x "), "\n\n")

md <- s@meta.data
for (cc in c("celltype_coarse", "celltype_fine",
             "sketch_celltype_coarse", "sketch_celltype_fine")) {
  cat(sprintf("%-24s non-NA %8d / %d\n", cc, sum(!is.na(md[[cc]])), nrow(md)))
}
stopifnot(!any(is.na(md$celltype_coarse)), !any(is.na(md$celltype_fine)))

cat("\ncoarse labels (all cells):\n")
print(sort(table(md$celltype_coarse), decreasing = TRUE))
cat("\nn distinct fine labels:", length(unique(md$celltype_fine)), "\n")
cat("smallest fine labels:\n")
print(head(sort(table(md$celltype_fine)), 8))

sk <- !is.na(md$sketch_celltype_coarse)
cat("\nsketch-vs-projected agreement on sketch cells: coarse",
    round(mean(as.character(md$celltype_coarse[sk]) == md$sketch_celltype_coarse[sk]) * 100, 1), "% | fine",
    round(mean(as.character(md$celltype_fine[sk]) == md$sketch_celltype_fine[sk]) * 100, 1), "%\n")
cat("\nprediction score (coarse, all cells):\n")
print(summary(md$celltype_coarse.score))
cat("\nOK\n")
