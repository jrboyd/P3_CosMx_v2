#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------------------
## step 22a -- extract the ligand-receptor genes for every analysable cell.
##
## WHY A SEPARATE EXTRACTION. Proximity co-expression is irreducibly a per-CELL statistic:
## it asks whether a ligand-positive cell of type S sits near a receptor-positive cell of
## type T. Pseudobulk cannot answer it, so the existing bin x celltype matrices are no use
## here. And the existing ~100-gene extract (output/16_regions/panel_counts_all.Rds) covers
## only 60 of the 1,125 genes the curated pair list needs.
##
## THE PAIR LIST. scDiffCom's curated LRI_human, restricted to SIMPLE 1:1 interactions whose
## ligand AND receptor are both on the CosMx 6K panel: 2,315 pairs over 1,125 genes. The
## multi-subunit complexes are dropped rather than approximated -- a complex scored on one
## subunit is a different claim from the one the database makes.
##
## WHICH CELLS. The symmetric `matched` QC tier, defect cells dropped per cell, and only
## cells carrying a step-19 cell-type label. Step 19's labels are the project's best cell
## identity and nothing downstream had yet consumed them (gap G1), so this step uses them.
## That also means every cell here sits in the step-19 high-confidence or rescued set, with
## `confidence` carried through so a statement can be restricted to `high` later.
##
## Written as its own job because it is the only expensive part of step 22: it touches the
## full pre-QC BPCells matrices, which is what 16_script_panel_counts.R needed 140 GB for.
## Everything after this runs off the output in minutes.
##
## Usage: sbatch --mem=160G --cpus-per-task=2 --time=3:00:00 \
##          scripts/sbatch/run_rscript.sh scripts/22_script_lr_extract.R
## ---------------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(data.table); setDTthreads(1)
})
options(width = 200)

SAMPLES <- c("IB1","RB1")
PRE_RDS <- "output/01_import_v1/seuratObject_%s_BPCells_full.Rds"
OUT     <- "output/22_lr"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
if (!file.exists("scripts/lib/cluster_labels.R")) stop("run from the project root", call. = FALSE)

PAIRS <- fread(file.path(OUT, "lr_pairs_on_panel.csv"))
GENES <- sort(unique(c(PAIRS$ligand, PAIRS$receptor)))
cat(sprintf("pairs %d over %d genes\n", nrow(PAIRS), length(GENES)))

## ---- the analysable cell set: matched tier, non-defect, step-19 labelled ----------------
md <- fread("output/16_regions/cells_meta_all.csv.gz",
            select = c("cell_id","samp","fov_id","x_slide_mm","y_slide_mm",
                       "nCount_RNA","is_defect","tier_matched"))
md <- md[tier_matched == TRUE & is_defect == FALSE]
lab <- fread("output/19_celltypes/all_cells_celltypes.csv.gz",
             select = c("cell_id","celltype_cluster","celltype_argmax","confidence"))
md <- merge(md, lab, by = "cell_id")
md <- md[!is.na(celltype_cluster)]
cat(sprintf("analysable cells: %s\n", format(nrow(md), big.mark = ",")))
print(md[, .(cells = .N), by = .(samp, confidence)][order(samp, -cells)])

## ---- pull the counts, one sample at a time ---------------------------------------------
## Column-subset immediately after the row-subset: the row-subset is what BPCells streams
## cheaply, and holding the full-width matrix is what costs the memory.
parts <- lapply(SAMPLES, function(SM) {
  f <- sprintf(PRE_RDS, SM)
  if (!file.exists(f)) stop("missing ", f, call. = FALSE)
  cat("\nreading ", SM, " ... ", sep = "")
  seu <- readRDS(f)
  cn  <- GetAssayData(seu, assay = "RNA", layer = "counts")
  g   <- intersect(GENES, rownames(cn))
  keep <- intersect(md[samp == SM, cell_id], colnames(cn))
  M <- as(cn[g, keep, drop = FALSE], "dgCMatrix")
  cat(sprintf("%d genes x %s cells\n", nrow(M), format(ncol(M), big.mark = ",")))
  rm(seu, cn); gc()
  M
})
g_common <- Reduce(intersect, lapply(parts, rownames))
M <- do.call(cbind, lapply(parts, function(x) x[g_common, , drop = FALSE]))
cat(sprintf("\ncombined: %d genes x %s cells\n", nrow(M), format(ncol(M), big.mark = ",")))
miss <- setdiff(GENES, g_common)
if (length(miss)) cat(sprintf("genes absent from the object (%d): %s\n",
                              length(miss), paste(head(miss, 20), collapse = ", ")))

md <- md[cell_id %in% colnames(M)]
setkey(md, cell_id); md <- md[colnames(M)]
stopifnot(identical(md$cell_id, colnames(M)))

saveRDS(M, file.path(OUT, "lr_counts.Rds"))
fwrite(md, file.path(OUT, "lr_cells.csv.gz"), compress = "gzip")
fwrite(PAIRS[ligand %in% g_common & receptor %in% g_common],
       file.path(OUT, "lr_pairs_used.csv"))
cat("wrote lr_counts.Rds, lr_cells.csv.gz, lr_pairs_used.csv\n")
cat(sprintf("pairs usable after gene availability: %d\n",
            nrow(PAIRS[ligand %in% g_common & receptor %in% g_common])))
writeLines(format(Sys.time()), file.path(OUT, "lr_extract.done"))
cat("\nOK\n")
