#
# Step 5: sample integration.
#
# Step 4 showed the joint clustering is almost completely sample-segregated: at
# res 0.4, 0 of 7 clusters fall within 1.5x of the expected IB1 share, and the
# joint UMAP is nearly sample-disjoint. Integration is worth trying ONLY if it
# buys confident B-cell identification in IB1 by co-clustering with RB1 B cells
# WITHOUT dissolving rare populations. Both are measured in 05b_evaluate.
#
# Method choice: RPCA is the primary. Seurat recommends it when a substantial
# fraction of cells in one dataset have no counterpart in the other -- exactly
# IB1, ~37% of whose retained cells are the low-information InSituType cluster i.
# CCA is run as the aggressive comparator: it forces more mixing, so it bounds
# how much rare-cell structure aggressive correction costs. (harmony is not
# installed on this cluster.)
#
# Sketch-based so 1.8M cells stay tractable: integrate a per-sample-balanced
# sketch, then project the correction and the cluster labels back to every cell.
#
# Usage: sbatch scripts/sbatch/05_integrate_wrapper.sh

suppressPackageStartupMessages({
  library(Seurat); library(BPCells); library(dplyr); library(data.table); setDTthreads(1)
  library(ggplot2)
})
options(future.globals.maxSize = 16 * 1024^3)

IN_RDS <- "output/04_merge/merged_processed.Rds"
OUT_DIR <- "output/05_integrate"
ASSAY <- "RNA"; SAMPLE_COL <- "merge_sample"
N_DIMS <- 30
N_SKETCH_PER_SAMPLE <- 50000     # per layer => balanced across samples
RES_LOW <- 0.4; RES_HIGH <- 0.8
set.seed(1)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

done_f <- function(x) file.path(OUT_DIR, paste0("merged_", x, ".done"))
is_done <- function(x) file.exists(done_f(x))
mark <- function(x) { writeLines(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), done_f(x))
                      message("marked complete: ", done_f(x)) }

SKETCH_RDS <- file.path(OUT_DIR, "merged_sketch_integrated.Rds")
FULL_RDS <- file.path(OUT_DIR, "merged_integrated.Rds")

## ------------------------------------------- 1. split + sketch --------------
if (!is_done("01_sketch")) {
  message("reading ", IN_RDS)
  obj <- readRDS(IN_RDS)
  DefaultAssay(obj) <- ASSAY
  message("cells: ", ncol(obj), " | counts class ",
          class(LayerData(obj[[ASSAY]], layer = "counts"))[1])

  # split counts by sample so IntegrateLayers has per-sample layers
  obj[[ASSAY]] <- split(obj[[ASSAY]], f = obj@meta.data[[SAMPLE_COL]])
  message("layers after split: ", paste(Layers(obj[[ASSAY]]), collapse = ", "))
  obj <- NormalizeData(obj)                    # lazy on BPCells
  obj <- FindVariableFeatures(obj)

  # SketchData samples per data layer -> equal cells per sample, which stops the
  # 22/78 cell-count imbalance from driving the correction
  obj <- SketchData(obj, ncells = N_SKETCH_PER_SAMPLE, method = "LeverageScore",
                    sketched.assay = "sketch")
  DefaultAssay(obj) <- "sketch"
  print(table(obj@meta.data[[SAMPLE_COL]][colnames(obj) %in% Cells(obj[["sketch"]])]))
  obj <- FindVariableFeatures(obj)
  obj <- ScaleData(obj)
  obj <- RunPCA(obj, npcs = N_DIMS, reduction.name = "pca.sketch")
  saveRDS(obj, SKETCH_RDS); mark("01_sketch")
} else message("skip sketch (done)")

## ------------------------------------------- 2. integrate -------------------
METHODS <- list(
  rpca = list(fn = RPCAIntegration, red = "integrated.rpca"),
  cca  = list(fn = CCAIntegration,  red = "integrated.cca")
)
for (m in names(METHODS)) {
  if (is_done(paste0("02_", m))) { message("skip ", m, " (done)"); next }
  obj <- readRDS(SKETCH_RDS)
  DefaultAssay(obj) <- "sketch"
  red <- METHODS[[m]]$red
  message("=== IntegrateLayers: ", m, " -> ", red, " ===")
  t0 <- Sys.time()
  obj <- IntegrateLayers(obj, method = METHODS[[m]]$fn, orig.reduction = "pca.sketch",
                         new.reduction = red, assay = "sketch", verbose = TRUE)
  message(m, " integration took ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min")
  obj <- FindNeighbors(obj, reduction = red, dims = 1:N_DIMS)
  obj <- FindClusters(obj, resolution = RES_LOW,
                      cluster.name = paste0("int_", m, "_clusters_low"))
  obj <- FindClusters(obj, resolution = RES_HIGH,
                      cluster.name = paste0("int_", m, "_clusters_high"))
  obj <- RunUMAP(obj, reduction = red, dims = 1:N_DIMS, return.model = TRUE,
                 reduction.name = paste0("umap.", m))
  saveRDS(obj, file.path(OUT_DIR, paste0("merged_sketch_", m, ".Rds")))
  mark(paste0("02_", m))
  rm(obj); invisible(gc())
}

## ------------------------------------------- 3. project to all cells --------
if (!is_done("03_projected")) {
  obj <- readRDS(file.path(OUT_DIR, "merged_sketch_rpca.Rds"))
  # carry the cca results onto the same object
  cca <- readRDS(file.path(OUT_DIR, "merged_sketch_cca.Rds"))
  obj[["integrated.cca"]] <- cca[["integrated.cca"]]
  obj[["umap.cca"]] <- cca[["umap.cca"]]
  for (cc in grep("^int_cca_clusters", colnames(cca@meta.data), value = TRUE)) {
    obj[[cc]] <- cca@meta.data[[cc]][match(colnames(obj), colnames(cca))]
  }
  rm(cca); invisible(gc())

  for (m in names(METHODS)) {
    red <- METHODS[[m]]$red
    message("=== ProjectIntegration + ProjectData: ", m, " ===")
    obj <- ProjectIntegration(obj, sketched.assay = "sketch", assay = ASSAY,
                              reduction = red, verbose = TRUE)
    full_red <- paste0(red, ".full")
    stopifnot(full_red %in% Reductions(obj))
    obj <- ProjectData(
      obj, assay = ASSAY, sketched.assay = "sketch",
      sketched.reduction = red, full.reduction = full_red, dims = 1:N_DIMS,
      refdata = setNames(
        list(paste0("int_", m, "_clusters_low"), paste0("int_", m, "_clusters_high")),
        paste0("int_", m, c("_low_full", "_high_full"))),
      verbose = TRUE)
  }
  DefaultAssay(obj) <- ASSAY
  # drop the dense sketch scale.data before saving
  LayerData(obj[["sketch"]], layer = "scale.data") <- NULL
  saveRDS(obj, FULL_RDS)
  message("wrote ", FULL_RDS, " (", format(file.size(FULL_RDS), big.mark = ","), " bytes)")
  print(sapply(grep("_full$", colnames(obj@meta.data), value = TRUE),
               function(k) sum(!is.na(obj@meta.data[[k]]))))
  mark("03_projected")
} else message("skip projection (done)")

message("Done. Outputs in: ", normalizePath(OUT_DIR))
