#
# use cosmx specific processing workflow
# https://github.com/Nanostring-Biostats/CosMx-Analysis-Scratch-Space/tree/Main/_code/scPearsonPCA

# Claude Opus 4.8
# =============================================================================
# CosMx 6K — B-cell processing & merge, MEMORY-BOUNDED (Seurat v5 sketch)
#
# Alternate to 10_bcell_process_and_merge.R. Strategy:
#   * Keep full-resolution counts DISK-BACKED (BPCells) — never densify them.
#   * Cluster/UMAP on a leverage-score SKETCH (~50k cells) that over-represents
#     rare states, so the dense ScaleData matrix stays tiny.
#   * ProjectData() propagates the sketch PCA/UMAP/cluster/celltype labels back
#     to every full-resolution cell.
#   * Cell-type = tonsil-atlas label transfer, done on the sketch then projected.
#   * Save DietSeurat objects (no scale.data, counts/data stay BPCells) => tiny.
#
# Two products:
#   (1) per-sample processed objects: cluster + UMAP + celltype   -> *_processed.Rds
#   (2) merged object: per-sample labels/UMAP preserved AS METADATA, plus a
#       freshly recomputed joint PCA/cluster/UMAP                  -> merged_processed.Rds
#
# NOTE ON "MINIMAL SIZE": the saved .Rds keeps counts/data as BPCells, i.e. it
# stores POINTERS to on-disk matrix directories (step-9 stores + the sketch/merge
# stores written under OUT_DIR). The .Rds itself is tiny but is NOT portable on
# its own — keep IN_DIR and OUT_DIR intact, or re-materialize before moving.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(BPCells)
  library(dplyr)
  library(data.table); setDTthreads(1)
  library(scPearsonPCA)
  library(ggplot2)
  library(ggrepel)
})
options(future.globals.maxSize = 8 * 1024^3)   # allow larger worker payloads

# ------------------------------------------------------------------ CONFIG ----
RDS_ROOT      <- "/netfiles/volaric_research/DLBCL/CosMx/export1/"
RDS_PATHS     <- dir(RDS_ROOT, pattern = "RDS$", full.names = TRUE)
IN_DIR        <- "09_bcell_selection_v3"                 # step-9 bcell rds live here

# v2
# added COARSE and FINE cell types
# clusters (not sketched) as default ident,
# sample_id instead of sample

# v3
# uses step 9 v2 to include potential HRS cells
# increased resolution

# v4
# uses selection v3 (marker based cell type)
# saves celltype probabilities for sketched cells
OUT_DIR       <- "output/03_filter_and_process"  # results go here
save_fig = function(filename, p, w = 8, h = 8, dpi = 150){
  fig_f = file.path(OUT_DIR, "figures", filename)
  dir.create(dirname(fig_f), showWarnings = FALSE, recursive = TRUE)
}
ASSAY         <- "RNA"
SAMPLE_COL    <- "Run_Tissue_name"
REF_ANNO_COL_COARSE  <- "annotation_level_1"      # tonsil-atlas label column (see below)
REF_ANNO_COL_FINE  <- "annotation_20230508"      # tonsil-atlas label column (see below)

N_SKETCH      <- 50000     # cells to sketch per sample (tune to available RAM)
N_SKETCH_MERGE<- 100000     # cells to sketch for the merged joint analysis
N_DIMS        <- 30        # PCs used everywhere
RES           <- 1       # clustering resolution (tune + inspect)
HIGH_RES      <- 2
REF_N         <- 50000     # tonsil-atlas cells to subsample for transfer
set.seed(1)

dir.create(OUT_DIR, showWarnings = FALSE)
BP_STORE <- file.path(OUT_DIR, "bpcells")   # stable on-disk home for BPCells mats
dir.create(BP_STORE, showWarnings = FALSE)

GENES = readRDS("output/00_tonsil_reference/cosmx_genes.Rds")

# ---------------------------------------------------------------- HELPERS -----


# Standardize the full-resolution UMAP reduction name produced by ProjectData
# (its exact name varies across Seurat point releases). Renames to `umap.full`.
.standardize_umap_full <- function(obj) {
  DefaultAssay(obj) <- ASSAY
  n_full <- ncol(obj)
  for (rn in Reductions(obj)) {
    emb <- Embeddings(obj, rn)
    if (grepl("umap", rn, ignore.case = TRUE) && nrow(emb) == n_full && rn != "umap.full") {
      obj[["umap.full"]] <- CreateDimReducObject(
        embeddings = emb, key = "umapfull_", assay = ASSAY)
      obj[[rn]] <- NULL
      break
    }
  }
  obj
}

# Minimal-size save: drop scale.data, drop the sketch assay + its reductions,
# drop graphs; keep RNA counts/data (BPCells) + full PCA/UMAP + metadata.
save_min <- function(obj, path, keep_reds) {
  DefaultAssay(obj) <- ASSAY
  keep_reds <- intersect(keep_reds, Reductions(obj))
  obj <- DietSeurat(obj, assays = ASSAY, dimreducs = keep_reds, scale.data = FALSE)
  saveRDS(obj, path)
  message("wrote ", path, "  (", format(file.size(path), big.mark = ","), " bytes)")
  invisible(obj)
}

# --------------------------------------------------- 1. TONSIL REFERENCE ------
# Load + process ONCE (label source for cell-type assignment). Cached.
REF_RDS <- "output/00_tonsil_reference/tonsil_ref_processed.Rds"
# message("load cached tonsil reference")
# ref <- readRDS(REF_RDS)
# stopifnot(all(c(REF_ANNO_COL_COARSE, REF_ANNO_COL_FINE) %in% colnames(ref[[]])))

qc_meta_df = data.table::fread("output/02_filter_QC/qc_spatial_cells.csv.gz")
rownames(qc_meta_df) = qc_meta_df$cell
qc_meta_df.by_sample = split(qc_meta_df, qc_meta_df$sample_id)
rownames(qc_meta_df) = qc_meta_df$cell_id_prefixed
# ------------------------------------------- 2. PER-SAMPLE SKETCH PIPELINE ----

filter_sample = function(seu, sample_id){
  #sync new metadata with seu object
  common = intersect(colnames(seu), qc_meta_df.by_sample[[sample_id]]$cell)
  stopifnot(length(common) == ncol(seu))
  new_meta = qc_meta_df.by_sample[[sample_id]]
  new_meta = as.data.frame(new_meta)
  rownames(new_meta) = new_meta$cell
  new_meta = new_meta[colnames(seu),]
  drop_col = intersect(colnames(seu@meta.data), colnames(new_meta))
  new_meta = new_meta[, !colnames(new_meta) %in% drop_col]
  #add
  seu = AddMetaData(seu, new_meta)
  table(seu@meta.data$qc_passing_all, seu@meta.data$qcCellsPassed)
  seu = UpdateSeuratObject(seu)
  seu = subset(seu, qc_passing_all & qcCellsPassed)
  seu
}

transfer_reference <- function(seu, sample_id) {

  # (c) cell-type = tonsil-atlas transfer on the sketch cells (coarse + fine).
  #     TransferData takes a named list -> both levels in one anchor set.
  anchors <- FindTransferAnchors(reference = ref, query = seu,
                                 normalization.method = "LogNormalize",
                                 reference.reduction = "pca", dims = 1:N_DIMS)
  pred <- TransferData(
    anchorset = anchors, dims = 1:N_DIMS,
    refdata = list(celltype_coarse = ref[[REF_ANNO_COL_COARSE]][, 1],
                   celltype_fine   = ref[[REF_ANNO_COL_FINE]][, 1]))
  seu$sketch_celltype_coarse <- pred$celltype_coarse$predicted.id
  seu$sketch_celltype_fine   <- pred$celltype_fine$predicted.id

  # Keep the FULL soft classification (predicted.id + per-class prediction scores)
  # for the sketch cells so it can be written to a sidecar file. Only predicted.id
  # goes into the object; the score matrix would bloat it. TransferData data.frames
  # are keyed by sketch-cell barcode (rownames), so the sidecar joins back by cell.
  seu@misc$celltype_HCATonsil <- list(sample_id = sample_id,
                                      coarse    = pred$celltype_coarse,
                                      fine      = pred$celltype_fine)
  seu
}

# ------------------------------------------------- 3. RUN PER-SAMPLE ----------
RDS_PATH = RDS_PATHS[1]
for (RDS_PATH in RDS_PATHS) {
  sample_id <- strsplit(basename(RDS_PATH), "[_\\.]")[[1]][2]
  message("begin filtering for ", sample_id)
  FILT_RDS = file.path(OUT_DIR, paste0(sample_id, "_01_filtered.Rds"))
  if (file.exists(FILT_RDS)) { message("cached processed ", sample_id); next }

  seu = readRDS(RDS_PATH)
  message("=== processing ", sample_id, " ===")
  seu = filter_sample(seu, sample_id)
  saveRDS(seu, FILT_RDS)
}

for (RDS_PATH in RDS_PATHS) {
  sample_id <- strsplit(basename(RDS_PATH), "[_\\.]")[[1]][2]
  message("begin processing for ", sample_id)
  FILT_RDS = file.path(OUT_DIR, paste0(sample_id, "_01_filtered.Rds"))
  PROC_RDS  <- file.path(OUT_DIR, paste0(sample_id, "_02_processed.Rds"))
  if(file.exists(PROC_RDS)){
    next
  }

  if (file.exists(PROC_RDS)) { message("cached processed ", sample_id); next }
  rm(seu)
  seu = readRDS(FILT_RDS)

  dim(seu)

  tc <- Matrix::colSums(seu[["RNA"]]@counts) ## total counts per cell (across all genes)
  genefreq <- scPearsonPCA::gene_frequency(seu[["RNA"]]@counts) ## gene frequency (across all cells)
  sum(genefreq)==1 # TRUE


  seu <- Seurat::FindVariableFeatures(seu, nfeatures = 900)
  hvgs <- seu@assays$RNA@var.features

  ### Returns a Seurat-style DimReduc object with
  ### hvgs x pcs feature loadings
  ### cells x pcs cell embeddings
  ### gene-length vector of the mean pearson of residuals
  ### gene-length vector of the standard deviation of pearson residuals
  pcaobj <-
    sparse_quasipoisson_pca_seurat(seu[["RNA"]]@counts[hvgs,]
                                   ,totalcounts = tc
                                   ,grate = genefreq[hvgs]
                                   ,scale.max = 10 ## PCs reflect clipping pearson residuals > 10 SDs above the mean pearson residual
                                   ,do.scale = TRUE ## PCs reflect as if pearson residuals for each gene were scaled to have standard deviation=1
                                   ,do.center = TRUE ## PCs reflect as if pearson residuals for each gene were centered to have mean=0
    )

  umapobj <- scPearsonPCA::make_umap(pcaobj)
  seu[["pearsonpca"]] <- pcaobj$reduction.data
  seu[["pearsonumap"]] <- umapobj$ump  ## umap
  seu[["pearsongraph"]] <- Seurat::as.Graph(umapobj$grph) ## nearest neighbors / adjacency matrix used for unsupervised clustering
  seu <- Seurat::FindClusters(seu, graph = "pearsongraph", cluster.name = "pearson_clusters_high", resolution = .8)
  seu <- Seurat::FindClusters(seu, graph = "pearsongraph", cluster.name = "pearson_clusters_low", resolution = .4)

  saveRDS(seu, PROC_RDS)

  umapplot_high <- scPearsonPCA::plot_umap(umapreduc = "pearsonumap", clustercol = "pearson_clusters_high", semuse = seu)
  umapplot_low <- scPearsonPCA::plot_umap(umapreduc = "pearsonumap", clustercol = "pearson_clusters_low", semuse = seu)
  save_fig(paste0(sample_id, "_umap_high.png"), umapplot_high)
  save_fig(paste0(sample_id, "_umap_low.png"), umapplot_low)
  print(umapplot)



  # seu = process_sample(seu, sample_id)

}

for (RDS_PATH in RDS_PATHS) {
  sample_id <- strsplit(basename(RDS_PATH), "[_\\.]")[[1]][2]
  message("begin reference transfer for ", sample_id)
  FILT_RDS = file.path(OUT_DIR, paste0(sample_id, "_01_filtered.Rds"))
  PROC_RDS  <- file.path(OUT_DIR, paste0(sample_id, "_02_processed.Rds"))
  PRED_RDS  <- file.path(OUT_DIR, paste0(sample_id, "_03_HCA_cell_typed.Rds"))
  DETAIL_RDS  <- file.path(OUT_DIR, paste0(sample_id, "_03_HCA_result_tables.Rds"))
  if(file.exists(PRED_RDS)){
    message("skip ", PRED_RDS)
    nexts
  }
  seu = transfer_reference(seu, sample_id)

  # Sidecar: full per-class transfer probabilities for the sketch cells, written
  # separately so the saved object stays tiny; then dropped from misc so it does
  # NOT ride along inside the object.
  saveRDS(seu@misc$celltype_pred, DETAIL_RDS)
  message("wrote ", DETAIL_RDS)
  seu@misc$celltype_pred <- NULL
  saveRDS(seu, PRED_RDS)
  message("wrote ", PRED_RDS)

}
