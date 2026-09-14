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
  library(data.table)
  setDTthreads(1)
  library(scPearsonPCA)
  library(ggplot2)
  library(ggrepel)
})
options(future.globals.maxSize = 8 * 1024^3)

# ------------------------------------------------------------------ CONFIG ----
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(2L, 3L)) {
  stop("Usage: Rscript scripts/03_script_filter_and_process.R <input_seurat_rds> <sample_id> (<output_prefix>)", call. = FALSE)
}
RDS_PATH <- normalizePath(args[1], winslash = "/", mustWork = TRUE)
sample_id <- args[2]
if(length(args) == 3){
  output_prefix = args[3]
}else{
  output_prefix = sample_id
}


OUT_DIR <- "output/03_filter_and_process"
ASSAY <- "RNA"
SAMPLE_COL <- "Run_Tissue_name"
REF_ANNO_COL_COARSE <- "annotation_level_1"
REF_ANNO_COL_FINE <- "annotation_20230508"
N_DIMS <- 30
set.seed(1)

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
BP_STORE <- file.path(OUT_DIR, "bpcells")
dir.create(BP_STORE, showWarnings = FALSE, recursive = TRUE)

GENES <- readRDS("output/00_tonsil_reference/cosmx_genes.Rds")

# =============================================================================
# HCA TONSIL REFERENCE TRANSFER -- DISABLED 2026-08-06
# =============================================================================
# The transfer runs without error but its labels are not usable, so the whole
# step (reference load, loadings graft, transfer_reference(), and the step-03
# execution block) is commented out below rather than deleted.
#
# WHY IT WAS DISABLED
#   * The calls are near-CONSTANT across clusters, so they cannot annotate a
#     grouping: 22 of 26 joint clusters came back "epithelial" at ~60%/~35%
#     NBC_MBC while real B content (strongB) ranged 9.9%-46% across those same
#     clusters. Two equally B-enriched InSituType clusters got opposite calls.
#   * The composition is impossible: `epithelial` is 0.11% of the reference
#     (493 of 462,352 cells) but took 37.6% of IB1 and 53.5% of RB1; `GCBC` is
#     19.95% of the reference and received 2 cells total.
#   * It is inverted at the cell level: 78% of RB1's strongB (B-marker-positive)
#     cells were called `epithelial`, and cells called `epithelial` detect MS4A1
#     (57%) and CD79A (41%) MORE often than cells called `NBC_MBC` (36%/20%),
#     while being keratin- and PanCK-negative. Not a dropout artifact -- dropout
#     cannot make marker-POSITIVE cells land elsewhere.
#
# ROOT CAUSE (diagnosed, not guessed)
#   The grafted loadings are CORRECT: reference cells pushed back through the
#   same path reproduce their own stored embeddings at rho = 1.000/0.999/0.997
#   with matching spread. The failure is cross-platform projection. CosMx query
#   cells land in a 5-6x COMPRESSED, DISPLACED cloud inside the scRNA reference
#   PCA space (per-PC sd 1.3-1.5 vs the reference's 7.6-8.2; PC4 mean -9.1 vs
#   +1.3, PC5 +8.1 vs -0.4). Every query cell is crammed into a small off-centre
#   volume whose nearest reference neighbours are whatever few cells sit there --
#   hence a 0.11% class absorbing half the query, and no anchors reaching GCBC.
#
# ALTERNATIVES TESTED (all still project onto a pre-built scRNA PCA, all fail):
#   S1 grafted pca.full : 38.5% epithelial, cluster-level rho +0.387
#   S2 sketch-consistent: 35.1% epithelial, rho +0.297
#   S3 de-novo PCA      : identical to S2
#   In every one, non-strongB cells were called B-lineage MORE often than
#   strongB cells (33.6% vs 15.7%) -- the labels anti-track marker evidence.
#   Joint-space anchors (CCA / RPCA) were never successfully tested.
#
# WHAT REPLACES IT
#   * cell-type annotation: marker `lineage` / `strongB` / `inSitu_clusters`,
#     which DO discriminate across clusters (see output/04_merge/figures/
#     fig_04_cluster_lineage_heatmap_*.png -- 18 of 26 clusters resolve as
#     B-lineage, with clean myeloid/fibro/endo/T clusters).
#   * B cells in IB1: step 5 CCA integration co-clusters 33,209 IB1 cells with
#     RB1's B cells at 7.25x strongB enrichment and BELOW-baseline cluster-i
#     contamination (15.5% vs IB1's 37%).
#
# CONSEQUENCE IF YOU RE-RUN STEP 3 FROM SCRATCH
#   Fresh *_03_HCA_cell_typed.Rds objects will NOT carry celltype_coarse /
#   celltype_fine. The objects already on disk still do. Two consumers assume
#   those columns exist and would need adjusting:
#     - scripts/04_script_merge_and_process.R  (stopifnot on both columns)
#     - scripts/04_merge.R                     ($sample_celltype_coarse/_fine)
#   To re-enable, uncomment this block and the two blocks marked
#   "DISABLED REFERENCE TRANSFER" below; the code is unchanged and still works
#   mechanically -- it is the biology of the output that is wrong.
# =============================================================================
# REF_RDS <- "output/00_tonsil_reference/tonsil_ref_processed.Rds"
# REF_RDS_EXISTS <- file.exists(REF_RDS)
# if (!REF_RDS_EXISTS) {
#   stop("reference RDS not found: ", REF_RDS, call. = FALSE)
# }
# ref <- readRDS(REF_RDS)
# stopifnot(all(c(REF_ANNO_COL_COARSE, REF_ANNO_COL_FINE) %in% colnames(ref[[]])))
#
# # ProjectData() (00_tonsil_reference.R) stores only cell embeddings on
# # `pca.full` -- no feature loadings. FindTransferAnchors() derives the
# # projection genes from those loadings, so with 0x0 loadings the query collapses
# # to zero features and ProjectCellEmbeddings() dies inside as.sparse() with
# # "no applicable method for 'as.sparse' applied to an object of class
# # 'c(double, numeric)'". `pca.full` IS the full-cell projection onto the sketch
# # PCA, so the sketch loadings are its loadings.
# if (nrow(Loadings(ref[["pca.full"]])) == 0) {
#   ref_loadings <- Loadings(ref[["pca"]])
#   colnames(ref_loadings) <- colnames(Embeddings(ref[["pca.full"]]))  # key must match
#   Loadings(ref[["pca.full"]], projected = FALSE) <- ref_loadings
#   message("grafted sketch PCA loadings onto reference pca.full (",
#           paste(dim(Loadings(ref[["pca.full"]])), collapse = " x "), ")")
# }
# # The reference's default assay is `sketch`, which holds only 100k of its
# # 462,352 cells; `pca.full` spans all of them, so anchor finding must use RNA.
# DefaultAssay(ref) <- ASSAY
# REF_FEATURES <- rownames(Loadings(ref[["pca.full"]]))
# stopifnot(length(REF_FEATURES) > 0)

qc_meta_df <- data.table::fread("output/02_filter_QC/qc_spatial_cells.csv.gz")
qc_meta_df$cell_id_prefixed <- qc_meta_df$cell
qc_meta_df.by_sample <- split(qc_meta_df, qc_meta_df$sample_id)

# ---------------------------------------------------------------- HELPERS -----

completion_file <- function(step_name) {
  file.path(OUT_DIR, paste0(output_prefix, "_", step_name, ".done"))
}

mark_complete <- function(step_name) {
  done_path <- completion_file(step_name)
  writeLines(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), done_path)
  message("marked complete: ", done_path)
}

is_complete <- function(step_name) {
  file.exists(completion_file(step_name))
}

save_fig <- function(filename, p, w = 8, h = 8, dpi = 150) {
  fig_f <- file.path(OUT_DIR, "figures", filename)
  dir.create(dirname(fig_f), showWarnings = FALSE, recursive = TRUE)
  ggsave(fig_f, plot = p, width = w, height = h, dpi = dpi)
}

# Standardize the full-resolution UMAP reduction name produced by ProjectData
# (its exact name varies across Seurat point releases). Renames to `umap.full`.
.standardize_umap_full <- function(obj) {
  DefaultAssay(obj) <- ASSAY
  n_full <- ncol(obj)
  for (rn in Reductions(obj)) {
    emb <- Embeddings(obj, rn)
    if (grepl("umap", rn, ignore.case = TRUE) && nrow(emb) == n_full && rn != "umap.full") {
      obj[["umap.full"]] <- CreateDimReducObject(
        embeddings = emb, key = "umapfull_", assay = ASSAY
      )
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

# ------------------------------------------- 2. PER-SAMPLE SKETCH PIPELINE ----

filter_sample <- function(seu, sample_id) {
  common <- intersect(colnames(seu), qc_meta_df.by_sample[[sample_id]]$cell)
  stopifnot(length(common) == ncol(seu))

  new_meta <- qc_meta_df.by_sample[[sample_id]]
  new_meta <- as.data.frame(new_meta)
  rownames(new_meta) <- new_meta$cell
  new_meta <- new_meta[colnames(seu), ]

  drop_col <- intersect(colnames(seu@meta.data), colnames(new_meta))
  new_meta <- new_meta[, !colnames(new_meta) %in% drop_col]

  seu <- AddMetaData(seu, new_meta)
  table(seu@meta.data$qc_passing_all, seu@meta.data$qcCellsPassed)
  seu <- UpdateSeuratObject(seu)
  seu <- subset(seu, qc_passing_all & qcCellsPassed)
  seu
}

process_sample <- function(seu, sample_id) {
  PROC_RDS <- file.path(OUT_DIR, paste0(sample_id, "_02_processed.Rds"))


  if(!is.null(seu[["RNA"]]$counts)){
    count_assay = seu[["RNA"]]$counts
  }else{
    count_assay = seu[["RNA"]]@counts
  }

  tc <- Matrix::colSums(count_assay)


  genefreq <- scPearsonPCA::gene_frequency(count_assay)
  stopifnot(near(sum(genefreq), 1))

  seu <- Seurat::FindVariableFeatures(seu, nfeatures = 900)
  hvgs <- VariableFeatures(seu)

  counts_hvgs = count_assay[hvgs, ]
  counts_hvgs <- as(counts_hvgs, "dgCMatrix")
  pcaobj <- sparse_quasipoisson_pca_seurat(
    counts_hvgs,
    totalcounts = tc,
    grate = genefreq[hvgs],
    scale.max = 10,
    do.scale = TRUE,
    do.center = TRUE
  )

  umapobj <- scPearsonPCA::make_umap(pcaobj)
  seu[["pearsonpca"]] <- pcaobj$reduction.data
  seu[["pearsonumap"]] <- umapobj$ump
  seu[["pearsongraph"]] <- Seurat::as.Graph(umapobj$grph)
  seu <- Seurat::FindClusters(
    seu, graph = "pearsongraph", cluster.name = "pearson_clusters_high", resolution = .8
  )
  seu <- Seurat::FindClusters(
    seu, graph = "pearsongraph", cluster.name = "pearson_clusters_low", resolution = .4
  )

  umapplot_high <- scPearsonPCA::plot_umap(
    umapreduc = "pearsonumap", clustercol = "pearson_clusters_high", semuse = seu
  )
  umapplot_low <- scPearsonPCA::plot_umap(
    umapreduc = "pearsonumap", clustercol = "pearson_clusters_low", semuse = seu
  )

  save_fig(filename = paste0(sample_id, "_umap_high.png"), p = umapplot_high)
  save_fig(filename = paste0(sample_id, "_umap_low.png"), p = umapplot_low)

  saveRDS(seu, PROC_RDS)
  mark_complete("02_processed")
  invisible(seu)
}

# ----------------------------------------------------------------------------
# DISABLED REFERENCE TRANSFER (see the note in the CONFIG section above).
# Kept verbatim for reference; the code runs, but its labels anti-track marker
# evidence because CosMx cells project into a compressed, displaced region of
# the scRNA reference PCA space. Uncomment together with the CONFIG block.
# ----------------------------------------------------------------------------
# transfer_reference <- function(seu, sample_id) {
#   N_SKETCH = 1e5
#
#   # SketchData() samples per *data* layer -- with a counts-only assay it finds no
#   # layers ("No layers found matching search pattern provided"), samples nothing,
#   # and the "sketch" assay silently ends up holding EVERY cell. Normalize first;
#   # on a BPCells-backed assay this stays lazy/on-disk.
#   DefaultAssay(seu) <- ASSAY
#   if (!"data" %in% Layers(seu[[ASSAY]])) {
#     seu <- NormalizeData(seu)
#   }
#   seu <- SketchData(seu, ncells = N_SKETCH, method = "LeverageScore",
#                     sketched.assay = "sketch")
#   n_sketch_cells <- length(Cells(seu[["sketch"]]))
#   message("sketched ", n_sketch_cells, " of ", ncol(seu), " cells")
#   # guard against the silent no-op above (also fires if ncells > ncol)
#   stopifnot(n_sketch_cells == min(N_SKETCH, ncol(seu)))
#
#   DefaultAssay(seu) <- "sketch"
#   seu <- NormalizeData(seu)
#   seu <- FindVariableFeatures(seu)
#   seu <- ScaleData(seu)
#   # explicit name: "pca" would collide with the pearsonpca key (PC_)
#   seu <- RunPCA(seu, npcs = N_DIMS, reduction.name = "pca.sketch")
#
#   # anchors are found on the sketch cells only (query.assay = "sketch")
#   anchors <- FindTransferAnchors(
#     reference = ref,
#     query = seu,
#     normalization.method = "LogNormalize",
#     reference.reduction = "pca.full",
#     features = REF_FEATURES,
#     dims = 1:N_DIMS
#   )
#
#   ref_md = ref@meta.data
#   rownames(ref_md) = paste0(colnames(ref), "_reference")
#   pred <- TransferData(
#     anchorset = anchors,
#     dims = 1:N_DIMS,
#     refdata = list(
#       celltype_coarse = ref_md[anchors@reference.cells, REF_ANNO_COL_COARSE],
#       celltype_fine = ref_md[anchors@reference.cells, REF_ANNO_COL_FINE]
#     )
#   )
#
#   # sketch-cell labels (NA elsewhere); pred rownames are the sketch cell names
#   stopifnot(setequal(rownames(pred$celltype_coarse), Cells(seu[["sketch"]])))
#   seu$sketch_celltype_coarse <- NA_character_
#   seu$sketch_celltype_fine   <- NA_character_
#   seu@meta.data[rownames(pred$celltype_coarse), "sketch_celltype_coarse"] <-
#     pred$celltype_coarse$predicted.id
#   seu@meta.data[rownames(pred$celltype_fine), "sketch_celltype_fine"] <-
#     pred$celltype_fine$predicted.id
#
#   # propagate sketch labels + PCA to every full-resolution cell. Produces
#   # celltype_coarse / celltype_fine (+ .score) on all cells, which 04_merge.R
#   # consumes, plus the `pca.full` reduction.
#   seu <- ProjectData(
#     seu,
#     assay = ASSAY,
#     sketched.assay = "sketch",
#     sketched.reduction = "pca.sketch",
#     full.reduction = "pca.full",
#     dims = 1:N_DIMS,
#     refdata = list(
#       celltype_coarse = "sketch_celltype_coarse",
#       celltype_fine = "sketch_celltype_fine"
#     )
#   )
#   stopifnot(!any(is.na(seu$celltype_coarse)), !any(is.na(seu$celltype_fine)))
#
#   # TransferData() inside ProjectData() rewrites '_' to '-' in the transferred
#   # labels ("Names of identity class contain underscores"), so the projected
#   # columns silently diverge from the reference vocabulary and from the
#   # sketch_celltype_* columns beside them (NBC_MBC -> NBC-MBC, CD4_T -> CD4-T,
#   # DZ_LZ transition -> DZ-LZ transition). Map them back by exact match against
#   # the reference's own label set, skipping any dashed form that is itself a
#   # real reference label so a legitimate dash is never rewritten.
#   restore_ref_labels <- function(x, ref_levels) {
#     ref_levels <- unique(as.character(ref_levels))
#     key <- setNames(ref_levels, gsub("_", "-", ref_levels, fixed = TRUE))
#     key <- key[!names(key) %in% ref_levels]          # no collisions with real labels
#     key <- key[!duplicated(names(key))]
#     x <- as.character(x)
#     hit <- x %in% names(key)
#     x[hit] <- unname(key[x[hit]])
#     x
#   }
#   seu$celltype_coarse <- restore_ref_labels(seu$celltype_coarse,
#                                             ref_md[[REF_ANNO_COL_COARSE]])
#   seu$celltype_fine   <- restore_ref_labels(seu$celltype_fine,
#                                             ref_md[[REF_ANNO_COL_FINE]])
#   stopifnot(
#     all(seu$celltype_coarse %in% as.character(ref_md[[REF_ANNO_COL_COARSE]])),
#     all(seu$celltype_fine   %in% as.character(ref_md[[REF_ANNO_COL_FINE]]))
#   )
#
#   DefaultAssay(seu) <- ASSAY
#   # dense scale.data over the sketch would dominate the saved .Rds
#   LayerData(seu[["sketch"]], layer = "scale.data") <- NULL
#
#   seu@misc$celltype_HCATonsil <- list(
#     sample_id = sample_id,
#     coarse = pred$celltype_coarse,
#     fine = pred$celltype_fine
#   )
#   seu
# }

# ------------------------------------------------- 3. RUN SAMPLE ----------
message("input RDS: ", RDS_PATH)
message("sample_id: ", sample_id)

FILT_RDS <- file.path(OUT_DIR, paste0(output_prefix, "_01_filtered.Rds"))
PROC_RDS <- file.path(OUT_DIR, paste0(output_prefix, "_02_processed.Rds"))
# unused while the reference transfer is disabled:
# PRED_RDS <- file.path(OUT_DIR, paste0(output_prefix, "_03_HCA_cell_typed.Rds"))
# DETAIL_RDS <- file.path(OUT_DIR, paste0(output_prefix, "_03_HCA_result_tables.Rds"))

if (!is_complete("01_filtered")) {
  message("begin filtering for ", output_prefix)
  seu <- readRDS(RDS_PATH)
  seu <- filter_sample(seu, sample_id)
  saveRDS(seu, FILT_RDS)
  mark_complete("01_filtered")
} else {
  message("skip filtering for ", sample_id, " (already completed)")
}

if (!is_complete("02_processed")) {
  message("begin processing for ", output_prefix)
  seu <- readRDS(FILT_RDS)
  seu <- process_sample(seu, sample_id)
  saveRDS(seu, PROC_RDS)
  mark_complete("02_processed")
} else {
  message("skip processing for ", sample_id, " (already completed)")
}

# ----------------------------------------------------------------------------
# DISABLED REFERENCE TRANSFER step (see the note in the CONFIG section above).
# With this commented out, step 3 produces *_02_processed.Rds as its final
# product; no *_03_HCA_cell_typed.Rds and no celltype_coarse/celltype_fine.
# ----------------------------------------------------------------------------
# if (!is_complete("03_HCA_cell_typed")) {
#   message("begin reference transfer for ", output_prefix)
#   seu <- readRDS(PROC_RDS)
#   seu <- transfer_reference(seu, sample_id)
#
#   saveRDS(seu@misc$celltype_HCATonsil, DETAIL_RDS)
#   message("wrote ", DETAIL_RDS)
#   seu@misc$celltype_HCATonsil <- NULL
#   saveRDS(seu, PRED_RDS)
#   message("wrote ", PRED_RDS)
#   mark_complete("03_HCA_cell_typed")
# } else {
#   message("skip reference transfer for ", sample_id, " (already completed)")
# }
