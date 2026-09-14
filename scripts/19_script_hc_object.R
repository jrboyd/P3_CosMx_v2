## ============================================================================
## 19_script_hc_object.R -- build and cluster the HIGH-CONFIDENCE object.
##
## The clean substrate every cell-type STATEMENT in this project will rest on: both samples,
## merged, with the over-fixation defect removed BEFORE clustering, jointly embedded and
## clustered once.
##
## WHY IT DOES NOT START FROM output/04_merge/merged_processed.Rds. Step 03 filters with
## `subset(seu, qc_passing_all & qcCellsPassed)` (03_script_filter_and_process.R:227), so
## every processed object in this project -- and therefore step 04's merged object -- is the
## qc_passing_all tier. That tier is ASYMMETRIC: IB1 lost 2,013,225 cells to a
## staining/section rule RB1 was never evaluated for (`qc_spatial_staining ==
## "not_evaluated"`). The damage is not subtle here:
##
##     high-confidence set                     IB1        RB1      IB1 share
##     qc_passing_all & non-defect         251,747  1,286,223        16.4%
##     tier_matched   & non-defect         767,395  1,286,223        37.4%
##
## RB1 is IDENTICAL in both rows -- it was never evaluated, so the two tiers coincide for it.
## The whole difference is IB1 losing two thirds of its usable cells to a rule the other
## sample never faced. Clustering on the existing object would make the joint embedding 83.6%
## RB1, and IB1 is the case: the composite node carrying the Hodgkin component. Cell types
## that exist only in IB1's Hodgkin tissue could fail to form a cluster at all.
##
## So this re-filters from the step-01 imports on the SYMMETRIC tier
## (`qc_spatial_fovpx == "ok" & qcCellsPassed`, i.e. `tier_matched`) and merges itself.
## It is step 03's filter + step 04's merge + step 03's process block, with one filter
## changed. The code is deliberately copied rather than shared: those scripts are the
## record of what produced the existing objects and must not change under them.
##
## DEFECT IS DROPPED PER CELL, BEFORE CLUSTERING, and that is the point of the whole design.
## `09c_script_cluster_consistency.R` showed that pooling LAUNDERS defect tissue: merged
## pearson_clusters_low cluster 1 was 93.6% IB1 cells that IB1's own audit calls
## defect_dominated at 70.1% defect, but pooling diluted the defect fraction to 42% -- under
## the DEFECT_MAX = 0.50 guard -- so it escaped the guard and acquired the label "Neutrophil".
## With no defect cells in the object there is nothing left to dilute. This is why
## CLAUDE.md's "do not use the merged labels" rule does not apply to this object, and the
## rule must be rewritten (not deleted) when this lands.
##
## Never via the FOV-level `qc_spatial_i_fraction` flag, which is far too aggressive and
## cuts IB1 to 196,316.
##
## A PER-CELL LOW-INFORMATION FLOOR IS APPLIED TOO, and it is not cosmetic. Dropping the
## defect leaves behind cells that are simply empty: the first build put 159 cells at
## UMAP1 ~ -18.3 while 99.99% of the object lay between -4.42 and +4.19, stretching the x
## axis to 4.5x the width the data needed and squashing two million cells into a fifth of
## every panel. Those cells carry a median of 24 counts over 17 genes against an object
## median of 523 over 357 -- they are not a cell type, they are noise that the embedding
## has nowhere to put. Filtering them at PLOT time would leave them in the object and in
## the clustering, where they form their own junk cluster and distort the neighbour graph,
## so the floor belongs here.
##
## The criteria are per-cell and defined on the cell's own evidence, never on its cluster,
## which is the rule the step-19b gate lays down for exactly this remedy.
MIN_COUNT   <- 50L   # removes 2.90% of the tier
MIN_FEATURE <- 30L   # removes 1.84%; the 159 outliers max out at 34 counts / 28 genes
##
## CLUSTERS ARE NAMED hc_clusters_low/high, never pearson_clusters_*. Those ids belong to
## step 03 and mean something else; a name collision here would silently join step-09b
## labels onto clusters that have nothing to do with them.
##
## TWO STAGES with .done markers -- the merge is expensive and must not be redone when the
## PCA fails.
##
## Usage: sbatch scripts/sbatch/19_hc_object_wrapper.sh
## ============================================================================
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(data.table); library(scPearsonPCA)
  setDTthreads(1)
})
if (!file.exists("output/16_regions/cells_meta_all.csv.gz")) stop("run from the project root", call.=FALSE)

IMPORT   <- "output/01_import_v1"
OUT      <- "output/19_celltypes"
META_F   <- "output/16_regions/cells_meta_all.csv.gz"
ASSAY    <- "RNA"
SAMPLES  <- c("IB1", "RB1")
RAW_RDS  <- file.path(OUT, "hc_01_merged.Rds")
CLU_RDS  <- file.path(OUT, "hc_02_clustered.Rds")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
done <- function(k) file.path(OUT, paste0(k, ".done"))
is_done <- function(k) file.exists(done(k))
mark <- function(k) writeLines(format(Sys.time()), done(k))

cat("============================================================\n")
cat("step 19a -- high-confidence object: merge, embed, cluster\n")
cat("============================================================\n")

## ---- the high-confidence cell list ------------------------------------------
md <- fread(META_F, select = c("cell_id","samp","tier_matched","is_defect","inSitu_clusters",
                               "nCount_RNA"))
## nFeature_RNA is not in cells_meta_all; it comes from the vendor metadata via
## scripts/18_script_centroids.sh.
nf <- rbindlist(lapply(SAMPLES, function(S) {
  f <- sprintf("output/18_leakage/%s_centroids.csv.gz", S)
  if (!file.exists(f)) stop("missing ", f, " -- run scripts/18_script_centroids.sh ", S, call.=FALSE)
  fread(f, select = c("cell_id","nFeature_RNA")) }))
md <- nf[md, on = .(cell_id)]
hc <- md[tier_matched == TRUE & is_defect == FALSE]
cat(sprintf("\ntier_matched & non-defect: %d\n", nrow(hc)))
lowinfo <- hc[nCount_RNA < MIN_COUNT | nFeature_RNA < MIN_FEATURE | is.na(nFeature_RNA)]
cat(sprintf("dropped as low-information (nCount < %d or nFeature < %d): %d (%.2f%%)\n",
            MIN_COUNT, MIN_FEATURE, nrow(lowinfo), 100*nrow(lowinfo)/nrow(hc)))
print(lowinfo[, .(cells = .N, med_nCount = as.numeric(median(nCount_RNA)),
                  med_nFeature = as.numeric(median(nFeature_RNA))), by = samp])
hc <- hc[nCount_RNA >= MIN_COUNT & nFeature_RNA >= MIN_FEATURE]
cat("(they are not lost -- they fall to the rescue tier and still receive argmax and\n")
cat(" projected labels for the spatial maps)\n")
cat("\n=== high-confidence set ===\n")
print(hc[, .(cells = .N), by = samp][order(samp)])
cat(sprintf("total %d of %d (%.1f%%)\n", nrow(hc), nrow(md), 100*nrow(hc)/nrow(md)))
if (hc[, sum(inSitu_clusters == "i")] > 0) stop("defect cells survived the filter", call.=FALSE)
KEEP <- split(hc$cell_id, hc$samp)

## ---- stage 1: subset each import and merge ----------------------------------
if (!is_done("hc_01_merged")) {
  objs <- list()
  for (S in SAMPLES) {
    f <- file.path(IMPORT, sprintf("seuratObject_%s_BPCells_full.Rds", S))
    if (!file.exists(f)) stop("missing ", f, call. = FALSE)
    cat("\nreading ", f, "\n", sep = "")
    o <- readRDS(f)
    cat(sprintf("  %s: %d cells in the import\n", S, ncol(o)))
    common <- intersect(colnames(o), KEEP[[S]])
    cat(sprintf("  high-confidence cells found: %d of %d wanted (%.2f%%)\n",
                length(common), length(KEEP[[S]]), 100*length(common)/length(KEEP[[S]])))
    if (length(common) < 0.99 * length(KEEP[[S]]))
      stop("cell ids do not line up between the import and cells_meta_all for ", S,
           " -- check the id convention before going further", call. = FALSE)
    o <- subset(o, cells = common)
    o$merge_sample <- S
    ## counts only: the lazy `data` layer and every reduction belong to a different
    ## filtering and would be carried in silently. process below uses counts.
    o <- DietSeurat(o, assays = ASSAY, layers = "counts", dimreducs = NULL, graphs = NULL)
    cat(sprintf("  kept %d cells | counts class %s\n", ncol(o),
                class(LayerData(o[[ASSAY]], layer = "counts"))[1]))
    objs[[S]] <- o; invisible(gc())
  }
  ## Cell ids are c_<slide>_<fov>_<cell> and the slide token differs between samples, so they
  ## do not collide -- but assert rather than assume, because add.cell.ids would rename every
  ## cell and break every join downstream.
  all_cells <- unlist(lapply(objs, colnames), use.names = FALSE)
  if (anyDuplicated(all_cells) > 0) stop("cell ids collide across samples", call. = FALSE)
  cat("\nmerging\n")
  merged <- merge(objs[[1]], y = objs[-1]); rm(objs); invisible(gc())
  DefaultAssay(merged) <- ASSAY
  merged <- JoinLayers(merged)      # one counts layer, so seu[["RNA"]]$counts works
  stopifnot("counts" %in% Layers(merged[[ASSAY]]))
  cat(sprintf("merged: %d genes x %d cells | counts class %s\n", nrow(merged), ncol(merged),
              class(LayerData(merged[[ASSAY]], layer = "counts"))[1]))
  print(table(merged$merge_sample))
  saveRDS(merged, RAW_RDS); mark("hc_01_merged")
  cat("wrote ", RAW_RDS, "\n", sep = "")
} else {
  cat("\nskip merge (already complete)\n")
}

## ---- stage 2: the pearson pipeline, verbatim from step 03 -------------------
if (!is_done("hc_02_clustered")) {
  merged <- readRDS(RAW_RDS)
  cat(sprintf("\njoint processing of %d cells\n", ncol(merged)))
  ca <- if (!is.null(merged[[ASSAY]]$counts)) merged[[ASSAY]]$counts else merged[[ASSAY]]@counts
  tc <- Matrix::colSums(ca)
  gf <- scPearsonPCA::gene_frequency(ca)
  stopifnot(abs(sum(gf) - 1) < 1e-8)
  merged <- Seurat::FindVariableFeatures(merged, nfeatures = 900)
  hv <- VariableFeatures(merged)
  ch <- as(ca[hv, ], "dgCMatrix")
  cat(sprintf("  %d HVGs x %d cells, %.2f%% nonzero\n", nrow(ch), ncol(ch),
              100*length(ch@x)/(as.numeric(nrow(ch))*ncol(ch))))
  po <- sparse_quasipoisson_pca_seurat(ch, totalcounts = tc, grate = gf[hv],
                                       scale.max = 10, do.scale = TRUE, do.center = TRUE)
  rm(ch); invisible(gc())
  uo <- scPearsonPCA::make_umap(po)
  merged[["pearsonpca"]]   <- po$reduction.data
  merged[["pearsonumap"]]  <- uo$ump
  merged[["pearsongraph"]] <- Seurat::as.Graph(uo$grph)
  merged <- Seurat::FindClusters(merged, graph = "pearsongraph",
                                 cluster.name = "hc_clusters_high", resolution = .8)
  merged <- Seurat::FindClusters(merged, graph = "pearsongraph",
                                 cluster.name = "hc_clusters_low",  resolution = .4)
  merged$seurat_clusters <- NULL
  cat("\n=== hc_clusters_low ===\n");  print(table(merged$hc_clusters_low))
  cat("\n=== hc_clusters_high ===\n"); print(table(merged$hc_clusters_high))
  cat("\n=== hc_clusters_low x sample ===\n")
  print(table(merged$hc_clusters_low, merged$merge_sample))
  saveRDS(merged, CLU_RDS); mark("hc_02_clustered")
  cat("wrote ", CLU_RDS, "\n", sep = "")

  ## The per-cell table everything downstream reads, so no other script has to load the object.
  emb <- as.data.table(Embeddings(merged, "pearsonumap"), keep.rownames = "cell_id")
  setnames(emb, 2:3, c("U1","U2"))
  out <- data.table(cell_id = colnames(merged),
                    merge_sample = as.character(merged$merge_sample),
                    hc_clusters_low  = as.character(merged$hc_clusters_low),
                    hc_clusters_high = as.character(merged$hc_clusters_high))
  out <- emb[out, on = .(cell_id)]
  fwrite(out, file.path(OUT, "hc_cells.csv.gz"))
  cat("wrote ", file.path(OUT, "hc_cells.csv.gz"), " (", nrow(out), " cells)\n", sep = "")
} else {
  cat("\nskip processing (already complete)\n")
}
cat("\ndone\n")
