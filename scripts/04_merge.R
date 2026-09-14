
# reload the per-sample processed objects (tiny) for the merge
proc_paths <- Sys.glob(file.path(OUT_DIR, "*_processed.Rds"))
proc_paths <- proc_paths[basename(proc_paths) != "merged_processed.Rds"]
proc_paths <- proc_paths[basename(proc_paths) != "tonsil_ref_processed.Rds"]
objs <- lapply(proc_paths, readRDS)
names(objs) <- sub("_processed.Rds$", "", basename(proc_paths))

# ---------------------------------------- 4. STASH PER-SAMPLE LABELS ----------
# Carry each sample's clusters/celltype/UMAP into metadata so they SURVIVE the
# merge (a merged object holds one joint reduction; per-sample UMAP goes to cols).
for (sid in names(objs)) {
  o <- objs[[sid]]
  o$sample_cluster         <- paste0(sid, "_", o$cluster)
  um <- Embeddings(o, "umap.full")
  o$sample_umap_1 <- um[, 1]
  o$sample_umap_2 <- um[, 2]
  objs[[sid]] <- o
}

# --------------------------------------------- 5. MERGE + JOINT RECALC --------
MERGED_RDS <- file.path(OUT_DIR, "merged_processed.Rds")
message("=== merging ", paste(names(objs), collapse = " + "), " ===")

merged <- merge(objs[[1]], y = objs[-1], add.cell.ids = names(objs))
DefaultAssay(merged) <- ASSAY
merged <- JoinLayers(merged)          # single counts/data layer (no integration)

# joint pipeline, sketched again on the combined set
merged <- NormalizeData(merged)
merged <- FindVariableFeatures(merged)
merged <- SketchData(merged, ncells = N_SKETCH_MERGE, method = "LeverageScore",
                     sketched.assay = "sketch")

DefaultAssay(merged) <- "sketch"
merged <- NormalizeData(merged)
merged <- FindVariableFeatures(merged)
merged <- ScaleData(merged)
merged <- RunPCA(merged, npcs = N_DIMS)
merged <- FindNeighbors(merged, dims = 1:N_DIMS)
merged <- FindClusters(merged, resolution = RES, cluster.name = "joint_sketch_clusters")
merged <- FindClusters(merged, resolution = HIGH_RES, cluster.name = "joint_sketch_clusters_high")
merged <- RunUMAP(merged, dims = 1:N_DIMS, return.model = TRUE,
                  reduction.name = "umap.sketch")

merged <- ProjectData(merged, assay = ASSAY,
                      full.reduction    = "pca.full",
                      sketched.assay    = "sketch",
                      sketched.reduction= "pca",
                      umap.model        = "umap.sketch",
                      dims              = 1:N_DIMS,
                      refdata = list(joint_cluster = "joint_sketch_clusters",
                                     joint_cluster_high = "joint_sketch_clusters_high"))

DefaultAssay(merged) <- ASSAY
merged <- .standardize_umap_full(merged)
# The HCA tonsil reference transfer was DISABLED in step 3 (2026-08-06): its
# labels anti-track marker evidence (78% of B-marker-positive cells were called
# `epithelial`; a 0.11% reference class absorbed half the query). See the block
# comment in scripts/03_script_filter_and_process.R. There are therefore no
# celltype_coarse / celltype_fine calls to carry forward -- use the marker
# `lineage` / `strongB` / `inSitu_clusters` columns for cell-type annotation.

# sanity: per-sample structure preserved (metadata) vs joint recompute (reduction)
message("per-sample clusters preserved in $sample_cluster; joint in $joint_cluster")
print(table(merged$sample_id, merged$joint_cluster))

Idents(merged) = merged@meta.data$joint_cluster
save_min(merged, MERGED_RDS, keep_reds = c("pca.full", "umap.full"))
message("Done. Outputs in: ", normalizePath(OUT_DIR))

# ------------------------------------------------------------- DECISION POINTS
# 1. Cell types transferred at two levels: REF_ANNO_COL_COARSE (annotation_level_1,
#    robust for a 6K panel) and REF_ANNO_COL_FINE (annotation_20230508, treat the
#    fine call as provisional — inspect the prediction scores before trusting it).
# 2. N_SKETCH / RES / N_DIMS are tunables — sweep and inspect UMAP/clusters.
# 3. Malignant (FL/HRS) calls are NOT made here — do those downstream on the
#    projected clusters (see step-9 GC/BCL2/light-chain diagnostics).
# 4. To move the .Rds off this filesystem, first re-materialize counts to memory
#    ( seu[["RNA"]]$counts <- as(seu[["RNA"]]$counts,"dgCMatrix") ) or write a
#    fresh BPCells store next to the destination — the saved object points at
#    the BPCells dirs under IN_DIR/OUT_DIR.

# # figure out cell types
# TonsilData_cellinfo_html("Naive")
# TonsilData_cellinfo_html("GC-Tfh-SAP" %>% gsub("-", " ", .))
# TonsilData_cellinfo_html("Tfh-LZ-GC Early")
# TonsilData_cellinfo_html("GC-commited NBC")
# ?TonsilData_cellinfo_html
#
# merged@meta.data %>% colnames %>% tail
# DimPlot(merged, group.by = "joint_sketch_clusters")
# DimPlot(merged, group.by = "joint_cluster")
#
# merged@meta.data$sketch_celltype %>% table %>% sort %>% rev
