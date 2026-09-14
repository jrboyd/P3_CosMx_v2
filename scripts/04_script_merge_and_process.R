#
# Step 4: merge the per-sample step-03 cell-typed objects, then run the SAME
# scPearsonPCA workflow used per sample (process_sample() from step 3) on the
# merged object to get a joint PCA / UMAP / clustering.
#
# Per-sample labels (cell types, pearson clusters, per-sample UMAP coords) are
# stashed into metadata as sample_* columns so they survive the merge; the joint
# reduction is recomputed from scratch.
#   NOTE: the HCA tonsil reference transfer was disabled in step 3 (2026-08-06)
#   because its labels anti-track marker evidence -- see the block comment in
#   scripts/03_script_filter_and_process.R. celltype_coarse / celltype_fine are
#   therefore NOT carried forward. Cell-type annotation comes from the marker
#   `lineage` / `strongB` / `inSitu_clusters` columns, which do discriminate
#   across clusters.
#
# Usage: Rscript scripts/04_script_merge_and_process.R [<prefix> ...]
#        (default prefixes: IB1 RB1)
#
# NOTE: counts stay BPCells (disk-backed) through the merge; only the 900 HVGs
# are materialized for the quasipoisson PCA. Keep IN_DIR/OUT_DIR intact -- the
# saved .Rds stores pointers into the BPCells stores.

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
prefixes <- if (length(args)) args else c("IB1", "RB1")

IN_DIR <- "output/03_filter_and_process"
OUT_DIR <- "output/04_merge"          # process_sample()/save_fig() read this global
ASSAY <- "RNA"
N_DIMS <- 30
set.seed(1)

# Pre-flight hook: MERGE_TEST_NCELLS=<n> subsets each sample to n cells and
# writes everything under the `merged_test` prefix, so the full production code
# path (BPCells merge -> JoinLayers -> quasipoisson PCA -> UMAP -> clusters) can
# be validated in minutes before committing to the ~1.8M-cell run. Unset = real run.
TEST_NCELLS <- suppressWarnings(as.integer(Sys.getenv("MERGE_TEST_NCELLS", "")))
IS_TEST <- !is.na(TEST_NCELLS) && TEST_NCELLS > 0
output_prefix <- if (IS_TEST) "merged_test" else "merged"  # used by completion_file()
if (IS_TEST) message("PRE-FLIGHT: subsetting each sample to ", TEST_NCELLS, " cells")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Reuse step 3's helpers verbatim rather than copying them: evaluate ONLY the
# top-level `name <- function(...)` definitions from that script, so none of its
# side effects run (no reference load, no per-sample pipeline).
STEP3 <- "scripts/03_script_filter_and_process.R"
local({
  n <- 0L
  for (e in parse(STEP3)) {
    if (is.call(e) && identical(as.character(e[[1]]), "<-") &&
        is.call(e[[3]]) && identical(as.character(e[[3]][[1]]), "function")) {
      eval(e, envir = globalenv())
      n <- n + 1L
    }
  }
  message("loaded ", n, " helper functions from ", STEP3)
})
stopifnot(is.function(process_sample), is.function(save_fig), is.function(is_complete))

MERGED_RAW_RDS <- file.path(OUT_DIR, paste0(output_prefix, "_01_merged.Rds"))
MERGED_RDS <- file.path(OUT_DIR, paste0(output_prefix, "_processed.Rds"))

# ------------------------------------------------------------------- MERGE ----
if (!is_complete("01_merged")) {
  objs <- list()
  for (sid in prefixes) {
    f <- file.path(IN_DIR, paste0(sid, "_03_HCA_cell_typed.Rds"))
    stopifnot(file.exists(f))
    message("reading ", f)
    o <- readRDS(f)

    if (!"sample_id" %in% colnames(o@meta.data)) o$sample_id <- sid
    o$merge_sample <- sid

    if (IS_TEST) {
      o <- subset(o, cells = sample(colnames(o), min(TEST_NCELLS, ncol(o))))
      message("  ", sid, ": subset to ", ncol(o), " cells (pre-flight)")
    }

    # stash per-sample results that a joint recompute would otherwise overwrite
    for (cl in c("pearson_clusters_high", "pearson_clusters_low")) {
      if (cl %in% colnames(o@meta.data)) {
        o[[paste0("sample_", cl)]] <- paste0(sid, "_", o@meta.data[[cl]])
      }
    }
    if ("pearsonumap" %in% Reductions(o)) {
      um <- Embeddings(o, "pearsonumap")
      o$sample_pearsonumap_1 <- um[, 1]
      o$sample_pearsonumap_2 <- um[, 2]
    }

    # objects produced BEFORE the transfer was disabled still carry these; drop
    # them so nothing downstream can pick up labels we do not trust
    drop_cols <- intersect(c("celltype_coarse", "celltype_fine",
                             "celltype_coarse.score", "celltype_fine.score",
                             "sketch_celltype_coarse", "sketch_celltype_fine"),
                           colnames(o@meta.data))
    if (length(drop_cols)) {
      for (dc in drop_cols) o[[dc]] <- NULL
      message("  ", sid, ": dropped disabled-transfer columns: ",
              paste(drop_cols, collapse = ", "))
    }

    # keep only what the joint pipeline needs: RNA counts + metadata.
    # (drops the sketch/falsecode/negprobes assays, all reductions and graphs;
    #  the lazy `data` layer is dropped too -- process_sample uses counts only.)
    o <- DietSeurat(o, assays = ASSAY, layers = "counts", dimreducs = NULL, graphs = NULL)
    message("  ", sid, ": ", ncol(o), " cells, counts class ",
            class(LayerData(o[[ASSAY]], layer = "counts"))[1])
    objs[[sid]] <- o
  }

  # only prefix cell names if they would actually collide
  all_cells <- unlist(lapply(objs, colnames), use.names = FALSE)
  collide <- anyDuplicated(all_cells) > 0
  message("cell name collision across samples: ", collide)

  merged <- if (collide) {
    merge(objs[[1]], y = objs[-1], add.cell.ids = names(objs))
  } else {
    merge(objs[[1]], y = objs[-1])
  }
  rm(objs); invisible(gc())

  DefaultAssay(merged) <- ASSAY
  merged <- JoinLayers(merged)   # one counts layer, so seu[["RNA"]]$counts works
  stopifnot("counts" %in% Layers(merged[[ASSAY]]))
  message("merged: ", nrow(merged), " x ", ncol(merged),
          " | counts class ", class(LayerData(merged[[ASSAY]], layer = "counts"))[1])
  print(table(merged$merge_sample))

  saveRDS(merged, MERGED_RAW_RDS)
  message("wrote ", MERGED_RAW_RDS, " (", format(file.size(MERGED_RAW_RDS), big.mark = ","), " bytes)")
  mark_complete("01_merged")
} else {
  message("skip merge (already completed)")
}

# ------------------------------------------- JOINT PROCESSING (step 3 code) ----
if (!is_complete("02_processed")) {
  merged <- readRDS(MERGED_RAW_RDS)
  message("joint processing of ", ncol(merged), " cells via process_sample()")
  # writes OUT_DIR/merged_02_processed.Rds + merged_umap_{high,low}.png,
  # and marks 02_processed complete
  merged <- process_sample(merged, output_prefix)

  # give process_sample's auto-named output the step-4 name (rename, no re-save)
  auto <- file.path(OUT_DIR, paste0(output_prefix, "_02_processed.Rds"))
  if (file.exists(auto)) {
    stopifnot(file.rename(auto, MERGED_RDS))
    message("wrote ", MERGED_RDS, " (", format(file.size(MERGED_RDS), big.mark = ","), " bytes)")
  }
} else {
  message("skip joint processing (already completed)")
}

message("Done. Outputs in: ", normalizePath(OUT_DIR))
