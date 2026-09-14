# FINAL attempt at the reference transfer. All five strategies, one job.
#
# Previous rounds failed on MY harness, not the strategies:
#   - TransferData with a single-element refdata returns the data.frame DIRECTLY,
#     so p$ct was NULL and every metric came back empty. Handled both shapes now,
#     and the row count is printed so it cannot fail silently again.
#   - CCA anchors carry no `pcaproject` weight.reduction; TransferData needs
#     weight.reduction = "cca".
#
# Diagnosis so far: the grafted pca.full loadings are CORRECT (reference cells
# re-project onto their own embeddings at rho = 1.000), but CosMx query cells land
# in a 5-6x compressed, displaced cloud (sd 1.3 vs 7.6; PC4 mean -9.1 vs +1.3), so
# projecting onto a pre-built scRNA PCA cannot work. S4/S5 build a joint space.
#
# Yardstick: does %predicted-B track %strongB ACROSS CLUSTERS (the requirement --
# useful for interpreting groupings, not per-cell perfection)?
suppressPackageStartupMessages(library(Seurat))
set.seed(1)
N_DIMS <- 30; NQ <- 20000; NREF <- 50000
COARSE <- "annotation_level_1"; BLAB <- c("NBC_MBC","GCBC","PC","preBC")

ref0 <- readRDS("output/00_tonsil_reference/tonsil_ref_processed.Rds")
if (nrow(Loadings(ref0[["pca.full"]])) == 0) {
  Lf <- Loadings(ref0[["pca"]]); colnames(Lf) <- colnames(Embeddings(ref0[["pca.full"]]))
  Loadings(ref0[["pca.full"]], projected = FALSE) <- Lf
}
q0 <- readRDS("output/03_filter_and_process/RB1_03_HCA_cell_typed.Rds")
DefaultAssay(q0) <- "sketch"
q0 <- subset(q0, cells = sample(Cells(q0[["sketch"]]), NQ))
qmd <- q0@meta.data
cat("query:", ncol(q0), "cells | strongB", sprintf("%.1f%%", 100*mean(qmd$strongB)), "\n")

get_pred <- function(p, cells) {
  d <- if (inherits(p, "data.frame")) p else if (!is.null(p$ct)) p$ct else p[[1]]
  ids <- as.character(d$predicted.id)
  names(ids) <- sub("_query$", "", rownames(d))
  cat("   [pred rows:", nrow(d), "| matched:", sum(cells %in% names(ids)), "]\n")
  stopifnot(sum(cells %in% names(ids)) > 0.9*length(cells))
  ids[cells]
}
evaluate <- function(tag, pred) {
  pred <- gsub("-", "_", as.character(pred), fixed = TRUE)
  isB <- pred %in% BLAB
  cat("\n##########", tag, "##########\n")
  print(round(100*sort(table(pred), decreasing = TRUE)/length(pred), 2))
  cat(sprintf("recall strongB -> B-lineage : %.1f%%   (non-strongB: %.1f%%)\n",
              100*mean(isB[qmd$strongB]), 100*mean(isB[!qmd$strongB])))
  cat(sprintf("predicted epithelial        : %.1f%%\n", 100*mean(pred == "epithelial")))
  g <- qmd$pearson_clusters_high; keep <- names(which(table(g) >= 100)); i <- g %in% keep
  a <- tapply(isB[i], g[i], mean); b <- tapply(qmd$strongB[i], g[i], mean)
  cat(sprintf("CLUSTER-LEVEL rho(%%pred-B, %%strongB) = %+.3f (%d clusters); spread %.1f%%-%.1f%%\n",
              cor(a, b, method = "spearman"), length(keep), 100*min(a), 100*max(a)))
}
run <- function(tag, ref, reduction, ..., weight.reduction = NULL) {
  a <- try(FindTransferAnchors(reference = ref, query = q0, normalization.method = "LogNormalize",
           reference.reduction = reduction, dims = 1:N_DIMS, verbose = FALSE, ...), silent = TRUE)
  if (inherits(a, "try-error")) { cat("\n", tag, " ANCHORS FAILED: ",
      conditionMessage(attr(a,"condition")), "\n", sep=""); return(invisible()) }
  rmd <- ref@meta.data; rownames(rmd) <- paste0(colnames(ref), "_reference")
  args <- list(anchorset = a, dims = 1:N_DIMS, verbose = FALSE,
               refdata = list(ct = rmd[a@reference.cells, COARSE]))
  if (!is.null(weight.reduction)) args$weight.reduction <- weight.reduction
  p <- try(do.call(TransferData, args), silent = TRUE)
  if (inherits(p, "try-error")) { cat("\n", tag, " TRANSFER FAILED: ",
      conditionMessage(attr(p,"condition")), "\n", sep=""); return(invisible()) }
  evaluate(tag, get_pred(p, colnames(q0)))
}

r1 <- ref0; DefaultAssay(r1) <- "RNA"
run("S1 current (RNA + grafted pca.full)", r1, "pca.full",
    features = rownames(Loadings(r1[["pca.full"]])))
rm(r1); invisible(gc())

r2 <- subset(ref0, cells = Cells(ref0[["sketch"]])); DefaultAssay(r2) <- "sketch"
run("S2 sketch-consistent (pca on sketch)", r2, "pca")
run("S3 de-novo PCA", r2, NULL, npcs = N_DIMS)
rm(r2); invisible(gc())

r3 <- subset(ref0, cells = sample(Cells(ref0[["sketch"]]), NREF)); DefaultAssay(r3) <- "sketch"
r3 <- ScaleData(FindVariableFeatures(r3, verbose = FALSE), verbose = FALSE)
q0 <- ScaleData(FindVariableFeatures(q0, verbose = FALSE), verbose = FALSE)
run("S4 cca (joint space)", r3, NULL, reduction = "cca", weight.reduction = "cca")
run("S5 rpca + project.query", r3, NULL, reduction = "rpca", project.query = TRUE)

cat("\nTARGET: rho clearly > 0, epithelial far below 53%, GCBC/NBC_MBC represented.\n\nOK\n")
