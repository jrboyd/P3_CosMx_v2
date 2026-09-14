# Is the reference transfer miscalibrated by the grafted pca.full loadings?
#
# pca.full embeddings were built by ProjectData() using the SKETCH's feature
# means/sds; FindTransferAnchors re-derives them from the FULL RNA assay. The
# sketch is leverage-weighted toward rare states, so those stats differ and query
# cells may be projected into a displaced region of PC space -- which would sort
# them by signal magnitude, not identity (epithelial 0.11% of ref -> 50% of calls,
# GCBC 19.95% of ref -> 2 cells).
#
# TEST A: self-transfer control. Reference cells used as their own query MUST get
#         their own labels back. If they land in `epithelial` too, the projection
#         path is broken independently of any CosMx/scRNA technology gap.
# TEST B: embedding distributions. Compare per-PC mean/sd of reference pca.full
#         embeddings vs freshly projected embeddings for ref / IB1 / RB1 cells.
#         A large offset or scale mismatch is direct evidence of miscalibration.
suppressPackageStartupMessages(library(Seurat))
set.seed(1)
N_DIMS <- 30; NQ <- 20000
COARSE <- "annotation_level_1"

ref <- readRDS("output/00_tonsil_reference/tonsil_ref_processed.Rds")
# same fix as the pipeline applies
if (nrow(Loadings(ref[["pca.full"]])) == 0) {
  L <- Loadings(ref[["pca"]]); colnames(L) <- colnames(Embeddings(ref[["pca.full"]]))
  Loadings(ref[["pca.full"]], projected = FALSE) <- L
}
DefaultAssay(ref) <- "RNA"
FEAT <- rownames(Loadings(ref[["pca.full"]]))
cat("features:", length(FEAT), "| ref cells:", ncol(ref), "\n")

emb_ref <- Embeddings(ref[["pca.full"]])
cat("\n=== TEST B: per-PC mean / sd ===\n")
cat(sprintf("%-22s %s\n", "source", paste(sprintf("PC%-2d", 1:6), collapse="   ")))
show <- function(tag, e) {
  cat(sprintf("%-22s %s\n", paste0(tag, " mean"),
              paste(sprintf("%+7.2f", colMeans(e[, 1:6, drop=FALSE])), collapse=" ")))
  cat(sprintf("%-22s %s\n", paste0(tag, " sd"),
              paste(sprintf("%7.2f", apply(e[, 1:6, drop=FALSE], 2, sd)), collapse=" ")))
}
show("REF pca.full", emb_ref)

qcells <- sample(colnames(ref), NQ)
qref <- subset(ref, cells = qcells)
DefaultAssay(qref) <- "RNA"
pe <- ProjectCellEmbeddings(query = qref, reference = ref, reduction = "pca.full",
                            normalization.method = "LogNormalize", dims = 1:N_DIMS,
                            verbose = FALSE)
show("ref-as-query proj", pe)
cat("  (should match REF pca.full closely -- same cells, same space)\n")
cat("  per-PC correlation with their own stored embeddings:",
    paste(sprintf("%.3f", sapply(1:6, function(i)
      cor(pe[qcells, i], emb_ref[qcells, i]))), collapse=" "), "\n")

for (s in c("IB1", "RB1")) {
  q <- readRDS(file.path("output/03_filter_and_process", paste0(s, "_03_HCA_cell_typed.Rds")))
  DefaultAssay(q) <- "sketch"
  q <- subset(q, cells = sample(Cells(q[["sketch"]]), NQ))
  pq <- ProjectCellEmbeddings(query = q, reference = ref, reduction = "pca.full",
                              normalization.method = "LogNormalize", dims = 1:N_DIMS,
                              verbose = FALSE)
  show(paste(s, "proj"), pq)
  rm(q, pq); invisible(gc())
}

cat("\n=== TEST A: self-transfer control (reference cells as query) ===\n")
a <- FindTransferAnchors(reference = ref, query = qref, normalization.method = "LogNormalize",
                         reference.reduction = "pca.full", features = FEAT,
                         dims = 1:N_DIMS, verbose = FALSE)
rmd <- ref@meta.data; rownames(rmd) <- paste0(colnames(ref), "_reference")
pred <- TransferData(anchorset = a, dims = 1:N_DIMS, verbose = FALSE,
                     refdata = list(ct = rmd[a@reference.cells, COARSE]))
truth <- as.character(ref@meta.data[rownames(pred$ct), COARSE])
got <- as.character(pred$ct$predicted.id)
got <- gsub("-", "_", got, fixed = TRUE)   # undo the underscore mangling
cat("self-transfer accuracy:", sprintf("%.1f%%\n", 100*mean(got == truth, na.rm=TRUE)))
cat("\npredicted composition vs truth (%):\n")
lv <- union(names(sort(table(truth), decreasing=TRUE)), unique(got))
print(data.frame(truth = round(100*table(factor(truth, lv))/length(truth), 2),
                 predicted = round(100*table(factor(got, lv))/length(got), 2))[, c(2,4)])
cat("\ntrue GCBC predicted as (%):\n")
print(round(100*sort(table(got[truth=="GCBC"]), decreasing=TRUE)/sum(truth=="GCBC"), 1))
cat("\nOK\n")
