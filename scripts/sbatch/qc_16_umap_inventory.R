## What is actually in the step-03 and step-04 objects, before building UMAP panels on them.
suppressPackageStartupMessages({library(Seurat); library(data.table); setDTthreads(1)})
F <- c(IB1="output/03_filter_and_process/IB1_02_processed.Rds",
       RB1="output/03_filter_and_process/RB1_02_processed.Rds",
       merged="output/04_merge/merged_processed.Rds")
WANT <- c("inSitu_clusters","pearson_clusters_low","pearson_clusters_high","lineage",
          "merge_sample","Run_Tissue_name","sample_pearson_clusters_low",
          "sample_pearson_clusters_high","nCount_RNA")
for (k in names(F)) {
  f <- F[[k]]
  if (!file.exists(f)) { cat("\n== ", k, ": MISSING ", f, "\n", sep=""); next }
  cat("\n========== ", k, "  (", f, ")\n", sep="")
  o <- readRDS(f)
  cat("cells: ", ncol(o), " | genes: ", nrow(o), "\n", sep="")
  cat("reductions: ", paste(names(o@reductions), collapse=", "), "\n", sep="")
  cat("graphs:     ", paste(names(o@graphs), collapse=", "), "\n", sep="")
  cat("assays:     ", paste(names(o@assays), collapse=", "), "\n", sep="")
  md <- o@meta.data
  have <- intersect(WANT, colnames(md))
  cat("wanted cols present: ", paste(have, collapse=", "), "\n", sep="")
  cat("wanted cols ABSENT : ", paste(setdiff(WANT, colnames(md)), collapse=", "), "\n", sep="")
  for (cc in have) {
    v <- md[[cc]]
    if (is.numeric(v)) next
    t <- sort(table(as.character(v)), decreasing=TRUE)
    cat(sprintf("  %-28s %2d levels: %s\n", cc, length(t),
                paste(sprintf("%s=%d", names(head(t,12)), head(t,12)), collapse=" ")))
  }
  rm(o); invisible(gc())
}
cat("\n========== 09b audited labels available for:\n")
print(list.dirs("output/09b_marker_relabel", recursive=FALSE, full.names=FALSE))
L <- fread("output/09b_marker_relabel/IB1/cluster_labels_final.csv")
cat("\nIB1 pearson_clusters_low audited labels (the panel-3 source):\n")
print(L[clustering=="pearson_clusters_low", .(cluster_id, label, status, cells)])
