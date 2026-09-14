## Does the object 09b audited (IB1_03_HCA_cell_typed.Rds) share a clustering with the object
## that carries the UMAP (_02_processed.Rds)? If the cluster ids do not align, joining 09b's
## audited labels onto the UMAP object by cluster_id silently mislabels every cell -- which is
## the exact failure mode scripts/lib/cluster_labels.R exists to prevent.
suppressPackageStartupMessages({library(Seurat); library(data.table); setDTthreads(1)})
L <- fread("output/09b_marker_relabel/IB1/cluster_labels_final.csv")[clustering=="pearson_clusters_low"]
for (S in c("IB1","RB1")) {
  for (tag in c("02_processed","03_HCA_cell_typed")) {
    f <- sprintf("output/03_filter_and_process/%s_%s.Rds", S, tag)
    if (!file.exists(f)) { cat(S, tag, ": MISSING\n"); next }
    o <- readRDS(f); md <- o@meta.data
    cat(sprintf("\n%s / %-18s cells=%7d  reductions=%s\n", S, tag, ncol(o),
                paste(names(o@reductions), collapse=",")))
    if ("pearson_clusters_low" %in% colnames(md)) {
      t <- sort(table(as.character(md$pearson_clusters_low)), decreasing=TRUE)
      cat("   pearson_clusters_low: ",
          paste(sprintf("%s=%d", names(t), as.integer(t)), collapse=" "), "\n", sep="")
    }
    if ("inSitu_clusters" %in% colnames(md))
      cat("   cluster i cells: ", sum(md$inSitu_clusters == "i"), "\n", sep="")
    rm(o); invisible(gc())
  }
}
cat("\n=== 09b's own per-cluster counts for IB1 pearson_clusters_low\n")
print(L[, .(cluster_id, label, cells)])
cat("\nIf these match ONE of the objects above exactly, that is the object to join on.\n")
