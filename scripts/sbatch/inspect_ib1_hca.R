suppressPackageStartupMessages({library(Seurat); library(data.table); setDTthreads(1)})
o <- readRDS("output/03_filter_and_process/IB1_03_HCA_cell_typed.Rds")
cat("cells:", ncol(o), " genes:", nrow(o), "\n")
cat("\nmeta.data columns:\n"); print(colnames(o@meta.data))
md <- o@meta.data
for (cl in c("pearson_clusters_low","pearson_clusters_high"))
  if (cl %in% colnames(md)) { cat("\n", cl, ":\n", sep=""); print(table(md[[cl]])) }
cat("\nassays:", names(o@assays), " layers:", paste(Layers(o[["RNA"]]), collapse=","), "\n")
cat("first cell ids:", paste(head(rownames(md),2), collapse=" | "), "\n")
l <- fread("output/09_hierannot/IB1/hierannot_cluster_labels.csv",
           select=c("cluster_id","annot_path","clustering","hierarchy"))
cat("\n=== pearson_clusters_low labels by hierarchy ===\n")
print(dcast(l[clustering=="pearson_clusters_low"], cluster_id ~ hierarchy, value.var="annot_path"))
cat("\n=== pearson_clusters_high labels by hierarchy ===\n")
print(dcast(l[clustering=="pearson_clusters_high"], cluster_id ~ hierarchy, value.var="annot_path"))
