library(Seurat)
library(tidyverse)
seu = readRDS("a")
dim(seu)

md = seu@meta.data %>% as.data.frame
md %>% colnames %>% tail(n = 30)
md$inSitu_clusters
DefaultDimReduc(seu) = "pearsonumap"
DimPlot(seu)
DimPlot(seu, group.by = "celltype_coarse")
DimPlot(seu, group.by = "lineage")
DimPlot(seu, group.by = "pearson_clusters_low")
DimPlot(seu, group.by = "pearson_clusters_high")
DimPlot(seu, group.by = "inSitu_clusters")
DimPlot(seu %>% subset(inSitu_clusters %in% c("a", "j")), group.by = "inSitu_clusters")
# DimPlot(seu, group.by = "celltype_fine")
seu@reductions
