library(Seurat)
library(SingleCellExperiment)
library(HCATonsilData)

flat_files = dir("/netfiles/volaric_research/DLBCL/CosMx/export1/flatFiles/IB1", full.names = TRUE)
names(flat_files) = basename(flat_files)
flat_files = as.list(flat_files)

seu_files = dir("/netfiles/volaric_research/DLBCL/CosMx/", pattern = "RDS$", full.names = TRUE)
names(seu_files) = basename(seu_files)
seu_files = as.list(seu_files)
seu = readRDS("output/01_import_v1/seuratObject_IB1_BPCells_full.Rds")
seu = readRDS("output/01_import_v1/seuratObject_RB1_BPCells_full.Rds")
GENES = rownames(seu)
# tx_df = data.table::fread(flat_files$IB1_tx_file.csv.gz)
# data.table::fread(flat_files$IB1_exprMat_file.csv.gz)
N_DIMS = 30
RES = 1
HIGH_RES = 2
ASSAY         <- "RNA"
# --------------------------------------------------- 1. TONSIL REFERENCE ------
# Load + process ONCE (label source for cell-type assignment). Cached.
REF_RDS <- file.path("output/00_tonsil_reference/", "tonsil_ref_processed.Rds")
dir.create(dirname(REF_RDS), recursive = TRUE, showWarnings = FALSE)
N_SKETCH<- 100000
if (file.exists(REF_RDS)) {
  message("load cached tonsil reference")
  ref <- readRDS(REF_RDS)
} else {
  message("building tonsil reference (HCATonsilData)")
  ref.sce <- HCATonsilData(assayType = "RNA", cellType = "All", version = "2.0")
  # print(names(colData(ref.sce)))   # <- inspect; prefer a COARSE lineage level
  # idx <- sample(ncol(ref.sce), min(REF_N, ncol(ref.sce)))

  # 1. sort the indices so column reads are sequential, faster than random access
  # idx = sort(idx)

  # 2. cut I/O volume: you only need panel-overlapping genes for anchor transfer.
  #    Subsetting rows too shrinks what must be read off disk.
  #    (genes = your CosMx panel rownames)
  keep_g <- intersect(rownames(ref.sce), GENES)



  # 3. realize the block ONCE into a sparse matrix, then build the object
  cnt <- as(counts(ref.sce)[keep_g, ], "CsparseMatrix")   # single materialization
  seu <- CreateSeuratObject(counts = cnt,
                            meta.data = as.data.frame(colData(ref.sce)))

  seu <- NormalizeData(seu)
  seu <- FindVariableFeatures(seu)
  seu <- SketchData(seu, ncells = N_SKETCH, method = "LeverageScore",
                    sketched.assay = "sketch")

  DefaultAssay(seu) <- "sketch"
  seu <- NormalizeData(seu)
  seu <- FindVariableFeatures(seu)
  seu <- ScaleData(seu)
  seu <- RunPCA(seu, npcs = N_DIMS)
  seu <- FindNeighbors(seu, dims = 1:N_DIMS)
  seu <- FindClusters(seu, resolution = RES, cluster.name = "sketch_clusters")
  seu <- FindClusters(seu, resolution = HIGH_RES, cluster.name = "sketch_clusters_high")
  seu <- RunUMAP(seu, dims = 1:N_DIMS, return.model = TRUE,
                 reduction.name = "umap.sketch")

  DimPlot(seu)
  DimPlot(seu, group.by = "annotation_level_1")

  colnames(seu@meta.data)

  seu <- ProjectData(seu, assay = ASSAY,
                     full.reduction    = "pca.full",
                     sketched.assay    = "sketch",
                     sketched.reduction= "pca",
                     umap.model        = "umap.sketch",
                     dims              = 1:N_DIMS,
                     refdata = list(
                       cluster         = "sketch_clusters",
                       cluster_high    = "sketch_clusters_high"
                     )
  )

  saveRDS(GENES, file.path(OUT_DIR, "output/00_tonsil_reference/cosmx_genes.Rds"))
  saveRDS(seu, REF_RDS)
}
