## ============================================================================
## 16_script_panel_counts.R -- the step-16 gene panel, all cells, both samples.
##
## Companion to 16_script_cell_meta.R: that script extracts metadata, this one extracts
## counts, and between them every later 16* script runs without touching a Seurat object.
##
## WHY A PANEL AND NOT THE WHOLE MATRIX. The full 6K x 4.4M pre-QC matrix is not something
## to materialise repeatedly, and only step 16d (genome-wide DE) actually needs it -- and
## 16d works at the region level, so it can stream. Everything else in step 16 -- the
## chemokine bins, the architecture axis, the cell-type subsets, the checkpoint readout --
## touches at most ~100 genes. Those are extracted once, here.
##
## The negprobes assay is saved alongside, because the checkpoint readout in 16e needs a
## per-cell background rate: CD274 is detected in 5.6% of IB1 cells and PDCD1LG2 in 5.4%,
## so "positive" is one or two transcripts and a negprobe-matched false-positive rate is
## the difference between a result and an artifact.
##
## Usage: sbatch --mem=140G --cpus-per-task=2 --time=2:00:00 \
##          scripts/sbatch/run_rscript.sh scripts/16_script_panel_counts.R
## ============================================================================

suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(data.table); setDTthreads(1)
})

SAMPLES <- c("IB1", "RB1")
PRE_RDS <- "output/01_import_v1/seuratObject_%s_BPCells_full.Rds"
OUT     <- "output/16_regions"

if (!file.exists("scripts/lib/cluster_labels.R")) stop("run from the project root", call. = FALSE)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

## ---- the panel ---------------------------------------------------------------
## Grouped by the role each gene plays in step 16, so that a later reader can tell a
## SELECTION gene (one that defines a region, and must therefore be excluded from any
## test on that region) from a READOUT gene. The groups are written out alongside the
## matrix for exactly that reason.
PANEL <- list(
  ## --- selection genes: these DEFINE regions and are circular within them ---
  ## "IGLC1/2" is one combined probe on this panel, not two genes -- 07_script_light
  ## chain_neighborhoods.R:67 uses that exact string. IGLC1/IGLC2 do not exist here.
  sel_lightchain = c("IGKC","IGLC1/2"),                          # define FL foci (step 07)
  sel_chemokine  = c("CCL17","CCL22"),                           # define cHL hotspots (step 11)
  sel_bcell      = c("MS4A1","CD79A","CD79B","BANK1"),           # define B_dense architecture
  ## --- checkpoint axis: the question ---
  checkpoint_primary = c("PDCD1","CD274","PDCD1LG2"),
  checkpoint_context = c("CTLA4","LAG3","TIGIT","HAVCR2","TOX","ENTPD1","TNFRSF9",
                         "BTLA","VSIR","LGALS9","IDO1","ICOS","CD80","CD86"),
  ## --- chemokine receptors and interferon chemokines ---
  chemokine_axis = c("CCR4","CXCR5","CCR7","CXCL9","CXCL10","CXCL11","CXCL13",
                     "CCL19","CCL21"),
  ## --- T cell and subsets ---
  t_core = c("CD3D","CD3E","CD3G","TRAC"),
  t_cd4  = c("CD4","IL7R"),
  t_cd8  = c("CD8A","GZMB","GZMK","NKG7","IFNG"),
  t_treg = c("FOXP3","IL2RA","IKZF2"),
  ## --- myeloid, DC, NK, plasma ---
  myeloid = c("CD68","CD163","C1QA","C1QB","MRC1","LYZ","ITGAX","CD14","CSF1R"),
  dc      = c("BATF3","CLEC9A","CD1C","FCER1A","LAMP3"),
  pdc     = c("LILRA4","IRF7","IL3RA"),
  nk      = c("GNLY","KLRD1","NCR1"),
  plasma  = c("MZB1","XBP1","PRDM1","SDC1","DERL3","TNFRSF17"),
  ## --- stroma, FDC, endothelium ---
  stroma  = c("COL1A1","COL1A2","PDGFRB","DCN","ACTA2","CR2","C3","VCAM1"),
  endo    = c("PECAM1","VWF","CLDN5","CD34"),
  ## --- the FL / cHL tumour programmes ---
  fl_prog = c("BCL2","BCL6","MME","AICDA","RGS13","LMO2"),        # MEF2B not on the panel
  chl_prog= c("CD40","TNFRSF8","CD40LG","IRF4","TRAF1","RELB","NFKB2","SOCS1","GATA3"),
  cycling = c("MKI67","TOP2A")
)
GENES <- unique(unlist(PANEL, use.names = FALSE))
cat(sprintf("panel: %d genes in %d groups\n", length(GENES), length(PANEL)))

## ---- extract -----------------------------------------------------------------
## One object at a time, dropped before the next is read, as step 11 does -- peak memory
## is one object rather than two.
read_one <- function(S) {
  f <- sprintf(PRE_RDS, S)
  if (!file.exists(f)) stop("missing ", f, call. = FALSE)
  cat("reading ", f, "\n", sep = "")
  o  <- readRDS(f)
  rn <- rownames(o@meta.data)
  cn <- LayerData(o[["RNA"]], layer = "counts")
  g  <- intersect(GENES, rownames(cn))
  M  <- as(cn[g, , drop = FALSE], "dgCMatrix"); colnames(M) <- rn
  ## negprobes: every probe, not a total -- 16e needs the per-probe detection rate to
  ## calibrate what "one transcript" means for a gene like CD274.
  np <- o@assays$negprobes@counts
  NP <- as(np[, , drop = FALSE], "dgCMatrix"); colnames(NP) <- rn
  rm(o, cn, np); invisible(gc())
  cat(sprintf("  %s: %d cells, %d/%d panel genes, %d negprobes\n",
              S, ncol(M), length(g), length(GENES), nrow(NP)))
  list(M = M, NP = NP, G = g)
}

parts <- lapply(SAMPLES, read_one)
G  <- Reduce(intersect, lapply(parts, `[[`, "G"))
NPG<- Reduce(intersect, lapply(parts, function(z) rownames(z$NP)))
absent <- setdiff(GENES, G)
if (length(absent)) message("absent from the panel and dropped: ",
                            paste(absent, collapse = ", "))

M  <- do.call(cbind, lapply(parts, function(z) z$M[G, , drop = FALSE]))
NP <- do.call(cbind, lapply(parts, function(z) z$NP[NPG, , drop = FALSE]))
rm(parts); invisible(gc())
stopifnot(identical(colnames(M), colnames(NP)))
cat(sprintf("\ncombined: %d genes x %d cells (%.1f%% nonzero), %d negprobes\n",
            nrow(M), ncol(M), 100 * length(M@x) / (as.numeric(nrow(M)) * ncol(M)), nrow(NP)))

## ---- write -------------------------------------------------------------------
saveRDS(M,  file.path(OUT, "panel_counts_all.Rds"))
saveRDS(NP, file.path(OUT, "negprobe_counts_all.Rds"))
grp <- rbindlist(lapply(names(PANEL), function(k)
  data.table(group = k, gene = PANEL[[k]])))[gene %in% G]
grp[, is_selection := grepl("^sel_", group)]
fwrite(grp, file.path(OUT, "panel_groups.csv"))
cat("wrote ", file.path(OUT, "panel_counts_all.Rds"), "\n", sep = "")
cat("wrote ", file.path(OUT, "negprobe_counts_all.Rds"), "\n", sep = "")
cat("wrote ", file.path(OUT, "panel_groups.csv"), "\n", sep = "")

## Detection rate per gene per sample, as a sanity table that can be read without loading
## the matrix -- and as a cross-check against the published numbers (CCL17 13.49%/8.46%,
## CD274 5.64%/4.11%, PDCD1LG2 5.37%/4.41% in output/10_hrs/chl_program_per_gene.csv,
## which were computed post-QC, so these pre-QC rates will differ and should not match).
samp <- sub("^c_1_.*$", "RB1", sub("^c_2_.*$", "IB1", colnames(M)))
det <- rbindlist(lapply(SAMPLES, function(S) {
  ii <- which(samp == S)
  data.table(samp = S, gene = rownames(M),
             pct_det = round(100 * Matrix::rowMeans(M[, ii, drop = FALSE] > 0), 3),
             mean_ct = round(Matrix::rowMeans(M[, ii, drop = FALSE]), 4))
}))
fwrite(dcast(det, gene ~ samp, value.var = c("pct_det","mean_ct")),
       file.path(OUT, "panel_detection_preqc.csv"))
cat("wrote ", file.path(OUT, "panel_detection_preqc.csv"), "\n", sep = "")
writeLines(format(Sys.time()), file.path(OUT, "panel_counts_all.done"))
cat("============================================================\nOK\n")
