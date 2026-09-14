## ============================================================================
## 16_script_celltype_pseudobulk.R -- pseudobulk stratified BY CELL TYPE.
##
## WHY. The cross-sample tissue comparison is not identifiable: IB1's FL foci are 17.6% B
## cells and RB1's follicles are 50.9%, and matching on composition matched 2 of 44 bins
## (16c). Conditioning on cell type removes the problem by construction instead of by
## matching -- a B cell is a B cell in either node, so comparing IB1 FL-focus B cells with
## RB1 follicular B cells has no composition term left to confound it. This is the same move
## step 14 made when it tested expression IN vs OUT *within each cluster* rather than pooled.
##
## TWO GRANULARITIES, because the panel forces it:
##
##   bin x COARSE type   B, T, myeloid, PC, NK, endo, fibro. A 100 um bin holds 6-70 B cells
##                       and 4-18 T cells, so these clear a 5-cell floor often enough to give
##                       a REGION-UNIT test -- the bin stays the unit of the p-value.
##
##   region x FINE type  Tfh, Treg, CD8_T, CD4_conv, Macrophage, cDC, Monocyte. These are
##                       0.03-1.5% of cells, i.e. under one cell per bin. They CANNOT support
##                       a region-unit test and are pooled per region class for effect sizes
##                       only. Said out loud here so no downstream table implies otherwise.
##
## Usage: Rscript scripts/16_script_celltype_pseudobulk.R [tier]
## ============================================================================

suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(data.table); setDTthreads(1)
})
if (!file.exists("scripts/lib/lineage.R")) stop("run from the project root", call. = FALSE)
args <- commandArgs(trailingOnly = TRUE)
TIER <- if (length(args) >= 1) args[1] else "strict"

SAMPLES <- c("IB1","RB1")
PRE_RDS <- "output/01_import_v1/seuratObject_%s_BPCells_full.Rds"
REG_F   <- file.path("output/16_regions", TIER, "cell_regions.csv.gz")
CT_F    <- "output/16_celltypes/cell_celltypes.csv.gz"
OUT     <- file.path("output/16_ct_pseudobulk", TIER)
CHUNK   <- 150000L
MIN_CT  <- 5L     # cells for a (bin, celltype) pseudobulk to exist -- step 14's MIN_CELLS_RC
COARSE  <- c("B","T","myeloid","PC","NK","endo","fibro")
for (f in c(REG_F, CT_F)) if (!file.exists(f)) stop("missing ", f, call. = FALSE)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

cat("============================================================\n")
cat("step 16j -- cell-type-stratified pseudobulk, tier = ", TIER, "\n", sep="")
cat("============================================================\n")

rg <- fread(REG_F, select = c("cell_id","samp","bx","by_","bin_region","arch","is_defect",
                              "unbinned","fov_id","nCount_RNA"))
ct <- fread(CT_F, select = c("cell_id","celltype","lineage16"))
rg <- merge(rg, ct, by = "cell_id", sort = FALSE)
n0 <- nrow(rg); rg <- rg[is_defect == FALSE & unbinned == FALSE]
cat(sprintf("cells %d -> %d after dropping defect and unbinned\n", n0, nrow(rg)))
rg[, bin := paste(samp, bx, by_, sep = "_")]

## ---- group definitions -------------------------------------------------------
rg[, coarse := fifelse(lineage16 %in% COARSE, lineage16, NA_character_)]
rg[, grp_bin := fifelse(is.na(coarse), NA_character_, paste(bin, coarse, sep = "|"))]
rg[, grp_reg := paste(samp, bin_region, celltype, sep = "|")]

GB <- rg[!is.na(grp_bin), .(cells = .N, samp = samp[1], bin = bin[1], coarse = coarse[1],
                            bin_region = bin_region[1], arch = as.character(arch)[1],
                            fov = fov_id[1], nCount = sum(nCount_RNA)), by = grp_bin]
GB <- GB[cells >= MIN_CT]
GR <- rg[, .(cells = .N, samp = samp[1], bin_region = bin_region[1], celltype = celltype[1],
             fovs = uniqueN(fov_id), bins = uniqueN(bin), nCount = sum(nCount_RNA)), by = grp_reg]
cat(sprintf("\n(bin x coarse type) groups with >= %d cells: %d\n", MIN_CT, nrow(GB)))
print(dcast(GB[, .N, by = .(samp, coarse)], coarse ~ samp, value.var = "N", fill = 0L))
cat("\nusable (bin x type) groups per region -- this is the n of every region-unit test:\n")
print(dcast(GB[, .N, by = .(samp, bin_region, coarse)][coarse %in% c("B","T","myeloid")],
            samp + bin_region ~ coarse, value.var = "N", fill = 0L))
cat(sprintf("\n(region x fine type) groups: %d\n", nrow(GR)))

rg[, gb_i := match(grp_bin, GB$grp_bin)]
rg[, gr_i := match(grp_reg, GR$grp_reg)]

## ---- accumulate --------------------------------------------------------------
accB_s <- accB_d <- accR_s <- accR_d <- NULL; GENES <- NULL
for (S in SAMPLES) {
  f <- sprintf(PRE_RDS, S); cat("\nreading ", f, "\n", sep = "")
  o <- readRDS(f); cn <- LayerData(o[["RNA"]], layer = "counts"); rn <- colnames(cn)
  if (is.null(GENES)) {
    GENES <- rownames(cn)
    mk <- function(k, nm) Matrix(0, nrow = length(GENES), ncol = k, sparse = TRUE,
                                 dimnames = list(GENES, nm))
    accB_s <- mk(nrow(GB), GB$grp_bin); accB_d <- accB_s
    accR_s <- mk(nrow(GR), GR$grp_reg); accR_d <- accR_s
  } else if (!identical(GENES, rownames(cn))) stop("gene order differs", call. = FALSE)
  sub <- rg[samp == S]; jj <- match(sub$cell_id, rn)
  if (anyNA(jj)) stop(sum(is.na(jj)), " ", S, " cells absent from counts", call. = FALSE)
  cat(sprintf("  %s: %d cells, %d chunks\n", S, length(jj), ceiling(length(jj)/CHUNK)))
  for (st in seq(1L, length(jj), by = CHUNK)) {
    en <- min(st + CHUNK - 1L, length(jj)); k <- st:en
    mc <- Matrix::drop0(as(cn[, jj[k], drop = FALSE], "dgCMatrix"))
    mb <- mc; mb@x <- rep(1, length(mb@x))
    ok <- !is.na(sub$gb_i[k])
    if (any(ok)) {
      I1 <- sparseMatrix(i = which(ok), j = sub$gb_i[k][ok], x = 1,
                         dims = c(length(k), nrow(GB)))
      accB_s <- accB_s + mc %*% I1; accB_d <- accB_d + mb %*% I1; rm(I1)
    }
    I2 <- sparseMatrix(i = seq_along(k), j = sub$gr_i[k], x = 1,
                       dims = c(length(k), nrow(GR)))
    accR_s <- accR_s + mc %*% I2; accR_d <- accR_d + mb %*% I2
    rm(mc, mb, I2); cat(sprintf("    %d-%d\n", st, en))
  }
  rm(o, cn); invisible(gc())
}

## ---- checks ------------------------------------------------------------------
## The region-level accumulator covers EVERY cell, so its totals must reconcile exactly.
## The bin-level one covers only cells in a qualifying (bin, coarse type) group, so it is
## checked against that subset rather than against everything.
csR <- Matrix::colSums(accR_s)
if (any(abs(csR - GR$nCount) > 1e-6*pmax(GR$nCount,1)))
  stop("region x celltype totals do not reconcile with nCount", call. = FALSE)
csB <- Matrix::colSums(accB_s)
if (any(abs(csB - GB$nCount) > 1e-6*pmax(GB$nCount,1)))
  stop("bin x celltype totals do not reconcile with nCount", call. = FALSE)
cat("\nboth accumulators reconcile with the per-cell depths exactly\n")

saveRDS(accB_s, file.path(OUT,"pb_binct_sum.Rds")); saveRDS(accB_d, file.path(OUT,"pb_binct_det.Rds"))
saveRDS(accR_s, file.path(OUT,"pb_regct_sum.Rds")); saveRDS(accR_d, file.path(OUT,"pb_regct_det.Rds"))
fwrite(GB, file.path(OUT,"binct_meta.csv.gz"), compress="gzip")
fwrite(GR, file.path(OUT,"regct_meta.csv.gz"), compress="gzip")
cat("wrote 4 matrices + 2 metadata tables to ", OUT, "\n", sep="")
writeLines(format(Sys.time()), file.path(OUT, "ct_pseudobulk.done"))
cat("============================================================\nOK\n")
