## ============================================================================
## 16_script_pseudobulk.R -- genome-wide gene x BIN pseudobulk, both samples, one tier.
##
## THE REGION UNIT IS THE 100 um BIN, and that is a deliberate choice with a cost.
##
## The obvious unit is the light-chain neighbourhood, as step 14 used. It cannot work here:
## only 16 of IB1's 154 neighbourhoods are kappa-restricted and non-artifact, which is
## barely above step 14's own MIN_REG = 12 floor and gives no power at all for a
## cross-sample contrast. The 100 um bin is the unit step 11 already established, it is
## plentiful (79 FL_focus / 29 hotspot / 2755 background bins in IB1 at the strict tier;
## 2605 B_dense / 2606 B_sparse in RB1), and a bin score is a share of the bin's
## transcriptome, so cellularity and depth are controlled by construction.
##
## THE COST, STATED UP FRONT: bins within one FOV are not independent. IB1's 79 FL_focus
## bins sit in only 11 FOVs. That is the same pseudo-replication trap that cost step 11 a
## Tfh result, one level up. So 16d tests every contrast TWICE -- once with the bin as the
## unit and once with bins collapsed to their FOV -- and reports both. A result that
## survives only the bin-level test is reported as unresolved.
##
## WHY GENOME-WIDE HERE AND PANEL-ONLY IN 16e. A gene x (bin x celltype) array would be
## ~6000 x 130,000 and is not worth materialising. Genome-wide work therefore happens at
## the bin level (this script), and anything stratified by cell type runs off the 103-gene
## panel matrix, which is cheap. Neither compromise is hidden: 16e says which genes it saw.
##
## Counts are accumulated in cell chunks so peak memory is one chunk, not one sample.
##
## Usage: Rscript scripts/16_script_pseudobulk.R [tier]      tier = strict | matched
##        sbatch --mem=200G --cpus-per-task=4 --time=6:00:00 \
##          scripts/sbatch/run_rscript.sh scripts/16_script_pseudobulk.R matched
## ============================================================================

suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(data.table); setDTthreads(1)
})

if (!file.exists("scripts/lib/lineage.R")) stop("run from the project root", call. = FALSE)
args <- commandArgs(trailingOnly = TRUE)
TIER <- if (length(args) >= 1) args[1] else "strict"
if (!TIER %in% c("strict","matched")) stop("tier must be strict or matched", call. = FALSE)

SAMPLES <- c("IB1","RB1")
PRE_RDS <- "output/01_import_v1/seuratObject_%s_BPCells_full.Rds"
REG_F   <- file.path("output/16_regions", TIER, "cell_regions.csv.gz")
OUT     <- file.path("output/16_pseudobulk", TIER)
CHUNK   <- 150000L
if (!file.exists(REG_F)) stop("missing ", REG_F, " -- run 16_script_regions.R ", TIER,
                              call. = FALSE)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

cat("============================================================\n")
cat("step 16d-pb -- gene x bin pseudobulk, tier = ", TIER, "\n", sep="")
cat("============================================================\n")

## ---- 1. the bins we are aggregating into ------------------------------------
## Defect cells and unbinned cells are both excluded here. Bins were SCORED on defect-free
## cells in 16a, so including cluster `i` now would make the pseudobulk describe a different
## cell set than the region label does.
## `bin_region`, not `region`: the per-cell `region` splits a bin that straddles an FL
## footprint boundary, and this script's own check below rejected 111 such bins. 16a now
## emits a bin-level label built from the SHARE of a bin's cells inside the footprint.
rg <- fread(REG_F, select = c("cell_id","samp","bx","by_","bin_region","arch","is_defect",
                              "unbinned","fov_id","nCount_RNA","Area.um2"))
n0 <- nrow(rg)
rg <- rg[is_defect == FALSE & unbinned == FALSE]
cat(sprintf("cells %d -> %d after dropping defect and unbinned\n", n0, nrow(rg)))
rg[, bin := paste(samp, bx, by_, sep = "_")]
BINS <- rg[, .(cells = .N, fov = fov_id[1], bin_region = bin_region[1],
               arch = as.character(arch)[1],
               nCount = sum(nCount_RNA), med_depth = as.numeric(median(nCount_RNA)),
               med_area = as.numeric(median(Area.um2)),
               bx = bx[1], by_ = by_[1]), by = .(samp, bin)]
## A bin must have ONE region and ONE arch label; if a bin spans two, the label is being
## assigned per cell somewhere upstream and every aggregate below would be silently mixed.
chk <- rg[, .(nr = uniqueN(bin_region), na = uniqueN(arch)), by = bin][nr > 1 | na > 1]
if (nrow(chk)) stop(nrow(chk), " bins carry more than one bin_region/arch label", call. = FALSE)
setkey(BINS, bin)
cat(sprintf("bins: %d\n", nrow(BINS)))
print(BINS[, .(bins = .N, cells = sum(cells), fovs = uniqueN(fov)), by = .(samp, bin_region)][order(samp, -bins)])

bin_lvl <- BINS$bin
rg[, bin_i := match(bin, bin_lvl)]

## ---- 2. accumulate, in chunks ----------------------------------------------
acc_sum <- NULL; acc_det <- NULL; GENES <- NULL
for (S in SAMPLES) {
  f <- sprintf(PRE_RDS, S)
  if (!file.exists(f)) stop("missing ", f, call. = FALSE)
  cat("\nreading ", f, "\n", sep = "")
  o  <- readRDS(f)
  cn <- LayerData(o[["RNA"]], layer = "counts")
  rn <- colnames(cn)
  if (is.null(GENES)) {
    GENES  <- rownames(cn)
    acc_sum <- Matrix(0, nrow = length(GENES), ncol = length(bin_lvl), sparse = TRUE,
                      dimnames = list(GENES, bin_lvl))
    acc_det <- acc_sum
  } else if (!identical(GENES, rownames(cn))) {
    stop("gene order differs between samples", call. = FALSE)
  }
  sub <- rg[samp == S]
  jj  <- match(sub$cell_id, rn)
  if (anyNA(jj)) stop(sum(is.na(jj)), " ", S, " cells absent from the counts matrix",
                      call. = FALSE)
  cat(sprintf("  %s: %d cells into %d bins, %d chunks of %d\n", S, length(jj),
              uniqueN(sub$bin_i), ceiling(length(jj)/CHUNK), CHUNK))
  for (st in seq(1L, length(jj), by = CHUNK)) {
    en <- min(st + CHUNK - 1L, length(jj))
    mc <- Matrix::drop0(as(cn[, jj[st:en], drop = FALSE], "dgCMatrix"))
    IND <- sparseMatrix(i = seq_len(en - st + 1L), j = sub$bin_i[st:en],
                        x = 1, dims = c(en - st + 1L, length(bin_lvl)))
    ## Detection indicator built by overwriting the values of the count matrix rather than
    ## via `mc > 0`, which yields an lgCMatrix whose %*% behaviour is version-dependent.
    ## drop0() above guarantees no explicit zeros, so x <- 1 is exactly "detected".
    mb <- mc; mb@x <- rep(1, length(mb@x))
    acc_sum <- acc_sum + mc %*% IND
    acc_det <- acc_det + mb %*% IND
    rm(mc, mb, IND)
    cat(sprintf("    %d-%d done\n", st, en))
  }
  rm(o, cn); invisible(gc())
}

## ---- 3. checks --------------------------------------------------------------
## The pseudobulk column sums must equal the per-cell depth sums the region table reports.
## If they do not, the cell-to-bin mapping is wrong and nothing downstream is meaningful.
cs  <- Matrix::colSums(acc_sum)
tgt <- BINS[match(bin_lvl, bin), nCount]
bad <- which(abs(cs - tgt) > 1e-6 * pmax(tgt, 1))
cat(sprintf("\n=== column-sum check: %d of %d bins disagree with the region table\n",
            length(bad), length(cs)))
if (length(bad)) {
  print(data.table(bin = bin_lvl[head(bad, 5)], pb = cs[head(bad,5)], meta = tgt[head(bad,5)]))
  stop("pseudobulk totals do not match nCount_RNA sums -- the cell-to-bin join is wrong",
       call. = FALSE)
}
cat("pseudobulk totals reconcile with the per-cell depths exactly\n")
cat(sprintf("genes %d | bins %d | nonzero %.1f%%\n", nrow(acc_sum), ncol(acc_sum),
            100*length(acc_sum@x)/(as.numeric(nrow(acc_sum))*ncol(acc_sum))))

## ---- 4. write ---------------------------------------------------------------
saveRDS(acc_sum, file.path(OUT, "pb_sum.Rds"))
saveRDS(acc_det, file.path(OUT, "pb_det.Rds"))
fwrite(BINS, file.path(OUT, "bin_meta.csv.gz"), compress = "gzip")
cat("wrote ", file.path(OUT, "pb_sum.Rds"), "\n", sep="")
cat("wrote ", file.path(OUT, "pb_det.Rds"), "\n", sep="")
cat("wrote ", file.path(OUT, "bin_meta.csv.gz"), "\n", sep="")
writeLines(format(Sys.time()), file.path(OUT, "pseudobulk.done"))
cat("============================================================\nOK\n")
