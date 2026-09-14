## ============================================================================
## qc_16_batf3_spillover.R -- why does BATF3 rise in EVERY cell type in a cHL hotspot?
##
## The region heatmap shows BATF3 elevated in the IB1|cHL_hotspot column of PC, NK, endo,
## fibro, T_unresolved and B -- cell types that have no business expressing a cDC
## transcription factor. Three explanations, and they are distinguishable:
##
##   1. small-n variance   those columns are 114-307 cells and BATF3 is a sparse gene
##                         (4.6% detection). Tested by the |z|-vs-n diagnostic: real but
##                         modest, median rho -0.25 for genes under 6% detection.
##   2. a real field       BATF3 is the #1 gene genome-wide correlating with the CCL17/CCL22
##                         field (rho 0.555, step 16g), so a paracrine/field effect is live.
##   3. SEGMENTATION SPILLOVER  hotspots are DC-rich; a BATF3-high cDC bleeding transcripts
##                         into its neighbours would raise BATF3 in every adjacent cell
##                         regardless of type, at roughly the 1-transcript-per-cell level the
##                         heatmap columns actually show.
##
## 2 and 3 make opposite predictions about DISTANCE TO THE NEAREST cDC. A field does not care;
## spillover decays within a cell diameter. So: BATF3 in non-cDC cells, in hotspot vs
## background, stratified by distance to the nearest cDC. If the hotspot elevation survives at
## distance, it is a field. If it vanishes, the heatmap row is an artifact of segmentation.
##
## Depth-matched, because BATF3 is sparse and hotspot bins are shallower than background.
##
## Usage: sbatch --mem=64G scripts/sbatch/run_rscript.sh scripts/sbatch/qc_16_batf3_spillover.R
## ============================================================================
suppressPackageStartupMessages({ library(Matrix); library(RANN); library(data.table); setDTthreads(1) })

TIER <- "matched"
rg <- fread(file.path("output/16_regions", TIER, "cell_regions.csv.gz"),
            select=c("cell_id","samp","bin_region","x_slide_mm","y_slide_mm",
                     "is_defect","unbinned","nCount_RNA","fov_id"))
ct <- fread("output/16_celltypes/cell_celltypes.csv.gz", select=c("cell_id","celltype"))
d  <- merge(rg, ct, by="cell_id", sort=FALSE)[is_defect==FALSE & unbinned==FALSE & samp=="IB1"]
M  <- readRDS("output/16_regions/panel_counts_all.Rds")
jj <- match(d$cell_id, colnames(M)); stopifnot(!anyNA(jj))
d[, batf3 := as.vector(M["BATF3", jj])]
rm(M); invisible(gc())
OUT <- "output/16_checkpoint"; dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

cdc <- d[celltype == "cDC"]
cat(sprintf("IB1 cells %d | cDC %d | BATF3 detected in %.2f%% overall\n",
            nrow(d), nrow(cdc), 100*d[, mean(batf3 > 0)]))
cat("\nBATF3 by cell type (all IB1 regions), cp10k:\n")
print(d[, .(cells=.N, cp10k=round(1e4*sum(batf3)/sum(nCount_RNA),2),
            pct_det=round(100*mean(batf3>0),2)), by=celltype][order(-cp10k)][1:6])

## distance from every cell to the nearest cDC (excluding itself)
xy  <- as.matrix(d[,   .(x_slide_mm*1000, y_slide_mm*1000)])
cxy <- as.matrix(cdc[, .(x_slide_mm*1000, y_slide_mm*1000)])
nn  <- RANN::nn2(cxy, xy, k=2)
d[, d_cdc := fifelse(celltype=="cDC", nn$nn.dists[,2], nn$nn.dists[,1])]
d[, dband := cut(d_cdc, c(0,15,30,60,120,Inf),
                 labels=c("0-15","15-30","30-60","60-120",">120"), right=FALSE)]

## Depth-match hotspot to background WITHIN each distance band, by depth decile, so a
## sparse-gene comparison is not reading a depth difference.
sub <- d[celltype != "cDC" & bin_region %in% c("cHL_hotspot","background")]
sub[, ddec := cut(nCount_RNA, quantile(nCount_RNA, seq(0,1,.1)), include.lowest=TRUE, labels=FALSE)]
res <- rbindlist(lapply(levels(sub$dband), function(b) {
  z <- sub[dband == b]
  if (z[bin_region=="cHL_hotspot", .N] < 50) return(NULL)
  ## per depth decile, pool hotspot and an equal-sized background sample, then combine
  set.seed(1)
  m <- rbindlist(lapply(split(seq_len(nrow(z)), z$ddec), function(k) {
    zz <- z[k]; nh <- zz[bin_region=="cHL_hotspot", .N]
    if (!nh) return(NULL)
    bg <- zz[bin_region=="background"]; if (!nrow(bg)) return(NULL)
    rbind(zz[bin_region=="cHL_hotspot"], bg[sample.int(.N, min(.N, nh*5))])
  }))
  if (!nrow(m)) return(NULL)
  m[, .(dband=b, cells=.N,
        cp10k=1e4*sum(batf3)/sum(nCount_RNA), pct_det=100*mean(batf3>0),
        med_depth=as.numeric(median(nCount_RNA))), by=bin_region]
}))
if (!nrow(res)) stop("no comparable strata", call.=FALSE)
w <- dcast(res, dband ~ bin_region, value.var=c("cells","cp10k","pct_det","med_depth"))
w[, ratio_cp10k := round(cp10k_cHL_hotspot/cp10k_background, 2)]
w[, ratio_det   := round(pct_det_cHL_hotspot/pct_det_background, 2)]
for (cc in grep("^cp10k|^pct_det|^med_depth", names(w), value=TRUE)) w[[cc]] <- round(w[[cc]], 2)
cat("\n=== BATF3 in NON-cDC cells: hotspot vs background, by distance to the nearest cDC\n")
cat("    depth-matched within decile; ratio > 1 means higher in the hotspot\n")
print(w[, .(dband, n_hot=cells_cHL_hotspot, n_bg=cells_background,
            cp10k_hot=cp10k_cHL_hotspot, cp10k_bg=cp10k_background, ratio_cp10k,
            det_hot=pct_det_cHL_hotspot, det_bg=pct_det_background, ratio_det,
            depth_hot=med_depth_cHL_hotspot, depth_bg=med_depth_background)])
fwrite(w, file.path(OUT, "batf3_spillover.csv"))
far <- w[dband == ">120"]
cat(sprintf("\nVERDICT: at >120 um from any cDC the hotspot:background BATF3 ratio is %s (cp10k) / %s (detection)\n",
            if (nrow(far)) far$ratio_cp10k else NA, if (nrow(far)) far$ratio_det else NA))
cat(if (nrow(far) && is.finite(far$ratio_cp10k) && far$ratio_cp10k > 1.3)
      "  -> elevation SURVIVES far from any cDC: a real regional field, not spillover.\n" else
      "  -> elevation does NOT survive far from cDC: consistent with segmentation spillover.\n")
cat("wrote ", file.path(OUT,"batf3_spillover.csv"), "\n", sep="")
