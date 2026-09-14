## Does the over-fixation FLATTENING confound the pre-QC chemokine result?
##
## THE PROBLEM WITH "0.97x". scripts/sbatch/qc_11_preqc_sensitivity.R reported
## CCL17+CCL22 at 9.98 cp10k in `qc_spatial_staining == "ok"` cells vs 10.27 in the 2.0M
## `dim_membrane` cells, and read that near-parity as "the spatial QC filter does not
## select the signal". Two things are wrong with resting on it:
##
##  1. It is a POINT ESTIMATE offered as evidence of ABSENCE, with no equivalence bound
##     and no uncertainty. At n = 527k vs 2.01M the SE is tiny, so 0.97 is likely
##     distinguishable from 1.0 -- which is not the same as it mattering.
##  2. C2 of the over-fixation diagnosis shows defect cells have profiles dragged toward
##     uniform: log2(detection ratio) falls with log2(baseline prevalence) at slope
##     -0.30. Every gene therefore has an EXPECTED ratio given its prevalence. A raw
##     ratio near 1 could be a real change cancelled by the flattening, or no change at
##     all -- the raw number cannot tell them apart.
##
## THE TEST. Reproduce C2's regression on IB1's own cells (dim_membrane vs ok, depth
## matched), then ask where CCL17 and CCL22 sit RELATIVE TO THE FITTED LINE. The
## residual is the part not explained by flattening.
##   residual ~ 0  -> the chemokines behave exactly as flattening predicts; the raw
##                    0.97x carries no information about selection either way
##   residual > 0  -> enriched in defect cells beyond flattening
##   residual < 0  -> depressed beyond flattening; pre-QC magnitudes are inflated
##
## Also reports where CCL17 sits in the panel's prevalence distribution, because the
## direction of the flattening effect depends on that and it was previously asserted
## rather than checked.
##
##   sbatch --mem=128G --cpus-per-task=2 --time=2:00:00 \
##     scripts/sbatch/run_rscript.sh scripts/sbatch/qc_11_flattening_vs_chemokine.R
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(data.table); library(ggplot2); setDTthreads(1)
})

IN_RDS <- "output/01_import_v1/seuratObject_IB1_BPCells_full.Rds"
## The step-02 output directory was renamed 02_filter_QC -> 02_spatial_QC partway
## through this project. Resolve whichever exists rather than hardcoding either, so this
## script works before, during and after that rename.
.qc_dir <- local({
  cand <- c("output/02_spatial_QC", "output/02_filter_QC")
  hit  <- cand[dir.exists(cand)]
  if (!length(hit)) stop("no step-02 QC directory found; looked for: ",
                         paste(cand, collapse=", "), call. = FALSE)
  hit[1]
})
QC_CSV <- file.path(.qc_dir, "qc_spatial_cells.csv.gz")
OUT    <- "output/11_hodgkin_niche"; FIG <- file.path(OUT, "figures")
CHEMO  <- c("CCL17","CCL22"); MARK <- c("TNFRSF8","MS4A1","CD3E","CD68","LAMP3","BATF3")
N_SUB  <- 60000L        # cells per group after depth matching
MIN_DET<- 0.002         # genes below this baseline prevalence are unstable; excluded

SURF<-"#fcfcfb"; INK<-"#0b0b0b"; INK2<-"#52514e"; MUTED<-"#898781"; GRID<-"#e1e0d9"
BLUE<-"#2a78d6"; ORANGE<-"#eb6834"
theme_viz <- function(base=11) theme_minimal(base_size=base) +
  theme(plot.background=element_rect(fill=SURF,colour=NA),
        panel.background=element_rect(fill=SURF,colour=NA),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_line(colour=GRID,linewidth=.25),
        plot.title=element_text(colour=INK,face="bold",size=base+2),
        plot.subtitle=element_text(colour=INK2,size=base-1),
        plot.caption=element_text(colour=MUTED,size=base-3,hjust=0),
        legend.position="top", legend.justification="left")
sf <- function(n,p,w,h){ggsave(file.path(FIG,n),p,width=w,height=h,dpi=150,bg=SURF)
                       cat("wrote",file.path(FIG,n),"\n")}
sc <- function(n,d){fwrite(d,file.path(OUT,n));cat("wrote",file.path(OUT,n),"\n")}

cat("============================================================\n")
cat("Is the pre-QC chemokine parity explained by over-fixation flattening?\n")
cat("============================================================\n")

obj <- readRDS(IN_RDS); rn <- rownames(obj@meta.data)
md <- as.data.table(obj@meta.data[, c("Run_Tissue_name","nCount_RNA"), drop = FALSE])
md[, cell_id := rn]
qc <- fread(QC_CSV, select = c("cell","sample_id","qc_spatial_staining"))
md <- merge(md, qc[sample_id == "IB1", .(cell_id = cell, stain = qc_spatial_staining)],
            by = "cell_id", all.x = FALSE, sort = FALSE)
md <- md[Run_Tissue_name == "IB1" & nCount_RNA > 0 & stain %in% c("ok","dim_membrane")]
cat("\ncells by staining flag:\n"); print(md[, .N, by = stain])

## ---- depth-matched subsample -------------------------------------------------
## Flattening and depth are entangled, so the two groups must be compared at matched
## library size or the slope below is partly a depth slope.
md[, dec := cut(nCount_RNA, unique(quantile(nCount_RNA, 0:20/20)),
                include.lowest = TRUE, labels = FALSE)]
set.seed(1)
per_dec <- ceiling(N_SUB / uniqueN(md$dec))
sel <- md[, {
  n <- min(.N_ok <- sum(stain == "ok"), .N_dm <- sum(stain == "dim_membrane"), per_dec)
  if (n < 20) list(i = integer(0)) else
    list(i = c(.I[stain == "ok"][sample.int(.N_ok, n)],
               .I[stain == "dim_membrane"][sample.int(.N_dm, n)]))
}, by = dec]$i
sub <- md[sel]
cat(sprintf("\ndepth-matched subsample: %d cells (%s)\n", nrow(sub),
            paste(sub[, .N, by=stain][, paste0(stain,"=",N)], collapse=" ")))
cat(sprintf("median depth: ok %.0f, dim_membrane %.0f\n",
            sub[stain=="ok", median(nCount_RNA)], sub[stain=="dim_membrane", median(nCount_RNA)]))

cnt <- LayerData(obj[["RNA"]], layer = "counts")
idx <- match(sub$cell_id, colnames(cnt)); stopifnot(!anyNA(idx))
M <- as(cnt[, idx, drop = FALSE], "dgCMatrix")
rm(obj, cnt); invisible(gc())
is_ok <- sub$stain == "ok"

## ---- per-gene detection rate in each group ------------------------------------
det_ok <- Matrix::rowSums(M[, is_ok,  drop=FALSE] > 0) / sum(is_ok)
det_dm <- Matrix::rowSums(M[, !is_ok, drop=FALSE] > 0) / sum(!is_ok)
G <- data.table(gene = rownames(M), det_ok = as.numeric(det_ok), det_dm = as.numeric(det_dm))
G <- G[det_ok >= MIN_DET & det_dm > 0]
G[, l2ratio := log2(det_dm / det_ok)]      # dim_membrane relative to ok
G[, l2prev  := log2(det_ok)]               # baseline prevalence, from the CLEAN cells
cat(sprintf("\ngenes retained (baseline detection >= %.1f%%): %d of %d\n",
            100*MIN_DET, nrow(G), nrow(M)))

## ---- reproduce C2's flattening slope -------------------------------------------
fit <- lm(l2ratio ~ l2prev, data = G)
sl  <- coef(fit)[2]; pv <- summary(fit)$coefficients[2,4]
cat(sprintf("\n=== flattening slope on IB1 (C2 reported -0.30)\n  slope %.3f (p = %.3g), R2 = %.3f\n",
            sl, pv, summary(fit)$r.squared))
G[, expected := predict(fit, G)]
G[, residual := l2ratio - expected]

## ---- where do the chemokines sit? ----------------------------------------------
G[, prev_pctile := 100 * frank(det_ok) / .N]
foc <- G[gene %in% c(CHEMO, MARK)][order(-det_ok)]
cat("\n=== the genes of interest, against what flattening predicts\n")
print(foc[, .(gene, det_ok = round(100*det_ok,2), det_dm = round(100*det_dm,2),
              prev_pctile = round(prev_pctile), obs_l2 = round(l2ratio,3),
              exp_l2 = round(expected,3), residual = round(residual,3),
              obs_ratio = round(2^l2ratio,3), resid_ratio = round(2^residual,3))])
sc("flattening_vs_chemokine.csv", foc)

rsd <- sd(G$residual)
cat(sprintf("\nresidual SD across %d genes: %.3f log2 units\n", nrow(G), rsd))
for (g in CHEMO) {
  r <- G[gene == g, residual]
  if (length(r) == 1)
    cat(sprintf("  %-6s residual %+.3f log2 (%.2f SD from the flattening expectation)\n",
                g, r, r/rsd))
}

## ---- verdict --------------------------------------------------------------------
cc <- G[gene %in% CHEMO]
cat("\n============================================================\n")
if (nrow(cc) == length(CHEMO)) {
  mz <- mean(cc$residual) / rsd
  cat(sprintf("CCL17/CCL22 mean residual %+.3f log2 (%.2f SD)\n", mean(cc$residual), mz))
  cat(sprintf("CCL17 sits at the %.0fth prevalence percentile of the panel -- so flattening\n  is expected to %s its detection in defect cells.\n",
              G[gene=="CCL17", prev_pctile],
              if (G[gene=="CCL17", prev_pctile] > 50) "DEPRESS" else "INFLATE"))
  if (abs(mz) < 0.5)
    cat("VERDICT: the chemokines behave as flattening predicts. The raw 0.97x is\n",
        "         consistent with no selection effect AND with the flattening model --\n",
        "         it does not by itself establish either. Treat the pre-QC magnitudes as\n",
        "         provisional but not demonstrably biased.\n", sep="")
  else if (mz > 0)
    cat("VERDICT: the chemokines are ENRICHED in defect cells beyond flattening. Pre-QC\n",
        "         bin magnitudes in defect-heavy regions are INFLATED and must be\n",
        "         corrected or reported as upper bounds.\n", sep="")
  else
    cat("VERDICT: the chemokines are DEPRESSED beyond flattening. The post-QC estimate is\n",
        "         conservative and pre-QC coverage would if anything understate the field.\n", sep="")
}

p1 <- ggplot(G, aes(l2prev, l2ratio)) +
  geom_point(colour = "#c9c8c2", size = .5, stroke = 0) +
  geom_smooth(method = "lm", se = TRUE, colour = INK, linewidth = .6) +
  geom_point(data = foc, colour = ORANGE, size = 2.4) +
  ggrepel::geom_text_repel(data = foc, aes(label = gene), size = 3, colour = INK2,
                           min.segment.length = 0, seed = 1) +
  labs(title = "Flattening explains most of what the raw ratio shows",
       subtitle = sprintf(paste0("Each point is a gene. Line = the over-fixation flattening ",
         "relationship (slope %.2f); C2 reported -0.30 independently.\nDistance from the ",
         "line is the part NOT explained by flattening."), sl),
       x = "log2 baseline detection rate in clean (`ok`) cells",
       y = "log2( detection in dim_membrane / detection in ok )",
       caption = sprintf("Data: %s/flattening_vs_chemokine.csv", OUT)) + theme_viz()
if (requireNamespace("ggrepel", quietly = TRUE)) sf("fig_11_flattening_vs_chemokine.png", p1, 8.5, 6) else
  sf("fig_11_flattening_vs_chemokine.png", p1 + geom_text(data=foc, aes(label=gene),
     nudge_y=.15, size=3, colour=INK2), 8.5, 6)
cat("============================================================\nOK\n")
