## ============================================================================
## qc_16_negprobe_floor.R -- is a checkpoint gene's detection above background at all?
##
## Every checkpoint gene on this panel detects in 2.5-10% of cells with a mean of 0.03-0.15
## counts (output/16_regions/panel_detection_preqc.csv). Step 10 killed per-cell HRS calling
## for exactly this shape of number -- TNFRSF8 maxing at 6 counts, 78% singletons, depth
## predicting detection at AUC 0.716. So before step 16e builds anything on "PDCD1-positive"
## or "CD274-positive" cells, this asks the only question that matters: does the gene detect
## MORE OFTEN than a negative control probe on the same cells?
##
## The 20 negprobes are the panel's own blank controls. Their per-probe detection rate is the
## false-positive floor at a given count threshold. A gene at the floor is not measurable per
## cell, whatever its biology.
##
## Usage: sbatch --mem=32G scripts/sbatch/run_rscript.sh scripts/sbatch/qc_16_negprobe_floor.R
## ============================================================================
suppressPackageStartupMessages({ library(Matrix); library(data.table); setDTthreads(1) })

M  <- readRDS("output/16_regions/panel_counts_all.Rds")
NP <- readRDS("output/16_regions/negprobe_counts_all.Rds")
md <- fread("output/16_regions/cells_meta_all.csv.gz",
            select = c("cell_id","samp","tier_strict","is_defect","nCount_RNA"))
stopifnot(identical(colnames(M), colnames(NP)))
ii <- match(md$cell_id, colnames(M)); stopifnot(!anyNA(ii))
M <- M[, ii, drop=FALSE]; NP <- NP[, ii, drop=FALSE]
OUT <- "output/16_checkpoint"; dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

GOI <- c("PDCD1","CD274","PDCD1LG2","CTLA4","LAG3","TIGIT","HAVCR2","TOX","IDO1","VSIR",
         "BTLA","ENTPD1","TNFRSF9","LGALS9","CD80","CD86","ICOS","CXCL13","FOXP3","CCR4",
         "CCL17","CCL22","TNFRSF8","CD68","CD163","MS4A1","CD3E","BATF3","IGKC")
GOI <- intersect(GOI, rownames(M))
THR <- c(1L, 2L, 3L)

res <- rbindlist(lapply(c("strict","all"), function(tier) {
  keep <- if (tier == "strict") md$tier_strict & !md$is_defect else rep(TRUE, nrow(md))
  rbindlist(lapply(c("IB1","RB1"), function(S) {
    k <- which(keep & md$samp == S); if (!length(k)) return(NULL)
    rbindlist(lapply(THR, function(t) {
      ## negprobe floor: detection rate of each of the 20 blank probes at this threshold
      nd <- as.numeric(Matrix::rowMeans(NP[, k, drop=FALSE] >= t))
      gd <- as.numeric(Matrix::rowMeans(M[GOI, k, drop=FALSE] >= t))
      data.table(tier = tier, samp = S, thr = t, gene = GOI,
                 pct_det = round(100*gd, 4),
                 neg_mean = round(100*mean(nd), 4),
                 neg_max  = round(100*max(nd), 4),
                 ## how many times the mean blank probe, and does it clear the WORST blank
                 ratio_vs_neg_mean = round(gd / pmax(mean(nd), 1e-9), 2),
                 above_worst_neg = gd > max(nd))
    }))
  }))
}))
fwrite(res, file.path(OUT, "negprobe_floor.csv"))
cat("wrote ", file.path(OUT, "negprobe_floor.csv"), "\n", sep="")

cat("\n=== detection vs the negative-probe floor, strict tier, threshold >= 1 count\n")
print(res[tier=="strict" & thr==1L][order(samp, -ratio_vs_neg_mean),
      .(samp, gene, pct_det, neg_mean, neg_max, ratio_vs_neg_mean, above_worst_neg)])
cat("\n=== threshold >= 2 counts (strict)\n")
print(res[tier=="strict" & thr==2L][order(samp, -ratio_vs_neg_mean),
      .(samp, gene, pct_det, neg_mean, ratio_vs_neg_mean, above_worst_neg)])
## ---- the recommended threshold per gene -------------------------------------
## A first version judged usability at >= 1 count only and concluded that CD274 (2.7x the
## blank rate) and PDCD1 (2.9x) were unusable per cell. That was the wrong reading. Raising
## the bar to >= 2 counts drops the blank rate to ~0.05% and the same genes become 25-31x
## the floor -- specific, just sparse (about 1% of cells). So the question is not "usable or
## not" but AT WHICH THRESHOLD, and that is what this table reports.
SPEC_MIN <- 10   # a call must detect >= 10x the mean blank probe to be worth making
w <- dcast(res[tier=="strict"], gene + thr ~ samp,
           value.var = c("pct_det","ratio_vs_neg_mean","above_worst_neg"))
w[, specific := above_worst_neg_IB1 & above_worst_neg_RB1 &
                ratio_vs_neg_mean_IB1 >= SPEC_MIN & ratio_vs_neg_mean_RB1 >= SPEC_MIN]
v <- w[specific == TRUE, .(rec_thr = min(thr)), by = gene]
v <- merge(w, v, by = "gene", all.x = TRUE)[thr == 1L]
v[, rec_thr := fifelse(is.na(rec_thr), NA_integer_, rec_thr)]
v <- merge(v[, .(gene, pct_det_thr1_IB1 = pct_det_IB1, pct_det_thr1_RB1 = pct_det_RB1,
                 ratio_thr1_IB1 = ratio_vs_neg_mean_IB1, rec_thr)],
           w[, .(gene, thr, pct_det_IB1, pct_det_RB1, ratio_vs_neg_mean_IB1,
                 ratio_vs_neg_mean_RB1)],
           by.x = c("gene","rec_thr"), by.y = c("gene","thr"), all.x = TRUE)
setnames(v, c("pct_det_IB1","pct_det_RB1","ratio_vs_neg_mean_IB1","ratio_vs_neg_mean_RB1"),
         c("pct_det_at_rec_IB1","pct_det_at_rec_RB1","ratio_at_rec_IB1","ratio_at_rec_RB1"))
cat(sprintf("\n=== recommended per-cell threshold (>= %gx the mean blank probe, both samples)\n",
            SPEC_MIN))
print(v[order(rec_thr, -ratio_at_rec_IB1)])
fwrite(v, file.path(OUT, "negprobe_verdict.csv"))
cat("\nwrote ", file.path(OUT, "negprobe_verdict.csv"), "\n", sep="")
for (t in sort(unique(na.omit(v$rec_thr))))
  cat(sprintf("  threshold >= %d counts: %s\n", t,
              paste(v[rec_thr == t, gene], collapse = ", ")))
if (v[is.na(rec_thr), .N])
  cat("  NEVER specific at any threshold -- region-level only: ",
      paste(v[is.na(rec_thr), gene], collapse=", "), "\n", sep="")

## ---- does DEPTH predict the call? -------------------------------------------
## Step 10 killed per-cell CD30 partly because sequencing depth predicted TNFRSF8 detection
## at AUC 0.716. A sparse >= 2-count call is even more exposed to that: a deeper cell simply
## has more chances to reach 2. Any gene whose call is strongly depth-predicted cannot be
## compared between regions of unequal depth without matching.
## as.numeric() on the counts is load-bearing. sum(<logical>) returns an INTEGER, and both
## n1*(n1+1) and n1*n0 overflow .Machine$integer.max once a group passes ~46,000 cells --
## R then returns NA rather than erroring. The first run of this script silently produced NA
## for every RB1 gene (n1 up to 860,399) and for IB1 at threshold 1, which looked like
## missing data rather than arithmetic. Doubles have no such limit here.
auc_of <- function(pos, v) {           # rank-based AUC, cf. 10_script_hrs_call.R:111
  r <- rank(v)
  n1 <- as.numeric(sum(pos)); n0 <- as.numeric(length(pos)) - n1
  if (n1 == 0 || n0 == 0) return(NA_real_)
  (sum(r[pos]) - n1*(n1+1)/2) / (n1*n0)
}
dep_chk <- rbindlist(lapply(c("IB1","RB1"), function(S) {
  k <- which(md$tier_strict & !md$is_defect & md$samp == S)
  dep <- md$nCount_RNA[k]
  rbindlist(lapply(GOI, function(g) rbindlist(lapply(THR, function(t) {
    pos <- as.vector(M[g, k] >= t)
    data.table(samp = S, gene = g, thr = t, n_pos = sum(pos),
               auc_depth = round(auc_of(pos, dep), 4))
  }))))
}))
bad_na <- dep_chk[is.na(auc_depth) & n_pos > 0]
if (nrow(bad_na)) stop(nrow(bad_na), " AUCs are NA despite n_pos > 0 -- arithmetic problem",
                       call. = FALSE)
fwrite(dep_chk, file.path(OUT, "depth_auc.csv"))
cat("\n=== does depth predict positivity? (AUC 0.5 = no, step 10's CD30 was 0.716)\n")
print(dcast(dep_chk[gene %in% c("PDCD1","CD274","PDCD1LG2","CTLA4","CXCL13","CD3E","CD68",
                                "TNFRSF8","IGKC")],
            gene + thr ~ samp, value.var = "auc_depth")[order(gene, thr)])
cat("\nworst depth-AUC at the recommended threshold, across all genes tested:\n")
mm <- merge(dep_chk, v[, .(gene, rec_thr)], by = "gene")[thr == rec_thr]
print(mm[order(-auc_depth)][1:10, .(samp, gene, thr, n_pos, auc_depth)])
cat("wrote ", file.path(OUT, "depth_auc.csv"), "\n", sep="")
