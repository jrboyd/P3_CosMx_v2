## ============================================================================
## 10_script_hrs_call.R -- identify HRS cells in IB1, with RB1 as the negative control.
##
## WHY A REBUILD. scripts/sbatch/label_hrs_candidates.R produced 4,089 IB1 candidates
## (2,657 tier-A) and fails four independent checks:
##   * RB1 -- the FL recurrence, NO Hodgkin component, so ~zero HRS expected -- yields
##     6,761 candidates (1,899 tier-A). IB1:RB1 tier-A rate ratio only 4.89x.
##   * tier-A median Area.um2 104.8 vs 63.9 background = 1.6x. HRS cells are giant and
##     multinucleated, 5-20x a lymphocyte. Nothing giant is being selected.
##   * tier-A REQUIRES cd30_pos, then reports TNFRSF8 at 14.29x. That is circular. The
##     non-circular corroborating markers run 1.85-2.17x, and the 0.00 ratios for the
##     B/T markers are the gate itself.
##   * HRS cells are CD45-NEGATIVE. PTPRC sits in the score's negative set but not in
##     the tier-A gate, and comes back 1.37x ENRICHED.
## Cross-checked against step 09: 77.6% of the old tier-A sits in pearson_clusters_low
## cluster 1, which HierAnnot annotates "Myeloid > Macrophage". The old call is largely
## selecting macrophages -- consistent with CD30 being expressed on activated
## lymphocytes and some myeloid cells.
##
## WHAT THIS DOES DIFFERENTLY.
##   1. TNFRSF8 is REMOVED from the score and held out as the independent validator.
##   2. Adjacent CD30+ segments are merged into objects (section A2). A giant
##      multinucleated cell that the segmenter split into fragments looks like several
##      normal-sized cells -- which is exactly the 1.6x area ratio observed. The unit of
##      analysis has to be the object, and the merge is tested against a
##      density-matched permutation null so "adjacent by chance" is excluded.
##   3. The threshold is calibrated on RB1 rather than fixed at a top-1% quantile.
##   4. Morphology is SCORED, never a hard gate -- a size filter would silently convert
##      a segmentation failure on a giant cell into a true negative.
##
## Usage: Rscript scripts/10_script_hrs_call.R
##        scripts/10_submit_hrs_call.sh
## ============================================================================

suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(RANN); library(igraph)
  library(ggplot2); library(data.table); setDTthreads(1)
})

IN_RDS  <- "output/04_merge/merged_processed.Rds"
HIER    <- "output/09_hierannot/IB1/hierannot_cluster_labels.csv"
OUT     <- "output/10_hrs"
FIG     <- file.path(OUT, "figures")

SAMPLES    <- c("IB1", "RB1")     # IB1 = composite cHL+FL; RB1 = FL recurrence, the control
SAMPLE_COL <- "merge_sample"
DEFECT_CL  <- "i"                 # over-fixation defect population
CD30       <- "TNFRSF8"
MERGE_UM   <- 12                  # centroid-to-centroid segment adjacency; swept below
MERGE_SWEEP<- c(8, 10, 12, 15, 20)
K_MERGE    <- 20                  # candidate neighbours before the radius filter
N_PERM     <- 100                 # permutations for the adjacency null
NN_UM      <- 30                  # rosette radius
K_NN       <- 80
MIN_NBR    <- 5
W_MORPH    <- 1.0                 # weight on the morphology rank; swept in the report
OLD_RATIO  <- 4.89                # the IB1:RB1 tier-A rate ratio this has to beat

if (!file.exists(IN_RDS)) stop("missing ", IN_RDS, call. = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

## ---- marker sets -----------------------------------------------------------
## All verified present in the IB1 BPCells row_names before writing this. Only
## classical HLA class I (HLA-A/B/C) is absent by name from the 6K panel.
##
## TNFRSF8 IS NOT IN ANY SCORED SET. It is the held-out validator: if the score is
## finding HRS cells, CD30 detection must rise with it, and that test is only
## meaningful because the score cannot see CD30.
POS <- c("CCL17","CCL22","IRF4","TRAF1","CD274","PDCD1LG2","FUT4","BIRC3","NFKB2",
         "RELB","CD40","TNFAIP3","GATA3","ID2","LGALS1","IL13","IL21","BATF3",
         "STAT6","JAK2","SOCS1")
## HRS cells shed the B programme and are not T cells; classic HRS are CD45-negative,
## so PTPRC belongs here and -- unlike the old script -- its behaviour is reported.
NEG <- c("MS4A1","CD19","CD79A","CD79B","CD22","FCER2","PAX5","POU2AF1",
         "CD3D","CD3E","CD3G","TRAC","CD2","PTPRC")
## Antigen-presentation loss is a cHL hallmark and was unused by the old script.
## Negative-signed: low APL expression supports HRS.
APL <- c("B2M","TAP1","TAP2","NLRC5","CIITA","HLA-DRA","CD58")
## Rosette / microenvironment readout, not scored.
TCELL <- c("CD3D","CD3E","CD3G","TRAC","CD2","CD4","PDCD1")
GENES <- unique(c(CD30, POS, NEG, APL, TCELL, "CD68","CD163"))

SURF<-"#fcfcfb"; INK<-"#0b0b0b"; INK2<-"#52514e"; MUTED<-"#898781"; GRID<-"#e1e0d9"; AXIS<-"#c3c2b7"
DIV_HI<-"#184f95"; DIV_LO<-"#d03b3b"
SAMP_COLS <- c(IB1 = "#4e79a7", RB1 = "#e15759")
theme_viz <- function(base=11, xgrid=FALSE) {
  th <- theme_minimal(base_size=base) +
    theme(plot.background=element_rect(fill=SURF,colour=NA),
          panel.background=element_rect(fill=SURF,colour=NA),
          panel.grid.major.y=element_line(colour=GRID,linewidth=.25),
          panel.grid.major.x=element_blank(), panel.grid.minor=element_blank(),
          axis.line.x=element_line(colour=AXIS,linewidth=.3), axis.ticks=element_blank(),
          axis.text=element_text(colour=MUTED,size=base-2),
          axis.title=element_text(colour=INK2,size=base-1),
          plot.title=element_text(colour=INK,face="bold",size=base+2),
          plot.subtitle=element_text(colour=INK2,size=base-1),
          plot.caption=element_text(colour=MUTED,size=base-3,hjust=0),
          strip.text=element_text(colour=INK2,face="bold",size=base-1),
          legend.text=element_text(colour=INK2,size=base-2),
          legend.title=element_text(colour=INK2,size=base-2),
          legend.position="top", legend.justification="left", legend.key.height=unit(9,"pt"))
  if (xgrid) th <- th + theme(panel.grid.major.x=element_line(colour=GRID,linewidth=.25))
  th
}
sf <- function(n,p,w,h){ggsave(file.path(FIG,n),p,width=w,height=h,dpi=150,bg=SURF)
                       cat("wrote",file.path(FIG,n),"\n")}
sc <- function(n,d){fwrite(d,file.path(OUT,n));cat("wrote",file.path(OUT,n),"\n")}
## rank AUC, as scripts/07_script_plasma_marker_check.R:126
## as.numeric on the counts: at this scale n1*n0 overflows integer arithmetic
## (17,650 x 234,097 = 4.1e9 > .Machine$integer.max) and silently returns NA.
auc_of <- function(v, lab) {
  n1 <- as.numeric(sum(lab)); n0 <- as.numeric(sum(!lab))
  if (!n1 || !n0) return(NA_real_)
  r <- rank(v); (sum(r[lab]) - n1*(n1+1)/2) / (n1*n0)
}

cat("============================================================\n")
cat("step 10 -- HRS call, IB1 with RB1 as the negative control\n")
cat("============================================================\n")

## ---- A1. object, markers, within-sample z-scores ---------------------------
cat("\nreading ", IN_RDS, "\n", sep = "")
obj <- readRDS(IN_RDS)
rn  <- rownames(obj@meta.data)
NEED <- c(SAMPLE_COL, "Area.um2", "NucArea", "nCount_RNA", "inSitu_clusters",
          "lineage", "x_slide_mm", "y_slide_mm",
          "pearson_clusters_low", "pearson_clusters_high")
miss <- setdiff(NEED, colnames(obj@meta.data))
if (length(miss)) stop("meta.data is missing: ", paste(miss, collapse = ", "), call. = FALSE)
md <- as.data.table(obj@meta.data[, NEED, drop = FALSE]); md[, cell_id := rn]

cnt <- LayerData(obj[["RNA"]], layer = "counts")
stopifnot(identical(colnames(cnt), rn))
HAVE <- intersect(GENES, rownames(cnt))
if (length(HAVE) < length(GENES))
  message("absent from the panel and dropped: ", paste(setdiff(GENES, HAVE), collapse = ", "))
M <- as(cnt[HAVE, , drop = FALSE], "dgCMatrix"); colnames(M) <- rn
rm(obj, cnt); invisible(gc())
POS <- intersect(POS, HAVE); NEG <- intersect(NEG, HAVE); APL <- intersect(APL, HAVE)
TCELL <- intersect(TCELL, HAVE)
if (!CD30 %in% HAVE) stop("TNFRSF8 absent -- the held-out validator is unavailable", call.=FALSE)

md <- md[get(SAMPLE_COL) %in% SAMPLES][!is.na(x_slide_mm) & !is.na(y_slide_mm) & nCount_RNA > 0]
md[, samp := factor(get(SAMPLE_COL), levels = SAMPLES)]
cat("\ncells with usable coordinates and counts:\n"); print(md[, .N, by = samp])

## Exclude the over-fixation defect BEFORE scoring. It is 37% of IB1 and was 37.7% of
## the old tier-A -- identical to background, so it contributed noise and no signal.
md <- md[inSitu_clusters != DEFECT_CL]
cat(sprintf("after excluding cluster \"%s\": %d cells\n", DEFECT_CL, nrow(md)))
print(md[, .N, by = samp])

gidx <- match(md$cell_id, colnames(M)); stopifnot(!anyNA(gidx))
Msub <- M[, gidx, drop = FALSE]
md[, cd30_ct := as.vector(Msub[CD30, ])]
md[, cd30_pos := cd30_ct > 0]

## log-normalise from raw counts, then z-score WITHIN SAMPLE. Rationale carried over
## from scripts/sbatch/label_hrs_candidates.R:83-98: IB1 and RB1 differ ~4x in
## detection, so a shared z would re-encode sample identity as signal.
tc <- pmax(md$nCount_RNA, 1)
## Only ~50 marker rows, so densify once and use sweep() over cells. Dividing a
## dgCMatrix by a vector silently densifies anyway, and doing it explicitly keeps the
## allocation count (and the peak) predictable.
ln <- log1p(sweep(as.matrix(Msub), 2, tc / median(tc), "/"))
zmat <- matrix(0, nrow = nrow(ln), ncol = ncol(ln), dimnames = list(rownames(ln), NULL))
for (s in levels(md$samp)) {
  i  <- which(md$samp == s)
  sub <- as.matrix(ln[, i, drop = FALSE])
  mu  <- rowMeans(sub); sdv <- apply(sub, 1, stats::sd); sdv[sdv == 0 | !is.finite(sdv)] <- 1
  zmat[, i] <- (sub - mu) / sdv
}
rm(ln); invisible(gc())

## morphology rank, within sample, on the geometric mean of area and nuclear area
md[, morph_rank := {
  v <- sqrt(pmax(Area.um2, 0) * pmax(NucArea, 0))
  frank(v, ties.method = "average") / .N }, by = samp]
md[, `:=`(z_pos = colMeans(zmat[POS, , drop = FALSE]),
          z_neg = colMeans(zmat[NEG, , drop = FALSE]),
          z_apl = colMeans(zmat[APL, , drop = FALSE]))]
md[, hrs_score := z_pos - z_neg - z_apl + W_MORPH * (morph_rank - 0.5)]
## Split by CD30 status: the z-scores are mean-0 within sample by construction, so a
## per-sample mean would print 0.000 for every component and tell us nothing. The
## informative comparison is CD30+ vs the rest -- and CD30 is not in the score, so any
## separation here is already the held-out signal.
cat("\n=== score components, CD30+ vs rest (z-scores are mean-0 within sample)\n")
print(md[, .(cells = .N, z_pos = round(mean(z_pos),3), z_neg = round(mean(z_neg),3),
             z_apl = round(mean(z_apl),3), morph = round(mean(morph_rank),3),
             score = round(mean(hrs_score),3)), by = .(samp, cd30_pos)][order(samp, -cd30_pos)])
cat(sprintf("CD30+ prevalence: IB1 %.2f%%, RB1 %.2f%%\n",
            100*md[samp=="IB1", mean(cd30_pos)], 100*md[samp=="RB1", mean(cd30_pos)]))

## ---- A2. merge adjacent CD30+ segments into objects ------------------------
## Connected components, NOT communities: the claim is physical contiguity.
merge_objects <- function(d, radius_um) {
  if (!nrow(d)) return(integer(0))
  xy <- as.matrix(d[, .(x_slide_mm, y_slide_mm)])
  k  <- min(K_MERGE, nrow(xy))
  nn <- RANN::nn2(xy, k = k)
  i  <- rep(seq_len(nrow(xy)), k - 1L)
  j  <- as.vector(nn$nn.idx[, -1, drop = FALSE])
  dd <- as.vector(nn$nn.dists[, -1, drop = FALSE]) * 1000
  keep <- dd <= radius_um & i != j
  gr <- igraph::make_empty_graph(n = nrow(xy), directed = FALSE)
  if (any(keep)) gr <- igraph::add_edges(gr, as.vector(rbind(i[keep], j[keep])))
  igraph::components(gr)$membership
}

cd30 <- md[cd30_pos == TRUE]
cat(sprintf("\nCD30+ segments: %d (IB1 %d, RB1 %d)\n", nrow(cd30),
            cd30[samp=="IB1", .N], cd30[samp=="RB1", .N]))

## radius sweep, reported so MERGE_UM is visibly not a cliff
sw <- rbindlist(lapply(MERGE_SWEEP, function(r) rbindlist(lapply(SAMPLES, function(s) {
  d <- cd30[samp == s]; m <- merge_objects(d, r)
  tb <- table(table(m))
  data.table(radius_um = r, samp = s, n_segments = nrow(d), n_objects = uniqueN(m),
             pct_multi = round(100 * sum(tb[names(tb) != "1"]) / uniqueN(m), 1),
             max_size = if (length(m)) max(table(m)) else 0L)
}))))
cat("\n=== merge-radius sweep\n"); print(sw)
sc("hrs_merge_radius_sweep.csv", sw)

## the chosen radius
cd30[, obj_id := paste0(samp, "_", merge_objects(.SD, MERGE_UM)), by = samp,
     .SDcols = c("x_slide_mm","y_slide_mm")]
md <- merge(md, cd30[, .(cell_id, obj_id)], by = "cell_id", all.x = TRUE, sort = FALSE)

## ADJACENCY NULL. CD30+ cells touch by chance wherever cells are dense, so the
## observed component-size distribution means nothing on its own. Permute the CD30+
## label among cells matched on LOCAL DENSITY (neighbour count within MERGE_UM) and
## recompute. What matters is the excess of multi-segment objects over that null.
cat(sprintf("\n=== adjacency null, %d permutations, density-matched\n", N_PERM))
null_rows <- list()
for (s in SAMPLES) {
  d  <- copy(md[samp == s])
  xy <- as.matrix(d[, .(x_slide_mm, y_slide_mm)])
  nn <- RANN::nn2(xy, k = min(K_MERGE, nrow(xy)))
  d[, dens := rowSums(nn$nn.dists * 1000 <= MERGE_UM) - 1L]   # -1 drops the self match
  d[, dens_bin := cut(dens, unique(quantile(dens, 0:10/10)), include.lowest = TRUE,
                      labels = FALSE)]
  obs <- table(table(d[cd30_pos == TRUE]$obj_id))
  obs_multi <- sum(obs[names(obs) != "1"])
  set.seed(1)
  ## Draw the same number of pseudo-CD30+ cells WITHIN each density bin, so the
  ## permuted set has the observed spatial density profile and only the labels move.
  ## `.I[sample.int(.N, k)]` rather than `sample(.I, k)`: sample() on a length-1
  ## vector silently means 1:n.
  perm_multi <- vapply(seq_len(N_PERM), function(p) {
    sel <- d[, .(i = .I[sample.int(.N, sum(cd30_pos))]), by = dens_bin]$i
    if (length(sel) < 2) return(0)
    m <- merge_objects(d[sel], MERGE_UM); tb <- table(table(m))
    sum(tb[names(tb) != "1"])
  }, numeric(1))
  null_rows[[s]] <- data.table(
    samp = s, cd30_segments = d[, sum(cd30_pos)],
    obs_objects = uniqueN(d[cd30_pos == TRUE]$obj_id), obs_multi_objects = obs_multi,
    null_multi_mean = round(mean(perm_multi), 1), null_multi_sd = round(sd(perm_multi), 1),
    enrichment = round(obs_multi / max(mean(perm_multi), 1e-9), 2),
    p_emp = (1 + sum(perm_multi >= obs_multi)) / (N_PERM + 1))
}
NULLT <- rbindlist(null_rows)
print(NULLT); sc("hrs_merge_null.csv", NULLT)

## ---- A3. object-level table -------------------------------------------------
OBJ <- md[!is.na(obj_id), .(
  samp = first(samp), n_segments = .N,
  area_sum = sum(Area.um2), nuc_sum = sum(NucArea),
  nCount_sum = sum(nCount_RNA), cd30_sum = sum(cd30_ct),
  hrs_score = mean(hrs_score), z_pos = mean(z_pos), z_neg = mean(z_neg),
  z_apl = mean(z_apl), morph_rank = mean(morph_rank),
  x_slide_mm = mean(x_slide_mm), y_slide_mm = mean(y_slide_mm),
  ## dist() is O(n^2); a runaway component would hang here, and its span is not
  ## informative anyway -- such an object is a dense patch, not a cell
  span_um = if (.N > 1 && .N <= 200) 1000*max(dist(cbind(x_slide_mm, y_slide_mm)))
            else if (.N > 200) NA_real_ else 0
), by = obj_id]
## the object score should reflect pooled size, not the mean of its fragments
OBJ[, area_rank := frank(area_sum, ties.method = "average")/.N, by = samp]
OBJ[, obj_score := z_pos - z_neg - z_apl + W_MORPH * (area_rank - 0.5)]
cat("\n=== objects by sample\n")
print(OBJ[, .(objects = .N, multi = sum(n_segments > 1),
              med_segments = as.numeric(median(n_segments)),
              med_area = round(median(area_sum), 1),
              max_area = round(max(area_sum), 1)), by = samp])

## ---- A4. threshold calibration on RB1 --------------------------------------
## RB1 has no Hodgkin component, so every RB1 call is a false positive. Sweep and
## report; do not pick a fixed quantile.
tot <- md[, .N, by = samp]; setnames(tot, "N", "cells")
## Quantiles taken over IB1's objects, since IB1 is where the call is being made, and
## the grid is widened at the strict end where the interesting operating points are.
IB1_OBJ_N <- OBJ[samp == "IB1", .N]
QGRID <- unique(sort(c(seq(0.50, 0.99, length.out = 20),
                       seq(0.99, 0.9995, length.out = 10))))
thr <- rbindlist(lapply(QGRID, function(q) {
  cut <- quantile(OBJ[samp == "IB1", obj_score], q)
  x <- OBJ[obj_score >= cut, .(called = .N, med_area = median(area_sum),
                               pct_multi = 100*mean(n_segments > 1),
                               cd30_per_obj = median(cd30_sum)), by = samp]
  x <- merge(tot, x, by = "samp", all.x = TRUE)
  x[is.na(called), called := 0L]
  ## +1 pseudocount on the RB1 numerator. Without it, a threshold that happens to
  ## yield zero RB1 calls returns an infinite ratio and always wins the argmax, which
  ## would pick the most extreme cut on the strength of one absent false positive.
  data.table(q = q, cut = cut,
             IB1 = x[samp=="IB1", called], RB1 = x[samp=="RB1", called],
             ratio = round((x[samp=="IB1", called] / x[samp=="IB1", cells]) /
                           ((x[samp=="RB1", called] + 1) / x[samp=="RB1", cells]), 2),
             IB1_med_area = round(x[samp=="IB1", med_area], 1),
             IB1_pct_multi = round(x[samp=="IB1", pct_multi], 1))
}))
cat(sprintf("\n=== threshold sweep (the old per-segment call achieved %.2fx)\n", OLD_RATIO))
print(thr)
sc("hrs_threshold_sweep.csv", thr)
## Require enough IB1 calls to be worth rendering, scaled to the object pool -- a fixed
## floor is wrong at both ends (too strict on a small pool, far too loose on a large one)
MIN_CALLS <- max(5L, ceiling(0.002 * IB1_OBJ_N))
cat(sprintf("minimum IB1 calls required of a candidate threshold: %d (0.2%% of %d objects)\n",
            MIN_CALLS, IB1_OBJ_N))
cand <- thr[is.finite(ratio) & IB1 >= MIN_CALLS]
if (!nrow(cand)) {
  message("no threshold yields >= ", MIN_CALLS, " IB1 objects; falling back to q=0.99")
  best <- thr[which.min(abs(q - 0.99))]
} else best <- cand[order(-ratio)][1]
CUT <- best$cut
cat(sprintf("\nchosen cut at q=%.4f -> IB1 %d objects, RB1 %d, ratio %.2fx\n",
            best$q, best$IB1, best$RB1, best$ratio))
OBJ[, is_hrs := obj_score >= CUT]
md[, is_hrs := obj_id %in% OBJ[is_hrs == TRUE, obj_id]]

## ---- A5. validation 1: CD30, held out ---------------------------------------
md[, score_dec := cut(hrs_score, unique(quantile(hrs_score, 0:10/10)),
                      include.lowest = TRUE, labels = FALSE), by = samp]
cd30v <- md[, .(cells = .N, pct_cd30_pos = round(100*mean(cd30_pos), 2),
                mean_cd30 = round(mean(cd30_ct), 3)), by = .(samp, score_dec)][order(samp, score_dec)]
cat("\n=== held-out CD30 detection by score decile (must rise monotonically)\n")
print(dcast(cd30v, score_dec ~ samp, value.var = "pct_cd30_pos"))
auc_tab <- md[, .(auc_score_vs_cd30 = round(auc_of(hrs_score, cd30_pos), 4),
                  auc_morph_only    = round(auc_of(morph_rank, cd30_pos), 4)), by = samp]
print(auc_tab)
## Strict monotonicity over 10 deciles is brittle -- one noisy bin at the detection
## floor fails it, which is a flaw in the test rather than evidence about the score
## (the same problem was hit and documented at
## scripts/07_script_plasma_marker_check.R:100-105). Spearman rho across the deciles is
## the robust version; strict monotonicity is kept alongside as the stricter reading.
mono <- cd30v[, .(monotonic = all(diff(pct_cd30_pos) >= 0),
                  rho = round(suppressWarnings(
                          cor(score_dec, pct_cd30_pos, method = "spearman")), 3)), by = samp]
print(mono); sc("hrs_cd30_validation.csv", merge(merge(cd30v, auc_tab, by="samp"), mono, by="samp"))

## ---- A6. validation 2: morphology -------------------------------------------
## the `level` label is added AFTER aggregating: a length-1 constant in `by` is not
## recycled, it errors
morph <- rbind(
  md[,  .(cells = .N, med_area = median(Area.um2), med_nuc = median(NucArea)),
      by = .(samp, group = fifelse(is_hrs, "called", "background"))][, level := "segment"],
  OBJ[, .(cells = .N, med_area = median(area_sum), med_nuc = median(nuc_sum)),
      by = .(samp, group = fifelse(is_hrs, "called", "background"))][, level := "object"])
mw <- dcast(morph, samp + level ~ group, value.var = "med_area")
mw[, fold := round(called / background, 2)]
cat("\n=== morphology: median area, called vs background (old per-segment call: 1.6x)\n")
print(mw); sc("hrs_morphology.csv", morph)

## ---- A7. validation 3: T-cell rosette ---------------------------------------
## microenv() carried over from scripts/08_script_ib1_hierannot_bcells.R with both of
## its corrections: mask the index cell out of its own neighbourhood, and leave
## fractions NA below MIN_NBR rather than dividing by pmax(tot, 1).
microenv <- function(idx_dt, pool_dt) {
  tree <- as.matrix(pool_dt[, .(x_slide_mm, y_slide_mm)])
  qry  <- as.matrix(idx_dt[, .(x_slide_mm, y_slide_mm)])
  nn   <- RANN::nn2(tree, qry, k = min(K_NN, nrow(tree)))
  keep <- nn$nn.dists <= NN_UM / 1000
  sm <- nn$nn.idx == match(idx_dt$cell_id, pool_dt$cell_id); sm[is.na(sm)] <- FALSE
  keep <- keep & !sm
  lin <- as.character(pool_dt$lineage)[nn$nn.idx]; dim(lin) <- dim(nn$nn.idx)
  lin[!keep] <- NA_character_
  tot <- rowSums(keep)
  list(f_T = ifelse(tot >= MIN_NBR, rowSums(lin == "T", na.rm = TRUE) / tot, NA_real_),
       f_myeloid = ifelse(tot >= MIN_NBR, rowSums(lin == "myeloid", na.rm=TRUE)/tot, NA_real_),
       f_B = ifelse(tot >= MIN_NBR, rowSums(lin == "B", na.rm = TRUE) / tot, NA_real_),
       tot = tot)
}
ros <- rbindlist(lapply(SAMPLES, function(s) {
  pool <- md[samp == s]
  me   <- microenv(pool, pool)
  d <- data.table(samp = s, is_hrs = pool$is_hrs, cd30_pos = pool$cd30_pos,
                  f_T = me$f_T, f_myeloid = me$f_myeloid, f_B = me$f_B, n_nbr = me$tot,
                  obj_id = pool$obj_id)
  ## density-matched controls: same neighbour-count decile, not HRS
  d[, dens_dec := cut(n_nbr, unique(quantile(n_nbr, 0:10/10)), include.lowest=TRUE,
                      labels = FALSE)]
  d[]
}))
## .N is integer; assigning the normalised double back into it truncates to NA and
## silently propagates NaN through every standardised fraction below
w <- ros[!is.na(f_T), .(w = as.numeric(.N)), by = .(samp, dens_dec)][, w := w/sum(w), by = samp]
rstd <- merge(ros[!is.na(f_T), .(v = mean(f_T), n = .N), by = .(samp, dens_dec, is_hrs)],
              w, by = c("samp","dens_dec"))
rstd <- rstd[, .(f_T_density_std = round(100*sum(v*w)/sum(w), 2), cells = sum(n)),
             by = .(samp, is_hrs)]
cat("\n=== T-cell fraction within ", NN_UM, " um, density-standardised\n", sep = "")
print(rstd[order(samp, -is_hrs)])
## object as the unit, not the segment
ros_obj <- ros[!is.na(obj_id) & !is.na(f_T), .(f_T = mean(f_T), samp = first(samp),
                                               is_hrs = first(is_hrs)), by = obj_id]
cat(sprintf("object-level, IB1 called vs not: %.3f vs %.3f (Wilcoxon p = %s)\n",
            ros_obj[samp=="IB1" & is_hrs, mean(f_T)], ros_obj[samp=="IB1" & !is_hrs, mean(f_T)],
            signif(tryCatch(wilcox.test(ros_obj[samp=="IB1" & is_hrs, f_T],
                                        ros_obj[samp=="IB1" & !is_hrs, f_T])$p.value,
                            error = function(e) NA_real_), 3)))
sc("hrs_rosette.csv", rstd)

## ---- A8. validation 4: cluster composition ----------------------------------
comp <- md[samp == "IB1" & is_hrs == TRUE, .N, by = .(lineage)][order(-N)]
comp[, pct := round(100*N/sum(N), 1)]
cat("\n=== lineage of IB1 called HRS cells\n"); print(comp)
pc <- md[samp == "IB1" & is_hrs == TRUE, .N, by = pearson_clusters_low][order(-N)]
pc[, pct := round(100*N/sum(N), 1)]
cat("\n=== pearson_clusters_low of IB1 called HRS (old tier-A was 77.6% cluster 1 =\n")
cat("    HierAnnot \"Myeloid > Macrophage\")\n")
if (file.exists(HIER)) {
  h <- fread(HIER)[hierarchy == "immune_core" & clustering == "pearson_clusters_low",
                   .(pearson_clusters_low = as.character(cluster_id),
                     hierannot = annot_export_label, path = annot_path)]
  pc[, pearson_clusters_low := as.character(pearson_clusters_low)]
  pc <- merge(pc, h, by = "pearson_clusters_low", all.x = TRUE, sort = FALSE)[order(-N)]
} else message("missing ", HIER, " -- HierAnnot cross-check skipped")
print(pc); sc("hrs_cluster_composition.csv", pc)

## ---- outputs and figures -----------------------------------------------------
fwrite(md[, .(cell_id, samp, obj_id, is_hrs, hrs_score, z_pos, z_neg, z_apl, morph_rank,
              cd30_ct, cd30_pos, Area.um2, NucArea, nCount_RNA, lineage,
              pearson_clusters_low, pearson_clusters_high, x_slide_mm, y_slide_mm)],
       file.path(OUT, "hrs_segments.csv.gz"))
cat("wrote ", file.path(OUT, "hrs_segments.csv.gz"), "\n", sep = "")
fwrite(OBJ, file.path(OUT, "hrs_objects.csv.gz"))
cat("wrote ", file.path(OUT, "hrs_objects.csv.gz"), "\n", sep = "")

p1 <- ggplot(thr[is.finite(ratio)], aes(100*(1-q), ratio)) +
  geom_hline(yintercept = OLD_RATIO, colour = DIV_LO, linetype = "22", linewidth = .5) +
  geom_line(colour = DIV_HI, linewidth = .7) + geom_point(colour = DIV_HI, size = 1.6) +
  scale_x_log10() +
  labs(title = "IB1:RB1 specificity across the threshold",
       subtitle = paste0("RB1 has no Hodgkin component, so every RB1 call is a false ",
                         "positive. Dashed line = the 4.89x\nthe existing per-segment ",
                         "call achieved; above it is an improvement."),
       x = "% of objects called (log scale)", y = "IB1:RB1 call-rate ratio",
       caption = sprintf("Data: %s/hrs_threshold_sweep.csv", OUT)) + theme_viz()
sf("fig_10_threshold_sweep.png", p1, 8, 5)

p2 <- ggplot(cd30v, aes(score_dec, pct_cd30_pos, colour = samp)) +
  geom_line(linewidth = .7) + geom_point(size = 1.8) +
  scale_colour_manual(values = SAMP_COLS, name = NULL) + scale_x_continuous(breaks = 1:10) +
  labs(title = "CD30 was held out of the score, so this is an independent test",
       subtitle = paste0("If the score finds HRS cells, TNFRSF8 detection must rise with ",
                         "it. A flat line means the\nscore is not selecting CD30+ cells."),
       x = "hrs_score decile", y = "% of cells with TNFRSF8 >= 1",
       caption = sprintf("Data: %s/hrs_cd30_validation.csv", OUT)) + theme_viz()
sf("fig_10_cd30_heldout.png", p2, 8, 5)

p3 <- ggplot(sw, aes(radius_um, pct_multi, colour = samp)) +
  geom_vline(xintercept = MERGE_UM, colour = INK, linetype = "22", linewidth = .4) +
  geom_line(linewidth = .7) + geom_point(size = 1.8) +
  scale_colour_manual(values = SAMP_COLS, name = NULL) +
  labs(title = "Merging adjacent CD30+ segments: radius sensitivity",
       subtitle = paste0("Dashed = the chosen ", MERGE_UM, " um. A cliff here would mean the ",
                         "merge is an artefact of the radius;\nsee hrs_merge_null.csv for ",
                         "whether the adjacency beats a density-matched permutation."),
       x = "merge radius (um)", y = "% of objects with >1 segment",
       caption = sprintf("Data: %s/hrs_merge_radius_sweep.csv", OUT)) + theme_viz()
sf("fig_10_merge_radius.png", p3, 8, 5)

hb <- OBJ[samp == "IB1"]
p4 <- ggplot(hb, aes(area_sum, fill = is_hrs)) +
  geom_histogram(bins = 60, position = "identity", alpha = .75) +
  scale_x_log10() + scale_fill_manual(values = c(`FALSE` = "#a9a7a1", `TRUE` = DIV_HI),
                                      name = NULL, labels = c("other CD30+ object","called HRS")) +
  labs(title = "Pooled object area, IB1",
       subtitle = paste0("HRS cells are giant. The existing per-segment call reached only ",
                         "1.6x background median\narea; merging fragments is what would ",
                         "recover the true size if segmentation split them."),
       x = "summed Area.um2 of the object (log scale)", y = "objects",
       caption = sprintf("Data: %s/hrs_objects.csv.gz", OUT)) + theme_viz()
sf("fig_10_object_area.png", p4, 8, 5)

p5 <- ggplot(md[samp == "IB1"], aes(x_slide_mm, y_slide_mm)) +
  geom_point(data = md[samp=="IB1" & is_hrs == FALSE], colour = "#e4e3dc", size = .12, stroke = 0) +
  geom_point(data = md[samp=="IB1" & is_hrs == TRUE], colour = DIV_HI, size = .9, stroke = 0) +
  coord_equal() +
  labs(title = "Where the called HRS cells sit in IB1",
       subtitle = "Blue = called HRS, pale grey = all other cells (cluster i excluded).",
       x = "x (mm)", y = "y (mm)",
       caption = sprintf("Data: %s/hrs_segments.csv.gz", OUT)) +
  theme_viz() + theme(panel.grid.major.y = element_blank())
sf("fig_10_hrs_spatial_IB1.png", p5, 8.5, 9)

## ---- verdict -----------------------------------------------------------------
.one <- function(v) if (length(v) == 1) as.numeric(v) else NA_real_
ib_fold <- .one(mw[samp=="IB1" & level=="object", fold])
SUM <- data.table(
  metric = c("IB1_objects_called","RB1_objects_called","IB1_RB1_ratio","old_ratio",
             "IB1_pct_multi_segment","merge_enrichment_vs_null_IB1","merge_p_emp_IB1",
             "merge_enrichment_vs_null_RB1",
             "auc_score_vs_cd30_IB1","cd30_decile_rho_IB1","cd30_monotonic_IB1",
             "object_area_fold_IB1","segment_area_fold_IB1",
             "fT_called_IB1","fT_background_IB1"),
  value = c(best$IB1, best$RB1, best$ratio, OLD_RATIO,
            round(100*OBJ[samp=="IB1" & is_hrs==TRUE, mean(n_segments>1)],1),
            .one(NULLT[samp=="IB1", enrichment]), .one(NULLT[samp=="IB1", p_emp]),
            .one(NULLT[samp=="RB1", enrichment]),
            .one(auc_tab[samp=="IB1", auc_score_vs_cd30]),
            .one(mono[samp=="IB1", rho]),
            as.numeric(.one(mono[samp=="IB1", monotonic])),
            ib_fold, .one(mw[samp=="IB1" & level=="segment", fold]),
            .one(rstd[samp=="IB1" & is_hrs==TRUE, f_T_density_std]),
            .one(rstd[samp=="IB1" & is_hrs==FALSE, f_T_density_std])))
sc("summary.csv", SUM)
cat("\n============================================================\nSUMMARY\n")
print(SUM)
cat(sprintf("\n[%s] specificity: IB1:RB1 %.2fx vs the old %.2fx\n",
            ifelse(is.finite(best$ratio) && best$ratio > OLD_RATIO, "ok", "NO"),
            best$ratio, OLD_RATIO))
## The gate needs an EFFECT SIZE, not just a p-value. With ~15k objects the
## permutation p goes below 0.05 on a 3% excess, which is not evidence that adjacent
## CD30+ segments are one cell. It also has to beat RB1 -- if the FL recurrence shows
## the same excess, the adjacency is a property of CD30+ density, not of Hodgkin.
enr_ib <- .one(NULLT[samp=="IB1", enrichment]); enr_rb <- .one(NULLT[samp=="RB1", enrichment])
merge_ok <- is.finite(enr_ib) && enr_ib >= 1.5 &&
            .one(NULLT[samp=="IB1", p_emp]) < 0.05 && enr_ib > 1.2 * enr_rb
cat(sprintf("[%s] merge vs its density-matched null: IB1 %.2fx (p=%.3f), RB1 %.2fx\n",
            ifelse(merge_ok, "ok", "NO"), enr_ib, .one(NULLT[samp=="IB1", p_emp]), enr_rb))
if (!merge_ok)
  cat("     Adjacent CD30+ segments are no more common than chance given local density,\n",
      "     and RB1 shows the same. Fragmentation of giant cells is NOT supported --\n",
      "     read the object-level results as segments, not as reconstructed cells.\n", sep="")
## AUC is the gate; rho is context. When CD30 positivity is only a few percent, most
## deciles sit on the detection floor and the rank correlation across them is noise --
## measured on the synthetic fixture, AUC 0.751 came with rho 0.30 purely for that
## reason. Reported so a low rho is not mistaken for a failing score.
cat(sprintf("[%s] held-out CD30 AUC %.3f  (decile rho %.2f, strictly monotonic: %s)\n",
            ifelse(.one(auc_tab[samp=="IB1", auc_score_vs_cd30]) > 0.6, "ok", "NO"),
            .one(auc_tab[samp=="IB1", auc_score_vs_cd30]),
            .one(mono[samp=="IB1", rho]),
            .one(mono[samp=="IB1", monotonic]) == 1))
cat(sprintf("[%s] object area %.2fx background (old per-segment call: 1.6x)\n",
            ifelse(is.finite(ib_fold) && ib_fold > 1.6, "ok", "NO"), ib_fold))
cat("\nNone of this is a diagnosis. Render the views next -- the project has always\n")
cat("held that HRS calls MUST be confirmed visually, and a giant multinucleated CD30+\n")
cat("cell in the image is the only thing that settles it:\n")
cat("  scripts/10_submit_hrs_views.sh\n")
cat("============================================================\nOK\n")
