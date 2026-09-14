## ============================================================================
## lineage.R -- the marker-ratio argmax lineage caller.
##
## PROVENANCE. `assign_lineage()` below is a VERBATIM lift of the function defined inline
## at scripts/01_import_cosmx_seurat.R:119-168, which produced the `lineage` column that
## every later step reads. It is copied rather than sourced because step 01 is a script
## with side effects (it reads 8 GB vendor objects and writes the BPCells stores), so it
## cannot be sourced for the function alone -- and refactoring a completed pipeline stage
## to export it would invalidate its outputs.
##
## If step 01's copy is ever changed, change this one too. `diff` them:
##   sed -n '110,168p' scripts/01_import_cosmx_seurat.R
##
## WHY STEP 16 NEEDS IT. The `lineage` column is the only cell-type call in this project
## that is comparable ACROSS samples -- `pearson_clusters_low/high` ids are per-sample and
## mean different things in IB1 and RB1, and the step-09b marker audit has only ever been
## run for IB1. But `lineage` has no T-cell or myeloid subsets, and the step-16 checkpoint
## question is specifically about Tfh (PD-1) and macrophages/cDC (PD-L1). So step 16b
## re-runs this same caller with a finer marker-set list, on both samples at once, so the
## fold-over-mean scaling is shared and the calls stay comparable.
## ============================================================================

suppressPackageStartupMessages({
  library(Matrix)
})

## Row/column extrema. Step 01 got `colMaxs`/`rowMax` from its attached packages; these
## are local equivalents, kept dependency-free.
##
## VECTORISED ON PURPOSE. The obvious apply(m, 1, max) is unusable here: step 16 calls this
## on all 4.4M cells at once (so the fold-over-mean scaling is shared between the two
## samples), and apply() over 4.4M rows builds 4.4M subvectors. Reducing with pmax over the
## columns touches each element once instead.
.colMaxs <- function(m) do.call(pmax, c(lapply(seq_len(nrow(m)), function(i) m[i, ]),
                                        list(na.rm = TRUE)))
.rowMax  <- function(m) do.call(pmax, c(lapply(seq_len(ncol(m)), function(j) m[, j]),
                                        list(na.rm = TRUE)))

# Coarse-lineage assignment (the k:l ratio logic, generalized to N marker sets).
# Per gene: fold-over-mean, winsorized at `cap`, so sets of differing brightness
# are comparable. Per set: mean of the scaled values (size-normalized). A cell is
# ASSIGNED its argmax set only if that set (a) clears background in RAW counts --
# else "unassigned" -- and (b) beats the runner-up by `margin` log2 -- else
# "ambiguous". The raw floor is per-cell (negative-probe rate), so spillover into
# an intrinsically-dim lineage cannot mint a call. `cts` = counts (BPCells ok),
# `neg_total` = per-cell negprobe totals, `n_neg` = number of negative probes.
assign_lineage <- function(cts, sets, neg_total, n_neg,
                           margin = 1.5, bg_mult = 1.5,
                           min_raw = .2, cap = 10, use_max_marker = TRUE) {
  sets <- lapply(sets, function(g) intersect(g, rownames(cts)))
  sets <- sets[lengths(sets) > 0]
  M    <- as.matrix(cts[unique(unlist(sets)), , drop = FALSE])   # markers x cells
  gmean <- rowSums(M) / pmax(rowSums(M > 0), 1)   # mean among detected cells
  gmean <- pmax(gmean, 1e-3)
  S    <- M / gmean
  S[S > cap] <- cap                           # fold-over-mean, winsorized

  marker_agg_FUN <- if (use_max_marker) .colMaxs else colMeans

  scaled <- sapply(sets, function(g) marker_agg_FUN(S[g, , drop = FALSE]))  # cells x sets
  raw    <- sapply(sets, function(g) marker_agg_FUN(M[g, , drop = FALSE]))
  bg     <- as.numeric(neg_total) / max(n_neg, 1)                # per-cell background

  ii     <- seq_len(nrow(scaled))
  top_i  <- max.col(scaled, ties.method = "random")               # argmax (vectorized)
  s2     <- scaled; s2[cbind(ii, top_i)] <- -Inf
  sec_i  <- max.col(s2, ties.method = "random")
  top_s  <- scaled[cbind(ii, top_i)]; sec_s <- s2[cbind(ii, sec_i)]
  top_raw<- raw[cbind(ii, top_i)]

  pass   <- top_raw > pmax(bg_mult * bg, min_raw)                # real signal vs background

  l2m    <- log2((top_s + 1e-6) / (sec_s + 1e-6))
  lab    <- colnames(scaled)[top_i]
  lab_sec<- colnames(scaled)[sec_i]

  lab[!pass]     <- "below_background"
  lab_sec[!pass] <- "below_background"
  all_sets_0 <- .rowMax(raw) == 0
  lab[all_sets_0]     <- "all_markers_0"
  lab_sec[all_sets_0] <- "all_markers_0"

  lab.raw <- lab
  lab[pass & l2m < log2(margin)] <- "ambiguous"
  list(label = lab, log2_margin = l2m, raw_label = lab.raw, sec_label = lab_sec)
}

## Step 01's parameters, kept here so 16b cannot silently drift from the `lineage` column
## it is meant to refine. 01_import_cosmx_seurat.R:221-227 calls with these plus
## use_max_marker = FALSE (note: NOT the function default).
LINEAGE_CALL_ARGS <- list(margin = 1.5, bg_mult = 1.5, min_raw = .2, cap = 10,
                          use_max_marker = FALSE)

## Step 01's coarse sets, verbatim from 01_import_cosmx_seurat.R:74-81, with its comments.
LINEAGE_SETS_COARSE <- list(
  B       = c("MS4A1", "CD79A", "CD79B", "BANK1"),          # drop CD19 (dropout)
  T       = c("CD3D", "CD3E", "CD3G", "TRAC"),              # drop CD2/CD5 (NK/B overlap)
  NK      = c("NKG7", "GNLY", "KLRD1", "NCR1"),             # drop NCAM1/KLRB1 (T overlap)
  myeloid = c("CD68", "CD163", "LYZ", "ITGAX", "C1QA", "C1QB", "CD14"),
  endo    = c("PECAM1", "VWF", "CLDN5", "CD34"),
  fibro   = c("COL1A1", "COL1A2", "PDGFRB", "DCN", "ACTA2"),
  PC      = c("MZB1", "XBP1", "PRDM1", "SDC1", "DERL3", "TNFRSF17"))
# FDC dropped: CR2/FCER2 are B markers; nothing FDC-exclusive survives on this panel
