#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------------------
## step 21c -- ARE THE TWO LESIONS ARRANGED DEPENDENTLY?  Replacing step 12's distance test.
##
## WHAT IS BEING REPLACED, AND WHY. Step 13's claim "no spatial relationship between the two
## diseases" rests entirely on one test in step 12: the distance from each Hodgkin hotspot
## bin to the nearest FL-like focus (median 1035 um) against the same distance for
## non-hotspot bins (1205 um), Wilcoxon p = 0.23 on 11 foci. Three problems, and they all
## push the same way -- toward a null that was never really tested:
##
##   1. THE UNIT IS THE BIN. Bins inside one FOV are not independent; both lesions are
##      spatially clustered, so the effective n is nearer the FOV count than the bin count.
##      This is the project's own pseudoreplication rule (region is the unit of every
##      p-value) applied to the one test that never got it. Step 11 lost a Tfh result this
##      way and step 17 lost the sign-consistency result this way (ICC 0.875).
##   2. THE COMPARISON SET IS UNMATCHED. "Non-hotspot bins" differ from hotspot bins in cell
##      count and depth, and a bin can only BE a hotspot if it has enough cells to compute
##      a stable cp10k -- so the reference carries a different propensity to be labelled.
##   3. TISSUE GEOMETRY IS IGNORED. Both lesions can only occur where tissue exists, and
##      IB1 is in three pieces (21a). A distance test that does not hold the tissue mask
##      fixed measures the shape of the section as much as the arrangement of the disease.
##
## THE TEST HERE. Two statistics, in both directions, with the FOV respected and the tissue
## mask held fixed by construction:
##
##      obs_h2f   median distance from each hotspot bin to the nearest FL bin
##      obs_f2h   median distance from each FL bin to the nearest hotspot bin
##
## against TWO INDEPENDENT NULLS that fail in opposite ways, so agreement between them means
## something that either alone would not:
##
##   NULL R (primary) -- FOCUS RELOCATION. Break the hotspot bins into connected foci, keep
##     each focus's size and internal shape exactly, and re-site each focus's centroid at a
##     random occupied bin. This PRESERVES the number of foci, their sizes and their internal
##     clustering, and the tissue mask; it randomises only WHERE each focus sits. That is the
##     null the question needs: it asks whether the two patterns are positioned dependently,
##     without assuming either is internally random. A label permutation cannot do this -- it
##     destroys the clustering it should hold fixed, which tightens the null and makes the
##     test anti-conservative.
##
##     A RIGID ROTATION of the whole pattern was tried first and IS NOT USABLE HERE: the
##     primary tissue piece fills only ~33% of its own bounding box, so a 23-point pattern
##     spanning the piece never lands back on tissue -- 0 of 5,000 draws were admissible.
##     Relocating foci one at a time fits because each focus is a few bins across.
##
##   NULL Rd (sensitivity) -- FOCUS RELOCATION, DEFECT-MATCHED. 21e measured that the
##     over-fixation defect deletes cells very unevenly across these regions: hotspot bins
##     lost 57.1% of their cells, background 47.4%, FL bins only 33.8%. Since the defect is
##     a peripheral rim, hotspots sit in more peripheral, more depleted tissue than the FL
##     foci do -- which would push the observed hotspot-to-FL distance UP and could MASK a
##     real co-location. Plain NULL R sites foci uniformly over occupied tissue, i.e. in
##     systematically less-depleted tissue than the real hotspots occupy. NULL Rd restricts
##     the siting pool to bins whose defect fraction falls inside the observed hotspots' own
##     interquartile range, so the null foci sit in tissue as depleted as the real ones.
##     If the verdict is the same under both, the defect geography is not carrying it.
##
##   NULL P (secondary) -- LABEL PERMUTATION among candidate bins matched on cell count and
##     depth deciles within the same tissue piece. This holds the labelling PROPENSITY fixed
##     and asks a narrower question: is the arrangement explained by where dense, deep tissue
##     happens to sit? It destroys clustering by construction, so it is reported as a
##     propensity check and never as the primary result.
##
## Restricted to the tissue piece carrying both lesions in quantity; the others are reported
## descriptively, because a rotation null inside a small irregular piece has almost no
## admissible draws.
##
## Bins only, so light enough for the login node.
## Usage: Rscript scripts/21_script_interface_arrangement.R [tier]
## ---------------------------------------------------------------------------------------
suppressMessages({library(data.table)})
setDTthreads(1); options(width = 200)

TIER <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "matched"
GEO  <- file.path("output/21_interface", TIER)
OUT  <- GEO
sc <- function(f, x) { fwrite(x, file.path(OUT, f)); cat("wrote", file.path(OUT, f), "\n") }
if (!file.exists(file.path(GEO, "bin_axes.csv.gz")))
  stop("run 21_script_interface_geometry.R first", call. = FALSE)

BIN_F <- file.path("output/16_regions", TIER, "bins_all_scales.csv.gz")
NROT  <- 5000L; NPERM <- 5000L; MIN_ON_TISSUE <- 0.90; SEED <- 20L
BIN_UM <- 100

cat("========================================================================\n")
cat("step 21c -- arrangement of the two lesions, tier =", TIER, "\n")
cat("========================================================================\n")

B <- fread(BIN_F)[bin_um == 100 & samp == "IB1" & bin_low_info == FALSE]
PZ <- fread(file.path(GEO, "tissue_pieces.csv"))
## Recover the piece labels the same way 21a did, so the two scripts cannot drift apart.
hc <- hclust(dist(as.matrix(B[, .(x_mm, y_mm)])), method = "single")
B[, piece := cutree(hc, h = 0.4)]
B[, `:=`(is_hot = bin_region == "cHL_hotspot", is_fl = bin_region == "FL_focus")]

cat("\n=== both lesions by tissue piece\n")
byp <- B[, .(bins = .N, fovs = uniqueN(fov), hotspots = sum(is_hot), fl = sum(is_fl)),
         by = piece][order(-bins)][bins >= 20]
print(byp)
PIECE <- byp[hotspots >= 10 & fl >= 20][order(-bins)][1, piece]
if (is.na(PIECE)) stop("no tissue piece carries both lesions in usable quantity", call. = FALSE)
cat(sprintf("\nprimary piece for the rotation null: %d (%d hotspots, %d FL bins)\n",
            PIECE, byp[piece == PIECE, hotspots], byp[piece == PIECE, fl]))
cat("other pieces are reported descriptively below.\n")

## ---- statistics -------------------------------------------------------------------------
## Integer bin-index space: a 100um lattice, so a rotation can be snapped back to the lattice
## and membership in the tissue mask is an exact set lookup rather than a tolerance.
nn_med <- function(from, to) {
  if (!nrow(from) || !nrow(to)) return(NA_real_)
  median(apply(from, 1, function(p) sqrt(min((to[,1]-p[1])^2 + (to[,2]-p[2])^2)))) * BIN_UM
}

Z  <- B[piece == PIECE]
occ <- unique(Z[, .(bx, by_)])
key <- function(a, b) paste(a, b, sep = "_")
OCCK <- key(occ$bx, occ$by_)
HOT <- as.matrix(Z[is_hot == TRUE, .(bx, by_)])
FL  <- as.matrix(Z[is_fl  == TRUE, .(bx, by_)])
obs_h2f <- nn_med(HOT, FL); obs_f2h <- nn_med(FL, HOT)
cat(sprintf("\nOBSERVED  hotspot -> nearest FL   %.0f um  (n = %d hotspots)\n", obs_h2f, nrow(HOT)))
cat(sprintf("OBSERVED  FL -> nearest hotspot   %.0f um  (n = %d FL bins)\n", obs_f2h, nrow(FL)))

## ---- NULL R: focus relocation ----------------------------------------------------------
set.seed(SEED)
## Foci = connected groups of hotspot bins. Single linkage at 200um on the 100um lattice, so
## bins that touch (including diagonally) are one focus.
if (nrow(HOT) >= 2) {
  fh <- cutree(hclust(dist(HOT), method = "single"), h = 2.0)
} else fh <- rep(1L, nrow(HOT))
FOCI <- lapply(split(seq_len(nrow(HOT)), fh), function(k) {
  m <- HOT[k, , drop = FALSE]; sweep(m, 2, round(colMeans(m))) })
cat(sprintf("\nNULL R: %d hotspot bins form %d connected foci (sizes %s)\n",
            nrow(HOT), length(FOCI), paste(vapply(FOCI, nrow, 0L), collapse = ",")))

OCCM <- as.matrix(occ)
## mnn_self = mean nearest-neighbour distance WITHIN the drawn hotspot pattern. It is the
## diagnostic that decides which null to believe: a null that preserves the observed
## clustering reproduces the observed mnn_self, and one that destroys it inflates it.
mnn_self <- function(m) { if (nrow(m) < 3) return(NA_real_)
  d <- as.matrix(dist(m)); diag(d) <- Inf; mean(apply(d, 1, min)) * BIN_UM }
OBS_MNN <- mnn_self(HOT)
rotR <- matrix(NA_real_, NROT, 3); n_acc <- 0L
for (i in seq_len(NROT)) {
  pts <- vector("list", length(FOCI)); okdraw <- TRUE
  for (j in seq_along(FOCI)) {
    placed <- NULL
    for (try in 1:25) {
      c0 <- OCCM[sample(nrow(OCCM), 1L), ]
      P  <- sweep(FOCI[[j]], 2, c0, "+")
      on <- key(P[,1], P[,2]) %in% OCCK
      if (mean(on) >= MIN_ON_TISSUE) { placed <- P[on, , drop = FALSE]; break }
    }
    if (is.null(placed)) { okdraw <- FALSE; break }
    pts[[j]] <- placed
  }
  if (!okdraw) next
  Pk <- unique(do.call(rbind, pts))
  if (nrow(Pk) < MIN_ON_TISSUE*nrow(HOT)) next
  n_acc <- n_acc + 1L
  rotR[n_acc, ] <- c(nn_med(Pk, FL), nn_med(FL, Pk), mnn_self(Pk))
}
rotR <- rotR[seq_len(n_acc), , drop = FALSE]
cat(sprintf("  %d of %d relocations admissible (>= %.0f%% of each focus back on tissue)\n",
            n_acc, NROT, 100*MIN_ON_TISSUE))
if (n_acc < 200L) cat("  WARNING: too few admissible draws for a stable null.\n")

pR <- function(obs, v) { v <- v[is.finite(v)]
  if (!length(v)) return(NA_real_)
  2*min((1+sum(v <= obs))/(1+length(v)), (1+sum(v >= obs))/(1+length(v))) }
cat(sprintf("  hotspot -> FL   obs %.0f um   null median %.0f  [%.0f, %.0f]   p = %.3f\n",
    obs_h2f, median(rotR[,1], na.rm=TRUE), quantile(rotR[,1],.025,na.rm=TRUE),
    quantile(rotR[,1],.975,na.rm=TRUE), pR(obs_h2f, rotR[,1])))
cat(sprintf("  FL -> hotspot   obs %.0f um   null median %.0f  [%.0f, %.0f]   p = %.3f\n",
    obs_f2h, median(rotR[,2], na.rm=TRUE), quantile(rotR[,2],.025,na.rm=TRUE),
    quantile(rotR[,2],.975,na.rm=TRUE), pR(obs_f2h, rotR[,2])))

## ---- NULL Rd: focus relocation restricted to comparably depleted tissue -----------------
DFB_F <- file.path(GEO, "defect_per_bin.csv.gz")
if (file.exists(DFB_F)) {
  DFB <- fread(DFB_F)
  Zd  <- merge(Z[, .(bx, by_, is_hot, is_fl)], DFB[, .(bx, by_, defect_frac)],
               by = c("bx","by_"), all.x = TRUE)
  hq  <- quantile(Zd[is_hot == TRUE, defect_frac], c(.25,.75), na.rm = TRUE)
  poolD <- as.matrix(Zd[!is.na(defect_frac) & defect_frac >= hq[1] & defect_frac <= hq[2],
                        .(bx, by_)])
  cat(sprintf("\nNULL Rd: hotspot defect fraction IQR [%.3f, %.3f]; siting pool %d of %d occupied bins\n",
              hq[1], hq[2], nrow(poolD), nrow(OCCM)))
  cat(sprintf("  (median defect fraction: hotspots %.3f, FL %.3f, all occupied %.3f)\n",
              median(Zd[is_hot==TRUE, defect_frac], na.rm=TRUE),
              median(Zd[is_fl==TRUE, defect_frac], na.rm=TRUE),
              median(Zd$defect_frac, na.rm=TRUE)))
  set.seed(SEED)
  rotD <- matrix(NA_real_, NROT, 3); nD <- 0L
  if (nrow(poolD) >= 50) for (i in seq_len(NROT)) {
    pts <- vector("list", length(FOCI)); okd <- TRUE
    for (j in seq_along(FOCI)) {
      placed <- NULL
      for (try in 1:25) {
        c0 <- poolD[sample(nrow(poolD), 1L), ]
        P  <- sweep(FOCI[[j]], 2, c0, "+")
        on <- key(P[,1], P[,2]) %in% OCCK
        if (mean(on) >= MIN_ON_TISSUE) { placed <- P[on, , drop = FALSE]; break }
      }
      if (is.null(placed)) { okd <- FALSE; break }
      pts[[j]] <- placed
    }
    if (!okd) next
    Pk <- unique(do.call(rbind, pts))
    if (nrow(Pk) < MIN_ON_TISSUE*nrow(HOT)) next
    nD <- nD + 1L; rotD[nD, ] <- c(nn_med(Pk, FL), nn_med(FL, Pk), mnn_self(Pk))
  }
  rotD <- rotD[seq_len(nD), , drop = FALSE]
  cat(sprintf("  %d of %d relocations admissible\n", nD, NROT))
  if (nD >= 200L) {
    cat(sprintf("  hotspot -> FL   obs %.0f um   null median %.0f  [%.0f, %.0f]   p = %.3f\n",
        obs_h2f, median(rotD[,1], na.rm=TRUE), quantile(rotD[,1],.025,na.rm=TRUE),
        quantile(rotD[,1],.975,na.rm=TRUE), pR(obs_h2f, rotD[,1])))
    cat(sprintf("  FL -> hotspot   obs %.0f um   null median %.0f  [%.0f, %.0f]   p = %.3f\n",
        obs_f2h, median(rotD[,2], na.rm=TRUE), quantile(rotD[,2],.025,na.rm=TRUE),
        quantile(rotD[,2],.975,na.rm=TRUE), pR(obs_f2h, rotD[,2])))
  } else cat("  too few admissible draws; NULL Rd not reported\n")
} else {
  cat("\nNULL Rd skipped: run 21_script_interface_defect_qc.R for defect_per_bin.csv.gz\n")
  rotD <- matrix(NA_real_, 0, 3); nD <- 0L
}

## ---- NULL P: propensity-matched label permutation ---------------------------------------
## Candidates are non-FL bins in the same piece, matched to the hotspots on cell-count and
## depth deciles. This asks whether the arrangement follows from where dense, deep tissue is.
Z[, `:=`(cell_d = cut(cells, quantile(cells, 0:5/5), include.lowest = TRUE, labels = FALSE),
         dep_d  = cut(med_depth, quantile(med_depth, 0:5/5), include.lowest = TRUE, labels = FALSE))]
strat <- Z[is_hot == TRUE, .N, by = .(cell_d, dep_d)]
pool  <- Z[is_fl == FALSE]
permP <- matrix(NA_real_, NPERM, 3)
for (i in seq_len(NPERM)) {
  pick <- pool[strat, on = .(cell_d, dep_d), {
    k <- min(.N, N[1]); if (!k || is.na(k)) NULL else .SD[sample(.N, k)] },
    by = .EACHI, .SDcols = c("bx","by_")]
  if (nrow(pick) < 0.8*nrow(HOT)) next
  Pk <- as.matrix(pick[, .(bx, by_)])
  permP[i, ] <- c(nn_med(Pk, FL), nn_med(FL, Pk), mnn_self(Pk))
}
cat(sprintf("\nNULL P: %d of %d permutations usable\n", sum(is.finite(permP[,1])), NPERM))
cat(sprintf("  hotspot -> FL   obs %.0f um   null median %.0f  [%.0f, %.0f]   p = %.3f\n",
    obs_h2f, median(permP[,1], na.rm=TRUE), quantile(permP[,1],.025,na.rm=TRUE),
    quantile(permP[,1],.975,na.rm=TRUE), pR(obs_h2f, permP[,1])))
cat(sprintf("  FL -> hotspot   obs %.0f um   null median %.0f  [%.0f, %.0f]   p = %.3f\n",
    obs_f2h, median(permP[,2], na.rm=TRUE), quantile(permP[,2],.025,na.rm=TRUE),
    quantile(permP[,2],.975,na.rm=TRUE), pR(obs_f2h, permP[,2])))

## ---- WHICH NULL TO BELIEVE, measured rather than asserted -------------------------------
## The two nulls disagree. The doc block above predicts why: NULL P destroys the clustering
## it should hold fixed, so its patterns are more spread out, which pushes every
## nearest-neighbour distance up and makes the observed value look unusually short. That
## prediction is checkable on the drawn patterns themselves.
cat("\n=== which null holds the observed clustering fixed?\n")
adj <- data.table(
  pattern = c("OBSERVED hotspots", "NULL R draws (uniform relocation)",
              "NULL Rd draws (defect-matched relocation)",
              "NULL P draws (label permutation)"),
  mnn_self_um = round(c(OBS_MNN, median(rotR[,3], na.rm = TRUE),
                        if (nD >= 200L) median(rotD[,3], na.rm = TRUE) else NA_real_,
                        median(permP[,3], na.rm = TRUE))))
adj[, ratio_to_observed := round(mnn_self_um / adj$mnn_self_um[1], 3)]
print(adj); sc("arrangement_null_adjudication.csv", adj)
## Index by pattern name, not by row number: adding NULL Rd shifted the rows and the first
## version of these lines then attributed NULL P's inflation to the wrong null.
.r <- function(pat) adj[grepl(pat, pattern, fixed = TRUE), ratio_to_observed][1]
cat(sprintf("\n  Both relocation nulls reproduce the observed within-pattern spacing:\n"))
cat(sprintf("    NULL R  %.0f%%   NULL Rd %.0f%%   -> both hold the clustering fixed\n",
            100*.r("NULL R draws"), 100*.r("NULL Rd draws")))
cat(sprintf("    NULL P  %.0f%%  -> inflates it by %.0f%%, i.e. it SPREADS the hotspots out\n",
            100*.r("NULL P draws"), 100*(.r("NULL P draws")-1)))
cat("  A spread-out null pattern sits farther from everything, the FL foci included, so\n")
cat(sprintf("  NULL P's nominal p = %.3f is an artifact of its null and is not evidence about\n",
            pR(obs_h2f, permP[,1])))
cat("  arrangement. The two relocation nulls are the admissible pair, and they differ only\n")
cat("  in whether they match the hotspots' own defect burden -- which is what brackets the\n")
cat("  answer rather than settling it.\n")

## ---- the FOV-level version of step 12's own statistic ----------------------------------
## Step 12 asked whether hotspot bins are closer to FL than other bins are, with the BIN as
## the unit. Recomputed with the FOV as the unit, which is the only change: each FOV
## contributes one hotspot value and one background value, and the test is paired where both
## exist. This isolates how much of p = 0.23 was the unit choice.
ALL <- B[piece == PIECE]
ALL[, d_fl_um := apply(as.matrix(.SD), 1, function(p)
      sqrt(min((FL[,1]-p[1])^2 + (FL[,2]-p[2])^2))) * BIN_UM, .SDcols = c("bx","by_")]
fv <- ALL[is_fl == FALSE, .(hot = median(d_fl_um[is_hot]), bg = median(d_fl_um[!is_hot]),
                            n_hot = sum(is_hot)), by = fov][n_hot > 0 & is.finite(bg)]
cat(sprintf("\n=== step 12's statistic, FOV as the unit: %d FOVs carry a hotspot\n", nrow(fv)))
if (nrow(fv) >= 5) {
  wt <- suppressWarnings(wilcox.test(fv$hot, fv$bg, paired = TRUE))
  cat(sprintf("  paired median: hotspot %.0f um vs other bins %.0f um in the same FOV\n",
              median(fv$hot), median(fv$bg)))
  cat(sprintf("  Wilcoxon signed rank p = %.3f  (step 12, bin as unit, reported p = 0.233)\n",
              wt$p.value))
} else cat("  too few FOVs to pair.\n")

## ---- descriptive: the other pieces ------------------------------------------------------
cat("\n=== the other tissue pieces, descriptively\n")
for (pz in setdiff(byp$piece, PIECE)) {
  z <- B[piece == pz]
  h <- as.matrix(z[is_hot == TRUE, .(bx, by_)]); f <- as.matrix(z[is_fl == TRUE, .(bx, by_)])
  cat(sprintf("  piece %d: %d hotspots, %d FL bins -> %s\n", pz, nrow(h), nrow(f),
      if (nrow(h) && nrow(f))
        sprintf("hotspot->FL %.0f um, FL->hotspot %.0f um (no test; too few)",
                nn_med(h,f), nn_med(f,h))
      else "one lesion absent; no within-piece distance is defined"))
}

sc("arrangement_verdict.csv", data.table(
  tier = TIER, piece = PIECE, n_hotspots = nrow(HOT), n_fl_bins = nrow(FL),
  obs_h2f_um = round(obs_h2f), obs_f2h_um = round(obs_f2h),
  nullR_foci = length(FOCI), nullR_draws = n_acc,
  nullR_h2f_med = round(median(rotR[,1], na.rm=TRUE)), nullR_h2f_p = round(pR(obs_h2f, rotR[,1]),4),
  nullR_f2h_med = round(median(rotR[,2], na.rm=TRUE)), nullR_f2h_p = round(pR(obs_f2h, rotR[,2]),4),
  nullP_h2f_med = round(median(permP[,1], na.rm=TRUE)), nullP_h2f_p = round(pR(obs_h2f, permP[,1]),4),
  nullP_f2h_med = round(median(permP[,2], na.rm=TRUE)), nullP_f2h_p = round(pR(obs_f2h, permP[,2]),4),
  nullRd_draws = nD,
  nullRd_h2f_med = if (nD >= 200L) round(median(rotD[,1], na.rm=TRUE)) else NA_real_,
  nullRd_h2f_p   = if (nD >= 200L) round(pR(obs_h2f, rotD[,1]), 4) else NA_real_,
  nullRd_f2h_med = if (nD >= 200L) round(median(rotD[,2], na.rm=TRUE)) else NA_real_,
  nullRd_f2h_p   = if (nD >= 200L) round(pR(obs_f2h, rotD[,2]), 4) else NA_real_,
  obs_mnn_self_um = round(OBS_MNN),
  nullR_mnn_self_um = round(median(rotR[,3], na.rm=TRUE)),
  nullP_mnn_self_um = round(median(permP[,3], na.rm=TRUE)),
  ## Two numeric columns rather than one "[lo, hi]" string: the embedded comma made the
  ## field need quoting, which every downstream consumer then has to handle.
  nullR_h2f_lo = round(quantile(rotR[,1], .025, na.rm = TRUE)),
  nullR_h2f_hi = round(quantile(rotR[,1], .975, na.rm = TRUE)),
  nullR_f2h_lo = round(quantile(rotR[,2], .025, na.rm = TRUE)),
  nullR_f2h_hi = round(quantile(rotR[,2], .975, na.rm = TRUE)),
  primary_null = "R (focus relocation)",
  step12_fov_paired_p = if (nrow(fv) >= 5)
      round(suppressWarnings(wilcox.test(fv$hot, fv$bg, paired = TRUE))$p.value, 4) else NA_real_,
  step12_bin_p_reported = 0.2333))
pk <- is.finite(permP[,1])
sc("arrangement_null_draws.csv",
   data.table(null = rep(c("R_relocation","Rd_relocation_defect_matched","P_permutation"),
                         c(nrow(rotR), nrow(rotD), sum(pk))),
              h2f_um      = c(rotR[,1], rotD[,1], permP[pk,1]),
              f2h_um      = c(rotR[,2], rotD[,2], permP[pk,2]),
              mnn_self_um = c(rotR[,3], rotD[,3], permP[pk,3])))
writeLines(format(Sys.time()), file.path(OUT, "arrangement.done"))
cat("\nOK\n")
