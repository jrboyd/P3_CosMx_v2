#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------------------
## step 21b -- DOES FL TISSUE RESPOND TO THE HODGKIN FIELD?  Within IB1, within cell type.
##
## THE QUESTION, AND WHY IT IS THE ONE WORTH ASKING. Step 13 concluded "no spatial
## relationship between the two diseases" on a single bin-level Wilcoxon over 11 foci
## (median hotspot-to-FL 1035um vs 1205um, p = 0.23 -- pseudoreplicated, since bins within
## an FOV are not independent, and the direction actually favoured proximity). Step 17 then
## named the design that would settle it: "a within-IB1 design comparing FL tissue near and
## far from the Hodgkin field, which needs no cross-sample subtraction". This is that design.
##
## IT IS THE ONLY REMAINING CROSS-LESION DESIGN THAT IS IDENTIFIABLE. Everything else is
## foreclosed: the IB1-vs-RB1 FL contrast is not identifiable (2 of 89 bins find a
## composition-matched partner), and for the exhaustion/chemokine genes the between-slide
## offset exceeds every within-slide regional effect (reference-zone ORs to 5.33 against a
## maximum within-slide 1.93). This contrast is within one sample, one section, one cell
## type at a time -- no slide term, no composition term, no cross-sample subtraction.
##
## THE AXIS IS SET BY 21a, NOT CHOSEN HERE. 21a tested both candidates and rejected
## distance-to-hotspot: its positive control fails (CCL17+CCL22 inside FL tissue is 23.3
## cp10k near vs 26.5 far) and field exposure is uncorrelated with it across the FL FOVs
## (rho = -0.126, p = 0.575). The cause is that the hotspot set is the top 1% of chemokine
## bins and does not enumerate the field's sources -- step 17 flagged exactly this -- and
## only 47 of the 80 hotspot bins survive the low-information filter. So the axis is the
## bin's OWN field level, CCL17+CCL22 cp10k, which spans 39x across the 22 FL FOVs.
##
## THE PRIMARY STATISTIC IS AN INTERACTION, NOT A SLOPE. A gene rising with field exposure
## inside FL tissue is not yet an FL-specific result: 1,723 of 6,175 genes are
## field-associated tissue-wide (16_gradient), so the generic field response has to be
## subtracted. So for each gene x cell type:
##
##     b_FL   slope of expression on log2(field exposure), among IB1 FL_focus bins
##     b_BG   the same slope among IB1 background bins
##     delta  b_FL - b_BG   <- THE RESULT: the field response specific to FL tissue
##
## Both arms are the same sample, same section, same cell type, and depth enters each fit
## as a covariate. delta = 0 means FL tissue responds to the field exactly as ordinary IB1
## tissue does, which is the null this step exists to test.
##
## CELL-TYPE CONDITIONING IS MANDATORY HERE, NOT A REFINEMENT. 21a measured field exposure
## and B share to be collinear inside FL tissue (fb_a SMD 1.41, p = 0.013), which is the
## already-known coupling in CLAUDE.md fact 4. A bin-level contrast on this axis would
## therefore re-measure B share -- precisely how step 16's C3 contrast failed, returning
## RB1's interfollicular markers as "up in IB1 FL". Conditioning on cell type removes the
## composition term by construction rather than by matching.
##
## THREE GUARDS, EACH FOR A FAILURE THIS PROJECT HAS ALREADY MADE.
##   1. FOV IS THE RESAMPLING UNIT (stratified cluster bootstrap, FL and background arms
##      separately). 83 FL bins are not 83 independent observations; they sit in 22 FOVs.
##      Step 11 lost a Tfh result to treating cells as independent and step 17 lost the
##      sign-consistency result to treating compartments as independent (ICC 0.875).
##   2. cp10k AND DETECTION RATE MUST AGREE IN DIRECTION. cp10k is compositional and fakes
##      depression in the panel's most abundant genes (CLAUDE.md fact 8); the axis gene
##      CCL17 sits at the 95th prevalence percentile, so its rise mechanically depresses
##      others. Discordant genes are reported and not interpreted.
##   3. THE MINIMUM DETECTABLE EFFECT IS REPORTED BESIDE EVERY NULL. The FL arm has 22 FOVs.
##      A null here may be a power failure, and step 17's 0-of-30 taught that "did not
##      survive" covers two different situations that must not be read the same way.
##
## TWO P-VALUES, BECAUSE A PERCENTILE BOOTSTRAP P HAS A FLOOR AND BH DOES NOT KNOW THAT.
## The first version of this script scored the percentile p as 2*min(mean(d<=0), mean(d>=0)),
## which returns EXACTLY ZERO when every resample falls one side of zero. BH then treats
## zero as infinitely significant, and the whole genome-wide hit list became precisely the
## set of genes at that floor -- 143 "hits" across three cell types, every one an artifact of
## the p-value rather than a result. So:
##
##   p_perc   percentile, with the standard (1 + count)/(1 + NBOOT) correction so it can
##            never be 0. Distribution-free but floored at 2/(NBOOT+1), which for a
##            genome-wide BH over ~6,000 genes can never reach q < 0.05. That is a real
##            resolution limit and is reported rather than hidden.
##   p_norm   Wald on the bootstrap SD, 2*pnorm(-|delta|/sd). Not floor-limited, so it can
##            support a genome-wide correction, at the cost of assuming the resampling
##            distribution is roughly normal.
##
## A hit must satisfy ALL THREE: q_norm < 0.05, the percentile interval excluding 0, and
## cp10k/detection concordance. Reporting either p alone would be misleading in one
## direction or the other.
##
## THE OVER-FIXATION DEFECT IS ADJUSTED FOR, NOT MERELY EXCLUDED. Defect cells are already
## dropped per cell upstream (16_script_regions.R:166, 16_script_celltype_pseudobulk.R:51),
## so no defect transcriptome enters a slope. But 21e measured that removing the cells does
## not leave the surviving AREA uniform, and both of the things that could distort this
## contrast are present:
##   - the two arms are differently depleted -- FL bins lost 33.8% of their cells to the
##     defect against background bins' 47.4% (per-bin fraction 0.330 vs 0.450, SMD -0.60,
##     p = 1.5e-8);
##   - the exposure axis partly tracks defect burden inside FL tissue (bin-level
##     rho = 0.412, p = 1.3e-4; FOV-level 0.379, p = 0.083; SMD 0.85 across the low/high
##     exposure split).
## The second is not surprising: qc_11_flattening_vs_chemokine.R already found the discarded
## cells to be chemokine-ENRICHED (CCL17 1.20x, CCL22 1.18x). So field exposure and defect
## burden are genuinely coupled in this tissue.
##
## Therefore every fit carries the bin's DEFECT FRACTION as a covariate alongside depth, and
## the unadjusted estimate is reported beside it so the size of the adjustment is visible
## rather than asserted. This follows the project's rule that a result which could be
## circular is reported twice, with and without, and only the second is evidence.
##
## CCL17 and CCL22 are excluded from every outcome: they define the axis.
##
## Usage: sbatch --mem=64G --time=2:00:00 scripts/sbatch/run_rscript.sh \
##          scripts/21_script_interface_response.R [tier]
## ---------------------------------------------------------------------------------------
suppressPackageStartupMessages({library(Matrix); library(data.table)})
setDTthreads(1); options(width = 200)

AG    <- commandArgs(TRUE)
TIER  <- if (length(AG) >= 1) AG[1] else "matched"
## SMOKE MODE. This job costs ~25 minutes and three of its early runs died on typos in the
## LAST few expressions -- a stale column name in a write -- discarding all of the compute.
## `smoke` runs the identical code path with a tiny bootstrap and a gene subsample, into a
## separate directory so it cannot clobber a real result. Run it before every real submit.
SMOKE <- length(AG) >= 2 && AG[2] == "smoke"
PB   <- file.path("output/16_ct_pseudobulk", TIER)
GEO  <- file.path("output/21_interface", TIER)
OUT  <- if (SMOKE) file.path(GEO, "smoke") else GEO
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
sc <- function(f, x) { fwrite(x, file.path(OUT, f)); cat("wrote", file.path(OUT, f), "\n") }

AXIS_GENES <- c("CCL17","CCL22")
CT_TESTED  <- c("B","T","myeloid")     # the three clearing step 16's own 5-cell bin floor
CT_DESCR   <- c("fibro","endo","NK","PC")
MIN_CELLS  <- 5L                       # per (bin x celltype) group, as 16_script_celltype_de.R
MIN_GROUPS <- 25L                      # per arm, per cell type
MIN_DET    <- 0.02                     # a gene detected in <2% of a type's cells cannot slope
NBOOT      <- 4000L
NBOOT_DESC <- 500L                     # descriptive-only types get a coarse CI, not a test
if (SMOKE) { NBOOT <- 30L; NBOOT_DESC <- 10L }
Q_CUT      <- 0.05
SEED       <- 20L

for (f in c("pb_binct_sum.Rds","pb_binct_det.Rds","binct_meta.csv.gz"))
  if (!file.exists(file.path(PB, f))) stop("missing ", file.path(PB, f), call. = FALSE)
if (!file.exists(file.path(GEO, "bin_axes.csv.gz")))
  stop("run 21_script_interface_geometry.R first (no bin_axes.csv.gz)", call. = FALSE)
GV <- fread(file.path(GEO, "gate_verdict.csv"))
if (!grepl("^B ", GV$axis_selected))
  stop("21a did not select axis B; refusing to proceed on a rejected axis: ",
       GV$axis_selected, call. = FALSE)

cat("========================================================================\n")
if (SMOKE) cat("*** SMOKE MODE: tiny bootstrap, gene subsample, output to ", OUT, " ***\n", sep="")
cat("step 21b -- FL tissue's response to the Hodgkin field, within cell type\n")
cat("            tier =", TIER, " axis =", GV$axis_selected, "\n")
cat("========================================================================\n")

## ---- the exposure axis, per bin --------------------------------------------------------
AX <- fread(file.path(GEO, "bin_axes.csv.gz"))
AX[, bin := paste0("IB1_", bx, "_", by_)]
## bin_region is carried by BOTH tables. Keep the pseudobulk's copy as authoritative and
## bring this one along under a different name purely to assert they agree -- a silent
## bin_region.x / bin_region.y from the merge is what broke the first run.
AX <- AX[bin_region %in% c("FL_focus","background"),
         .(bin, fov, bin_region_ax = bin_region, exposure = chemo_cp10k, fb_a, f_assigned)]
cat(sprintf("\nbins on the axis: %d FL_focus, %d background\n",
            AX[bin_region_ax=="FL_focus", .N], AX[bin_region_ax=="background", .N]))

## ---- per-bin over-fixation burden, from 21e ---------------------------------------------
DF_F <- file.path(GEO, "defect_per_bin.csv.gz")
if (!file.exists(DF_F))
  stop("run 21_script_interface_defect_qc.R first (no defect_per_bin.csv.gz)", call. = FALSE)
DFB <- fread(DF_F)[, .(bin = paste0("IB1_", bx, "_", by_), defect_frac)]
cat(sprintf("per-bin defect burden available for %d bins\n", nrow(DFB)))

## ---- pseudobulk ------------------------------------------------------------------------
S  <- readRDS(file.path(PB, "pb_binct_sum.Rds"))
Dt <- readRDS(file.path(PB, "pb_binct_det.Rds"))
GB <- fread(file.path(PB, "binct_meta.csv.gz"))
setkey(GB, grp_bin); GB <- GB[colnames(S)]
stopifnot(identical(GB$grp_bin, colnames(S)), identical(colnames(S), colnames(Dt)))
GENES <- rownames(S)
cat(sprintf("pseudobulk: %d genes x %d (bin x celltype) groups\n", length(GENES), ncol(S)))

GB[, row := .I]
GB <- merge(GB[samp == "IB1" & cells >= MIN_CELLS], AX, by = c("bin","fov"), all = FALSE)
GB <- merge(GB, DFB, by = "bin", all.x = TRUE)
if (GB[is.na(defect_frac), .N]) {
  cat(sprintf("  %d groups have no defect burden and are dropped\n", GB[is.na(defect_frac), .N]))
  GB <- GB[!is.na(defect_frac)]
}
if (GB[bin_region != bin_region_ax, .N])
  stop("bin_region disagrees between the pseudobulk metadata and the axis table for ",
       GB[bin_region != bin_region_ax, .N], " groups", call. = FALSE)
GB[, bin_region_ax := NULL]
cat(sprintf("IB1 groups on the axis with >=%d cells: %d (bin_region agrees in both tables)\n",
            MIN_CELLS, nrow(GB)))
print(GB[, .(groups = .N, bins = uniqueN(bin), fovs = uniqueN(fov)),
         by = .(coarse, bin_region)][order(coarse, bin_region)])

## ---- least squares for every gene at once ----------------------------------------------
## beta = (X'X)^-1 X'Y' with X = [1, log2(exposure+1), log depth]. One 3x3 solve and one
## matrix product per arm, so a 1000-rep bootstrap over 6,175 genes is cheap. The exposure
## coefficient is column 2.
## `w` is the optional over-fixation covariate; pass NULL for the unadjusted fit. The
## exposure coefficient is column 2 either way.
slopes <- function(Y, x, z, w = NULL) {
  X <- if (is.null(w)) cbind(1, x, z) else cbind(1, x, z, w)
  ok <- apply(X, 1, function(r) all(is.finite(r)))
  if (sum(ok) < 8L) return(rep(NA_real_, nrow(Y)))
  X <- X[ok, , drop = FALSE]; Y <- Y[, ok, drop = FALSE]
  if (qr(X)$rank < ncol(X)) return(rep(NA_real_, nrow(Y)))
  as.numeric(tcrossprod(Y, solve(crossprod(X), t(X)))[, 2])
}

res_all <- list(); comp_all <- list()
set.seed(SEED)

for (CT in c(CT_TESTED, CT_DESCR)) {
  g <- GB[coarse == CT]
  nFL <- g[bin_region == "FL_focus", .N]; nBG <- g[bin_region == "background", .N]
  fFL <- g[bin_region == "FL_focus", uniqueN(fov)]
  testable <- CT %in% CT_TESTED && nFL >= MIN_GROUPS && nBG >= MIN_GROUPS
  cat(sprintf("\n--- %s: %d FL groups (%d FOVs), %d background groups (%d FOVs) -> %s\n",
              CT, nFL, fFL, nBG, g[bin_region=="background", uniqueN(fov)],
              if (testable) "TESTED" else "descriptive only"))
  if (nFL < 8L || nBG < 8L) { cat("    too few groups for even a descriptive slope\n"); next }

  idx <- g$row
  cp  <- as.matrix(S[, idx, drop = FALSE] %*% Diagonal(x = 1e4/pmax(g$nCount, 1)))
  det <- as.matrix(Dt[, idx, drop = FALSE] %*% Diagonal(x = 1/pmax(g$cells, 1)))
  keep <- rowMeans(det) >= MIN_DET & !(GENES %in% AXIS_GENES)
  if (SMOKE) {  # keep every Tier-1 gene plus a sample of the rest, so both tiers are exercised
    t1g <- unique(fread("output/16_regions/panel_groups.csv")[
             group %in% c("checkpoint_primary","checkpoint_context","chemokine_axis",
                          "t_treg","fl_prog","chl_prog"), gene])
    set.seed(1); sub <- keep & (GENES %in% t1g | seq_along(GENES) %in% sample(length(GENES), 200))
    keep <- sub
  }
  cat(sprintf("    genes with pooled detection >= %.0f%%: %d of %d\n",
              100*MIN_DET, sum(keep), length(GENES)))
  Ycp <- log2(cp[keep, , drop = FALSE] + 1); Ydt <- det[keep, , drop = FALSE]
  gn  <- GENES[keep]
  x <- log2(g$exposure + 1); z <- log(g$nCount/g$cells); w <- g$defect_frac
  isFL <- g$bin_region == "FL_focus"

  ## PRIMARY: adjusted for the bin's over-fixation burden.
  b_fl_cp <- slopes(Ycp[, isFL, drop=FALSE],  x[isFL],  z[isFL],  w[isFL])
  b_bg_cp <- slopes(Ycp[, !isFL, drop=FALSE], x[!isFL], z[!isFL], w[!isFL])
  b_fl_dt <- slopes(Ydt[, isFL, drop=FALSE],  x[isFL],  z[isFL],  w[isFL])
  b_bg_dt <- slopes(Ydt[, !isFL, drop=FALSE], x[!isFL], z[!isFL], w[!isFL])
  d_cp <- b_fl_cp - b_bg_cp; d_dt <- b_fl_dt - b_bg_dt

  ## COMPARISON: the same thing with no defect covariate, so the size of the adjustment is
  ## visible. No interval -- it exists only to be differenced against the primary.
  d_cp_unadj <- slopes(Ycp[, isFL, drop=FALSE],  x[isFL],  z[isFL]) -
                slopes(Ycp[, !isFL, drop=FALSE], x[!isFL], z[!isFL])

  ## Stratified FOV cluster bootstrap: resample FL FOVs and background FOVs separately, so
  ## each arm keeps its own size and the interval reflects between-FOV variation rather
  ## than between-bin variation.
  ## Precompute FOV -> row indices once. The first version called which() per FOV per rep,
  ## which is 553 scans of a 4,700-row vector x NBOOT and dominated the runtime.
  ifl <- split(which(isFL),  g$fov[isFL])
  ibg <- split(which(!isFL), g$fov[!isFL])
  fov_fl <- names(ifl); fov_bg <- names(ibg)
  ## Descriptive-only types (below the 5-cell/25-group floor) get no q-value, so a full
  ## 4,000-rep bootstrap on them is wasted: fibro alone carries 4,141 background groups and
  ## dominated the runtime of the first full run. They get a coarse interval instead.
  NB <- if (testable) NBOOT else NBOOT_DESC
  BS_d <- matrix(NA_real_, nrow(Ycp), NB); BS_fl <- matrix(NA_real_, nrow(Ycp), NB)
  for (b in seq_len(NB)) {
    kf <- unlist(ifl[sample(fov_fl, length(fov_fl), TRUE)], use.names = FALSE)
    kb <- unlist(ibg[sample(fov_bg, length(fov_bg), TRUE)], use.names = FALSE)
    if (length(kf) < 8L || length(kb) < 8L) next
    bf <- slopes(Ycp[, kf, drop=FALSE], x[kf], z[kf], w[kf])
    bb <- slopes(Ycp[, kb, drop=FALSE], x[kb], z[kb], w[kb])
    BS_fl[, b] <- bf; BS_d[, b] <- bf - bb
  }
  qs   <- t(apply(BS_d,  1, quantile, c(.025,.975), na.rm = TRUE))
  sdd  <- apply(BS_d,  1, sd, na.rm = TRUE)
  ## (1 + count)/(1 + n) so the percentile p can never be exactly 0 -- the same form as
  ## 21c's pR(). Without it BH reads the floor as infinite significance.
  ppc  <- function(v) { v <- v[is.finite(v)]; n <- length(v)
    if (!n) return(NA_real_)
    2*min((1 + sum(v <= 0))/(1 + n), (1 + sum(v >= 0))/(1 + n)) }
  pd   <- apply(BS_d,  1, ppc)
  pfl  <- apply(BS_fl, 1, ppc)
  sdfl <- apply(BS_fl, 1, sd, na.rm = TRUE)

  res_all[[CT]] <- data.table(
    celltype = CT, gene = gn, testable = testable,
    n_fl_groups = nFL, n_fl_fovs = fFL, n_bg_groups = nBG,
    b_fl = b_fl_cp, b_bg = b_bg_cp, delta = d_cp,
    ci_lo = qs[,1], ci_hi = qs[,2], boot_sd = sdd,
    p_perc = pd, p_norm = 2*pnorm(-abs(d_cp)/sdd),
    n_boot = NB, p_perc_floor = 2/(NB+1),
    ci_excludes_0 = is.finite(qs[,1]) & is.finite(qs[,2]) & (qs[,1] > 0 | qs[,2] < 0),
    mde_log2 = 1.96*sdd, mde_fold = 2^(1.96*sdd),
    b_fl_det = b_fl_dt, b_bg_det = b_bg_dt, delta_det = d_dt,
    delta_unadj = d_cp_unadj,
    defect_shift = d_cp - d_cp_unadj,
    p_fl_perc = pfl, p_fl_norm = 2*pnorm(-abs(b_fl_cp)/sdfl),
    concordant = is.finite(d_cp) & is.finite(d_dt) & sign(d_cp) == sign(d_dt))

  ## Composition on the same axis, same unit. Reported as CONFIRMATORY: the field/B-density
  ## coupling in IB1 is already recorded (CLAUDE.md fact 4), and 21a measured it here at
  ## SMD 1.41. Included so the compositional and cell-state legs sit side by side, which is
  ## the distinction step 16 established as the important one for this tissue.
  cf <- g[, .(cells = sum(cells)), by = .(fov, bin, bin_region, exposure)]
  comp_all[[CT]] <- cf[, .(celltype = CT, bins = .N, fovs = uniqueN(fov),
        rho_bin = suppressWarnings(cor(cells, exposure, method = "spearman")),
        rho_fov = { v <- .SD[, .(c = mean(cells), e = mean(exposure)), by = fov]
                    if (nrow(v) < 6) NA_real_ else
                      suppressWarnings(cor(v$c, v$e, method = "spearman")) }),
     by = bin_region]
}

R <- rbindlist(res_all, fill = TRUE)
CMP <- rbindlist(comp_all, fill = TRUE)

## ---- gene tiers and multiplicity -------------------------------------------------------
## Tier 1 is pre-declared and small, so it carries the power. Tier 2 is the genome-wide
## scan, corrected across everything and expected to be null at 22 FL FOVs -- it is run so
## that the FL-restricted version of 16_gradient's scan exists, since that scan covered all
## non-hotspot bins and never the FL tissue alone.
PG <- fread("output/16_regions/panel_groups.csv")
T1_SETS <- c("checkpoint_primary","checkpoint_context","chemokine_axis",
             "t_treg","fl_prog","chl_prog")
T1 <- unique(PG[group %in% T1_SETS, .(gene, set = group)])
R <- merge(R, T1, by = "gene", all.x = TRUE)
R[, tier := ifelse(!is.na(set), "T1_declared", "T2_genomewide")]
R[testable == TRUE, q_delta := p.adjust(p_norm,    "BH"), by = .(tier, celltype)]
R[testable == TRUE, q_perc  := p.adjust(p_perc,    "BH"), by = .(tier, celltype)]
R[testable == TRUE, q_fl    := p.adjust(p_fl_norm, "BH"), by = .(tier, celltype)]
R[, verdict := ifelse(!testable, "descriptive",
              ifelse(!is.finite(q_delta), "unestimable",
              ifelse(q_delta < Q_CUT & ci_excludes_0 & concordant, "FL-SPECIFIC FIELD RESPONSE",
              ifelse(q_delta < Q_CUT & ci_excludes_0 & !concordant, "discordant, not interpreted",
              ifelse(q_delta < Q_CUT, "q only; percentile interval spans 0",
                     "not significant")))))]

cat("\n\n========================================================================\n")
cat("TIER 1 -- pre-declared sets, FL-specific field response (delta)\n")
cat("========================================================================\n")
for (CT in CT_TESTED) {
  cat(sprintf("\n--- %s ---\n", CT))
  print(R[celltype == CT & tier == "T1_declared"][order(p_norm)][1:12,
    .(gene, set, b_fl = round(b_fl,3), b_bg = round(b_bg,3), delta = round(delta,3),
      ci = sprintf("[%.2f, %.2f]", ci_lo, ci_hi),
      d_unadj = round(delta_unadj,3),
      p_norm = signif(p_norm,3), p_perc = signif(p_perc,3), q = signif(q_delta,3),
      mde_fold = round(mde_fold,2), conc = concordant, verdict)])
}
cat("\n=== TIER 1 hits at q < 0.05 and concordant, all cell types\n")
h1 <- R[tier == "T1_declared" & verdict == "FL-SPECIFIC FIELD RESPONSE"]
if (nrow(h1)) print(h1[order(p_norm), .(celltype, gene, set, delta = round(delta,3),
      ci = sprintf("[%.2f, %.2f]", ci_lo, ci_hi), q_norm = signif(q_delta,3),
      p_perc = signif(p_perc,3), q_perc = signif(q_perc,3))]) else
  cat("  none.\n")

cat("\n=== TIER 2 genome-wide hits at q < 0.05 and concordant\n")
h2 <- R[tier == "T2_genomewide" & verdict == "FL-SPECIFIC FIELD RESPONSE"]
if (nrow(h2)) print(h2[order(p_norm)][1:min(25,.N),
      .(celltype, gene, delta = round(delta,3), q_norm = signif(q_delta,3),
        p_perc = signif(p_perc,3))]) else
  cat("  none.\n")

cat("\n=== POWER: what a null here does and does not exclude\n")
pw <- R[testable == TRUE, .(genes = .N,
        med_mde_fold = round(median(mde_fold, na.rm = TRUE), 2),
        med_mde_log2 = round(median(mde_log2, na.rm = TRUE), 3),
        max_abs_delta = round(max(abs(delta), na.rm = TRUE), 2)),
      by = .(celltype, tier)][order(celltype, tier)]
print(pw)
cat("\n  Read this beside every null above: a |delta| smaller than med_mde_log2 could not\n")
cat("  have been resolved, so its absence is not evidence of absence.\n")
cat(sprintf("\n  The percentile bootstrap p is floored at %.4f (2/(NBOOT+1), NBOOT = %d). Under a\n",
            2/(NBOOT+1), NBOOT))
cat("  genome-wide BH over ~6,000 genes that floor alone gives q ~ 12, so TIER 2 CANNOT\n")
cat("  REACH q < 0.05 ON THE PERCENTILE P AT ANY EFFECT SIZE. Tier 2 is therefore judged on\n")
cat("  the Wald p, with the percentile interval required to exclude 0 as the check on it.\n")
cat("  Tier 1, at 44 tests, is not floor-limited in the same way.\n")

cat("\n=== the FL field response on its own (b_fl vs 0) -- the FL-restricted version of\n")
cat("    16_gradient's genome-wide scan, which was only ever run on all non-hotspot bins\n")
for (CT in CT_TESTED) {
  s <- R[celltype == CT & testable == TRUE & is.finite(q_fl)]
  cat(sprintf("  %-8s %d genes, %d with q_fl < 0.05 (%d up, %d down)\n", CT, nrow(s),
              s[q_fl < Q_CUT, .N], s[q_fl < Q_CUT & b_fl > 0, .N],
              s[q_fl < Q_CUT & b_fl < 0, .N]))
}

cat("\n=== how much did adjusting for over-fixation burden move delta?\n")
ds <- R[testable == TRUE, .(genes = .N,
        med_abs_shift = round(median(abs(defect_shift), na.rm = TRUE), 4),
        p90_abs_shift = round(quantile(abs(defect_shift), .9, na.rm = TRUE), 4),
        med_abs_delta = round(median(abs(delta), na.rm = TRUE), 4),
        pct_sign_flip = round(100*mean(sign(delta) != sign(delta_unadj), na.rm = TRUE), 1)),
      by = .(celltype, tier)][order(celltype, tier)]
print(ds); sc("defect_adjustment_shift.csv", ds)
cat("  A shift small against med_abs_delta means the defect was not driving the contrast.\n")

cat("\n=== COMPOSITION on the same axis (confirmatory -- see CLAUDE.md fact 4)\n")
print(CMP[order(bin_region, -abs(rho_fov))][, .(celltype, bin_region, bins, fovs,
      rho_bin = round(rho_bin,3), rho_fov = round(rho_fov,3))])

sc("interface_response.csv.gz", R)
## SUPERSEDED: this block regresses per-bin cell COUNT on exposure, which measures the bin's
## total rather than its composition (total cells per bin falls with exposure inside FL tissue
## at rho -0.34). Replaced by 21_script_interface_composition.R, which uses the SHARE and
## adjusts for defect burden. Kept under a name that cannot be mistaken for the live result,
## following the project's *_SUPERSEDED convention.
sc("interface_composition_COUNT_SUPERSEDED.csv", CMP)
sc("interface_power.csv", pw)
sc("interface_tier1.csv", R[tier == "T1_declared", ][order(celltype, p_norm)])
writeLines(format(Sys.time()), file.path(OUT, "interface_response.done"))
cat("\nOK\n")
