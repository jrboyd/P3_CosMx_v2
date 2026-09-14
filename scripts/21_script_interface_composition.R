#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------------------
## step 21f -- COMPOSITION against field exposure, done as a SHARE and defect-adjusted.
##
## WHY THIS EXISTS AS A SEPARATE SCRIPT. 21b carried a composition block that regressed each
## cell type's per-bin CELL COUNT on field exposure. That block is WITHDRAWN and replaced
## here, because the count is confounded by the bin's total: inside FL tissue total cells per
## bin falls steeply with exposure (Spearman -0.35 to -0.62 depending on the cell type's FOV
## set) against -0.07 to -0.23 in background tissue, so every count correlated negatively
## with exposure whether or not its SHARE moved. The diagnostic that settles it is endothelium,
## whose association REVERSES between the two definitions: rho -0.378 on counts and +0.236 on
## shares. A measure that flips sign when you change the denominator was measuring the
## denominator.
##
## WHAT THIS DOES INSTEAD.
##   outcome   share = the cell type's cells / all assigned cells in the same bin
##   model     share ~ log2(exposure) + log depth + DEFECT FRACTION, fitted separately in
##             IB1 FL bins and IB1 background bins
##   result    delta = b_FL - b_BG, the exposure-composition slope specific to FL tissue
##   inference stratified FOV cluster bootstrap, both arms resampled at FOV level
##
## The defect covariate is not optional here. 21e measured field exposure and over-fixation
## burden to be coupled inside FL tissue (bin-level rho 0.412), and the defect deletes cells,
## so it moves shares directly whenever it deletes cell types unevenly. Without the covariate
## this analysis would partly re-measure the fixation gradient.
##
## The unadjusted estimate is reported beside the adjusted one, as in 21b.
##
## Light: reads only binct_meta.csv.gz (581 KB) and the two small step-21 tables.
## Usage: Rscript scripts/21_script_interface_composition.R [tier]
## ---------------------------------------------------------------------------------------
suppressMessages({library(data.table)})
setDTthreads(1); options(width = 200)

TIER <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "matched"
PB   <- file.path("output/16_ct_pseudobulk", TIER)
OUT  <- file.path("output/21_interface", TIER)
sc <- function(f, x) { fwrite(x, file.path(OUT, f)); cat("wrote", file.path(OUT, f), "\n") }
NBOOT <- 4000L; MIN_CELLS <- 5L; MIN_FOV <- 8L; SEED <- 20L; Q_CUT <- 0.05

cat("========================================================================\n")
cat("step 21f -- composition vs field exposure, share-based, tier =", TIER, "\n")
cat("========================================================================\n")

GB <- fread(file.path(PB, "binct_meta.csv.gz"))[samp == "IB1" & cells >= MIN_CELLS]
AX <- fread(file.path(OUT, "bin_axes.csv.gz"))
AX[, bin := paste0("IB1_", bx, "_", by_)]
AX <- AX[bin_region %in% c("FL_focus","background"), .(bin, fov, exposure = chemo_cp10k)]
DFB <- fread(file.path(OUT, "defect_per_bin.csv.gz"))[, .(bin = paste0("IB1_", bx, "_", by_),
                                                          defect_frac)]
G <- merge(merge(GB, AX, by = c("bin","fov")), DFB, by = "bin")
G[, tot := sum(cells), by = bin][, share := cells/tot]
G[, `:=`(x = log2(exposure + 1), z = log(nCount/cells), w = defect_frac)]
cat(sprintf("\ngroups %d in %d bins, %d FOVs\n", nrow(G), uniqueN(G$bin), uniqueN(G$fov)))

## The confound, stated as a measurement so the withdrawal above is checkable.
cat("\n=== total assigned cells per bin against exposure (why counts could not be used)\n")
TT <- unique(G[, .(bin, fov, bin_region, tot, x)])
print(TT[, .(bins = .N, fovs = uniqueN(fov),
             rho_tot_vs_exposure = round(cor(tot, x, method = "spearman"), 3)),
         by = bin_region])

fitb <- function(d, adj = TRUE) {
  X <- if (adj) cbind(1, d$x, d$z, d$w) else cbind(1, d$x, d$z)
  ok <- apply(X, 1, function(r) all(is.finite(r))) & is.finite(d$share)
  if (sum(ok) < 8L) return(NA_real_)
  X <- X[ok, , drop = FALSE]; y <- d$share[ok]
  if (qr(X)$rank < ncol(X)) return(NA_real_)
  unname(solve(crossprod(X), crossprod(X, y))[2])
}

set.seed(SEED)
res <- rbindlist(lapply(sort(unique(G$coarse)), function(ct) {
  fl <- G[coarse == ct & bin_region == "FL_focus"]
  bg <- G[coarse == ct & bin_region == "background"]
  if (uniqueN(fl$fov) < MIN_FOV || uniqueN(bg$fov) < MIN_FOV) return(NULL)
  b_fl <- fitb(fl); b_bg <- fitb(bg); d <- b_fl - b_bg
  d_un <- fitb(fl, FALSE) - fitb(bg, FALSE)
  ifl <- split(seq_len(nrow(fl)), fl$fov); ibg <- split(seq_len(nrow(bg)), bg$fov)
  bs <- vapply(seq_len(NBOOT), function(i) {
    kf <- unlist(ifl[sample(names(ifl), length(ifl), TRUE)], use.names = FALSE)
    kb <- unlist(ibg[sample(names(ibg), length(ibg), TRUE)], use.names = FALSE)
    fitb(fl[kf]) - fitb(bg[kb])
  }, 0)
  bs <- bs[is.finite(bs)]
  data.table(celltype = ct, fl_bins = nrow(fl), fl_fovs = uniqueN(fl$fov),
             bg_bins = nrow(bg), bg_fovs = uniqueN(bg$fov),
             b_fl = b_fl, b_bg = b_bg, delta = d, delta_unadj = d_un,
             defect_shift = d - d_un,
             ci_lo = quantile(bs, .025), ci_hi = quantile(bs, .975),
             boot_sd = sd(bs),
             p_perc = 2*min((1+sum(bs <= 0))/(1+length(bs)), (1+sum(bs >= 0))/(1+length(bs))),
             p_norm = 2*pnorm(-abs(d)/sd(bs)),
             mde = 1.96*sd(bs))
}))
res[, q_norm := p.adjust(p_norm, "BH")]
res[, ci_excludes_0 := is.finite(ci_lo) & is.finite(ci_hi) & (ci_lo > 0 | ci_hi < 0)]
res[, verdict := ifelse(q_norm < Q_CUT & ci_excludes_0,
                        "FL-SPECIFIC COMPOSITION SHIFT", "not significant")]

cat("\n=== share of each cell type against field exposure: FL slope minus background slope\n")
print(res[order(p_norm), .(celltype, fl_fovs, bg_fovs,
      b_fl = signif(b_fl,3), b_bg = signif(b_bg,3), delta = signif(delta,3),
      d_unadj = signif(delta_unadj,3),
      ci = sprintf("[%.4f, %.4f]", ci_lo, ci_hi), p_norm = signif(p_norm,3),
      p_perc = signif(p_perc,3), q = signif(q_norm,3), mde = signif(mde,3), verdict)])
cat("\n  b_fl and b_bg are change in SHARE per doubling of CCL17+CCL22, depth and\n")
cat("  over-fixation burden held fixed. delta is the FL-specific part.\n")
sc("composition_share.csv", res)
writeLines(format(Sys.time()), file.path(OUT, "composition_share.done"))
cat("\nOK\n")
