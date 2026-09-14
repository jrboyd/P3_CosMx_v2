#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------------------
## Does the Hodgkin-field negative control survive a between-slide offset?
##
## THE GAP THIS FILLS. The project's headline Hodgkin result rests on a cross-sample
## negative control: IB1's chemokine hotspots are defined as bins clearing IB1's own q99 of
## CCL17+CCL22, and **zero RB1 bins clear that cut** (step 11, restated in steps 13 and 16).
## Step 17 then showed the between-slide offset for other genes reaches OR 5.33 in tissue
## that shares no disease-region membership. If CCL17/CCL22 carry a comparable offset, then
## IB1's q99 is an inflated bar and RB1's failure to clear it is guaranteed by the slide, not
## by the absence of Hodgkin. That check has never been run: every CCL17/CCL22 contrast in
## the project is within-IB1 (verified against 16_de/de_contrasts.csv.gz and
## 16_ct_de/ct_de.csv.gz -- no cross-sample reference-zone contrast exists for these genes).
##
## TWO TESTS, and the second is the one that matters.
##   1. Is there an offset at all? Compare the two samples' REFERENCE tissue (IB1 background
##      vs RB1 interfollicular) on chemo_cp10k, FOV as the unit.
##   2. Is the hotspot enrichment offset-free? Compare each sample's own top-1% bins to its
##      OWN background as a fold enrichment, then compare those folds across samples. A
##      within-sample fold cancels any multiplicative slide offset exactly, so this is the
##      offset-free version of "hotspots are Hodgkin-specific and not just the top 1% of
##      anything".
##
## Runs off bins_all_scales.csv.gz only.
## ---------------------------------------------------------------------------------------
suppressMessages(library(data.table))
options(width = 200)

TIER <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "matched"
OUT  <- file.path("output/17_pd1", TIER)
B <- fread(file.path("output/16_regions", TIER, "bins_all_scales.csv.gz"))
B <- B[bin_um == 100 & bin_low_info == FALSE]
sc <- function(f, x) { fwrite(x, file.path(OUT, f)); cat("wrote", file.path(OUT, f), "\n") }

cat("bins at 100um, informative:\n"); print(B[, .N, by = samp])
cat("\nbin_region levels:\n"); print(B[, .N, by = .(samp, bin_region)][order(samp, -N)])

## Reference tissue: neither hotspot nor FL in either sample. Use each sample's own
## non-disease region, which is what "background" and "interfollicular" already mean.
REF <- c(IB1 = "background", RB1 = "interfollicular")
B[, is_ref := bin_region == REF[samp]]
cat(sprintf("\nreference bins: IB1 %d, RB1 %d\n",
            B[samp == "IB1" & is_ref, .N], B[samp == "RB1" & is_ref, .N]))

## ---- TEST 1: is there an offset in the chemokine readout itself? ----------------------
## FOV as the unit, because bins within a field share its staining batch -- the same rule
## the rest of the project uses.
fov <- B[is_ref == TRUE, .(chemo = sum(chemo) / sum(nCount) * 1e4, bins = .N),
         by = .(samp, fov)][bins >= 3]
w  <- wilcox.test(chemo ~ samp, data = fov)
md <- fov[, .(n_fov = .N, median_chemo_cp10k = median(chemo),
              mean_chemo_cp10k = mean(chemo)), by = samp][order(samp)]
ratio <- md[samp == "IB1", median_chemo_cp10k] / md[samp == "RB1", median_chemo_cp10k]
cat("\n============================================================\n")
cat("TEST 1: CCL17+CCL22 in REFERENCE tissue -- is there a slide offset?\n")
cat("============================================================\n")
print(md)
cat(sprintf("\nIB1 : RB1 ratio in reference tissue = %.3f   (Mann-Whitney p = %.3g, FOV unit)\n",
            ratio, w$p.value))

## ---- TEST 2: within-sample fold enrichment of each sample's OWN top 1% ----------------
## A multiplicative slide offset cancels exactly in a within-sample ratio, so comparing the
## two folds is the offset-free test.
tops <- rbindlist(lapply(c("IB1", "RB1"), function(s) {
  z   <- B[samp == s]
  cut <- quantile(z$chemo_cp10k, .99, na.rm = TRUE)
  bg  <- z[is_ref == TRUE, sum(chemo) / sum(nCount) * 1e4]
  tp  <- z[chemo_cp10k >= cut, sum(chemo) / sum(nCount) * 1e4]
  ## FOV-level bootstrap of the fold, so the interval respects the batch structure.
  fv  <- unique(z[chemo_cp10k >= cut | is_ref == TRUE, fov])
  bs <- vapply(seq_len(2000), function(i) {
    f <- sample(fv, length(fv), TRUE)
    zz <- z[data.table(fov = f), on = "fov", allow.cartesian = TRUE, nomatch = 0L]
    b <- zz[is_ref == TRUE, sum(chemo) / sum(nCount)]
    t <- zz[chemo_cp10k >= cut, sum(chemo) / sum(nCount)]
    if (!is.finite(b) || b <= 0 || !is.finite(t)) return(NA_real_)
    log(t / b)
  }, 0)
  bs <- bs[is.finite(bs)]
  data.table(samp = s, own_q99_cut = cut, n_top_bins = z[chemo_cp10k >= cut, .N],
             bg_cp10k = bg, top_cp10k = tp, fold = tp / bg,
             lo = exp(quantile(bs, .025)), hi = exp(quantile(bs, .975)), boot_ok = length(bs))
}))
cat("\n============================================================\n")
cat("TEST 2: each sample's own top 1% vs its OWN background (offset-free)\n")
cat("============================================================\n")
print(tops[, .(samp, own_q99_cut = round(own_q99_cut, 2), n_top_bins,
               bg_cp10k = round(bg_cp10k, 2), top_cp10k = round(top_cp10k, 2),
               fold = round(fold, 2), ci = sprintf("[%.2f, %.2f]", lo, hi))])
fi <- tops[samp == "IB1", fold]; fr <- tops[samp == "RB1", fold]
cat(sprintf("\nfold ratio IB1 / RB1 = %.2f\n", fi / fr))

## Does IB1's q99 cut sit above RB1's entire range, and how much of that is the offset?
ib_cut <- tops[samp == "IB1", own_q99_cut]
cat(sprintf("\nRB1 bins clearing IB1's q99 cut (%.2f cp10k): %d of %d\n",
            ib_cut, B[samp == "RB1" & chemo_cp10k >= ib_cut, .N], B[samp == "RB1", .N]))
cat(sprintf("RB1 bins clearing IB1's cut DEFLATED by the reference ratio (%.2f -> %.2f): %d\n",
            ib_cut, ib_cut / ratio, B[samp == "RB1" & chemo_cp10k >= ib_cut / ratio, .N]))

sc("chemo_offset_check.csv",
   merge(tops[, .(samp, own_q99_cut = round(own_q99_cut, 3), n_top_bins,
                  bg_cp10k = round(bg_cp10k, 3), top_cp10k = round(top_cp10k, 3),
                  fold_over_own_bg = round(fold, 3),
                  ci = sprintf("[%.2f, %.2f]", lo, hi))],
         md[, .(samp, n_ref_fov = n_fov, ref_median_cp10k = round(median_chemo_cp10k, 3))],
         by = "samp"))

cat("\n=== verdict\n")
if (w$p.value < .05 && (ratio > 1.5 || ratio < 1/1.5)) {
  cat(sprintf("  An offset EXISTS in the chemokine readout (%.2fx, p = %.3g).\n", ratio, w$p.value))
} else {
  cat(sprintf("  No large offset in the chemokine readout (%.2fx, p = %.3g).\n", ratio, w$p.value))
}
cat(sprintf("  Offset-free comparison: IB1's top 1%% is %.2fx its own background, RB1's is\n", fi))
cat(sprintf("  %.2fx its own -- a ratio of %.2f. This number is what the Hodgkin-field\n", fr, fi/fr))
cat("  specificity claim should rest on, because it cancels any multiplicative slide effect.\n")
writeLines(format(Sys.time()), file.path(OUT, "chemo_offset_check.done"))
cat("\nOK\n")
