#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------------------
## Is IB1's elevated background chemokine a HODGKIN FIELD or a SLIDE OFFSET?
##
## THE CONFLICT THIS RESOLVES. Two reports read the same tissue incompatibly:
##   step 13  -- ~27% of IB1 bins clear RB1's q99, called "a tissue-wide field, not a niche",
##               i.e. IB1's BACKGROUND carries Hodgkin chemokine. Biology.
##   step 17  -- IB1 background vs RB1 interfollicular differ 3.56x on CCL17+CCL22
##               (chemo_offset_check.csv), read as a between-slide offset. Technical.
## Both cannot be the primary explanation, and no report flags the clash.
##
## THE DISCRIMINATOR. A diffusible field must DECAY with distance from its source. A slide
## offset cannot -- it is flat in space by construction. So:
##   - chemokine declining with distance from the nearest hotspot  => field
##   - chemokine flat in distance                                  => offset
##   - declining but plateauing ABOVE RB1's level                  => both, and the plateau
##     is the offset while the decaying part is the field. This is the quantity that
##     actually settles it, so it is computed explicitly as the FAR-FIELD PLATEAU RATIO.
##
## THE CONTROL THAT MAKES IT MEAN ANYTHING. chemo_cp10k is spatially autocorrelated, so the
## top 1% of ANY such field has lower values around it -- a decline would appear in RB1 too,
## where there is no Hodgkin. So the identical construction is run in RB1 using RB1's OWN
## top-1% bins as pseudo-peaks. Only an IB1 gradient that exceeds RB1's is evidence.
##
## Entirely within-sample for the gradient, so immune to the offset it is testing.
## Runs off bins_all_scales.csv.gz. Region is the unit; FOV is the conservative second unit.
## ---------------------------------------------------------------------------------------
suppressMessages({library(data.table); library(ggplot2)})
options(width = 200)

TIER <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "matched"
OUT  <- file.path("output/17_pd1", TIER)
FIG  <- file.path(OUT, "figures"); dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
sc <- function(f, x) { fwrite(x, file.path(OUT, f)); cat("wrote", file.path(OUT, f), "\n") }

B <- fread(file.path("output/16_regions", TIER, "bins_all_scales.csv.gz"))
B <- B[bin_um == 100 & bin_low_info == FALSE]
stopifnot(nrow(B) > 0, all(c("x_mm","y_mm","chemo","nCount","chemo_cp10k","fov","samp") %in% names(B)))

## Peaks: each sample's OWN top 1%, so the construction is identical in both.
B[, q99 := quantile(chemo_cp10k, .99, na.rm = TRUE), by = samp]
B[, is_peak := chemo_cp10k >= q99]
cat("peaks (own top 1%):\n"); print(B[, .(bins = .N, peaks = sum(is_peak),
      q99_cut = round(q99[1], 2)), by = samp])

## Distance from every non-peak bin to the nearest peak, in um. A few thousand bins against
## a few dozen peaks, so the direct distance matrix is cheap and exact -- no need for nn2.
D <- rbindlist(lapply(unique(B$samp), function(s) {
  z <- B[samp == s]
  pk <- as.matrix(z[is_peak == TRUE, .(x_mm, y_mm)])
  np <- z[is_peak == FALSE]
  if (!nrow(pk) || !nrow(np)) return(NULL)
  xy <- as.matrix(np[, .(x_mm, y_mm)])
  dmin <- apply(xy, 1, function(p) sqrt(min((pk[,1]-p[1])^2 + (pk[,2]-p[2])^2)))
  cbind(np, d_um = dmin * 1000)
}))

## Rings. The first is dropped from the headline gradient: bins touching a peak share its
## diffusion halo and its segmentation, so an elevation there is not independent evidence.
BR <- c(0, 100, 200, 400, 800, 1600, 3200, Inf)
D[, ring := cut(d_um, BR, labels = c("0-100","100-200","200-400","400-800","800-1600",
                                     "1600-3200",">3200"), right = FALSE)]
prof <- D[, .(bins = .N, fovs = uniqueN(fov), med_d_um = round(median(d_um)),
              chemo_cp10k = sum(chemo) / sum(nCount) * 1e4), by = .(samp, ring)][order(samp, ring)]
cat("\n=== chemokine profile away from each sample's own top-1% peaks\n")
print(prof)

## ---- the gradient, per sample, excluding the adjacent ring -----------------------------
res <- rbindlist(lapply(unique(D$samp), function(s) {
  z <- D[samp == s & d_um >= 100]
  ct <- suppressWarnings(cor.test(z$chemo_cp10k, z$d_um, method = "spearman"))
  ## FOV-level means as the conservative unit: bins within a field share its batch.
  fv <- z[, .(chemo = sum(chemo)/sum(nCount)*1e4, d = median(d_um)), by = fov][, .N, by = fov][, .N]
  fvz <- z[, .(chemo = sum(chemo)/sum(nCount)*1e4, d = median(d_um)), by = fov]
  cf <- suppressWarnings(cor.test(fvz$chemo, fvz$d, method = "spearman"))
  data.table(samp = s, n_bins = nrow(z), rho_bin = ct$estimate, p_bin = ct$p.value,
             n_fov = nrow(fvz), rho_fov = cf$estimate, p_fov = cf$p.value)
}))
cat("\n=== gradient (Spearman, >=100um from a peak so the adjacent ring cannot carry it)\n")
print(res[, .(samp, n_bins, rho_bin = round(rho_bin, 3), p_bin = signif(p_bin, 3),
              n_fov, rho_fov = round(rho_fov, 3), p_fov = signif(p_fov, 3))])

## ---- THE DECIDING NUMBER: the far-field plateau ----------------------------------------
## If the whole 3.56x reference difference is a decaying field, IB1 should converge on RB1
## far from any peak. If a flat multiplicative offset remains out there, it is not the field.
far <- D[d_um >= 1600, .(bins = .N, fovs = uniqueN(fov),
                         chemo_cp10k = sum(chemo)/sum(nCount)*1e4), by = samp]
near <- D[d_um >= 100 & d_um < 400, .(chemo_cp10k = sum(chemo)/sum(nCount)*1e4), by = samp]
cat("\n=== near field (100-400um) vs far field (>1600um), each sample\n")
print(merge(near[, .(samp, near_cp10k = round(chemo_cp10k, 2))],
            far[, .(samp, far_cp10k = round(chemo_cp10k, 2), far_bins = bins, far_fovs = fovs)],
            by = "samp"))
r_near <- near[samp=="IB1", chemo_cp10k] / near[samp=="RB1", chemo_cp10k]
r_far  <- far[samp=="IB1", chemo_cp10k]  / far[samp=="RB1", chemo_cp10k]
decay_i <- D[samp=="IB1" & d_um>=100 & d_um<400, sum(chemo)/sum(nCount)] /
           D[samp=="IB1" & d_um>=1600, sum(chemo)/sum(nCount)]
decay_r <- D[samp=="RB1" & d_um>=100 & d_um<400, sum(chemo)/sum(nCount)] /
           D[samp=="RB1" & d_um>=1600, sum(chemo)/sum(nCount)]
cat(sprintf("\nIB1:RB1 ratio near field (100-400um)  = %.2f\n", r_near))
cat(sprintf("IB1:RB1 ratio FAR field  (>1600um)   = %.2f   <-- the plateau offset\n", r_far))
cat(sprintf("\nwithin-sample decay near/far: IB1 %.2fx, RB1 %.2fx  (ratio %.2f)\n",
            decay_i, decay_r, decay_i/decay_r))

## ---- is the decay DIFFERENCE real, or two noisy point estimates? ----------------------
## Bootstrap FOVs within each sample and recompute near/far decay, so the interval respects
## the batch structure. The quantity of interest is the RATIO of the two samples' decays:
## if its interval excludes 1, IB1 has a spatial decay that RB1's top-1% construction does not.
set.seed(11)
bs <- replicate(2000, {
  v <- sapply(c("IB1","RB1"), function(s) {
    z  <- D[samp == s]
    fv <- unique(z$fov); f <- sample(fv, length(fv), TRUE)
    zz <- z[data.table(fov = f), on = "fov", allow.cartesian = TRUE, nomatch = 0L]
    n <- zz[d_um >= 100 & d_um < 400, sum(chemo)/sum(nCount)]
    f2 <- zz[d_um >= 1600, sum(chemo)/sum(nCount)]
    if (!is.finite(n) || !is.finite(f2) || f2 <= 0) return(NA_real_)
    n / f2
  })
  c(v[["IB1"]], v[["RB1"]], v[["IB1"]]/v[["RB1"]])
})
dq <- function(i) quantile(bs[i, is.finite(bs[i,])], c(.025, .975))
cat("\n=== decay (near/far), FOV cluster bootstrap, 2000 resamples\n")
cat(sprintf("  IB1 decay        %.2fx  [%.2f, %.2f]\n", decay_i, dq(1)[1], dq(1)[2]))
cat(sprintf("  RB1 decay        %.2fx  [%.2f, %.2f]\n", decay_r, dq(2)[1], dq(2)[2]))
cat(sprintf("  IB1/RB1 decay    %.2fx  [%.2f, %.2f]  %s\n", decay_i/decay_r, dq(3)[1], dq(3)[2],
            if (dq(3)[1] > 1) "-> IB1 decays MORE; a real field" else "-> not separable from RB1"))

## ---- is the FAR-FIELD residual spatially structured, or uniform? ----------------------
## A plateau flat in distance-to-peak is NOT proof of a technical offset: diffuse HRS cells
## below the hotspot threshold would also be invisible to distance-to-peak. The two differ in
## another way -- biology stays spatially CLUSTERED, a slide effect does not. So within the far
## field only, ask whether the top-quartile chemokine bins sit closer together than chance,
## using the project's nearest-neighbour permutation idiom.
mnn <- function(xy) { if (nrow(xy) < 3) return(NA_real_)
  mean(apply(as.matrix(dist(xy)) + diag(Inf, nrow(xy)), 1, min)) }
coh <- rbindlist(lapply(c("IB1","RB1"), function(s) {
  z <- D[samp == s & d_um >= 1600]
  if (nrow(z) < 40) return(NULL)
  cut4 <- quantile(z$chemo_cp10k, .75, na.rm = TRUE)
  hi <- which(z$chemo_cp10k >= cut4); k <- length(hi)
  xy <- as.matrix(z[, .(x_mm, y_mm)])
  obs <- mnn(xy[hi, , drop = FALSE])
  nul <- replicate(1000, mnn(xy[sample(nrow(xy), k), , drop = FALSE]))
  data.table(samp = s, far_bins = nrow(z), hi_bins = k,
             obs_mnn_um = obs * 1000, null_mnn_um = mean(nul, na.rm = TRUE) * 1000,
             ratio = obs / mean(nul, na.rm = TRUE),
             p_perm = (1 + sum(nul <= obs, na.rm = TRUE)) / (1 + sum(is.finite(nul))))
}))
cat("\n=== far-field (>1600um) spatial coherence of the top-quartile chemokine bins\n")
cat("    clustered => the residual is structured tissue, not a flat slide effect\n")
print(coh[, .(samp, far_bins, hi_bins, obs_mnn_um = round(obs_mnn_um, 1),
              null_mnn_um = round(null_mnn_um, 1), ratio = round(ratio, 3),
              p_perm = signif(p_perm, 3))])
sc("field_extent_coherence.csv", coh)

sc("field_extent_profile.csv", prof)
sc("field_extent_gradient.csv",
   res[, .(samp, n_bins, rho_bin = round(rho_bin,4), p_bin, n_fov,
           rho_fov = round(rho_fov,4), p_fov)])
sc("field_extent_verdict.csv", data.table(
   ratio_near_100_400um = round(r_near, 3), ratio_far_gt1600um = round(r_far, 3),
   decay_IB1 = round(decay_i, 3), decay_RB1 = round(decay_r, 3),
   decay_ratio = round(decay_i/decay_r, 3),
   rho_IB1 = round(res[samp=="IB1", rho_bin], 4), rho_RB1 = round(res[samp=="RB1", rho_bin], 4),
   decay_ratio_lo = round(dq(3)[1], 3), decay_ratio_hi = round(dq(3)[2], 3)))

## ---- figure: profile, both samples, log y ---------------------------------------------
SAMP_COLS <- c(IB1 = "#2a78d6", RB1 = "#eb6834")
p <- ggplot(prof[is.finite(med_d_um)], aes(med_d_um, chemo_cp10k, colour = samp)) +
  geom_line(linewidth = .8) + geom_point(size = 2.4) +
  scale_colour_manual(values = SAMP_COLS, name = NULL) +
  scale_x_continuous(trans = "log10") + scale_y_continuous(trans = "log10") +
  labs(title = "CCL17+CCL22 away from each sample's own top-1% chemokine bins",
       subtitle = sprintf("A field decays with distance; a slide offset cannot. Far-field IB1:RB1 = %.2fx", r_far),
       x = "distance to nearest peak bin (um, log)", y = "pooled cp10k (log)",
       caption = "Peaks are each sample's OWN q99, so the construction is identical in both. RB1 is the no-Hodgkin control.") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top", legend.justification = "left",
        panel.grid.minor = element_blank())
ggsave(file.path(FIG, "fig_17_field_extent.png"), p, width = 7.4, height = 4.6, dpi = 150)
cat("wrote", file.path(FIG, "fig_17_field_extent.png"), "\n")

writeLines(format(Sys.time()), file.path(OUT, "field_extent.done"))
cat("\nOK\n")
