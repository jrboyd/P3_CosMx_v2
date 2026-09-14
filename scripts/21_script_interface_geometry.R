#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------------------
## step 21a -- AXIS SELECTION AND GATE for the FL/Hodgkin interface design.
##
## THE DESIGN. Within IB1 only: does the FL tissue change with exposure to the Hodgkin
## chemokine field? Step 17 names this as the highest-value remaining analysis because it
## needs no cross-sample subtraction, so it is immune to the between-slide offset
## (reference-zone ORs to 5.33) that forecloses every IB1-vs-RB1 regional contrast, and
## immune to the unattributed 4.07x far-field residual.
##
## WHY A GATE RUNS BEFORE ANY CONTRAST. The project has already lost a design to
## non-identifiability discovered AFTER the contrast was computed -- step 16 C3/C4, where
## 2 of 89 FL bins found a composition-matched partner and the difference-in-differences
## meant to rescue it had to be withdrawn. The lesson is to ask whether the contrast EXISTS,
## and on what axis, before asking what it says.
##
## TWO CANDIDATE AXES, AND THIS SCRIPT CHOOSES BETWEEN THEM.
##
##   AXIS A -- distance to the nearest cHL hotspot bin. The intuitive choice, and the one
##     step 12 used. Tested here and REJECTED; the reasons are the point of this script.
##
##   AXIS B -- the bin's OWN local field level, CCL17+CCL22 cp10k. Already used successfully
##     by 16_script_field_gradient.R for eight pre-declared genes (CCR4 rho 0.415 inside FL
##     foci). Tested here and ACCEPTED.
##
## FL IS DEFINED BY LIGHT CHAIN, NOT BY THE ARCHITECTURE AXIS. bin_region == "FL_focus"
## comes from proximity to step-07 kappa-restricted neighbourhoods. CLAUDE.md fact 4 records
## that IB1's arch B_dense tissue carries the Hodgkin chemokines and NOT the FL
## germinal-centre programme, so an arch-based FL call would define the "FL" tissue by the
## very field being tested. That would be circular by construction.
##
## Reads only the 100um bin table (3.8 MB gz). Light enough for the login node.
## Usage: Rscript scripts/21_script_interface_geometry.R [tier]
## ---------------------------------------------------------------------------------------
suppressMessages({library(data.table)})
setDTthreads(1); options(width = 200)

TIER  <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "matched"
BIN_F <- file.path("output/16_regions", TIER, "bins_all_scales.csv.gz")
OUT   <- file.path("output/21_interface", TIER)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
sc <- function(f, x) { fwrite(x, file.path(OUT, f)); cat("wrote", file.path(OUT, f), "\n") }
if (!file.exists(BIN_F)) stop("missing ", BIN_F, call. = FALSE)

## Gate constants. MIN_MEDIAN_SEP_UM is one field length: the field decays 2.71x over
## 100-400um -> >1600um and is mostly spent by 400-800um, so strata whose medians differ by
## less than ~800um are two points on the same plateau. Read from the measurement, not
## invented here.
MIN_FOV_PER_STRATUM <- 5L
MIN_MEDIAN_SEP_UM   <- 800
PIECE_GAP_MM        <- 0.4   # bins are 100um; a >400um centre gap means no contiguous tissue
SMD_BAD             <- 0.5

cat("========================================================================\n")
cat("step 21a -- FL/Hodgkin interface: which axis, and is it identifiable?\n")
cat("            tier =", TIER, "\n")
cat("========================================================================\n")

B <- fread(BIN_F)[bin_um == 100 & samp == "IB1" & bin_low_info == FALSE]
stopifnot(nrow(B) > 0,
          all(c("x_mm","y_mm","fov","bin_region","chemo","nCount","cells",
                "fb_a","f_assigned","seam_um","med_depth") %in% names(B)))
cat(sprintf("\nIB1 100um bins, non-low-info: %d in %d FOVs\n", nrow(B), uniqueN(B$fov)))
print(B[, .(bins = .N, fovs = uniqueN(fov), cells = sum(cells)), by = bin_region][order(-bins)])

smd <- function(a, b) { s <- sqrt((var(a) + var(b))/2)
                        if (!is.finite(s) || s == 0) return(NA_real_); (mean(a)-mean(b))/s }
bal <- function(D, grp, hi, lo, vars) rbindlist(lapply(vars, function(v) {
  a <- D[get(grp) == hi, get(v)]; b <- D[get(grp) == lo, get(v)]
  data.table(variable = v, hi = round(mean(a), 3), lo = round(mean(b), 3),
             smd = round(smd(a, b), 3),
             p_wilcox = signif(suppressWarnings(wilcox.test(a, b)$p.value), 3))
}))[, flag := ifelse(is.na(smd), "", ifelse(abs(smd) >= SMD_BAD, "IMBALANCED",
                                     ifelse(abs(smd) >= 0.25, "watch", "ok")))][]

## =======================================================================================
## Tissue pieces. Euclidean distance is only a diffusion distance along contiguous tissue,
## so a section in several pieces makes AXIS A partly meaningless before any statistics.
## =======================================================================================
hc <- hclust(dist(as.matrix(B[, .(x_mm, y_mm)])), method = "single")
B[, piece := cutree(hc, h = PIECE_GAP_MM)]
pc <- B[, .(bins = .N, fovs = uniqueN(fov),
            x_mm = paste0(round(min(x_mm),1), "-", round(max(x_mm),1)),
            y_mm = paste0(round(min(y_mm),1), "-", round(max(y_mm),1)),
            hotspots = sum(bin_region == "cHL_hotspot"),
            fl_bins  = sum(bin_region == "FL_focus"),
            chemo_cp10k = round(sum(chemo)/sum(nCount)*1e4, 1)),
        by = piece][order(-bins)][bins >= 20]
cat(sprintf("\n=== IB1 is in %d tissue pieces of >=20 bins (single linkage, %g um gap)\n",
            nrow(pc), PIECE_GAP_MM*1000))
print(pc); sc("tissue_pieces.csv", pc)
cat("\n  NOTE: a piece with hotspots but no FL, or FL but no hotspots, contributes\n")
cat("  cross-piece distances to AXIS A that are not diffusion paths.\n")

## =======================================================================================
## AXIS A -- distance to the nearest hotspot bin.
## =======================================================================================
## The hotspot SOURCE SET is smaller than the 80 the region summary reports: the low-info
## filter applied above removes some. Stated because AXIS A's source set being incomplete is
## one of the two reasons it fails, and an unstated 80 -> 47 would hide it.
HOT <- as.matrix(B[bin_region == "cHL_hotspot", .(x_mm, y_mm)])
n_hot_all <- fread(file.path("output/16_regions", TIER, "bin_region_summary.csv"))[
                samp == "IB1" & bin_region == "cHL_hotspot", bins]
cat(sprintf("\n\n### AXIS A -- distance to the nearest of %d hotspot bins\n", nrow(HOT)))
cat(sprintf("    (of %d cHL_hotspot bins at this tier; %d are low-information and dropped)\n",
            n_hot_all, n_hot_all - nrow(HOT)))
if (nrow(HOT) < 10) stop("too few hotspots to define a distance axis", call. = FALSE)
NH <- B[bin_region != "cHL_hotspot"]
NH[, d_hot_um := apply(as.matrix(.SD), 1, function(p)
      sqrt(min((HOT[,1]-p[1])^2 + (HOT[,2]-p[2])^2))) * 1000, .SDcols = c("x_mm","y_mm")]
## Same-piece distance, so the cost of the cross-piece problem is visible rather than
## assumed. Computed per piece rather than per row: the row-wise version is easy to write
## with a self-referential `piece` on both sides of the comparison, which silently returns
## nonsense (it reported a median of 11,687 um -- wider than the whole section).
HP <- B[bin_region == "cHL_hotspot", .(x_mm, y_mm, piece)]
NH[, d_hot_same_piece_um := NA_real_]
for (pz in unique(NH$piece)) {
  h <- HP[piece == pz]
  idx <- which(NH$piece == pz)
  if (!nrow(h) || !length(idx)) next          # no hotspot in this piece -> stays NA
  NH[idx, d_hot_same_piece_um := apply(as.matrix(.SD), 1, function(p)
       sqrt(min((h$x_mm - p[1])^2 + (h$y_mm - p[2])^2))) * 1000, .SDcols = c("x_mm","y_mm")]
}

dprof <- NH[, .(bins = .N, fovs = uniqueN(fov),
                q1 = round(quantile(d_hot_um,.25)), med = round(median(d_hot_um)),
                q3 = round(quantile(d_hot_um,.75)), max = round(max(d_hot_um)),
                med_same_piece = round(median(d_hot_same_piece_um, na.rm = TRUE)),
                pct_no_hotspot_in_piece = round(100*mean(is.na(d_hot_same_piece_um)), 1),
                pct_cross_piece = round(100*mean(
                   is.na(d_hot_same_piece_um) |
                   d_hot_same_piece_um > d_hot_um + 1e-6), 1)),
            by = bin_region][order(med)]
print(dprof); sc("axisA_distance_profile.csv", dprof)

FL <- NH[bin_region == "FL_focus"]
FLA <- FL[, .(fl_bins = .N, d_hot_um = median(d_hot_um),
              chemo_cp10k = sum(chemo)/sum(nCount)*1e4, cells = sum(cells),
              depth = median(med_depth), fb_a = median(fb_a),
              f_assigned = median(f_assigned), seam_um = median(seam_um)),
          by = fov][order(d_hot_um)]
SPLIT_A <- median(FLA$d_hot_um)
FLA[, stratum := ifelse(d_hot_um <= SPLIT_A, "near", "far")]
sepA <- FLA[, .(fovs = .N, fl_bins = sum(fl_bins), cells = sum(cells),
                med_d_um = round(median(d_hot_um)),
                chemo_cp10k = round(sum(chemo_cp10k*cells)/sum(cells), 2)),
            by = stratum][order(med_d_um)]
cat(sprintf("\nsplit at the FL FOVs' own median distance, %.0f um:\n", SPLIT_A))
print(sepA)
sepA_um <- abs(diff(sepA$med_d_um))
A_g1 <- sepA_um >= MIN_MEDIAN_SEP_UM && min(sepA$fovs) >= MIN_FOV_PER_STRATUM
cat(sprintf("\n  G1 separation: %.0f um over %d/%d FOVs -> %s\n", sepA_um,
            sepA[1, fovs], sepA[2, fovs], if (A_g1) "PASS" else "FAIL"))

## THE POSITIVE CONTROL, and the reason AXIS A dies. If distance to a hotspot orders field
## exposure inside FL tissue, CCL17+CCL22 must be higher in the near stratum. It is not.
## These are the genes that DEFINE the field, so this is a check on the axis, never a result.
near_chemo <- sepA[stratum == "near", chemo_cp10k]; far_chemo <- sepA[stratum == "far", chemo_cp10k]
A_pc <- near_chemo > far_chemo
cat(sprintf("  POSITIVE CONTROL: CCL17+CCL22 near %.2f vs far %.2f cp10k -> %s\n",
            near_chemo, far_chemo,
            if (A_pc) "PASS" else "FAIL -- distance does not order field exposure"))
rA <- suppressWarnings(cor.test(FLA$chemo_cp10k, FLA$d_hot_um, method = "spearman"))
cat(sprintf("  rho(field exposure, distance) across the %d FL FOVs = %.3f, p = %.3f\n",
            nrow(FLA), rA$estimate, rA$p.value))
cat("\n=== AXIS A balance (FOV means)\n"); balA <- bal(FLA, "stratum", "near", "far",
  c("depth","cells","fb_a","f_assigned","seam_um","fl_bins")); print(balA)
sc("axisA_balance.csv", balA); sc("axisA_fl_fovs.csv", FLA)

## =======================================================================================
## AXIS B -- the bin's own local field level.
## =======================================================================================
cat("\n\n### AXIS B -- the bin's own CCL17+CCL22 cp10k\n")
FLB <- copy(FLA)
SPLIT_B <- median(FLB$chemo_cp10k)
FLB[, stratum := ifelse(chemo_cp10k <= SPLIT_B, "low", "high")]
sepB <- FLB[, .(fovs = .N, fl_bins = sum(fl_bins), cells = sum(cells),
                med_chemo = round(median(chemo_cp10k), 1),
                med_d_hot_um = round(median(d_hot_um))), by = stratum][order(med_chemo)]
cat(sprintf("\nexposure range across FL FOVs: %.1f - %.1f cp10k (%.0fx spread)\n",
            min(FLB$chemo_cp10k), max(FLB$chemo_cp10k),
            max(FLB$chemo_cp10k)/min(FLB$chemo_cp10k)))
cat(sprintf("hotspot cut at this tier is %.1f cp10k, field cut %.1f -- so the FL tissue spans\n",
            176.07, 11.38))
cat("from below the field cut to well above it, without reaching hotspot level.\n")
cat(sprintf("\nsplit at the FL FOVs' own median exposure, %.1f cp10k:\n", SPLIT_B))
print(sepB)
B_g1 <- min(sepB$fovs) >= MIN_FOV_PER_STRATUM &&
        sepB[stratum == "high", med_chemo] / sepB[stratum == "low", med_chemo] >= 2
cat(sprintf("  G1: %d/%d FOVs, %.1fx exposure contrast -> %s\n",
            sepB[1, fovs], sepB[2, fovs],
            sepB[stratum=="high", med_chemo]/sepB[stratum=="low", med_chemo],
            if (B_g1) "PASS" else "FAIL"))
cat("\n=== AXIS B balance (FOV means)\n"); balB <- bal(FLB, "stratum", "high", "low",
  c("depth","cells","fb_a","f_assigned","seam_um","fl_bins","d_hot_um")); print(balB)
sc("axisB_balance.csv", balB)

## The collinearity that shapes the whole downstream design.
fb_smd <- balB[variable == "fb_a", smd]
cat(sprintf("\n  fb_a (B share among assigned cells) SMD = %.2f.\n", fb_smd))
if (abs(fb_smd) >= SMD_BAD) {
  cat("  Field exposure and B density are COLLINEAR inside FL tissue. This is not news --\n")
  cat("  CLAUDE.md fact 4 already records that IB1's B-dense tissue is where the Hodgkin\n")
  cat("  chemokines concentrate (CCL17 21.8 -> 7.3 across the B_dense boundary). The\n")
  cat("  consequence for step 21b is that CELL-TYPE CONDITIONING IS MANDATORY, not\n")
  cat("  optional: a bin-level contrast on this axis would re-measure B share, which is\n")
  cat("  exactly how step 16's C3 contrast failed.\n")
}

## ---- verdict ---------------------------------------------------------------------------
axis <- if (B_g1 && A_pc) "A (distance)" else if (B_g1) "B (own field level)" else "NONE"
cat("\n========================================================================\n")
cat("VERDICT\n")
cat(sprintf("  AXIS A (distance to hotspot): G1 %s, positive control %s -> %s\n",
            if (A_g1) "PASS" else "FAIL", if (A_pc) "PASS" else "FAIL",
            if (A_g1 && A_pc) "USABLE" else "REJECTED"))
cat(sprintf("  AXIS B (own field level):     G1 %s -> %s\n",
            if (B_g1) "PASS" else "FAIL", if (B_g1) "USABLE" else "REJECTED"))
cat(sprintf("  step 21b proceeds on AXIS %s\n", axis))
cat("========================================================================\n")

sc("gate_verdict.csv", data.table(
  tier = TIER, fl_bins = nrow(FL), fl_fovs = uniqueN(FL$fov), hotspot_bins = nrow(HOT),
  tissue_pieces = nrow(pc),
  axisA_split_um = round(SPLIT_A), axisA_sep_um = round(sepA_um),
  axisA_g1 = A_g1, axisA_poscontrol = A_pc,
  axisA_chemo_near = near_chemo, axisA_chemo_far = far_chemo,
  axisA_rho_exposure_vs_dist = round(unname(rA$estimate), 3),
  axisA_p = signif(rA$p.value, 3),
  axisB_split_cp10k = round(SPLIT_B, 2),
  axisB_spread_fold = round(max(FLB$chemo_cp10k)/min(FLB$chemo_cp10k), 1),
  axisB_g1 = B_g1, axisB_fb_a_smd = fb_smd,
  axis_selected = axis))

fwrite(NH[, .(fov, bx, by_, bin_region, piece,
              d_hot_um = round(d_hot_um, 1),
              d_hot_same_piece_um = round(d_hot_same_piece_um, 1),
              cells, nCount, chemo, chemo_cp10k, fb_a, f_assigned, med_depth,
              seam_um, seam_flag)],
       file.path(OUT, "bin_axes.csv.gz"), compress = "gzip")
cat("wrote", file.path(OUT, "bin_axes.csv.gz"), "\n\nOK\n")
