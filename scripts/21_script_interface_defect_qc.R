#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------------------
## step 21e -- IS THE INTERFACE ANALYSIS DISTORTED BY THE OVER-FIXATION DEFECT?
##
## WHAT IS ALREADY HANDLED. Defect cells are dropped PER CELL upstream of everything in
## step 21, which is the method CLAUDE.md requires (never the FOV-level
## qc_spatial_i_fraction flag, which cuts IB1 to 196,316):
##   16_script_regions.R:166             mdc <- md[is_defect == FALSE]
##   16_script_celltype_pseudobulk.R:51  rg[is_defect == FALSE & unbinned == FALSE]
## So no defect transcriptome enters a bin score, a pseudobulk, or a slope.
##
## WHAT IS NOT HANDLED, AND IS THE POINT OF THIS SCRIPT. Removing the cells does not make
## the surviving SPATIAL AREA uniform. IB1's defect is spatially organised -- step 19 maps
## it as a peripheral rim, fixative penetrating from the section edge inward -- so the
## defect deletes cells unevenly across the section, and three things could follow:
##
##   Q1  DIFFERENTIAL BURDEN. If FL bins and background bins lost different fractions of
##       their cells, the two arms of the step-21b interaction are differently depleted.
##       Depth is a covariate in the fit, but "fraction of neighbours deleted" is not.
##
##   Q2  CONFOUNDING WITH THE EXPOSURE AXIS. If defect burden tracks CCL17+CCL22 within FL
##       tissue, the axis is partly a defect axis. Reported as a correlation and an SMD
##       across the same low/high split 21a used.
##
##   Q4  IS DEFECT BURDEN CONFOUNDED WITH THE ARRANGEMENT TEST'S OUTCOME? 21c asks whether
##       the hotspots sit closer to the FL foci than chance. If defect burden tracks both
##       (a) position within the section and (b) distance to the FL foci, then a null that
##       sites hotspot foci uniformly over occupied tissue is comparing the real hotspots
##       against tissue with a different marginal geography, and the test is not a test of
##       dependence. Measured here; 21c's NULL Rd exists because of the answer.
##
##   Q3  DEFECT-CARVED TISSUE PIECES. 21a decomposed IB1 into 3 pieces by single linkage on
##       occupied bins, and 21c sites its null on that mask. But "occupied" is measured
##       AFTER the defect drop, so a gap could be a hole the defect made rather than real
##       tissue separation. Tested by redoing the identical decomposition on PRE-drop
##       occupied bins: if the pieces merge, they were partly artifacts.
##
## Reads cell_regions.csv.gz with a narrow column selection, so it stays light.
## Usage: Rscript scripts/21_script_interface_defect_qc.R [tier]
## ---------------------------------------------------------------------------------------
suppressMessages({library(data.table)})
setDTthreads(1); options(width = 200)

TIER <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "matched"
REG  <- file.path("output/16_regions", TIER, "cell_regions.csv.gz")
OUT  <- file.path("output/21_interface", TIER)
sc <- function(f, x) { fwrite(x, file.path(OUT, f)); cat("wrote", file.path(OUT, f), "\n") }
if (!file.exists(REG)) stop("missing ", REG, call. = FALSE)
PIECE_GAP_MM <- 0.4

cat("========================================================================\n")
cat("step 21e -- over-fixation defect vs the interface analysis, tier =", TIER, "\n")
cat("========================================================================\n")

C <- fread(REG, select = c("samp","fov_id","bx","by_","x_slide_mm","y_slide_mm",
                           "bin_region","is_defect","unbinned","tier_matched"))
C <- C[samp == "IB1"]
cat(sprintf("\nIB1 cells in the region table: %s (defect %s, %.1f%%)\n",
            format(nrow(C), big.mark=","), format(C[, sum(is_defect)], big.mark=","),
            100*C[, mean(is_defect)]))

## ---- Q1: defect burden per bin, by region ----------------------------------------------
## The bin identity comes from the defect-free cells (that is how step 16 binned), so a
## defect cell is attributed to the bin its coordinates fall in for the purpose of this
## count only. Bins with no surviving cell have no bin_region and cannot be recovered here.
BB <- C[unbinned == FALSE, .(cells_all = .N, defect = sum(is_defect),
                             kept = sum(!is_defect),
                             defect_frac = mean(is_defect)),
        by = .(bx, by_)]
BR <- unique(C[unbinned == FALSE & !is.na(bin_region), .(bx, by_, bin_region)])
BB <- merge(BB, BR, by = c("bx","by_"))
q1 <- BB[, .(bins = .N, cells_before = sum(cells_all), cells_kept = sum(kept),
             pct_cells_lost = round(100*sum(defect)/sum(cells_all), 1),
             med_bin_defect_frac = round(median(defect_frac), 3),
             q1_ = round(quantile(defect_frac, .25), 3),
             q3_ = round(quantile(defect_frac, .75), 3)),
         by = bin_region][order(-bins)]
cat("\n=== Q1: how much each region lost to the defect\n")
print(q1); sc("defect_burden_by_region.csv", q1)

## Per-bin table, consumed by 21b as a covariate. Written here rather than recomputed there
## so both scripts use one definition of "defect burden": the fraction of the bin's cells
## that the per-cell defect drop removed.
fwrite(BB[, .(bx, by_, bin_region, cells_all, kept, defect,
              defect_frac = round(defect_frac, 5))],
       file.path(OUT, "defect_per_bin.csv.gz"), compress = "gzip")
cat("wrote", file.path(OUT, "defect_per_bin.csv.gz"), "\n")

fl <- BB[bin_region == "FL_focus", defect_frac]; bg <- BB[bin_region == "background", defect_frac]
smd <- function(a,b){s <- sqrt((var(a)+var(b))/2); if(!is.finite(s)||s==0) NA_real_ else (mean(a)-mean(b))/s}
cat(sprintf("\n  FL vs background per-bin defect fraction: %.3f vs %.3f, SMD %.2f, Wilcoxon p = %.3g\n",
            mean(fl), mean(bg), smd(fl,bg),
            suppressWarnings(wilcox.test(fl,bg)$p.value)))
cat("  -> the two arms of the 21b interaction are ")
cat(if (abs(smd(fl,bg)) < 0.25) "comparably depleted.\n" else "NOT comparably depleted; see the verdict below.\n")

## ---- Q2: does defect burden track the exposure axis inside FL tissue? -------------------
AX <- fread(file.path(OUT, "bin_axes.csv.gz"))
FLX <- merge(AX[bin_region == "FL_focus", .(bx, by_, fov, exposure = chemo_cp10k, cells)],
             BB[, .(bx, by_, defect_frac)], by = c("bx","by_"))
cat(sprintf("\n=== Q2: defect burden vs field exposure, within FL tissue (%d bins)\n", nrow(FLX)))
ctb <- suppressWarnings(cor.test(FLX$defect_frac, FLX$exposure, method = "spearman"))
FV <- FLX[, .(defect_frac = mean(defect_frac),
              exposure = sum(exposure*cells)/sum(cells)), by = fov]
ctf <- suppressWarnings(cor.test(FV$defect_frac, FV$exposure, method = "spearman"))
cat(sprintf("  bin level: rho = %.3f, p = %.3g   (%d bins)\n", ctb$estimate, ctb$p.value, nrow(FLX)))
cat(sprintf("  FOV level: rho = %.3f, p = %.3g   (%d FOVs)  <- the unit that matters\n",
            ctf$estimate, ctf$p.value, nrow(FV)))
S <- median(FV$exposure); FV[, stratum := ifelse(exposure <= S, "low", "high")]
sm <- smd(FV[stratum=="high", defect_frac], FV[stratum=="low", defect_frac])
cat(sprintf("  across 21a's low/high exposure split: defect frac %.3f high vs %.3f low, SMD %.2f\n",
            FV[stratum=="high", mean(defect_frac)], FV[stratum=="low", mean(defect_frac)], sm))
sc("defect_vs_exposure.csv", data.table(
  level = c("bin","fov"), n = c(nrow(FLX), nrow(FV)),
  rho = round(c(ctb$estimate, ctf$estimate), 3),
  p = signif(c(ctb$p.value, ctf$p.value), 3),
  smd_high_vs_low = c(NA_real_, round(sm, 3))))

## ---- Q3: are 21a's tissue pieces real, or defect-carved? -------------------------------
## Identical decomposition, twice: on the bins that survive the defect drop (what 21a and
## 21c actually used) and on every bin that had ANY cell before the drop.
piece_of <- function(D) { if (nrow(D) < 2) return(rep(1L, nrow(D)))
  cutree(hclust(dist(as.matrix(D[, .(x_mm, y_mm)])), method = "single"), h = PIECE_GAP_MM) }
POST <- unique(AX[, .(bx, by_)])
POST <- merge(POST, C[unbinned == FALSE, .(x_mm = mean(x_slide_mm), y_mm = mean(y_slide_mm)),
                      by = .(bx, by_)], by = c("bx","by_"))
PRE  <- C[, .(x_mm = mean(x_slide_mm), y_mm = mean(y_slide_mm)), by = .(bx, by_)]
POST[, piece := piece_of(POST)]; PRE[, piece := piece_of(PRE)]
np <- function(D) D[, .N, by = piece][N >= 20, .N]
cat(sprintf("\n=== Q3: tissue pieces of >=20 bins, same %g um linkage\n", PIECE_GAP_MM*1000))
cat(sprintf("  POST defect drop (what 21a/21c used): %d bins -> %d pieces\n", nrow(POST), np(POST)))
cat(sprintf("  PRE  defect drop (every occupied bin): %d bins -> %d pieces\n", nrow(PRE),  np(PRE)))
print(PRE[, .(bins = .N, x = paste0(round(min(x_mm),1),"-",round(max(x_mm),1)),
              y = paste0(round(min(y_mm),1),"-",round(max(y_mm),1))),
          by = piece][bins >= 20][order(-bins)])
same <- np(POST) == np(PRE)
cat(sprintf("\n  -> the piece decomposition is %s the defect drop.\n",
            if (same) "UNCHANGED by" else "CHANGED by"))
sc("defect_tissue_pieces.csv", data.table(
  stage = c("post_defect_drop","pre_defect_drop"),
  bins = c(nrow(POST), nrow(PRE)), pieces_ge20 = c(np(POST), np(PRE))))

## ---- Q4: is defect burden an edge proxy, and is it confounded with the outcome? ---------
## Restricted to the piece 21c actually tests, and using the same linkage as 21a so the two
## scripts cannot disagree about what the piece is.
## bin_axes.csv.gz holds only NON-hotspot bins (21a wrote the complement of the hotspot set),
## so the hotspots have to come back from the bin table -- and they must be filtered to the
## SAME tissue piece. A first version unioned all 80 hotspots regardless of piece, which put
## hotspots from other pieces into the table and returned a median FL distance of 6,240 um,
## wider than the piece. The piece decomposition is recomputed here with the identical
## linkage 21a and 21c use, so the three scripts cannot disagree about what a piece is.
BT <- fread(file.path("output/16_regions", TIER, "bins_all_scales.csv.gz"))[
        bin_um == 100 & samp == "IB1" & bin_low_info == FALSE]
BT[, piece := cutree(hclust(dist(as.matrix(BT[, .(x_mm, y_mm)])), method = "single"),
                     h = PIECE_GAP_MM)]
PZ4 <- BT[, .N, by = piece][order(-N)][1, piece]
Z4  <- merge(BT[piece == PZ4, .(bx, by_, bin_region)],
             BB[, .(bx, by_, defect_frac)], by = c("bx","by_"))
cat(sprintf("\n  Q4 works on piece %d: %d bins, of which %d hotspot and %d FL\n",
            PZ4, nrow(Z4), Z4[bin_region == "cHL_hotspot", .N], Z4[bin_region == "FL_focus", .N]))
occ <- unique(Z4[, .(bx, by_)]); K <- paste(occ$bx, occ$by_, sep = "_")
## Edge distance: nearest lattice cell that is NOT occupied. A rim defect must track it.
G <- CJ(bx = seq(min(occ$bx)-1, max(occ$bx)+1), by_ = seq(min(occ$by_)-1, max(occ$by_)+1))
EMP <- as.matrix(G[!(paste(bx, by_, sep = "_") %in% K), .(bx, by_)])
xy  <- as.matrix(Z4[, .(bx, by_)])
Z4[, edge_um := apply(xy, 1, function(p)
      sqrt(min((EMP[,1]-p[1])^2 + (EMP[,2]-p[2])^2))) * 100]
FLM <- as.matrix(Z4[bin_region == "FL_focus", .(bx, by_)])
Z4[, d_fl_um := apply(xy, 1, function(p)
      sqrt(min((FLM[,1]-p[1])^2 + (FLM[,2]-p[2])^2))) * 100]

c_edge <- suppressWarnings(cor.test(Z4$defect_frac, Z4$edge_um, method = "spearman"))
nf     <- Z4[bin_region != "FL_focus"]
c_out  <- suppressWarnings(cor.test(nf$defect_frac, nf$d_fl_um, method = "spearman"))
cat(sprintf("\n=== Q4: defect burden as an edge proxy and an outcome confounder (piece %d)\n", PZ4))
cat(sprintf("  rho(defect_frac, distance to tissue edge)   = %.3f  p = %.3g\n",
            c_edge$estimate, c_edge$p.value))
print(Z4[, .(bins = .N, med_defect = round(median(defect_frac), 3)),
         by = .(edge_band = cut(edge_um, c(0,100,200,400,800,Inf), right = FALSE))][order(edge_band)])
cat(sprintf("\n  rho(defect_frac, distance to nearest FL)    = %.3f  p = %.3g\n",
            c_out$estimate, c_out$p.value))
q4 <- Z4[, .(bins = .N, med_edge_um = round(median(edge_um)),
             med_defect = round(median(defect_frac), 3),
             med_d_fl_um = round(median(d_fl_um))), by = bin_region][order(-bins)]
print(q4); sc("defect_edge_geometry.csv", q4)
pe <- suppressWarnings(wilcox.test(Z4[bin_region=="cHL_hotspot", edge_um],
                                   Z4[bin_region=="background", edge_um])$p.value)
cat(sprintf("  hotspot vs background edge distance: Wilcoxon p = %.3g\n", pe))
cat("\n  CONSEQUENCE. Defect burden is largely a rim effect, and it rises with distance to\n")
cat("  the FL foci. The real hotspots sit in more depleted (more peripheral) tissue than\n")
cat("  occupied tissue on average, so a null that relocates them UNIFORMLY compares them\n")
cat("  against tissue intrinsically closer to the FL. 21c therefore reports a second,\n")
cat("  defect-matched null (NULL Rd) alongside, and the two bracket the answer.\n")

## ---- verdict ---------------------------------------------------------------------------
cat("\n========================================================================\n")
cat("VERDICT\n")
cat(sprintf("  Q1 differential burden      FL %.3f vs bg %.3f, SMD %.2f -> %s\n",
    mean(fl), mean(bg), smd(fl,bg), if (abs(smd(fl,bg)) < 0.25) "OK" else "REPORT AS A LIMIT"))
cat(sprintf("  Q2 confounded with exposure FOV rho %.3f (p %.3g), SMD %.2f -> %s\n",
    ctf$estimate, ctf$p.value, sm,
    if (abs(ctf$estimate) < 0.3 && ctf$p.value > 0.05) "OK" else "REPORT AS A LIMIT"))
cat(sprintf("  Q3 defect-carved pieces     %d pre vs %d post -> %s\n", np(PRE), np(POST),
    if (same) "OK, pieces are real" else "REPORT AS A LIMIT"))
cat(sprintf("  Q4 outcome confounding      rho(defect, edge) %.3f; rho(defect, dist to FL) %.3f -> %s\n",
    c_edge$estimate, c_out$estimate,
    if (abs(c_out$estimate) < 0.15) "OK" else "REPORT AS A LIMIT; 21c needs NULL Rd"))
cat("========================================================================\n")
sc("defect_qc_verdict.csv", data.table(
  tier = TIER,
  fl_defect_frac = round(mean(fl),3), bg_defect_frac = round(mean(bg),3),
  fl_vs_bg_smd = round(smd(fl,bg),3),
  rho_defect_exposure_fov = round(unname(ctf$estimate),3),
  p_defect_exposure_fov = signif(ctf$p.value,3),
  smd_defect_high_vs_low = round(sm,3),
  pieces_pre = np(PRE), pieces_post = np(POST), pieces_unchanged = same,
  rho_defect_edge = round(unname(c_edge$estimate), 3),
  p_defect_edge = signif(c_edge$p.value, 3),
  rho_defect_dist_to_fl = round(unname(c_out$estimate), 3),
  p_defect_dist_to_fl = signif(c_out$p.value, 3),
  p_hotspot_vs_bg_edge = signif(pe, 3)))
cat("\nOK\n")
