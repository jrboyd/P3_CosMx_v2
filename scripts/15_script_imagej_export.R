## ============================================================================
## 15_script_imagej_export.R -- STAGE 1 of the ImageJ/Fiji export.
##
## Turns the two spatial readouts into ONE table of canvas-pixel boxes:
##   hodgkin_hotspot  the CCL17/CCL22 focal hotspot bins   (step 11, 100 um bins)
##   fl_focus         kappa-restricted FL foci             (step 07, 0.15 mm grid)
##   fl_nb            FL-like neighbourhoods               (step 07, r40um, span_mm)
##
## Output: output/15_imagej/<SAMPLE>_regions_canvas_px.csv
## Stage 2 (scripts/15_imagej_tiles.py) reads that file and writes the TIFFs and
## the ImageJ ROI sets. Nothing here touches the 87 GB OME-TIFF.
##
## WHY A FITTED AFFINE. The exact transform (cosmx_local_to_px) needs FOV-local
## px, which only CELLS have. Bin corners and focus centres are recorded in slide
## mm, so they need the fitted mm->px affine. The fit block below is copied from
## scripts/12_script_fl_vs_hodgkin_overlay.R:69-89 -- including both assertions,
## because a positive slope means the fit was done on FOV origins instead of
## cells and every box lands mirrored. See scripts/lib/cosmx_tiff.R:357-414 and
## its lines 44-61 "DO NOT SIMPLIFY" list.
##
##   sbatch --mem=64G --time=1:00:00 scripts/sbatch/run_rscript.sh \
##     scripts/15_script_imagej_export.R IB1
##
## Usage: Rscript scripts/15_script_imagej_export.R [sample] [radius_um] [qc_mode]
## ============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

if (!file.exists("scripts/lib/cosmx_tiff.R"))
  stop("run from the project root: cd /gpfs1/pi/avolaric/files_jrboyd/P3_CosMx_v2",
       call. = FALSE)
source("scripts/lib/cosmx_tiff.R")

args    <- commandArgs(trailingOnly = TRUE)
SAMPLE  <- if (length(args) >= 1) args[1] else "IB1"
RADIUS  <- if (length(args) >= 2) as.integer(args[2]) else 40L
QC_MODE <- if (length(args) >= 3 && args[3] %in% c("postqc", "preqc")) args[3] else "postqc"

LC_DIR  <- file.path("output/07_lightchain", SAMPLE)
NI_DIR  <- if (QC_MODE == "preqc") "output/11_hodgkin_niche_preqc" else "output/11_hodgkin_niche"
OUT_DIR <- "output/15_imagej"
FLAT    <- "/netfiles/volaric_research/DLBCL/CosMx/export1/flatFiles/%s/%s_metadata_file.csv.gz"
SEG_F   <- "output/10_hrs/hrs_segments.csv.gz"
HOT_F   <- file.path(NI_DIR, "niche_bins_called.csv")
FOCI_F  <- file.path(LC_DIR, "lightchain_restricted_foci_persistence.csv")
NB_F    <- file.path(LC_DIR, sprintf("lightchain_neighborhoods_r%dum.csv", RADIUS))
for (f in c(SEG_F, HOT_F, FOCI_F, NB_F))
  if (!file.exists(f)) stop("missing ", f, call. = FALSE)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

BIN_UM   <- 100      # step 11's FOCUS_BIN
FOCI_MM  <- 0.15     # step 07's foci grid pitch
NB_MIN_R <- 0.030    # mm; floor on the FL-neighbourhood radius so a 0.05 mm span
                     # is still a visible circle rather than a dot

## Canvas extent of the assembled OME-TIFF, from its OME-XML SizeX/SizeY. Used
## only to report boxes that fall outside the image -- stage 2 clips them.
CANVAS <- list(IB1 = c(x = 124852L, y = 170561L),
               RB1 = c(x = 116538L, y = 161773L))[[SAMPLE]]
if (is.null(CANVAS)) stop("no canvas extent recorded for ", SAMPLE, call. = FALSE)

cat("============================================================\n")
cat("step 15 -- ImageJ region export, ", SAMPLE, " [", QC_MODE, ", r", RADIUS, "um]\n", sep = "")
cat("============================================================\n")

## ---- mm -> canvas px affine ----------------------------------------------------
flat <- fread(sprintf(FLAT, SAMPLE, SAMPLE),
              select = c("cell", "fov", "CenterX_local_px", "CenterY_local_px"))
px <- cosmx_local_to_px(SAMPLE, flat$fov, flat$CenterX_local_px, flat$CenterY_local_px)
flat[, `:=`(key = cosmx_seu_key(cell), X_px = px$X_px, Y_px = px$Y_px)]
SEG <- fread(SEG_F, select = c("cell_id", "samp", "x_slide_mm", "y_slide_mm"))
SEG[, key := cosmx_seu_key(cell_id)]
d <- merge(SEG[samp == SAMPLE, .(key, x_slide_mm, y_slide_mm)],
           flat[, .(key, X_px, Y_px)], by = "key")
fx <- lm(X_px ~ x_slide_mm, data = d); fy <- lm(Y_px ~ y_slide_mm, data = d)
AFF <- list(ax = unname(coef(fx)[2]), bx = unname(coef(fx)[1]),
            ay = unname(coef(fy)[2]), by = unname(coef(fy)[1]))
if (AFF$ax > 0 || AFF$ay > 0) stop("mm->px slope must be negative (mirrored x axis)")
if (max(abs(abs(c(AFF$ax, AFF$ay)) - 1 / COSMX_MM_PER_PX) * COSMX_MM_PER_PX) > 0.01)
  stop("mm->px slope is >1% from the expected ", round(1 / COSMX_MM_PER_PX))
cat(sprintf("n=%d cells; mm->px affine ax %.3f bx %.1f ay %.3f by %.1f, resid %.2f/%.2f px\n",
            nrow(d), AFF$ax, AFF$bx, AFF$ay, AFF$by, sd(resid(fx)), sd(resid(fy))))
mm2px <- function(x, y) list(X = AFF$ax * x + AFF$bx, Y = AFF$ay * y + AFF$by)

## Both slopes are negative, so an mm box's min corner maps to the px MAX corner.
## Normalise with pmin/pmax rather than assuming an order -- the same reason
## cosmx_px_rect() exists (scripts/lib/cosmx_tiff.R:426-436).
mm_box <- function(x_lo, y_lo, x_hi, y_hi) {
  a <- mm2px(x_lo, y_lo); b <- mm2px(x_hi, y_hi)
  data.table(x1 = pmin(a$X, b$X), x2 = pmax(a$X, b$X),
             y1 = pmin(a$Y, b$Y), y2 = pmax(a$Y, b$Y))
}

## ---- 1. Hodgkin chemokine hotspots ---------------------------------------------
## Geometric bin bounds from (bx, by_), NOT the x_mm/y_mm member centroid: the
## centroid drifts inside the bin, so centroid-derived boxes do not abut. Same
## choice as 11_script_hodgkin_niche_views.R:213-220 and
## 12_script_fl_vs_hodgkin_overlay.R:110-113.
HOT <- fread(HOT_F)[samp == SAMPLE & bin_um == BIN_UM]
hb <- mm_box(HOT$bx * BIN_UM / 1000, HOT$by_ * BIN_UM / 1000,
             (HOT$bx + 1) * BIN_UM / 1000, (HOT$by_ + 1) * BIN_UM / 1000)
hot <- data.table(
  group = "hodgkin_hotspot",
  id    = sprintf("%s_%d_%d", HOT$samp, HOT$bx, HOT$by_),
  shape = "rect",
  hb,
  metric      = round(HOT$chemo_cp10k, 1),
  metric_name = "chemo_cp10k",
  n_cells     = HOT$cells,
  flag        = ifelse(HOT$seam_flag, "seam", ""),
  note        = sprintf("CCL17=%d CCL22=%d CD30=%d seam_um=%.0f",
                        HOT$g_CCL17, HOT$g_CCL22, HOT$g_TNFRSF8, HOT$seam_um))
cat(sprintf("hodgkin_hotspot: %d bins (%d seam-flagged), %.0f px per %d um bin\n",
            nrow(hot), sum(hot$flag == "seam"),
            median(hot$x2 - hot$x1), BIN_UM))

## ---- 2. FL foci (persistence across radii, 0.15 mm grid) -----------------------
## gx = round(x/0.15), so the cell the index names spans (gx +/- 0.5)*0.15.
FOC <- fread(FOCI_F)
fb <- mm_box((FOC$gx - 0.5) * FOCI_MM, (FOC$gy - 0.5) * FOCI_MM,
             (FOC$gx + 0.5) * FOCI_MM, (FOC$gy + 0.5) * FOCI_MM)
foc <- data.table(
  group = "fl_focus",
  id    = sprintf("%s_%d_%d", SAMPLE, FOC$gx, FOC$gy),
  shape = "rect",
  fb,
  metric      = FOC$kappa,
  metric_name = "kappa_frac",
  n_cells     = FOC$cells,
  flag        = ifelse(FOC$n_radii >= 3, "persistent", ""),
  note        = sprintf("radii=%s n_radii=%d", FOC$radii, FOC$n_radii))
cat(sprintf("fl_focus: %d foci, %.0f px per %.2f mm cell\n",
            nrow(foc), median(foc$x2 - foc$x1), FOCI_MM))

## ---- 3. FL-like neighbourhoods (r40um, real extent) ---------------------------
## span_mm is the max radial extent from the centroid, so it is a radius. Drop
## artifact_suspect -- step 12 does the same, these are the tiling/edge artefacts.
NB <- fread(NB_F)[class == "kappa-restricted (FL-like)" & artifact_suspect == FALSE]
rad <- pmax(NB$span_mm, NB_MIN_R)
nb_b <- mm_box(NB$x - rad, NB$y - rad, NB$x + rad, NB$y + rad)
nbr <- data.table(
  group = "fl_nb",
  id    = sprintf("%s_nb%d", SAMPLE, NB$nb),
  shape = "oval",
  nb_b,
  metric      = NB$kappa_frac,
  metric_name = "kappa_frac",
  n_cells     = NB$n_cells,
  flag        = ifelse(NB$span_mm < NB_MIN_R, "radius_floored", ""),
  note        = sprintf("span_mm=%.3f q=%.2g n_kappa=%d n_lambda=%d",
                        NB$span_mm, NB$q, NB$n_kappa, NB$n_lambda))
cat(sprintf("fl_nb: %d FL-like neighbourhoods (%d radius-floored to %.0f um)\n",
            nrow(nbr), sum(nbr$flag == "radius_floored"), NB_MIN_R * 1000))

## ---- combine, report out-of-canvas, write --------------------------------------
R <- rbindlist(list(hot, foc, nbr), use.names = TRUE)
R[, `:=`(cx = (x1 + x2) / 2, cy = (y1 + y2) / 2,
         w = x2 - x1, h = y2 - y1)]
R[, in_canvas := x1 >= 0 & y1 >= 0 & x2 <= CANVAS[["x"]] & y2 <= CANVAS[["y"]]]
R[, `:=`(x1 = round(x1, 1), x2 = round(x2, 1),
         y1 = round(y1, 1), y2 = round(y2, 1),
         cx = round(cx, 1), cy = round(cy, 1),
         w = round(w, 1), h = round(h, 1))]
setorder(R, group, -metric)

cat("\ncanvas is ", CANVAS[["x"]], " x ", CANVAS[["y"]], " px\n", sep = "")
oob <- R[in_canvas == FALSE]
if (nrow(oob)) {
  cat(sprintf("WARNING: %d box(es) fall partly outside the canvas -- stage 2 clips them:\n",
              nrow(oob)))
  print(oob[, .(group, id, x1, y1, x2, y2)])
} else {
  cat("all boxes are inside the canvas\n")
}
cat("\nper group:\n"); print(R[, .(n = .N, med_w_px = median(w)), by = group])

out_f <- file.path(OUT_DIR, sprintf("%s_regions_canvas_px.csv", SAMPLE))
fwrite(R[, .(group, id, shape, x1, y1, x2, y2, cx, cy, w, h,
             metric, metric_name, n_cells, flag, in_canvas, note)], out_f)
cat("\nwrote ", out_f, " (", nrow(R), " regions)\n", sep = "")

aff_f <- file.path(OUT_DIR, sprintf("%s_mm_to_px_affine.csv", SAMPLE))
fwrite(data.table(sample_id = SAMPLE, ax = AFF$ax, bx = AFF$bx,
                  ay = AFF$ay, by = AFF$by,
                  resid_x_px = sd(resid(fx)), resid_y_px = sd(resid(fy)),
                  n_cells = nrow(d),
                  canvas_x_px = CANVAS[["x"]], canvas_y_px = CANVAS[["y"]],
                  mm_per_px = COSMX_MM_PER_PX), aff_f)
cat("wrote ", aff_f, "\n", sep = "")
writeLines(format(Sys.time()), file.path(OUT_DIR, sprintf("%s_regions.done", SAMPLE)))
