## ============================================================================
## 11_script_hodgkin_niche_global.R -- whole-tissue view of the Hodgkin chemokine
## signal, the counterpart to scripts/07_script_lightchain_global.R.
##
## Four panels, matching the step-07 global layout:
##   0. slide locator
##   1. image
##   2. CALLED HOTSPOTS -- categorical: hotspot / elevated field / background
##   3. BIN SCORE       -- continuous CCL17+CCL22 cp10k per bin
##
## Two thresholds, because they answer different questions (see
## scripts/11_script_hodgkin_niche.R):
##   field   >= RB1's 99th percentile. 29% of IB1 bins clear it -- IB1's whole
##           distribution is shifted, so this is a tissue-wide field, not a niche.
##   hotspot >= IB1's OWN 99th percentile. 36 bins; ZERO RB1 bins clear the same cut.
##
## WINDOW. Unlike step 07's fixed 10 x 10 mm box, the window here is derived from the
## data extent: IB1 is two tissue pieces and the top-scoring hotspot sits at y = 6.2 mm,
## outside step 07's y 9.2-19.4 mm window. Copying that box would have silently cropped
## the strongest hotspot out of the figure.
##
## Usage: Rscript scripts/11_script_hodgkin_niche_global.R [sample] [bin_um]
##        scripts/11_submit_hodgkin_niche_global.sh
## ============================================================================

if (!file.exists("scripts/lib/cosmx_tiff.R"))
  stop("run from the project root", call. = FALSE)
source("scripts/lib/cosmx_tiff.R")
suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(data.table); setDTthreads(1)
})

args   <- commandArgs(trailingOnly = TRUE)
SAMPLE <- if (length(args) >= 1) args[1] else "IB1"
BIN_UM  <- if (length(args) >= 2) as.numeric(args[2]) else 100
QC_MODE <- if (length(args) >= 3 && args[3] %in% c("postqc","preqc")) args[3] else "postqc"
PREQC   <- QC_MODE == "preqc"

IN_DIR  <- if (PREQC) "output/11_hodgkin_niche_preqc" else "output/11_hodgkin_niche"
OUT_DIR <- "output/11_hodgkin_niche_global"
FLAT    <- "/netfiles/volaric_research/DLBCL/CosMx/export1/flatFiles/%s/%s_metadata_file.csv.gz"
BINS_F  <- file.path(IN_DIR, "niche_bins_all_scales.csv.gz")
SEG_F   <- "output/10_hrs/hrs_segments.csv.gz"
if (!file.exists(BINS_F))
  stop("missing ", BINS_F, " -- run scripts/11_submit_hodgkin_niche.sh first", call. = FALSE)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

NULL_Q  <- 0.99      # RB1 quantile -> the field
NICHE_Q <- 0.99      # IB1's own quantile -> the hotspots
PAD_PX  <- 3000L
BG_N    <- 250000L   # tissue-outline centroids
PT_BG   <- 0.10

SURF<-"#fcfcfb"; INK<-"#0b0b0b"; INK2<-"#52514e"; MUTED<-"#898781"; AXIS<-"#c3c2b7"
## Categorical panel 2 uses validated slot 1 at two lightnesses plus neutral: hotspot
## and field are two LEVELS of one thing, not unrelated categories.
CAT_COLS <- c(hotspot = "#0b3d7a", field = "#7ba3d0", background = "#dedcd4")
theme_v <- theme_minimal(base_size = 11) +
  theme(plot.background = element_rect(fill = SURF, colour = NA),
        panel.background = element_rect(fill = SURF, colour = NA),
        panel.grid = element_blank(),
        axis.text = element_text(colour = MUTED, size = 8),
        axis.title = element_text(colour = INK2, size = 9),
        plot.title = element_text(colour = INK2, face = "bold", size = 11),
        legend.text = element_text(colour = INK2, size = 9),
        legend.title = element_text(colour = INK2, size = 9),
        legend.position = "bottom")

cat("============================================================\n")
cat("step 11 global -- ", SAMPLE, ", ", BIN_UM, " um bins\n", sep = "")
cat("============================================================\n")

## ---- coordinates --------------------------------------------------------------
flat <- fread(sprintf(FLAT, SAMPLE, SAMPLE),
              select = c("cell","fov","CenterX_local_px","CenterY_local_px"))
px <- cosmx_local_to_px(SAMPLE, flat$fov, flat$CenterX_local_px, flat$CenterY_local_px)
flat[, `:=`(key = cosmx_seu_key(cell), X_px = px$X_px, Y_px = px$Y_px)]

## Affine on CELLS: the canvas x axis is mirrored inside each FOV, so the mm->px slope
## must be NEGATIVE (~ -8314 px/mm). Fitting on FOV origins gives a positive slope and
## puts everything off the slide -- asserted below rather than trusted.
SEG <- fread(SEG_F, select = c("cell_id","samp","x_slide_mm","y_slide_mm"))
SEG[, key := cosmx_seu_key(cell_id)]
d  <- merge(SEG[samp == SAMPLE, .(key, x_slide_mm, y_slide_mm)],
            flat[, .(key, X_px, Y_px)], by = "key")
stopifnot(nrow(d) > 1000)
fx <- lm(X_px ~ x_slide_mm, data = d); fy <- lm(Y_px ~ y_slide_mm, data = d)
AFF <- list(ax = unname(coef(fx)[2]), bx = unname(coef(fx)[1]),
            ay = unname(coef(fy)[2]), by = unname(coef(fy)[1]))
exp_slope <- 1 / COSMX_MM_PER_PX
if (AFF$ax > 0 || AFF$ay > 0) stop("mm->px slope must be negative; got ", round(AFF$ax))
if (max(abs(abs(c(AFF$ax, AFF$ay)) - exp_slope) / exp_slope) > 0.01)
  stop("mm->px slope ", round(AFF$ax), " is >1% from the expected ", round(exp_slope))
cat(sprintf("mm->px affine: ax %.1f ay %.1f (expected |%.1f|), resid %.2f/%.2f px, n=%d\n",
            AFF$ax, AFF$ay, exp_slope, sd(resid(fx)), sd(resid(fy)), nrow(d)))
mm2px <- function(x, y) list(X = AFF$ax * x + AFF$bx, Y = AFF$ay * y + AFF$by)
UMPX  <- COSMX_MM_PER_PX * 1000
.um   <- function(v) round(v * UMPX)

## ---- bins, and the two thresholds ---------------------------------------------
ALL <- fread(BINS_F)[bin_um == BIN_UM]
ib  <- ALL[samp == SAMPLE]; rb <- ALL[samp != SAMPLE]
FIELD_CUT <- quantile(rb$chemo_cp10k, NULL_Q)
HOT_CUT   <- quantile(ib$chemo_cp10k, NICHE_Q)
ib[, cat := fifelse(chemo_cp10k >= HOT_CUT, "hotspot",
            fifelse(chemo_cp10k >= FIELD_CUT, "field", "background"))]
ib[, cat := factor(cat, levels = c("hotspot","field","background"))]
cat(sprintf("\nfield cut (RB1 q%.0f) %.1f | hotspot cut (IB1 q%.0f) %.1f\n",
            100*NULL_Q, FIELD_CUT, 100*NICHE_Q, HOT_CUT))
print(ib[, .(bins = .N, med_cp10k = round(median(chemo_cp10k), 1)), by = cat][order(cat)])
cat(sprintf("RB1 bins clearing the hotspot cut: %d of %d\n",
            sum(rb$chemo_cp10k >= HOT_CUT), nrow(rb)))

c1 <- mm2px( ib$bx      * BIN_UM/1000,  ib$by_      * BIN_UM/1000)
c2 <- mm2px((ib$bx + 1) * BIN_UM/1000, (ib$by_ + 1) * BIN_UM/1000)
ib[, `:=`(x1 = pmin(c1$X, c2$X), x2 = pmax(c1$X, c2$X),
          y1 = pmin(c1$Y, c2$Y), y2 = pmax(c1$Y, c2$Y))]

## ---- window from the DATA, not a fixed box -------------------------------------
pl  <- cosmx_fov_placement(SAMPLE)
CAN <- c(x = max(pl$x0) + COSMX_FOV_PX - 1L, y = max(pl$y0) + COSMX_FOV_PX - 1L)
WIN <- list(x = c(max(0, min(ib$x1) - PAD_PX), min(CAN[["x"]], max(ib$x2) + PAD_PX)),
            y = c(max(0, min(ib$y1) - PAD_PX), min(CAN[["y"]], max(ib$y2) + PAD_PX)))
rect <- cosmx_px_rect(WIN$x[1], WIN$x[2], WIN$y[1], WIN$y[2], name = "global")
mmx <- sort((WIN$x - AFF$bx) / AFF$ax); mmy <- sort((WIN$y - AFF$by) / AFF$ay)
cat(sprintf("window %.0f-%.0f x %.0f-%.0f px  (%.1f-%.1f x %.1f-%.1f mm, %d FOVs)\n",
            WIN$x[1], WIN$x[2], WIN$y[1], WIN$y[2], mmx[1], mmx[2], mmy[1], mmy[2],
            nrow(cosmx_fovs_in_rect(SAMPLE, rect))))

.coord <- function() list(scale_x_reverse(labels = .um),
                          coord_fixed(xlim = WIN$x, ylim = WIN$y))
set.seed(1)
bg <- flat[sample.int(.N, min(.N, BG_N)), .(X_px, Y_px)]

## ---- panels ---------------------------------------------------------------------
p_loc <- {
  b <- copy(bg); b[, `:=`(x_mm = (X_px - AFF$bx)/AFF$ax, y_mm = (Y_px - AFF$by)/AFF$ay)]
  ggplot() +
    geom_point(data = b, aes(x_mm, y_mm), colour = "#cfcec9", size = 0.04) +
    annotate("rect", xmin = mmx[1], xmax = mmx[2], ymin = mmy[1], ymax = mmy[2],
             colour = "#008300", fill = NA, linewidth = 0.6) +
    coord_fixed() +
    labs(title = sprintf("0. slide position, mm\n%.1f-%.1f x %.1f-%.1f mm",
                         mmx[1], mmx[2], mmy[1], mmy[2])) +
    theme_void() +
    theme(plot.title = element_text(size = 9),
          plot.background = element_rect(fill = "white", colour = "#b5b3ac", linewidth = .3))
}

img <- tryCatch(cosmx_image(SAMPLE, rect) %>% add_um_scale_bar(), error = function(e) NULL)
p_img <- if (is.null(img)) ggplot() + theme_v + labs(title = "1. image (unavailable)") else
  img@plots$colorized + .coord() + scale_y_reverse(labels = .um) +
  labs(x = "x (um)", y = "y (um)", title = "1. image") + theme_v

## Panel 2 -- called hotspots. Background bins are drawn too, so the reader sees the
## denominator: the hotspots are 1% of a tissue that is 29% "field".
p_cat <- ggplot() +
  geom_point(data = bg, aes(X_px, Y_px), colour = "#eceae4", size = PT_BG, stroke = 0) +
  geom_rect(data = ib[cat != "hotspot"],
            aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2, fill = cat), colour = NA) +
  geom_rect(data = ib[cat == "hotspot" & !seam_flag],
            aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2, fill = cat),
            colour = "#0b3d7a", linewidth = 0.5) +
  ## seam-flagged hotspots: same fill, amber outline. Not a fourth fill colour --
  ## panel 2 is categorical and only three fills clear the all-pairs separation bar.
  geom_rect(data = ib[cat == "hotspot" & seam_flag],
            aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2, fill = cat),
            colour = "#eb6834", linewidth = 0.5) +
  scale_fill_manual(values = CAT_COLS, name = NULL, drop = FALSE) +
  .coord() + scale_y_reverse(labels = .um) +
  labs(x = "x (um)", y = "y (um)",
       ## short title: the long form collided with panel 3's at this width
       title = sprintf("2. %d hotspots -- %d seam-free (blue), %d seam-flagged (amber)",
                       ib[cat=="hotspot", .N], ib[cat=="hotspot" & !seam_flag, .N],
                       ib[cat=="hotspot" & seam_flag, .N])) + theme_v

## Panel 3 -- the score itself. Sequential = one hue, light to dark; log10 because the
## range is 0 to ~280 and a linear ramp flattens everything below the peaks.
## CLAMPED to a percentile window rather than the data range. 80% of post-QC bins sit
## between 1 and 20 cp10k while the max reaches 272, so a ramp stretched to the maximum
## renders four-fifths of the map in one pale band and hides the field/hotspot structure
## the panel exists to show. The window is stated in the legend so a squished bin is not
## mistaken for an in-range one.
CLAMP_Q <- c(0.05, 0.99)
sc_lim <- unname(quantile(ib$chemo_cp10k + 0.1, CLAMP_Q))
sc_lim[1] <- max(sc_lim[1], 0.1)          # log10 scale cannot take a zero lower bound
p_score <- ggplot() +
  geom_point(data = bg, aes(X_px, Y_px), colour = "#eceae4", size = PT_BG, stroke = 0) +
  geom_rect(data = ib, aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2,
                           fill = chemo_cp10k + 0.1), colour = NA) +
  scale_fill_gradient(low = "#eef3fa", high = "#0b3d7a", trans = "log10",
                      name = sprintf("CCL17+CCL22 cp10k (clamped %.1f-%.0f)",
                                     sc_lim[1], sc_lim[2]),
                      limits = sc_lim, oob = scales::squish,
                      breaks = c(1, 10, 100), labels = c("1","10","100")) +
  .coord() + scale_y_reverse(labels = .um) +
  labs(x = "x (um)", y = "y (um)",
       title = sprintf("3. CCL17+CCL22 per bin (median %.1f, max %.0f)",
                       median(ib$chemo_cp10k), max(ib$chemo_cp10k))) + theme_v

## 2x2 rather than a single row of four. At nrow=1 the figure is 24 x 8.5in with an
## aspect near 2.8, and embedded at 100% page width every panel renders too small to read.
pg <- p_loc + p_img + p_cat + p_score +
  plot_layout(nrow = 2, ncol = 2) +
  plot_annotation(
    title = sprintf("%s global Hodgkin-chemokine summary, %g um bins [%s]",
                    SAMPLE, BIN_UM, toupper(QC_MODE)),
    ## Wrapped for the 2x2 width; at 14in the single-line form ran off the edge.
    subtitle = sprintf(paste0("window %.1f-%.1f x %.1f-%.1f mm | bins need >=20 cells; ",
      "blanks are bins that did not\nhotspot = IB1 q%.0f (%.1f cp10k), field = RB1 q%.0f ",
      "(%.1f) on a %d-bin field of %d | RB1 bins clearing the hotspot cut: %d%s"),
      mmx[1], mmx[2], mmy[1], mmy[2], 100*NICHE_Q, HOT_CUT, 100*NULL_Q, FIELD_CUT,
      ib[cat=="field", .N], nrow(ib), sum(rb$chemo_cp10k >= HOT_CUT),
      if (PREQC) "  [RB1 NOT COMPARABLE pre-QC]" else ""))
out_f <- file.path(OUT_DIR, sprintf("%s_global_bin%gum_%s.png", SAMPLE, BIN_UM, QC_MODE))
ggsave(out_f, pg, width = 14, height = 16, dpi = 300, limitsize = FALSE, bg = SURF)
cat("wrote ", out_f, "\n", sep = "")
writeLines(format(Sys.time()),
           file.path(OUT_DIR, sprintf("%s_bin%gum_%s.done", SAMPLE, BIN_UM, QC_MODE)))
cat("============================================================\nOK\n")
