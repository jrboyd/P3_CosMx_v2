## ============================================================================
## 12_script_fl_vs_hodgkin_overlay.R -- the two diseases on one map.
##
## IB1 is a COMPOSITE node: classic Hodgkin plus follicular lymphoma. Two independent
## readouts now exist for it, built from disjoint gene sets and different units:
##   FL       step 07 light-chain neighbourhoods (IGKC vs IGLC1/2), r40um, per
##            NEIGHBOURHOOD -- 17 kappa-restricted (FL-like) of 154 scored
##   Hodgkin  step 11 CCL17+CCL22 hotspots, per 100 um BIN
## Neither uses any gene the other uses, so where they sit relative to each other is a
## real question rather than a restatement.
##
## WINDOW is step 07's global box verbatim -- WIN_X 15,000-100,000, WIN_Y 5,000-90,000
## canvas px (3.0-13.2 x 9.2-19.4 mm) -- so this figure is directly comparable to
## output/07_lightchain_global/IB1/IB1_global_r40um.png. NOTE that box
## excludes IB1's lower tissue piece and therefore drops some step-11 hotspots; the
## count dropped is printed and put in the subtitle rather than left implicit.
##
## Panels: locator | image | FL neighbourhoods | Hodgkin bins | OVERLAY
##
## Usage: Rscript scripts/12_script_fl_vs_hodgkin_overlay.R [sample] [radius_um] [qc_mode]
## ============================================================================

if (!file.exists("scripts/lib/cosmx_tiff.R"))
  stop("run from the project root", call. = FALSE)
source("scripts/lib/cosmx_tiff.R")
suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(data.table); setDTthreads(1)
})

args    <- commandArgs(trailingOnly = TRUE)
SAMPLE  <- if (length(args) >= 1) args[1] else "IB1"
RADIUS  <- if (length(args) >= 2) as.integer(args[2]) else 40L
QC_MODE <- if (length(args) >= 3 && args[3] %in% c("postqc","preqc")) args[3] else "postqc"

LC_DIR  <- file.path("output/07_lightchain", SAMPLE)
NI_DIR  <- if (QC_MODE == "preqc") "output/11_hodgkin_niche_preqc" else "output/11_hodgkin_niche"
OUT_DIR <- "output/12_fl_vs_hodgkin"
FLAT    <- "/netfiles/volaric_research/DLBCL/CosMx/export1/flatFiles/%s/%s_metadata_file.csv.gz"
SEG_F   <- "output/10_hrs/hrs_segments.csv.gz"
NB_F    <- file.path(LC_DIR, sprintf("lightchain_neighborhoods_r%dum.csv", RADIUS))
BIN_F   <- file.path(NI_DIR, "niche_bins_all_scales.csv.gz")
for (f in c(NB_F, BIN_F)) if (!file.exists(f)) stop("missing ", f, call. = FALSE)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

## step 07's global window, verbatim
WIN_X <- c(15000, 100000); WIN_Y <- c(5000, 90000)
BIN_UM <- 100; NICHE_Q <- 0.99; NULL_Q <- 0.99; BG_N <- 250000L

## Two diseases, two hues, from the validated first slots. FL orange, Hodgkin blue --
## and nothing else on the overlay competes with them.
FL_COL  <- "#eb6834"; FL_DIM <- "#f6c4ac"
HK_COL  <- "#0b3d7a"; HK_MID <- "#7ba3d0"; HK_LOW <- "#eef3fa"
SURF<-"#fcfcfb"; INK<-"#0b0b0b"; INK2<-"#52514e"; MUTED<-"#898781"; AXIS<-"#c3c2b7"
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
cat("step 12 -- FL neighbourhoods vs Hodgkin niches, ", SAMPLE, " [", QC_MODE, "]\n", sep="")
cat("============================================================\n")

## ---- coordinates ---------------------------------------------------------------
flat <- fread(sprintf(FLAT, SAMPLE, SAMPLE),
              select = c("cell","fov","CenterX_local_px","CenterY_local_px"))
px <- cosmx_local_to_px(SAMPLE, flat$fov, flat$CenterX_local_px, flat$CenterY_local_px)
flat[, `:=`(key = cosmx_seu_key(cell), X_px = px$X_px, Y_px = px$Y_px)]
SEG <- fread(SEG_F, select = c("cell_id","samp","x_slide_mm","y_slide_mm"))
SEG[, key := cosmx_seu_key(cell_id)]
d <- merge(SEG[samp == SAMPLE, .(key, x_slide_mm, y_slide_mm)],
           flat[, .(key, X_px, Y_px)], by = "key")
fx <- lm(X_px ~ x_slide_mm, data = d); fy <- lm(Y_px ~ y_slide_mm, data = d)
AFF <- list(ax=unname(coef(fx)[2]), bx=unname(coef(fx)[1]),
            ay=unname(coef(fy)[2]), by=unname(coef(fy)[1]))
if (AFF$ax > 0 || AFF$ay > 0) stop("mm->px slope must be negative (mirrored x axis)")
if (max(abs(abs(c(AFF$ax,AFF$ay)) - 1/COSMX_MM_PER_PX)*COSMX_MM_PER_PX) > 0.01)
  stop("mm->px slope is >1% from the expected ", round(1/COSMX_MM_PER_PX))
cat(sprintf("mm->px affine ax %.1f ay %.1f, resid %.2f/%.2f px\n",
            AFF$ax, AFF$ay, sd(resid(fx)), sd(resid(fy))))
mm2px <- function(x,y) list(X = AFF$ax*x + AFF$bx, Y = AFF$ay*y + AFF$by)
UMPX <- COSMX_MM_PER_PX*1000; .um <- function(v) round(v*UMPX)
rect <- cosmx_px_rect(WIN_X[1], WIN_X[2], WIN_Y[1], WIN_Y[2], name = "global")
mmx <- sort((WIN_X - AFF$bx)/AFF$ax); mmy <- sort((WIN_Y - AFF$by)/AFF$ay)
.coord <- function() list(scale_x_reverse(labels = .um),
                          coord_fixed(xlim = WIN_X, ylim = WIN_Y))
inwin <- function(X, Y) X >= WIN_X[1] & X <= WIN_X[2] & Y >= WIN_Y[1] & Y <= WIN_Y[2]

## ---- FL neighbourhoods (step 07) -------------------------------------------------
NB <- fread(NB_F)[artifact_suspect == FALSE]
p_nb <- mm2px(NB$x, NB$y); NB[, `:=`(X = p_nb$X, Y = p_nb$Y)]
NB[, fl := class == "kappa-restricted (FL-like)"]
NB[, inw := inwin(X, Y)]
cat(sprintf("\nFL neighbourhoods r%dum (artifact-free): %d, of which %d FL-like\n",
            RADIUS, nrow(NB), NB[, sum(fl)]))
cat(sprintf("  inside the step-07 window: %d (%d FL-like); OUTSIDE and therefore not drawn: %d (%d FL-like)\n",
            NB[, sum(inw)], NB[inw==TRUE, sum(fl)], NB[, sum(!inw)], NB[inw==FALSE, sum(fl)]))

## ---- Hodgkin bins (step 11) -------------------------------------------------------
ALL <- fread(BIN_F)[bin_um == BIN_UM]
ib <- ALL[samp == SAMPLE]; rb <- ALL[samp != SAMPLE]
HOT_CUT <- quantile(ib$chemo_cp10k, NICHE_Q); FIELD_CUT <- quantile(rb$chemo_cp10k, NULL_Q)
ib[, hot := chemo_cp10k >= HOT_CUT]
c1 <- mm2px(ib$bx*BIN_UM/1000, ib$by_*BIN_UM/1000)
c2 <- mm2px((ib$bx+1)*BIN_UM/1000, (ib$by_+1)*BIN_UM/1000)
ib[, `:=`(x1=pmin(c1$X,c2$X), x2=pmax(c1$X,c2$X), y1=pmin(c1$Y,c2$Y), y2=pmax(c1$Y,c2$Y))]
ib[, cx := (x1+x2)/2][, cy := (y1+y2)/2]
ib[, inw := inwin(cx, cy)]
cat(sprintf("Hodgkin bins: %d (%d hotspots)\n", nrow(ib), ib[, sum(hot)]))
cat(sprintf("  inside the window: %d bins, %d hotspots; OUTSIDE and not drawn: %d bins, %d hotspots\n",
            ib[, sum(inw)], ib[inw==TRUE, sum(hot)], ib[, sum(!inw)], ib[inw==FALSE, sum(hot)]))
ibw <- ib[inw == TRUE]; nbw <- NB[inw == TRUE]

## ---- do they avoid each other? ---------------------------------------------------
## Distance from each Hodgkin hotspot to the nearest FL-like focus, against the same
## distance for non-hotspot bins. Both restricted to the window so the geometry is
## shared. Neither gene set informs the other, so this is a real test.
fl_xy <- as.matrix(nbw[fl == TRUE, .(X, Y)])
if (nrow(fl_xy) >= 2 && ibw[, sum(hot)] >= 3) {
  d_hot <- RANN::nn2(fl_xy, as.matrix(ibw[hot==TRUE, .(cx,cy)]), k=1)$nn.dists[,1]*UMPX
  d_bg  <- RANN::nn2(fl_xy, as.matrix(ibw[hot==FALSE,.(cx,cy)]), k=1)$nn.dists[,1]*UMPX
  pw <- suppressWarnings(tryCatch(wilcox.test(d_hot, d_bg)$p.value, error=function(e) NA_real_))
  cat(sprintf("\n=== distance to the nearest FL-like focus, inside the window\n"))
  cat(sprintf("  hotspot bins    median %.0f um (n=%d)\n  other bins      median %.0f um (n=%d)\n",
              median(d_hot), length(d_hot), median(d_bg), length(d_bg)))
  cat(sprintf("  Wilcoxon p = %s -> %s\n", signif(pw,3),
              if (!is.finite(pw) || pw >= 0.05) "no spatial relationship either way"
              else if (median(d_hot) > median(d_bg)) "Hodgkin hotspots AVOID the FL foci"
              else "Hodgkin hotspots sit CLOSER to the FL foci"))
  fwrite(data.table(n_fl_foci=nrow(fl_xy), med_hot_um=round(median(d_hot),1),
                    med_bg_um=round(median(d_bg),1), p=pw),
         file.path(OUT_DIR, sprintf("%s_fl_vs_hodgkin_distance_%s.csv", SAMPLE, QC_MODE)))
} else { pw <- NA_real_; d_hot <- d_bg <- numeric(0)
         message("too few FL foci or hotspots in the window for the distance test") }

## ---- panels -------------------------------------------------------------------
set.seed(1); bg <- flat[sample.int(.N, min(.N, BG_N)), .(X_px, Y_px)]
p_loc <- {
  b <- copy(bg); b[, `:=`(x_mm=(X_px-AFF$bx)/AFF$ax, y_mm=(Y_px-AFF$by)/AFF$ay)]
  ggplot() + geom_point(data=b, aes(x_mm,y_mm), colour="#cfcec9", size=0.04) +
    annotate("rect", xmin=mmx[1], xmax=mmx[2], ymin=mmy[1], ymax=mmy[2],
             colour="#008300", fill=NA, linewidth=0.6) + coord_fixed() +
    labs(title=sprintf("0. slide position, mm\n%.1f-%.1f x %.1f-%.1f", mmx[1],mmx[2],mmy[1],mmy[2])) +
    theme_void() + theme(plot.title=element_text(size=9),
      plot.background=element_rect(fill="white", colour="#b5b3ac", linewidth=.3))
}
img <- tryCatch(cosmx_image(SAMPLE, rect) %>% add_um_scale_bar(), error=function(e) NULL)
p_img <- if (is.null(img)) ggplot()+theme_v+labs(title="1. image (unavailable)") else
  img@plots$colorized + .coord() + scale_y_reverse(labels=.um) +
  labs(x="x (um)", y="y (um)", title="1. image") + theme_v

## Panel 2 -- FL, in step 07's own idiom: circles sized by neighbourhood size
p_fl <- ggplot() +
  geom_point(data=bg, aes(X_px,Y_px), colour="#eceae4", size=0.10, stroke=0) +
  geom_point(data=nbw[fl==FALSE], aes(X,Y,size=n_cells), colour=FL_DIM, alpha=.8, stroke=0) +
  geom_point(data=nbw[fl==TRUE],  aes(X,Y,size=n_cells), colour=FL_COL, alpha=.95, stroke=0) +
  scale_size_area(max_size=7, name="cells in nb") +
  .coord() + scale_y_reverse(labels=.um) +
  labs(x="x (um)", y="y (um)",
       title=sprintf("2. FL: %d light-chain neighbourhoods, %d kappa-restricted (dark)",
                     nrow(nbw), nbw[, sum(fl)])) + theme_v

## Panel 3 -- Hodgkin, bins by score
## CLAMPED to a percentile window rather than the data range. 80% of post-QC bins sit
## between 1 and 20 cp10k while the max reaches 272, so a ramp stretched to the maximum
## renders four-fifths of the map in one pale band and hides the field/hotspot structure
## the panel exists to show. The window is stated in the legend so a squished bin is not
## mistaken for an in-range one.
CLAMP_Q <- c(0.05, 0.99)
sc_lim <- unname(quantile(ibw$chemo_cp10k + 0.1, CLAMP_Q))
sc_lim[1] <- max(sc_lim[1], 0.1)
p_hk <- ggplot() +
  geom_point(data=bg, aes(X_px,Y_px), colour="#eceae4", size=0.10, stroke=0) +
  geom_rect(data=ibw, aes(xmin=x1,xmax=x2,ymin=y1,ymax=y2, fill=chemo_cp10k+0.1), colour=NA) +
  scale_fill_gradient(low=HK_LOW, high=HK_COL, trans="log10",
                      name=sprintf("CCL17+CCL22 cp10k (clamped %.1f-%.0f)",
                                   sc_lim[1], sc_lim[2]),
                      limits=sc_lim, oob=scales::squish,
                      breaks=c(1,10,100), labels=c("1","10","100")) +
  .coord() + scale_y_reverse(labels=.um) +
  labs(x="x (um)", y="y (um)",
       title=sprintf("3. Hodgkin: CCL17+CCL22 per %g um bin (%d hotspots)",
                     BIN_UM, ibw[, sum(hot)])) + theme_v

## Panel 4 -- THE OVERLAY. Hodgkin as filled tiles underneath, FL as rings on top, so
## neither hides the other and co-location is judged by eye rather than asserted.
p_ov <- ggplot() +
  geom_point(data=bg, aes(X_px,Y_px), colour="#eceae4", size=0.10, stroke=0) +
  geom_rect(data=ibw[hot==FALSE & chemo_cp10k >= FIELD_CUT],
            aes(xmin=x1,xmax=x2,ymin=y1,ymax=y2), fill=HK_MID, alpha=.55, colour=NA) +
  geom_rect(data=ibw[hot==TRUE],
            aes(xmin=x1,xmax=x2,ymin=y1,ymax=y2), fill=HK_COL, colour=HK_COL, linewidth=.4) +
  geom_point(data=nbw[fl==TRUE], aes(X,Y,size=n_cells),
             shape=21, colour=FL_COL, fill=NA, stroke=1.1) +
  scale_size_area(max_size=8, name="FL nb size") +
  .coord() + scale_y_reverse(labels=.um) +
  labs(x="x (um)", y="y (um)",
       title=sprintf("4. OVERLAY -- Hodgkin hotspots (filled blue) vs FL foci (orange rings)")) +
  theme_v

## Wrapped: the figure is now 14in wide (2x2), not 29in (1x5), and this subtitle was
## running off the right edge.
sub <- sprintf(paste0("step-07 window %.1f-%.1f x %.1f-%.1f mm | FL r%dum, ",
  "Hodgkin %g um bins [%s]\nNOT DRAWN, outside this box: %d bins incl. %d hotspots, ",
  "%d neighbourhoods incl. %d FL-like%s"),
  mmx[1],mmx[2],mmy[1],mmy[2], RADIUS, BIN_UM, toupper(QC_MODE),
  ib[, sum(!inw)], ib[inw==FALSE, sum(hot)], NB[, sum(!inw)], NB[inw==FALSE, sum(fl)],
  if (is.finite(pw)) sprintf(" | hotspot-to-FL distance p = %.3g", pw) else "")
## Five panels do not fit a 2x2, and a single row of five renders every panel too small
## to read at page width. The locator is not a data panel -- it only says where on the
## slide this window sits -- so it becomes an inset on panel 1 and the four content
## panels take the 2x2.
p_img_i <- p_img + inset_element(p_loc + theme(plot.title = element_blank()),
                                 left = 0.62, bottom = 0.62, right = 0.995, top = 0.995,
                                 align_to = "panel")
pg <- p_img_i + p_fl + p_hk + p_ov +
  plot_layout(nrow=2, ncol=2) +
  plot_annotation(title=sprintf("%s: follicular lymphoma and Hodgkin signal on one map",
                                SAMPLE), subtitle=sub)
out_f <- file.path(OUT_DIR, sprintf("%s_fl_vs_hodgkin_r%dum_%s.png", SAMPLE, RADIUS, QC_MODE))
ggsave(out_f, pg, width=14, height=16, dpi=300, limitsize=FALSE, bg=SURF)
cat("wrote ", out_f, "\n", sep="")
cat("============================================================\nOK\n")
