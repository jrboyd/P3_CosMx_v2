## ============================================================================
## 11_script_hodgkin_niche_views.R -- image views of the CCL17/CCL22 hotspots.
##
## WHY. scripts/01_import_cosmx_seurat.R:25-26 -- Hodgkin claims get confirmed
## visually. Step 11 located 36 focal hotspots (top 1% of IB1's own 100 um bins, cut
## 124.4 cp10k against an IB1 median of 4.0) with ZERO RB1 bins clearing the same cut,
## spatially clustered (444 um vs 664 um permuted, p = 0.005), T-enriched (11.7% vs
## 6.4%) and CD30-elevated (2.19 vs 1.53 cp10k). None of that has been looked at.
##
## CONTROLS, because a plausible picture proves nothing on its own:
##   hotspot     top IB1 bins by CCL17+CCL22
##   IB1bg       IB1 bins at the median -- same tissue, no chemokine
##   RB1top      RB1's OWN highest bins. Zero RB1 bins clear the IB1 cut, so this shows
##               what the most chemokine-rich place in a Hodgkin-free node looks like.
##
## Panels follow scripts/07_script_lightchain_views.R and scripts/10_script_hrs_views.R.
## Scale is REGIONAL, not per-cell: a 100 um bin plus context, because CCL17/CCL22 are
## secreted and mark a region rather than a cell.
##
## Usage: Rscript scripts/11_script_hodgkin_niche_views.R [n_hotspot] [n_control]
##        scripts/11_submit_hodgkin_niche_views.sh
## ============================================================================

if (!file.exists("scripts/lib/cosmx_tiff.R"))
  stop("run from the project root", call. = FALSE)
source("scripts/lib/cosmx_tiff.R")
suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(data.table); setDTthreads(1)
})

args    <- commandArgs(trailingOnly = TRUE)
N_HOT   <- if (length(args) >= 1) as.integer(args[1]) else 12L
N_CTRL  <- if (length(args) >= 2) as.integer(args[2]) else 4L

IN_DIR  <- "output/11_hodgkin_niche"
OUT_DIR <- "output/11_hodgkin_niche_views"
FLAT    <- "/netfiles/volaric_research/DLBCL/CosMx/export1/flatFiles/%s/%s_metadata_file.csv.gz"
BINS_F  <- file.path(IN_DIR, "niche_bins_all_scales.csv.gz")
HOT_F   <- file.path(IN_DIR, "niche_bins_called.csv")
for (f in c(BINS_F, HOT_F)) if (!file.exists(f))
  stop("missing ", f, " -- run scripts/11_submit_hodgkin_niche.sh first", call. = FALSE)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

TIGHT_PX <- 4200L                    # ~505 um -- the bin plus a ring of context
ZOOMS    <- c(tight = 1, context = 3)
LOC_N    <- 60000L

## Transcripts: three families on the validated first three palette slots. Panel 3 is
## a scatter, so every pair is adjacent at once and only three categories clear the
## all-pairs bar -- the same constraint that forced the step-10 redesign.
FAM <- c(CCL17 = "chemo", CCL22 = "chemo",
         TNFRSF8 = "CD30",
         CD3E = "Tcell", CD4 = "Tcell")
GENES     <- names(FAM)
FAM_LEVELS<- c("chemo", "CD30", "Tcell")
FAM_COLS  <- c(chemo = "#2a78d6", CD30 = "#eb6834", Tcell = "#1baf7a")
FAM_LABS  <- c(chemo = "CCL17/CCL22", CD30 = "CD30 (TNFRSF8)", Tcell = "T cell (CD3E/CD4)")
FAM_SIZE  <- c(chemo = 1.5, CD30 = 2.6, Tcell = 1.0)   # CD30 is rare: enlarge it
## Panel 2 cell fills, same three slots.
LIN_COLS  <- c(B = "#2a78d6", T = "#1baf7a", other = "#dedcd4")

SURF<-"#fcfcfb"; INK<-"#0b0b0b"; INK2<-"#52514e"; MUTED<-"#898781"; AXIS<-"#c3c2b7"
theme_v <- theme_minimal(base_size = 10) +
  theme(plot.background = element_rect(fill = SURF, colour = NA),
        panel.background = element_rect(fill = SURF, colour = NA),
        panel.grid = element_blank(),
        axis.text = element_text(colour = MUTED, size = 7),
        axis.title = element_text(colour = INK2, size = 8),
        plot.title = element_text(colour = INK2, face = "bold", size = 9),
        legend.text = element_text(colour = INK2, size = 8),
        legend.title = element_blank(), legend.position = "bottom")

## Axis units. add_um_scale_bar()'s `unit_per_pixel` defaults to NULL because the
## TiffData object carries the value; pull it from there when the slot exists and fall
## back to the documented header constant (COSMX_MM_PER_PX, the ImPixelSize_nm field).
## Which source was used is printed once, so a silent wrong scale is not possible.
UMPX <- NA_real_
.um_per_px <- function(img) {
  v <- tryCatch({
    sl  <- slotNames(class(img))
    hit <- intersect(c("unit_per_pixel","um_per_pixel","pixel_size","resolution"), sl)
    if (length(hit)) as.numeric(slot(img, hit[1])[1]) else NA_real_
  }, error = function(e) NA_real_)
  if (length(v) && is.finite(v) && v > 0) {
    if (is.na(UMPX)) { UMPX <<- v; cat(sprintf("um/px from the TiffData object: %.6f\n", v)) }
  } else if (is.na(UMPX)) {
    UMPX <<- COSMX_MM_PER_PX * 1000
    cat(sprintf("um/px from COSMX_MM_PER_PX (object slot absent): %.6f\n", UMPX))
  }
  UMPX
}
.um <- function(v) round(v * (if (is.na(UMPX)) COSMX_MM_PER_PX*1000 else UMPX))

## scale_x_reverse(), never reversed coord limits -- ggplot sorts limits passed to
## coord_fixed, so coord_fixed(xlim = c(max, min)) silently does nothing.
.coord_px <- function(rect) {
  q <- rect@coords[1L, ]
  list(scale_x_reverse(labels = .um),
       coord_fixed(xlim = c(q$xmin, q$xmax), ylim = c(q$ymin, q$ymax)))
}
.gp_subset <- function(gp, rect) {
  if (is.null(gp) || !nrow(gp)) return(NULL)
  q <- rect@coords[1L, ]
  out <- gp[gp$X_px >= q$xmin & gp$X_px <= q$xmax &
            gp$Y_px >= q$ymin & gp$Y_px <= q$ymax, , drop = FALSE]
  if (!nrow(out)) NULL else out
}

cat("============================================================\n")
cat("step 11 views -- CCL17/CCL22 hotspots\n")
cat("============================================================\n")

## ---- 1. selection ------------------------------------------------------------
HOT <- fread(HOT_F); ALL <- fread(BINS_F)
BIN_UM <- HOT$bin_um[1]
B100 <- ALL[bin_um == BIN_UM]
setorder(HOT, -chemo_cp10k)
sel <- rbind(
  HOT[seq_len(min(N_HOT, .N))][, grp := "hotspot"],
  B100[samp == "IB1"][order(abs(chemo_cp10k - median(B100[samp=="IB1", chemo_cp10k])))][
       seq_len(min(N_CTRL, .N))][, grp := "IB1bg"],
  B100[samp == "RB1"][order(-chemo_cp10k)][seq_len(min(N_CTRL, .N))][, grp := "RB1top"],
  fill = TRUE)
sel[, id := sprintf("%s_%d_%d", samp, bx, by_)]
cat(sprintf("bin size %g um | rendering %d windows\n", BIN_UM, nrow(sel)))
print(sel[, .(windows = .N, med_chemo = round(median(chemo_cp10k), 1)), by = grp])
fwrite(sel, file.path(OUT_DIR, "selected_windows.csv"))

## ---- 2. exact canvas px, and the mm->px inverse for the bin outline -----------
SAMPS <- unique(sel$samp)
FLATS <- lapply(SAMPS, function(S) {
  f <- fread(sprintf(FLAT, S, S),
             select = c("cell","fov","CenterX_local_px","CenterY_local_px"))
  px <- cosmx_local_to_px(S, f$fov, f$CenterX_local_px, f$CenterY_local_px)
  f[, `:=`(key = cosmx_seu_key(cell), X_px = px$X_px, Y_px = px$Y_px)]
  f
}); names(FLATS) <- SAMPS

## Affine fitted on CELLS, which carry both coordinate systems.
##
## DO NOT fit this on the FOV origins. The first version did, and it is wrong: a FOV
## origin satisfies x0 = (X_mm - min(X_mm)) / COSMX_MM_PER_PX, a POSITIVE slope, while
## a CELL's canvas px is X_px = x0 + (COSMX_FOV_PX - 1) - x_local -- the x axis is
## MIRRORED inside each FOV, so the cell-level slope is NEGATIVE (~ -8314 px/mm, see
## scripts/lib/cosmx_tiff.R). Using the origin fit for cell coordinates put every
## window off the slide and rendered 40 blank panels. Worse, its residual was
## 0.0000 mm -- the FOV origins are exactly linear in mm by construction, so that
## check confirmed something true and irrelevant. The assertion below is the one that
## would actually have caught it.
SEG <- fread("output/10_hrs/hrs_segments.csv.gz",
             select = c("cell_id","samp","lineage","x_slide_mm","y_slide_mm"))
SEG[, key := cosmx_seu_key(cell_id)]
FIT <- lapply(SAMPS, function(S) {
  d <- merge(SEG[samp == S, .(key, x_slide_mm, y_slide_mm)],
             FLATS[[S]][, .(key, X_px, Y_px)], by = "key")
  stopifnot(nrow(d) > 1000)
  fx <- lm(X_px ~ x_slide_mm, data = d); fy <- lm(Y_px ~ y_slide_mm, data = d)
  a <- list(ax = unname(coef(fx)[2]), bx = unname(coef(fx)[1]),
            ay = unname(coef(fy)[2]), by = unname(coef(fy)[1]),
            res = c(sd(resid(fx)), sd(resid(fy))), n = nrow(d))
  exp_slope <- 1 / COSMX_MM_PER_PX          # ~8314 px per mm
  if (a$ax > 0 || a$ay > 0)
    stop(S, ": mm->px slope is positive; the canvas x axis is mirrored so it must be ",
         "negative. Got ax=", round(a$ax), " ay=", round(a$ay), call. = FALSE)
  if (max(abs(abs(c(a$ax, a$ay)) - exp_slope) / exp_slope) > 0.01)
    stop(S, ": mm->px slope ", round(a$ax), "/", round(a$ay), " is more than 1% from the ",
         "expected ", round(exp_slope), " px/mm", call. = FALSE)
  a
}); names(FIT) <- SAMPS
for (S in SAMPS) cat(sprintf("%s cell mm->px affine: ax %.1f ay %.1f | resid %.1f / %.1f px | n=%d\n",
                             S, FIT[[S]]$ax, FIT[[S]]$ay, FIT[[S]]$res[1], FIT[[S]]$res[2],
                             FIT[[S]]$n))
## FIT is parameterised px = a*mm + b, so mm->px is the forward direction and px->mm
## is the inverse. Getting these the wrong way round is what produced the blank panels.
mm2px <- function(S, x_mm, y_mm) {
  a <- FIT[[S]]
  list(X_px = a$ax * x_mm + a$bx, Y_px = a$ay * y_mm + a$by)
}
px2mm <- function(S, X_px, Y_px) {
  a <- FIT[[S]]
  list(x = (X_px - a$bx) / a$ax, y = (Y_px - a$by) / a$ay)
}

## Canvas extent per sample, so an edge window is clamped instead of producing a rect
## with negative bounds (idiom from 07_script_lightchain_views.R:268-272).
CANVAS <- lapply(SAMPS, function(S) {
  pl <- cosmx_fov_placement(S)
  c(x = max(pl$x0) + COSMX_FOV_PX - 1L, y = max(pl$y0) + COSMX_FOV_PX - 1L)
}); names(CANVAS) <- SAMPS

set.seed(1)
LOCBG <- lapply(SAMPS, function(S) {
  f <- FLATS[[S]][sample.int(.N, min(.N, LOC_N)), .(X_px, Y_px)]
  mm <- px2mm(S, f$X_px, f$Y_px); f[, `:=`(x_mm = mm$x, y_mm = mm$y)]; f
}); names(LOCBG) <- SAMPS

.locator <- function(S, rect) {
  q <- rect@coords[1L, ]
  wx <- range(px2mm(S, c(q$xmin, q$xmax), c(q$ymin, q$ymax))$x)
  wy <- range(px2mm(S, c(q$xmin, q$xmax), c(q$ymin, q$ymax))$y)
  ggplot() +
    geom_point(data = LOCBG[[S]], aes(x_mm, y_mm), colour = "#cfcec9", size = 0.04) +
    annotate("rect", xmin = wx[1], xmax = wx[2], ymin = wy[1], ymax = wy[2],
             colour = "#008300", fill = NA, linewidth = 0.6) +
    coord_fixed() +
    labs(title = sprintf("0. slide position, mm\n(%.2f, %.2f)", mean(wx), mean(wy))) +
    theme_void() +
    theme(plot.title = element_text(size = 8),
          plot.background = element_rect(fill = "white", colour = "#b5b3ac", linewidth = .3))
}

## Bin geometry in canvas px, once per sample, for the new heatmap panel. Uses the
## bin's GEOMETRIC bounds (bx*BIN_UM .. (bx+1)*BIN_UM), not the centroid of its member
## cells, so the tiles abut exactly and the grid reads as a grid.
BINPX <- lapply(SAMPS, function(S) {
  b  <- B100[samp == S]
  c1 <- mm2px(S,  b$bx      * BIN_UM/1000,  b$by_      * BIN_UM/1000)
  c2 <- mm2px(S, (b$bx + 1) * BIN_UM/1000, (b$by_ + 1) * BIN_UM/1000)
  data.table(x1 = pmin(c1$X_px, c2$X_px), x2 = pmax(c1$X_px, c2$X_px),
             y1 = pmin(c1$Y_px, c2$Y_px), y2 = pmax(c1$Y_px, c2$Y_px),
             v  = b$chemo_cp10k, cells = b$cells)
}); names(BINPX) <- SAMPS

## ---- 4. per-cell lineage, for panel 2 -------------------------------------------
## Read once from the step-10 segment table, which already carries lineage per cell.
## SEG is already in memory from the affine fit above.
SEG[, lin3 := fifelse(lineage == "B", "B", fifelse(lineage == "T", "T", "other"))]
LINMAP <- lapply(SAMPS, function(S) SEG[samp == S, .(key, lin3)])
names(LINMAP) <- SAMPS

## ---- 5. panels ------------------------------------------------------------------
.cells_plot <- function(S, rect, binbox) {
  poly <- tryCatch(cosmx_cells_in_rect(S, rect), error = function(e) NULL)
  if (is.null(poly)) return(list(p = ggplot() + theme_v + labs(title = "2. cells (none)"),
                                 df = NULL))
  df <- poly@coords %>% dplyr::group_by(name) %>%
        dplyr::reframe(x = unlist(x), y = unlist(y))
  ## Do the match OUTSIDE `[.data.table`. Inside it, `key` resolves to LINMAP's own
  ## `key` COLUMN rather than this local variable, so match(key, LINMAP$key) matches
  ## the column against itself, returns 1:nrow, and hands back the entire table.
  pk  <- cosmx_poly_key(df$name)
  lm  <- LINMAP[[S]]
  lin <- if (nrow(lm)) lm$lin3[match(pk, lm$key)] else rep(NA_character_, length(pk))
  stopifnot(length(lin) == nrow(df))
  df$fill <- unname(LIN_COLS[ifelse(is.na(lin), "other", lin)])
  p <- ggplot(df, aes(x, y, group = name)) +
    geom_polygon(aes(fill = fill), colour = "#5a5a5a", linewidth = 0.18) +
    scale_fill_identity() +
    annotate("rect", xmin = binbox$x1, xmax = binbox$x2,
             ymin = binbox$y1, ymax = binbox$y2,
             colour = "#c1121f", fill = NA, linewidth = 0.8) +
    .coord_px(rect) + scale_y_reverse(labels = .um) +
    labs(x = "x (um)", y = "y (um)",
         title = sprintf("2. cells -- B (blue), T (green); red box = the %g um bin", BIN_UM)) +
    theme_v
  list(p = p, df = df)
}

## Panel 4: the bin score itself, as a spatial heatmap. Sequential = ONE hue,
## light to dark (never a rainbow). log10 because the window spans 4 to ~280 cp10k and
## a linear ramp would flatten everything below the peak into the same pale tile.
## Bins dropped for <20 cells are left blank rather than imputed.
.heat_plot <- function(S, rect, binbox) {
  q <- rect@coords[1L, ]
  d <- BINPX[[S]][x2 >= q$xmin & x1 <= q$xmax & y2 >= q$ymin & y1 <= q$ymax]
  p <- ggplot()
  if (nrow(d))
    p <- p + geom_rect(data = d, aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2,
                                     fill = v + 0.1),
                       colour = "#ffffff", linewidth = 0.12) +
      scale_fill_gradient(low = "#eef3fa", high = "#0b3d7a", trans = "log10",
                          name = "cp10k", na.value = "#f2f1ec")
  p + annotate("rect", xmin = binbox$x1, xmax = binbox$x2,
               ymin = binbox$y1, ymax = binbox$y2,
               colour = "#c1121f", fill = NA, linewidth = 0.9) +
    .coord_px(rect) + scale_y_reverse(labels = .um) +
    labs(x = "x (um)", y = "y (um)",
         title = sprintf("4. CCL17+CCL22 per %g um bin (%d bins in view)",
                         BIN_UM, nrow(d))) + theme_v
}

.genes_plot <- function(gp, df, rect, binbox, scale = 1) {
  p <- ggplot()
  if (!is.null(df) && nrow(df))
    p <- p + geom_polygon(data = df, aes(x, y, group = name), fill = "#f2f1ec", colour = NA)
  n <- setNames(integer(length(FAM_LEVELS)), FAM_LEVELS)
  if (!is.null(gp) && nrow(gp)) {
    d <- as.data.table(gp)
    d[, family := factor(FAM[target], levels = FAM_LEVELS)]
    d <- d[!is.na(family)]
    n[FAM_LEVELS] <- as.integer(table(factor(d$family, levels = FAM_LEVELS)))
    for (fm in FAM_LEVELS) if (nrow(d[family == fm]))
      p <- p + geom_point(data = d[family == fm], aes(X_px, Y_px, colour = family),
                          size = FAM_SIZE[[fm]] * scale, stroke = 0)
    p <- p + scale_colour_manual(values = FAM_COLS, labels = FAM_LABS, name = NULL,
                                 drop = FALSE,
                                 guide = guide_legend(override.aes = list(size = 3)))
  }
  if (!is.null(df) && nrow(df))
    p <- p + geom_polygon(data = df, aes(x, y, group = name), fill = NA,
                          colour = "#33333355", linewidth = 0.12)
  p + annotate("rect", xmin = binbox$x1, xmax = binbox$x2,
               ymin = binbox$y1, ymax = binbox$y2,
               colour = "#c1121f", fill = NA, linewidth = 0.8) +
    .coord_px(rect) + scale_y_reverse(labels = .um) +
    labs(x = "x (um)", y = "y (um)",
         title = sprintf("3. transcripts -- CCL17/22 %d, CD30 %d, T cell %d",
                         n[["chemo"]], n[["CD30"]], n[["Tcell"]])) + theme_v
}

## Transcripts for EVERY FOV a window touches, cached per FOV.
##
## The first version read only the window's modal FOV, which left a hard vertical edge
## across panel 3 wherever the window straddled a seam -- a FOV is 512 um and the
## window is 505 um placed arbitrarily, so it almost always spans 2-4 of them. The
## transcript counts in the title were then counts of a partly-empty window.
## Caching is cheap because cosmx_gene_points returns only these 5 genes' points; the
## ~120 MB is the csv read, which is what the cache avoids repeating across windows.
GPCACHE <- new.env(parent = emptyenv())
gp_for_rect <- function(S, rect) {
  fv <- cosmx_fovs_in_rect(S, rect)
  if (is.null(fv) || !nrow(fv)) return(NULL)
  parts <- lapply(fv$fov, function(f) {
    k <- paste0(S, "_", f)
    if (is.null(GPCACHE[[k]])) {
      row <- fv[fv$fov == f, , drop = FALSE]
      fr  <- cosmx_px_rect(row$x0, row$x0 + COSMX_FOV_PX,
                           row$y0, row$y0 + COSMX_FOV_PX, name = paste0("fov", f))
      g <- tryCatch(cosmx_gene_points(S, fr, GENES), error = function(e) NULL)
      GPCACHE[[k]] <- if (is.null(g)) data.frame() else as.data.frame(g)
    }
    GPCACHE[[k]]
  })
  parts <- parts[vapply(parts, nrow, 0L) > 0L]
  if (!length(parts)) return(NULL)
  do.call(rbind, parts)
}

## ---- 6. render ---------------------------------------------------------------------
sel[, fov := as.integer(fov)]
setorder(sel, samp, fov, -chemo_cp10k)
n_done <- 0L
for (S in SAMPS) {
  for (FV in sort(unique(sel[samp == S, fov]))) {
    if (is.na(FV)) next
    ws <- sel[samp == S & fov == FV]
    pl <- cosmx_fov_placement(S); row <- pl[pl$fov == FV, ]
    if (!nrow(row)) { warning("FOV ", FV, " absent for ", S); next }
    for (r in seq_len(nrow(ws))) {
      w <- ws[r]
      ct <- mm2px(S, w$x_mm, w$y_mm)
      ## the bin itself, as a box in canvas px
      c1 <- mm2px(S, w$bx * BIN_UM/1000, w$by_ * BIN_UM/1000)
      c2 <- mm2px(S, (w$bx+1) * BIN_UM/1000, (w$by_+1) * BIN_UM/1000)
      binbox <- list(x1 = min(c1$X_px, c2$X_px), x2 = max(c1$X_px, c2$X_px),
                     y1 = min(c1$Y_px, c2$Y_px), y2 = max(c1$Y_px, c2$Y_px))
      for (z in names(ZOOMS)) {
        s2 <- TIGHT_PX * ZOOMS[[z]] / 2
        rect <- cosmx_px_rect(max(0, ct$X_px - s2), min(CANVAS[[S]][["x"]], ct$X_px + s2),
                              max(0, ct$Y_px - s2), min(CANVAS[[S]][["y"]], ct$Y_px + s2),
                              name = sprintf("%s_%s", w$id, z))
        nfov <- nrow(cosmx_fovs_in_rect(S, rect))
        img <- tryCatch(cosmx_image(S, rect) %>% add_um_scale_bar(),
                        error = function(e) NULL)
        if (!is.null(img)) invisible(.um_per_px(img)) else invisible(.um_per_px(NULL))
        ## TiffPlotR's colorized plot ALREADY carries scale_y_reverse internally, so
        ## relabelling to um must re-apply scale_y_reverse -- adding scale_y_continuous
        ## replaces the built-in reverse and flips the image against panels 2-4. That
        ## bug shipped once here and once in step 12; the axis directions of panels 1
        ## and 2 must be compared on any regenerated figure.
        p_img <- if (is.null(img)) ggplot() + theme_v + labs(title = "1. image (unavailable)")
                 else img@plots$colorized + .coord_px(rect) +
                      scale_y_reverse(labels = .um) +
                      labs(x = "x (um)", y = "y (um)", title = "1. image") + theme_v
        cw <- .cells_plot(S, rect, binbox)
        ## Transcripts on BOTH zooms. The context view is where the field's structure
        ## beyond the bin is visible, which is the whole point of having it. Points are
        ## shrunk at 3x because the same window carries ~9x the area and would otherwise
        ## solidify into a mat. Cost is bounded by the per-FOV cache: a context window
        ## touches ~9 FOVs, but hotspots cluster so most are already read.
        pg3 <- .genes_plot(.gp_subset(gp_for_rect(S, rect), rect), cw$df, rect, binbox,
                           scale = if (z == "tight") 1 else 0.42)
        lab <- sprintf("%s %s | %s | CCL17+CCL22 %.0f cp10k | CD30 %.2f | T %.0f%% | %d cells",
                       S, w$grp, w$id, w$chemo_cp10k, w$cd30_cp10k, 100*w$f_T, w$cells)
        p_heat <- .heat_plot(S, rect, binbox)
        ## 2x2: image | cells over transcripts | heatmap, with the whole-slide locator
        ## inset into the top-left of the image panel rather than taking its own column.
        pgg <- (((p_img + inset_element(.locator(S, rect), left = 0.0, bottom = 0.60,
                                        right = 0.32, top = 1.0, align_to = "panel")) |
                 cw$p) / (pg3 | p_heat)) +
          plot_annotation(title = lab,
            subtitle = sprintf(paste0("%s -- %.0f um window over a %g um bin, FOV %d (%d in ",
                       "window). Red box = the bin. IB1 median bin is 4.0 cp10k."),
                       z, 2*s2*COSMX_MM_PER_PX*1000, BIN_UM, FV, nfov))
        out_f <- file.path(OUT_DIR, sprintf("%s_%s_%s.png", w$grp, w$id, z))
        ggsave(out_f, plot = pgg, width = 14, height = 12, dpi = 140, limitsize = FALSE)
        cat("wrote ", out_f, "\n", sep = ""); n_done <- n_done + 1L
      }
    }
    invisible(gc())
  }
}
writeLines(format(Sys.time()), file.path(OUT_DIR, "views.done"))
cat("\n============================================================\n")
cat(sprintf("%d panels written to %s\n", n_done, OUT_DIR))
cat("A real hotspot should show a dense CCL17/CCL22 cluster inside the red box, with\n")
cat("T cells around it. The IB1bg panels are the same tissue without the chemokine,\n")
cat("and the RB1top panels are the most chemokine-rich places in a Hodgkin-free node --\n")
cat("if those look the same as the hotspots, the call is not specific.\n")
cat("============================================================\nOK\n")
