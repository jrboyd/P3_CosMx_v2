## ============================================================================
## 07_script_lightchain_global.R -- one whole-tissue summary view, following the
## same four-panel pattern as scripts/07_script_lightchain_views.R.
##
## Differences from the per-neighbourhood views:
##   * ONE fixed window covering the tissue, not a neighbourhood-derived rect and
##     not a zoom series: canvas px x 15,000-100,000 by y 5,000-90,000.
##   * Panel 3 is NOT transcript points. At this scale a per-FOV transcript read is
##     ~400 FOVs (~48 GB of csv) and the dots would be a solid smear. Instead each
##     neighbourhood is coloured by its AGGREGATE kappa:lambda ratio -- the quantity the
##     analysis actually tests -- on a divergent scale about the polytypic expectation.
##
## Panels, all four sharing one orientation:
##   0  slide locator, in mm, with the window boxed
##   1  colorized OME-TIFF image
##   2  cells by neighbourhood assignment, adjacent neighbourhoods distinctly
##      coloured, with the per-neighbourhood convex hulls
##   3  neighbourhoods by their AGGREGATE LOG2 kappa:lambda ratio, red-grey-blue diverging
##      about the polytypic expectation log2(0.62/0.38) = 0.71. Hulls are filled with the
##      ratio and outlined black; cells in no scored neighbourhood are grey. Aggregate,
##      not per cell: one cell's ratio is nearly binary (47.8% pure kappa, 21.3% pure
##      lambda here) so a per-cell map shows dot density rather than clonality. log2
##      rather than a fraction, so a 2x kappa excess and a 2x lambda excess sit equally
##      far from the expectation.
##
## The window spans ~85,000 px (~10 mm), which is far past the point where cell
## polygons and transcript csvs are affordable, so panels 2 and 3 draw cell CENTROIDS.
## That is the same fallback the per-neighbourhood views use above MAX_DETAIL_FOVS.
##
## Usage: Rscript scripts/07_script_lightchain_global.R <sample_id> [radius_um]
##   COSMX_LC_DIR=output/07_lightchain sbatch --mem=96G --cpus-per-task=2 --time=2:00:00 \
##     scripts/sbatch/run_rscript.sh scripts/07_script_lightchain_global.R IB1 20
## ============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (!length(args))
  stop("Usage: Rscript scripts/07_script_lightchain_global.R <sample_id> [radius_um]",
       call. = FALSE)
SAMPLE_ID <- args[1]
RADIUS    <- if (length(args) >= 2) as.integer(args[2]) else 20L

if (!file.exists("scripts/lib/cosmx_tiff.R"))
  stop("run from the project root: cd /gpfs1/pi/avolaric/files_jrboyd/P3_CosMx_v2",
       call. = FALSE)
source("scripts/lib/cosmx_tiff.R")
suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(data.table); library(RANN)
})
setDTthreads(1)

## ------------------ CONFIG ----
IN_DIR  <- Sys.getenv("COSMX_LC_DIR", file.path("output/07_lightchain", SAMPLE_ID))
OUT_DIR <- file.path("output/07_lightchain_global", basename(IN_DIR))
FLAT    <- "/netfiles/volaric_research/DLBCL/CosMx/export1/flatFiles/%s/%s_metadata_file.csv.gz"

## the fixed global window, canvas px
WIN_X <- c(15000, 100000)
WIN_Y <- c(5000, 90000)

## Value scales run blue -> red (CLAUDE.md): red is the high end. Here the value is
## log2(kappa:lambda), so kappa -- the high end -- is red and lambda is blue. This also
## settles a standing inconsistency: kappa-restricted means FL-like, and FL is orange/warm
## everywhere else in the project, but this scale used to render it blue.
DIV_HI <- "#d03b3b"   # kappa, the HIGH end of the ratio
DIV_MID <- "#f0efec"
DIV_LO <- "#184f95"   # lambda, the LOW end
UNINV    <- "#d9d8d3" # not light-chain-high
UNSCORED <- "#a9a7a1" # light-chain-high, neighbourhood too small to score
PLASMA_COL <- "#3b1f47"
HALO_COL   <- "#9b7fa6"
NULL_K   <- 0.62      # polytypic kappa fraction (the ~2:1 norm)
EXP_L2   <- log2(NULL_K / (1 - NULL_K))   # 0.706 -- expected log2(kappa/lambda)
L2_SPAN  <- 3         # divergent scale half-width, symmetric about EXP_L2. Covers the
                      # 1st-99th percentile of the observed deviation (-1.9 to +3.0)
LOC_N    <- 60000L
ADJ_MULT <- 1.5
NB_PAL   <- c("#4e79a7", "#f28e2b", "#59a14f", "#e15759",
              "#b07aa1", "#76b7b2", "#edc948", "#9c755f")
HULL_INFLATE <- 1.03
PT_SIZE  <- 0.16      # one dot per cell across ~10 mm

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
theme_v <- theme(panel.background = element_blank(), panel.grid = element_blank(),
                 legend.position = "bottom", plot.title = element_text(size = 10))

## ------------------ 1. tables ----
nb_f <- file.path(IN_DIR, sprintf("lightchain_neighborhoods_r%dum.csv", RADIUS))
lc_f <- file.path(IN_DIR, sprintf("lightchain_cells_r%dum.csv", RADIUS))
for (f in c(nb_f, lc_f)) if (!file.exists(f)) stop("missing ", f, call. = FALSE)
nbs <- fread(nb_f)
lc  <- fread(lc_f)
if (!"tier" %in% names(lc)) lc[, tier := "B_clean"]
lc[, scored := !is.na(nb) & nb %in% nbs$nb]
lc[, key := cosmx_seu_key(cell_id)]
## per-cell kappa fraction -- this is what panel 3 shows
lc[, k_frac := ifelse(lc_total > 0, k_ct / lc_total, NA_real_)]

## ------------------ 2. exact canvas px for every cell ----
flat <- fread(sprintf(FLAT, SAMPLE_ID, SAMPLE_ID),
              select = c("cell", "fov", "CenterX_local_px", "CenterY_local_px"))
px <- cosmx_local_to_px(SAMPLE_ID, flat$fov, flat$CenterX_local_px, flat$CenterY_local_px)
flat[, `:=`(key = cosmx_seu_key(cell), X_px = px$X_px, Y_px = px$Y_px)]
lc <- merge(lc, flat[, .(key, X_px, Y_px)], by = "key", all.x = TRUE)
cat(nrow(lc), " pool cells, ", sum(is.na(lc$X_px)), " without canvas coords\n", sep = "")

rect <- cosmx_px_rect(WIN_X[1], WIN_X[2], WIN_Y[1], WIN_Y[2], name = "global")
q <- rect@coords[1L, ]
in_win <- function(d) d[X_px >= q$xmin & X_px <= q$xmax & Y_px >= q$ymin & Y_px <= q$ymax]
allc <- in_win(flat[, .(key, X_px, Y_px)])
pool <- in_win(lc[!is.na(X_px)])
cat(sprintf("window %.0f x %.0f px (%.1f x %.1f mm): %s cells, %s in the light-chain pool\n",
            diff(WIN_X), diff(WIN_Y), diff(WIN_X)*COSMX_MM_PER_PX, diff(WIN_Y)*COSMX_MM_PER_PX,
            format(nrow(allc), big.mark = ","), format(nrow(pool), big.mark = ",")))
cat("tiers in window: ",
    paste(names(table(pool$tier)), table(pool$tier), sep = "=", collapse = "  "), "\n", sep = "")

## ------------------ 3. shared orientation and the locator ----
## Same convention as the per-neighbourhood views: canvas-px panels get their x LIMITS
## reversed so they read the same way round as the mm locator, since x_slide_mm and
## canvas X_px run in opposite directions. Reversing limits rather than the scale keeps
## geom_raster (the image panel) happy.
## scale_x_reverse, not reversed coord limits: ggplot2 sorts coord limits, so passing
## them backwards does nothing. Same convention as the per-neighbourhood views.
.coord_px <- function() list(scale_x_reverse(),
                             coord_fixed(xlim = c(q$xmin, q$xmax), ylim = c(q$ymin, q$ymax)))
.coord_win <- function() list(.coord_px(), scale_y_reverse(),
                              labs(x = "pixel", y = "pixel"), theme_v)

PX2MM <- local({
  fx <- lm(x_slide_mm ~ X_px, data = lc[!is.na(X_px)])
  fy <- lm(y_slide_mm ~ Y_px, data = lc[!is.na(X_px)])
  list(ax = unname(coef(fx)[2]), bx = unname(coef(fx)[1]),
       ay = unname(coef(fy)[2]), by = unname(coef(fy)[1]))
})
.px2mm <- function(x, y) list(x = PX2MM$ax * x + PX2MM$bx, y = PX2MM$ay * y + PX2MM$by)

set.seed(1)
loc_bg <- flat[sample.int(.N, min(.N, LOC_N)), .(X_px, Y_px)]
loc_bg[, c("x_mm", "y_mm") := .px2mm(X_px, Y_px)]
wx <- range(.px2mm(c(q$xmin, q$xmax), c(q$ymin, q$ymax))$x)
wy <- range(.px2mm(c(q$xmin, q$xmax), c(q$ymin, q$ymax))$y)
p_loc <- ggplot() +
  geom_point(data = loc_bg, aes(x_mm, y_mm), colour = "#cfcec9", size = 0.04) +
  annotate("rect", xmin = wx[1], xmax = wx[2], ymin = wy[1], ymax = wy[2],
           colour = "#c1121f", fill = NA, linewidth = 0.4) +
  coord_fixed() +
  labs(title = sprintf("0. slide position, mm\n%.1f-%.1f x %.1f-%.1f mm",
                       wx[1], wx[2], wy[1], wy[2])) +
  theme_void() +
  theme(plot.title = element_text(size = 9),
        plot.background = element_rect(fill = "white", colour = "#b5b3ac", linewidth = 0.3))

## ------------------ 4. panel 1: image ----
img   <- cosmx_image(SAMPLE_ID, rect) %>% add_um_scale_bar()
p_img <- img@plots$colorized + .coord_px() + theme_v + labs(title = "1. image")

## ------------------ 5. panel 2: neighbourhood assignment ----
## Same greedy adjacency colouring as the per-neighbourhood views: nominal colours,
## chosen so no two spatially adjacent neighbourhoods share one, and spread across the
## whole palette rather than collapsing onto the first few.
.nb_adjacent_colors <- function(nb, x, y, min_dist_px) {
  lev <- sort(unique(nb)); if (!length(lev)) return(character(0))
  if (length(lev) == 1L) return(stats::setNames(NB_PAL[1], as.character(lev)))
  g <- match(nb, lev)
  k <- min(48L, length(nb))
  nn <- RANN::nn2(cbind(x, y), k = k, searchtype = "radius", radius = min_dist_px)
  idx <- nn$nn.idx
  i <- rep.int(seq_len(nrow(idx)), times = ncol(idx)); j <- as.vector(idx)
  pad <- j > 0L; i <- i[pad]; j <- j[pad]
  ok <- g[i] != g[j]; a <- g[i][ok]; b <- g[j][ok]
  adj <- lapply(split(c(b, a), c(a, b)), unique)
  ord <- order(-vapply(seq_along(lev), function(z) length(adj[[as.character(z)]]), 1L))
  greedy <- function(spread) {
    cols <- integer(length(lev))
    for (z in ord) {
      used <- cols[adj[[as.character(z)]]]
      avail <- setdiff(seq_along(NB_PAL), used[used > 0L])
      if (!length(avail)) return(NULL)
      cols[z] <- if (spread) {
        tab <- tabulate(cols[cols > 0L], nbins = length(NB_PAL)); avail[which.min(tab[avail])]
      } else avail[1L]
    }
    cols
  }
  cols <- greedy(TRUE); if (is.null(cols)) cols <- greedy(FALSE)
  if (is.null(cols))
    return(stats::setNames(grDevices::hcl.colors(length(lev), "Dark 3"), as.character(lev)))
  bad <- sum(vapply(seq_along(lev),
                    function(z) sum(cols[adj[[as.character(z)]]] == cols[z]), 1L))
  if (bad) stop("greedy colouring failed: ", bad, " adjacent pairs share a colour",
                call. = FALSE)
  cat(sprintf("  %d neighbourhoods in window, %d colours\n", length(lev),
              length(unique(cols))))
  stats::setNames(NB_PAL[cols], as.character(lev))
}

.nb_hulls <- function(x, y, nb, col) {
  d <- data.frame(x = x, y = y, nb = nb, col = col, stringsAsFactors = FALSE)
  d <- d[!is.na(d$nb), , drop = FALSE]; if (!nrow(d)) return(NULL)
  parts <- lapply(split(d, d$nb), function(z) {
    if (nrow(z) < 3L) return(NULL)
    h <- grDevices::chull(z$x, z$y); if (length(h) < 3L) return(NULL)
    hx <- z$x[h]; hy <- z$y[h]; cx <- mean(hx); cy <- mean(hy)
    data.frame(x = cx + (hx - cx) * HULL_INFLATE, y = cy + (hy - cy) * HULL_INFLATE,
               nb = z$nb[1], col = z$col[1], stringsAsFactors = FALSE)
  })
  parts <- parts[!vapply(parts, is.null, TRUE)]
  if (!length(parts)) return(NULL)
  do.call(rbind, parts)
}

.hull_layers <- function(hulls, linewidth = 0.25) {
  if (is.null(hulls) || !nrow(hulls)) return(NULL)
  lapply(split(hulls, hulls$col), function(z)
    geom_polygon(data = z, aes(x = x, y = y, group = nb), fill = NA,
                 colour = z$col[1], linewidth = linewidth, inherit.aes = FALSE))
}

sc  <- pool[scored == TRUE]
nb_cols <- .nb_adjacent_colors(sc$nb, sc$X_px, sc$Y_px,
                               ADJ_MULT * (RADIUS / 1000) / COSMX_MM_PER_PX)
sc[, col := unname(nb_cols[as.character(nb)])]
hulls <- .nb_hulls(sc$X_px, sc$Y_px, sc$nb, sc$col)

p_nb <- ggplot() +
  geom_point(data = allc[!key %in% pool$key], aes(X_px, Y_px),
             colour = UNINV, size = PT_SIZE * 0.7) +
  geom_point(data = pool[tier == "B_clean" & scored == FALSE], aes(X_px, Y_px),
             colour = UNSCORED, size = PT_SIZE) +
  geom_point(data = pool[tier == "halo"],   aes(X_px, Y_px), colour = HALO_COL,   size = PT_SIZE) +
  geom_point(data = pool[tier == "plasma"], aes(X_px, Y_px), colour = PLASMA_COL, size = PT_SIZE) +
  geom_point(data = sc, aes(X_px, Y_px, colour = col), size = PT_SIZE * 1.3) +
  scale_colour_identity() +
  .hull_layers(hulls) + .coord_win() +
  labs(title = sprintf("2. %s cells by neighbourhood (%d scored nb; plum = plasma/halo, excluded)",
                       format(nrow(pool), big.mark = ","), length(nb_cols)))

## ------------------ 6. panel 3: kappa fraction per cell ----
## Every cell that carries light chain, coloured by its own k/(k+l). Diverging about
## NULL_K, the polytypic expectation, so red = kappa-skewed, blue = lambda-skewed and
## the reactive background sits at the pale midpoint -- the same colour logic as the
## step-07 whole-slide maps.
## The ratio is aggregated PER NEIGHBOURHOOD, not per cell. A single cell's k/(k+l) is
## nearly binary -- measured over this window, 47.8% of pool cells are pure kappa and
## 21.3% pure lambda -- so a per-cell map saturates at both ends and shows dot density
## rather than clonality. The neighbourhood kappa_frac is the quantity the analysis
## actually tests, and it is what distinguishes a kappa-restricted focus from reactive
## background.
## LOG2 RATIO, not kappa fraction. A fraction compresses exactly where the interesting
## values live: 0.85 vs 0.95 is 0.10 of the scale, while 0.50 vs 0.60 is also 0.10, even
## though the first is a 3x change in the underlying ratio and the second is 1.5x. log2
## makes a 2x kappa excess and a 2x lambda excess equidistant from the expectation.
## +0.5 on both counts (Haldane): n_lambda == 0 does occur (one r20um neighbourhood), and
## a raw ratio would be infinite there.
nbk  <- nbs[, .(nb, nb_l2 = log2((n_kappa + 0.5) / (n_lambda + 0.5)))]
sc3  <- merge(pool[scored == TRUE], nbk, by = "nb")
## Hulls FILLED with the neighbourhood's ratio, outlined black: at 10 mm across a 20-50
## cell neighbourhood is only a few pixels, so the filled hull is what carries the ratio
## and the outline is what makes it findable regardless of where it sits on the scale.
h3 <- .nb_hulls(sc3$X_px, sc3$Y_px, sc3$nb, rep("x", nrow(sc3)))
if (!is.null(h3)) h3 <- merge(h3, nbk, by = "nb")
w3 <- nbk[nb %in% sc3$nb]
cat(sprintf("panel 3: %d scored neighbourhoods, %s member cells\n",
            uniqueN(sc3$nb), format(nrow(sc3), big.mark = ",")))
cat(sprintf("           log2(k/l): median %.2f (%.2f:1), range %.2f to %.2f | expected %.2f (%.2f:1)\n",
            median(w3$nb_l2), 2^median(w3$nb_l2), min(w3$nb_l2), max(w3$nb_l2),
            EXP_L2, 2^EXP_L2))
cat(sprintf("           above expectation: %d | below: %d | beyond the +/-%g span (squished): %d\n",
            sum(w3$nb_l2 > EXP_L2), sum(w3$nb_l2 < EXP_L2),
            L2_SPAN, sum(abs(w3$nb_l2 - EXP_L2) > L2_SPAN)))

p_kf <- ggplot() +
  geom_point(data = allc[!key %in% pool$key], aes(X_px, Y_px),
             colour = UNINV, size = PT_SIZE * 0.7) +
  geom_point(data = pool[scored == FALSE], aes(X_px, Y_px),
             colour = UNSCORED, size = PT_SIZE * 0.7) +
  ## Black outline on the hulls: a polytypic neighbourhood sits at the pale midpoint of
  ## the divergent scale and is otherwise nearly invisible against the grey background,
  ## so the outline is what makes every neighbourhood findable regardless of its ratio.
  {if (!is.null(h3)) geom_polygon(data = h3, aes(x = x, y = y, group = nb, fill = nb_l2),
                                  colour = "black", linewidth = 0.25,
                                  inherit.aes = FALSE)} +
  ## Limits SYMMETRIC about the expectation, or the divergent scale would give one
  ## direction more colour range than the other. Breaks are labelled as ratios rather
  ## than log2 values, since that is what the numbers mean.
  scale_fill_gradient2(low = DIV_LO, mid = DIV_MID, high = DIV_HI,
                       midpoint = EXP_L2,
                       limits = EXP_L2 + c(-L2_SPAN, L2_SPAN),
                       oob = scales::squish,
                       breaks = c(-2, 0, EXP_L2, 2, 3),
                       labels = c("1:4", "1:1", "1.6:1\nexpected", "4:1", "8:1"),
                       name = "kappa:lambda per neighbourhood") +
  ## No member-cell points. They carried the same value as the hull they sit inside, so
  ## they added nothing, and at this scale they overplotted the outlines that make the
  ## neighbourhoods findable. The grey context points below stay -- they are drawn before
  ## the hulls, so they sit underneath and do not compete with the outlines.
  .coord_win() +
  labs(title = sprintf("3. %d neighbourhoods by aggregate log2 kappa:lambda (red = kappa, blue = lambda)",
                       uniqueN(sc3$nb)))

## ------------------ 7. assemble ----
pg <- p_loc + p_img + p_nb + p_kf +
  ## 2x2, not a row of four -- see the note in 11_script_hodgkin_niche_global.R.
  plot_layout(nrow = 2, ncol = 2) +
  plot_annotation(
    title = sprintf("%s global light-chain summary, r%dum neighbourhoods", SAMPLE_ID, RADIUS),
    subtitle = sprintf(paste0("window x %s-%s, y %s-%s canvas px (%.1f x %.1f mm) | ",
                              "cell centroids, not polygons, at this scale | ",
                              "panel 3 ratio is per NEIGHBOURHOOD | kappa %s, lambda %s, midpoint %.2f"),
                       format(WIN_X[1], big.mark=","), format(WIN_X[2], big.mark=","),
                       format(WIN_Y[1], big.mark=","), format(WIN_Y[2], big.mark=","),
                       diff(WIN_X)*COSMX_MM_PER_PX, diff(WIN_Y)*COSMX_MM_PER_PX,
                       DIV_HI, DIV_LO, NULL_K))
out_f <- file.path(OUT_DIR, sprintf("%s_global_r%dum.png", SAMPLE_ID, RADIUS))
## 300 dpi: the field is ~85,000 px wide rendered into a ~7 inch panel, so thin strokes
## (the hull outlines) were being resampled away at 150. Costs ~4x the file size.
ggsave(out_f, plot = pg, width = 14, height = 16, dpi = 300, limitsize = FALSE)
cat("wrote ", out_f, "\n", sep = "")
writeLines(format(Sys.time()), file.path(OUT_DIR, sprintf("%s_r%dum.done", SAMPLE_ID, RADIUS)))
