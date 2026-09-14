## ============================================================================
## 07_script_lightchain_views.R -- image views of the light-chain neighbourhoods
## called in step 04, for comparing putative FL-enriched kappa-dominant foci
## against "normal" (polytypic) light-chain neighbourhoods.
##
## sbatch entry point -- submit it with scripts/07_submit_lightchain_views.sh,
## which records the arguments that produced the current figures.
##
## Four panels per neighbourhood:
## All three canvas-px panels have their x axis REVERSED (via the coord limits, see
## .coord_px) so they read the same way round as the mm-space locator. x_slide_mm and
## canvas X_px run in opposite directions, so without this the locator is mirrored
## against the panels beside it.
##
##   0  slide locator -- the whole slide with this window boxed, so every figure can
##      be placed without cross-referencing anything. Drawn in SLIDE MM to match the
##      orientation of step 07's whole-slide maps; x_slide_mm runs opposite to canvas
##      X_px, so a px-space locator would be mirrored against them.
##   1  colorized OME-TIFF image
##   2  cell boundaries filled by NEIGHBOURHOOD assignment, plus a convex hull per
##      neighbourhood (inflated 3% about its centre, same colour as its cells) so a
##      scattered set of cells reads as one object. Colours are nominal, chosen by a
##      greedy graph colouring so adjacent neighbourhoods never share one and the
##      whole 8-colour palette gets used; the focal neighbourhood is outlined in
##      black. Cells in a neighbourhood too small to be scored are MID grey; cells
##      that are not light-chain-high at all are light grey. PLASMA cells and their
##      spillover HALO -- both dropped before the graph was built -- get their own plum
##      fills, so an excluded cell never reads as background. Every cell edge gets a
##      black underlay under a thinner white line, so boundaries stay visible against
##      the light fills.
##   3  IGKC / IGLC1/2 transcript points in PALER kappa/lambda tints. Layered bottom
##      to top: light grey FILL of the scored-neighbourhood cells, then the
##      transcripts, then every cell's HOLLOW outline, then the neighbourhood hulls --
##      so the fill never masks a dot and the edges stay crisp over them. No image
##      behind it; that made the dots unreadable.
##
## The neighbourhood hulls from panel 2 are overplotted on panels 1 and 3 as well,
## so all three share one frame of reference. They are added as plain layers with
## literal colours rather than a colour aesthetic, which is what lets them coexist
## with the image's fill scale and the gene colour scale.
##
## Each at three zooms: tight to the neighbourhood, 3x out, and 10x out (the tight
## window is outlined on both zoomed-out views).
##
## At 10x the window can span far more tissue than polygons and per-FOV transcript
## csvs can cover -- r60um at 10x is 164 FOVs, ~19 GB of csv and ~656k polygons. So
## panels 2 and 3 fall back to CENTROID POINTS above MAX_DETAIL_FOVS, keeping the
## same colours and the same question at cell rather than transcript level. Every
## figure states in its subtitle which rendering it used.
##
## One radius per invocation, written to its own subdirectory, so the radii can be
## compared side by side:
##   output/07_lightchain_views/r<R>um/<group>_nb<id>_<zoom>.png
##
##   for R in 20 30 40 60; do
##     sbatch --mem=64G --cpus-per-task=2 --time=2:00:00 --job-name=lcv-$R \
##       scripts/sbatch/run_rscript.sh scripts/07_script_lightchain_views.R IB1 $R
##   done
##
## Usage: Rscript scripts/07_script_lightchain_views.R <sample_id> <radius_um> [n_each] [focus_xy]
##
## Two selection modes:
##   default          the most kappa-skewed neighbourhood at this radius, plus a
##                    polytypic control sitting on the expected ~2:1 ratio.
##   focus_xy="x,y"   the neighbourhood nearest that slide-mm point, so the SAME
##                    physical focus can be followed across radii.
##
## The two modes answer different questions, and the difference matters: the most
## kappa-skewed neighbourhood is a DIFFERENT object at every radius (1.0-3.0 mm
## apart, against spans of 0.09-0.30 mm), so the default mode compares scales, not
## the same tissue. To judge whether a focus is real or a scale artefact, pin it
## with focus_xy and see whether it survives:
##   for R in 20 30 40 60; do
##     sbatch ... scripts/07_script_lightchain_views.R IB1 $R 1 7.424,12.953
##   done
##
## Coordinates: neighbourhood windows are built from the EXACT per-cell route
## (flat-file CenterX/Y_local_px -> cosmx_local_to_px), not from the fitted mm
## affine, so no Seurat object is loaded and no fit residual enters the geometry.
## The neighbourhood table's own x/y/span_mm are used only as a cross-check.
## ============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L)
  stop("Usage: Rscript scripts/07_script_lightchain_views.R <sample_id> <radius_um> [n_each] [focus_xy]",
       call. = FALSE)
SAMPLE_ID <- args[1]
RADIUS    <- as.integer(args[2])
N_EACH    <- if (length(args) >= 3) as.integer(args[3]) else 1L
FOCUS_XY  <- if (length(args) >= 4) as.numeric(strsplit(args[4], ",")[[1]]) else NULL
if (!is.null(FOCUS_XY) && length(FOCUS_XY) != 2L)
  stop("focus_xy must be \"x,y\" in slide mm", call. = FALSE)

if (!file.exists("scripts/lib/cosmx_tiff.R"))
  stop("run from the project root: cd /gpfs1/pi/avolaric/files_jrboyd/P3_CosMx_v2",
       call. = FALSE)
source("scripts/lib/cosmx_tiff.R")
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(data.table)
})
setDTthreads(1)

## ------------------ CONFIG ----
## Step-04 output to read. Override with COSMX_LC_DIR to point at a re-run of
## lightchain_neighborhoods.R (e.g. one built from a pre-QC object), which is the
## only way to change which cells carry a call/neighbourhood -- the boundary
## polygons come from the raw per-FOV CellBoundaries csvs and are unaffected.
IN_DIR  <- Sys.getenv("COSMX_LC_DIR", file.path("output/07_lightchain", SAMPLE_ID))
## Optional parallel set: restrict the candidate neighbourhoods to those sitting in
## FOVs that pass step 02's spatial verdict. COSMX_FOV_QC is the minimum fraction of a
## neighbourhood's member cells that must lie in passing FOVs:
##   COSMX_FOV_QC=1     every member cell in a passing FOV (strict; the default sense
##                      of "inside a QC-passing FOV", and what the shipped set used)
##   COSMX_FOV_QC=0.9   relax to 90%, letting a neighbourhood straddle one bad FOV
## Any other non-empty value enables the gate at 0.9. Output lands in a sibling
## "<in_dir>_fovqc" tree, so the unrestricted set is preserved for comparison.
FOV_QC_ENV <- Sys.getenv("COSMX_FOV_QC", "")
FOV_QC     <- nzchar(FOV_QC_ENV)
FOVQC_MIN_FRAC <- if (FOV_QC && !is.na(suppressWarnings(as.numeric(FOV_QC_ENV))) &&
                      as.numeric(FOV_QC_ENV) <= 1) as.numeric(FOV_QC_ENV) else 0.9
FOV_QC_CSV <- "output/02_filter_QC/qc_spatial_fov.csv"
## Keyed by which step-04-style input was used, so figures from the QC-filtered
## tables (04_merge) and the pre-QC-filter tables (07_lightchain) coexist.
OUT_DIR <- file.path("output/07_lightchain_views",
                     paste0(basename(IN_DIR), if (FOV_QC) "_fovqc" else ""),
                     sprintf("r%dum", RADIUS))
FLAT    <- "/netfiles/volaric_research/DLBCL/CosMx/export1/flatFiles/%s/%s_metadata_file.csv.gz"

## Palette and thresholds carried over from scripts/sbatch/lightchain_neighborhoods.R
## so these figures match the step-04 plots.
DIV_HI <- "#184f95"   # kappa
DIV_LO <- "#d03b3b"   # lambda
INK2   <- "#52514e"   # tie -- dark grey
UNINV  <- "#d9d8d3"   # not light-chain-high -- light grey
UNSCORED <- "#a9a7a1" # light-chain-high but in a neighbourhood too small to score
## Cells EXCLUDED from the graph because they are plasma cells or sit in the spillover
## halo of one. A plum family, deliberately outside NB_PAL and outside the greys, so
## "excluded as plasma" never reads as "a neighbourhood" or as "background".
PLASMA_COL <- "#3b1f47"
HALO_COL   <- "#9b7fa6"
GENE_CELL_FILL <- "#d9d8d3"   # panel 3: fill for cells in a scored neighbourhood
GENE_CELL_LINE <- "#33333355"  # panel 3: cell outline drawn over the transcripts
LOC_N    <- 60000L    # cells subsampled for the slide-locator panel
NULL_K    <- 0.62     # polytypic kappa fraction (normal ~2:1 kappa:lambda)
## NO size floor here. The producing script already applied one (MIN_CELLS, recorded in
## lightchain_call_params.csv) when it decided which communities to score, so re-applying
## a second copy of that threshold here only creates a constant that can drift -- and it
## did: lowering the calling floor to 20 left this at 30, silently hiding every
## neighbourhood of 20-29 cells from view selection, including the top focus at 24.

## Paler tints of the step-04 kappa/lambda hues, for the transcript panel only.
## Full-strength DIV_HI/DIV_LO competed with the neighbourhood hull outlines drawn
## over the same panel; these keep the hue identity but sit back visually.
GENE_HI <- "#7fa6cf"   # paler kappa / IGKC
GENE_LO <- "#e59a9a"   # paler lambda / IGLC1/2
CALL_COLS <- c(kappa = GENE_HI, lambda = GENE_LO, tie = "#a8a6a1")
GENE_COLS <- c(IGKC = GENE_HI, "IGLC1/2" = GENE_LO)

PAD        <- 1.30    # tight window = neighbourhood bbox * this
MIN_SIDE   <- 700     # px, so a very tight focus is still legible
ZOOMS      <- c(tight = 1, zoom3x = 3, zoom10x = 10)

## Detail budget. Per-FOV transcript csvs are ~120 MB each and cell boundaries are
## drawn as polygons, so both scale with window AREA. Measured for IB1: r20um at
## 10x spans 9 FOVs (fine), but r60um at 10x spans 164 -- ~19 GB of transcript csv
## and ~656k polygons, which is neither renderable nor readable.
## Above this many FOVs the two cell panels fall back to CENTROID POINTS drawn from
## data already in memory: panel 2 keeps neighbourhood colours, panel 3 switches
## from transcripts to light-chain-high cells coloured by call. Same question, same
## palette, aggregated one level up -- at 10x a cell is sub-pixel anyway.
## Observed on IB1: every tight window and every 3x window up to r40um keeps full
## detail; 10x keeps it up to r40um (12-19 FOVs). Only r60um crosses the line --
## its 3x is 29-40 FOVs and its 10x is 165-216, depending on how large the selected
## neighbourhood is. The subtitle of every figure states which rendering it used, so
## this is never silent.
MAX_DETAIL_FOVS <- 30

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

theme_v <- theme(panel.background = element_blank(),
                 panel.grid       = element_blank(),
                 legend.position  = "bottom",
                 plot.title       = element_text(size = 10))

## ------------------ 1. pick the neighbourhoods ----
nb_f <- file.path(IN_DIR, sprintf("lightchain_neighborhoods_r%dum.csv", RADIUS))
lc_f <- file.path(IN_DIR, sprintf("lightchain_cells_r%dum.csv", RADIUS))
for (f in c(nb_f, lc_f))
  if (!file.exists(f)) stop("missing ", f, call. = FALSE)

nbs <- fread(nb_f)
elig <- nbs[artifact_suspect == FALSE]
if (!nrow(elig)) stop("no eligible neighbourhoods at r", RADIUS, "um", call. = FALSE)

if (FOV_QC) {
  ## step 02's own FOV-level spatial verdict (02_spatial_QC.R:366): a whole FOV is
  ## rejected by fov_flag, and dim staining fails too. For IB1 only 145 of 758 FOVs
  ## pass, so this is a heavy restriction and is reported rather than assumed.
  if (!file.exists(FOV_QC_CSV)) stop("missing ", FOV_QC_CSV, call. = FALSE)
  fq <- fread(FOV_QC_CSV)[sample_id == SAMPLE_ID]
  fq[, fov_pass := fov_flag == "ok" & qc_spatial_staining %in% c("ok", "not_evaluated")]
  ok_fov <- fq[fov_pass == TRUE]$fov
  cat(sprintf("FOV QC: %d of %d %s FOVs pass step 02's spatial verdict\n",
              length(ok_fov), nrow(fq), SAMPLE_ID))

  lc_all <- fread(lc_f)
  lc_all[, fov := cosmx_fov_from_id(cell_id)]
  frac <- lc_all[nb %in% elig$nb, .(frac_fov_ok = mean(fov %in% ok_fov)), by = nb]
  elig <- merge(elig, frac, by = "nb")
  n_before <- nrow(elig)
  elig <- elig[frac_fov_ok >= FOVQC_MIN_FRAC]
  cat(sprintf("FOV QC: %d of %d eligible neighbourhoods have >= %.0f%% of members in passing FOVs\n",
              nrow(elig), n_before, 100 * FOVQC_MIN_FRAC))
  if (!nrow(elig))
    stop("no neighbourhoods survive the FOV QC gate at r", RADIUS, "um", call. = FALSE)
  rm(lc_all); invisible(gc())
}

MAX_FOCUS_MM <- 0.25   # how far from focus_xy still counts as the same focus

if (!is.null(FOCUS_XY)) {
  ## Follow one physical focus across radii. Failing to find it is a RESULT, not
  ## an error: it means this radius does not recover the focus at all.
  elig[, dist_mm := sqrt((x - FOCUS_XY[1])^2 + (y - FOCUS_XY[2])^2)]
  near <- elig[dist_mm <= MAX_FOCUS_MM][order(dist_mm)]
  if (!nrow(near)) {
    cat(sprintf("\n=== r%dum: NO neighbourhood within %.2f mm of (%.3f, %.3f).\n",
                RADIUS, MAX_FOCUS_MM, FOCUS_XY[1], FOCUS_XY[2]))
    cat("    The focus is not recovered at this radius -- nothing to render.\n")
    nearest <- elig[order(dist_mm)][1]
    cat(sprintf("    nearest eligible neighbourhood is nb %d at %.3f mm (kappa_frac %.3f, %s)\n",
                nearest$nb, nearest$dist_mm, nearest$kappa_frac, nearest$class))
    fwrite(data.table(radius_um = RADIUS, focus_x = FOCUS_XY[1], focus_y = FOCUS_XY[2],
                      recovered = FALSE, nearest_nb = nearest$nb,
                      nearest_dist_mm = nearest$dist_mm),
           file.path(OUT_DIR, "focus_not_recovered.csv"))
    quit(save = "no", status = 0)
  }
  sel <- near[seq_len(min(N_EACH, .N))]
  sel[, group := "focus"]
} else {
  ## kappa-dominant: the most kappa-skewed foci, whether or not they cleared the
  ## FL-like threshold at this radius -- the class is reported on every figure so
  ## a radius that loses the call is visible rather than silently absent.
  kappa_sel <- elig[order(-kappa_frac)][seq_len(min(N_EACH, .N))]
  kappa_sel[, group := "kappa"]

  ## normal: polytypic neighbourhoods sitting on the expected ~2:1 ratio
  poly <- elig[class == "polytypic (reactive-like)"]
  if (!nrow(poly)) poly <- elig[class == "intermediate"]
  normal_sel <- poly[order(abs(kappa_frac - NULL_K))][seq_len(min(N_EACH, .N))]
  normal_sel[, group := "normal"]

  sel <- rbind(kappa_sel, normal_sel)
}
cat("=== r", RADIUS, "um: ", nrow(elig), " eligible neighbourhoods, ",
    nrow(sel), " selected\n", sep = "")
print(sel[, .(group, nb, n_cells, n_kappa, n_lambda, n_tie, kappa_frac,
              q = signif(q, 3), span_mm, class)])
## keep the two modes' records separate so neither clobbers the other
fwrite(sel, file.path(OUT_DIR, if (is.null(FOCUS_XY)) "selected_neighborhoods.csv"
                               else "selected_focus.csv"))

## Canvas bounds. A 10x window is often wider than the slide, and an unclamped
## rect spends most of the figure on blank space outside the image.
CANVAS <- local({
  pl <- cosmx_fov_placement(SAMPLE_ID)
  c(x = max(pl$x0) + COSMX_FOV_PX - 1L, y = max(pl$y0) + COSMX_FOV_PX - 1L)
})
cat(sprintf("canvas: 0..%d x 0..%d px\n", CANVAS[["x"]], CANVAS[["y"]]))

## ------------------ 2. exact canvas px for every light-chain-high cell ----
lc <- fread(lc_f)
if (!all(lc$call %in% names(CALL_COLS)))
  stop("unexpected call values: ",
       paste(setdiff(unique(lc$call), names(CALL_COLS)), collapse = ", "), call. = FALSE)
lc$key <- cosmx_seu_key(lc$cell_id)
## Only SCORED neighbourhoods get a colour and a hull. Measured on the pre-QC pool
## at r20um: 57,043 neighbourhoods exist but the median size is ONE cell and 90.5%
## have fewer than 5, so colouring every one buried the real structure under tens of
## thousands of one-cell hulls. Cells in an unscored neighbourhood are drawn mid grey
## -- still visibly light-chain-high, but not pretending to be a neighbourhood.
## Older table sets (output/04_merge) predate the plasma tiering; treat them as all-B.
if (!"tier" %in% names(lc)) lc[, tier := "B_clean"]
lc[, scored := !is.na(nb) & nb %in% nbs$nb]
## Four tiers. plasma/halo were dropped before the graph, so they can never be scored;
## giving them their own palette keys keeps that visible instead of showing them as
## background.
lc[, nb_chr := ifelse(tier == "plasma", "plasma",
              ifelse(tier == "halo",   "halo",
              ifelse(scored, as.character(nb), "unscored")))]
cat(sprintf("tiers in this table: %s\n",
            paste(names(table(lc$tier)), table(lc$tier), sep="=", collapse="  ")))
cat(sprintf("%d of %d light-chain-high cells sit in one of the %d scored neighbourhoods\n",
            sum(lc$scored), nrow(lc), length(unique(nbs$nb))))

flat <- fread(sprintf(FLAT, SAMPLE_ID, SAMPLE_ID),
              select = c("cell", "fov", "CenterX_local_px", "CenterY_local_px"))
px <- cosmx_local_to_px(SAMPLE_ID, flat$fov,
                        flat$CenterX_local_px, flat$CenterY_local_px)
flat[, `:=`(key = cosmx_seu_key(cell), X_px = px$X_px, Y_px = px$Y_px)]

lc <- merge(lc, flat[, .(key, X_px, Y_px)], by = "key", all.x = TRUE)
n_lost <- sum(is.na(lc$X_px))
cat(nrow(lc), " light-chain-high cells, ", n_lost, " without canvas coords\n", sep = "")
if (n_lost) warning(n_lost, " light-chain-high cells did not join the flat file")

## Slide-locator background: a thin subsample of all cells, so every figure can show
## where its window sits on the whole slide.
##
## Drawn in SLIDE MM, not canvas px. x_slide_mm and X_px run in OPPOSITE directions
## (the mm->px affine has ax ~ -8314), so a px-space whole-slide view comes out
## mirrored in x against every other whole-slide figure in the project -- including
## step 07's own fig_07_lc_nb_map / fig_07_lc_cells_map. Panels 1-3 stay in px because
## they are windows and their axes say "pixel"; panel 0 is mm and says so.
##
## The px->mm map is fitted from the light-chain cells, which carry both coordinate
## systems, so no Seurat object or separate affine is needed.
PX2MM <- local({
  fx <- lm(x_slide_mm ~ X_px, data = lc)
  fy <- lm(y_slide_mm ~ Y_px, data = lc)
  out <- list(ax = unname(coef(fx)[2]), bx = unname(coef(fx)[1]),
              ay = unname(coef(fy)[2]), by = unname(coef(fy)[1]),
              resid_mm = c(sd(resid(fx)), sd(resid(fy))))
  cat(sprintf("px->mm for the locator: resid SD %.4f / %.4f mm over %d cells\n",
              out$resid_mm[1], out$resid_mm[2], nrow(lc)))
  if (max(out$resid_mm) > 0.01)
    warning("px->mm fit residual exceeds 10 um; locator positions may be off")
  out
})
.px2mm <- function(x, y) list(x = PX2MM$ax * x + PX2MM$bx,
                              y = PX2MM$ay * y + PX2MM$by)

set.seed(1)
loc_bg <- flat[sample.int(.N, min(.N, LOC_N)), .(X_px, Y_px)]
loc_bg[, c("x_mm", "y_mm") := .px2mm(X_px, Y_px)]

## ------------------ helpers ----

## Subset gene points, keeping the palette attribute that cosmx_add_genes needs.
.gp_subset <- function(gp, rect) {
  if (is.null(gp)) return(NULL)
  pal <- attr(gp, "palette")
  q <- rect@coords[1L, ]
  out <- gp[gp$X_px >= q$xmin & gp$X_px <= q$xmax &
            gp$Y_px >= q$ymin & gp$Y_px <= q$ymax, , drop = FALSE]
  if (!nrow(out)) return(NULL)
  attr(out, "palette") <- pal
  out
}

## Neighbourhood colours are nominal, so they only have to separate neighbours
## from each other. NB_PAL is cycled by a greedy graph colouring; two
## neighbourhoods whose cells come within ADJ_MULT * radius of each other are
## forbidden the same colour. Greedy colouring on this sort of near-planar graph
## normally settles in 4-5 colours, so the palette is comfortably large.
## 8 well-separated hues. The greedy pass below spreads assignment across all of
## them rather than collapsing onto the first few, so panel 2 uses 8 colours where
## the graph allows it -- minimising colour count made every panel look like 3.
NB_PAL   <- c("#4e79a7", "#f28e2b", "#59a14f", "#e15759",
              "#b07aa1", "#76b7b2", "#edc948", "#9c755f")
ADJ_MULT <- 1.5

## nb -> colour, such that spatially adjacent neighbourhoods differ.
.nb_adjacent_colors <- function(nb, x, y, min_dist_px) {
  lev <- sort(unique(nb))
  if (!length(lev)) return(character(0))
  if (length(lev) == 1L) return(stats::setNames(NB_PAL[1], as.character(lev)))
  g <- match(nb, lev)

  ## Adjacency: any two cells from different neighbourhoods within min_dist_px.
  ## A radius search, not k-nearest: with a fixed k, all k neighbours of a cell in
  ## a dense region can belong to its own neighbourhood, so genuine cross-
  ## neighbourhood contacts go undetected.
  k   <- min(48L, length(nb))
  nn  <- RANN::nn2(cbind(x, y), k = k, searchtype = "radius", radius = min_dist_px)
  idx <- nn$nn.idx
  ## index off the returned matrix, not off length(nb), so i and j cannot end up
  ## different lengths and silently recycle
  stopifnot(nrow(idx) == length(nb))
  i  <- rep.int(seq_len(nrow(idx)), times = ncol(idx))
  j  <- as.vector(idx)                # 0 padding where there is no further hit
  stopifnot(length(i) == length(j))
  ## Drop the padding BEFORE comparing. Indexing with 0 drops elements rather than
  ## yielding NA, so g[j] would come back shorter than g[i] and recycle silently.
  pad <- j > 0L
  i <- i[pad]; j <- j[pad]
  ok <- g[i] != g[j]
  a <- g[i][ok]; b <- g[j][ok]
  ## Symmetrise. radius/kNN hit lists are directional -- b can appear in a's list
  ## without a appearing in b's -- and a one-sided list lets the greedy pass colour
  ## b without ever seeing a. That is what the assertion below used to catch.
  adj <- lapply(split(c(b, a), c(a, b)), unique)

  ## Greedy, highest-degree first. `spread` decides which of the allowed colours to
  ## take: the least-used so far, which pushes usage out across the whole palette,
  ## or the lowest index, which minimises how many colours get used at all.
  ord <- order(-vapply(seq_along(lev),
                       function(z) length(adj[[as.character(z)]]), 1L))
  greedy <- function(spread) {
    cols <- integer(length(lev))
    for (z in ord) {
      used  <- cols[adj[[as.character(z)]]]
      avail <- setdiff(seq_along(NB_PAL), used[used > 0L])
      if (!length(avail)) return(NULL)          # palette exhausted for this node
      cols[z] <- if (spread) {
        tab <- tabulate(cols[cols > 0L], nbins = length(NB_PAL))
        avail[which.min(tab[avail])]
      } else avail[1L]
    }
    cols
  }
  ## Spread first, so the panel actually uses the whole palette. Minimising is the
  ## fallback because it needs fewer colours and so survives denser graphs.
  cols <- greedy(spread = TRUE)
  if (is.null(cols)) cols <- greedy(spread = FALSE)
  ## Dense windows can genuinely need more colours than NB_PAL holds. Rather than
  ## fail, give every neighbourhood its own colour -- less distinguishable, but the
  ## adjacent-differ invariant still holds.
  if (is.null(cols)) {
    message("    NB_PAL exhausted; one colour per neighbourhood instead")
    return(stats::setNames(grDevices::hcl.colors(length(lev), "Dark 3"),
                           as.character(lev)))
  }
  ## the whole point of the colouring, so check it rather than assume it
  bad <- sum(vapply(seq_along(lev), function(z) {
    sum(cols[adj[[as.character(z)]]] == cols[z])
  }, 1L))
  if (bad) stop("greedy colouring failed: ", bad,
                " adjacent neighbourhood pairs share a colour", call. = FALSE)

  n_used <- length(unique(cols))
  cat(sprintf("    %d neighbourhoods in window, %d colours, %d adjacent pairs (<= %.0f px)\n",
              length(lev), n_used, length(unlist(adj)), min_dist_px))
  stats::setNames(NB_PAL[cols], as.character(lev))
}

## Convex hull per neighbourhood, inflated about its own centre and drawn in that
## neighbourhood's colour. Scattered same-coloured cells do not read as one object;
## the hull is what makes the neighbourhood visible as a unit. Grouping is by nb, not
## by colour -- non-adjacent neighbourhoods deliberately share colours.
HULL_INFLATE <- 1.03

.nb_hulls <- function(x, y, nb, col, inflate = HULL_INFLATE) {
  d <- data.frame(x = x, y = y, nb = nb, col = col, stringsAsFactors = FALSE)
  d <- d[!is.na(d$nb), , drop = FALSE]
  if (!nrow(d)) return(NULL)
  parts <- lapply(split(d, d$nb), function(z) {
    if (nrow(z) < 3L) return(NULL)
    h <- grDevices::chull(z$x, z$y)
    if (length(h) < 3L) return(NULL)          # collinear or degenerate
    hx <- z$x[h]; hy <- z$y[h]
    cx <- mean(hx); cy <- mean(hy)            # centre of the hull, not of the cells
    data.frame(x = cx + (hx - cx) * inflate,
               y = cy + (hy - cy) * inflate,
               nb = z$nb[1], col = z$col[1], stringsAsFactors = FALSE)
  })
  parts <- parts[!vapply(parts, is.null, TRUE)]
  if (!length(parts)) return(NULL)
  do.call(rbind, parts)
}

## Where this window sits on the slide. Drawn as its own slim panel rather than an
## inset, because a tight window is ~2% of the slide width and vanishes as an inset.
## The red box is the window, and is the only marker -- a centre ring alongside it
## was redundant.
.locator_plot <- function(rect) {
  q  <- rect@coords[1L, ]
  ## rect corners in mm; ax/ay are negative so min/max swap, hence range()
  wx <- range(.px2mm(c(q$xmin, q$xmax), c(q$ymin, q$ymax))$x)
  wy <- range(.px2mm(c(q$xmin, q$xmax), c(q$ymin, q$ymax))$y)
  cx <- mean(wx); cy <- mean(wy)
  ggplot() +
    geom_point(data = loc_bg, aes(x_mm, y_mm), colour = "#cfcec9", size = 0.04) +
    annotate("rect", xmin = wx[1], xmax = wx[2], ymin = wy[1], ymax = wy[2],
             colour = "#c1121f", fill = NA, linewidth = 0.4) +
    coord_fixed() +
    labs(title = sprintf("0. slide position, mm\n(%.1f, %.1f) mm", cx, cy)) +
    theme_void() +
    theme(plot.title = element_text(size = 9),
          plot.background = element_rect(fill = "white", colour = "#b5b3ac",
                                         linewidth = 0.3))
}

## Hull outlines as plain layers, one per distinct colour with the colour passed as
## a literal parameter rather than an aesthetic. That keeps them free of any colour
## scale, so the same hulls can be overplotted on the image panel (which owns a fill
## scale) and on the transcript panel (which owns a colour scale for the genes).
.hull_layers <- function(hulls, linewidth = 0.6) {
  if (is.null(hulls) || !nrow(hulls)) return(NULL)
  lapply(split(hulls, hulls$col), function(z)
    geom_polygon(data = z, aes(x = x, y = y, group = nb),
                 fill = NA, colour = z$col[1], linewidth = linewidth,
                 inherit.aes = FALSE))
}

## Panel 2 as a standalone polygon plot rather than an overlay on the image: a
## solid fill over the image washes both out, and shape_annotate() cannot surface
## the value -> colour mapping anyway. Panel 1 shows the same window.
## No legend -- neighbourhood ids are arbitrary and there are dozens per window;
## the focal one is outlined instead.
.nb_fill_plot <- function(poly_nb, rect, focal_names) {
  df <- poly_nb@coords %>% group_by(name) %>% reframe(x = unlist(x), y = unlist(y))
  df <- merge(df, poly_nb@meta[, c("name", "fill")], by = "name")
  ## nb per polygon, for hull grouping. Cannot group on fill: non-adjacent
  ## neighbourhoods share colours by design.
  m <- match(cosmx_poly_key(df$name), lc$key)
  df$nb <- lc$nb[m]
  df$nb[!lc$scored[m] %in% TRUE] <- NA          # unscored: no hull, mid grey fill
  hulls <- .nb_hulls(df$x, df$y, df$nb, df$fill)

  q <- rect@coords[1L, ]
  ## A black underlay slightly thicker than the white line drawn over it, so every
  ## cell edge reads as a dark band. White alone was near-invisible against the light
  ## grey fills that most cells carry.
  p <- ggplot(df, aes(x = x, y = y, group = name)) +
    geom_polygon(aes(fill = fill), colour = "black", linewidth = 0.35) +
    geom_polygon(fill = NA, colour = "white", linewidth = 0.15) +
    scale_fill_identity() +
    .hull_layers(hulls, linewidth = 0.6)
  if (length(focal_names)) {
    fo <- df[df$name %in% focal_names, , drop = FALSE]
    if (nrow(fo))
      p <- p + geom_polygon(data = fo, fill = NA, colour = "black", linewidth = 0.5)
  }
  list(p = p + .coord_px(rect) + scale_y_reverse() +
         labs(x = "pixel", y = "pixel") +
         theme_v + theme(legend.position = "none"),
       hulls = hulls)
}

## Panel 3 without the image behind it -- the transcript colours are hard to read
## over a busy multichannel background. Faint cell outlines keep the spatial
## reference instead.
.gene_plot <- function(gp, poly, rect, size, hulls = NULL) {
  q <- rect@coords[1L, ]
  ## Explicit z-order, bottom to top:
  ##   1  filled cells   -- scored-neighbourhood cells, light grey, no outline
  ##   2  transcripts    -- on top of the fill, so nothing masks them
  ##   3  hollow cells   -- every cell's outline, over the dots, to keep edges crisp
  ##   4  neighbourhood hulls
  ## The fill can be fully opaque now that it sits UNDER the transcripts rather than
  ## over them, which is why GENE_CELL_ALPHA is gone.
  ## NOTE the grey means "in a neighbourhood" HERE, whereas the same light grey in
  ## panel 2 means "not light-chain-high" -- the panels answer different questions.
  df <- NULL
  p  <- ggplot()
  if (!is.null(poly)) {
    df <- poly@coords %>% group_by(name) %>% reframe(x = unlist(x), y = unlist(y))
    m  <- match(cosmx_poly_key(df$name), lc$key)
    df$in_nb <- lc$scored[m] %in% TRUE
    if (any(df$in_nb))
      p <- p + geom_polygon(data = df[df$in_nb, ], aes(x = x, y = y, group = name),
                            fill = GENE_CELL_FILL, colour = NA)
  }
  p <- cosmx_add_genes(p, gp, size = size)
  if (!is.null(df))
    p <- p + geom_polygon(data = df, aes(x = x, y = y, group = name),
                          fill = NA, colour = GENE_CELL_LINE, linewidth = 0.15)
  p <- p + .hull_layers(hulls, linewidth = 0.6)
  p + .coord_px(rect) + scale_y_reverse() +
    labs(x = "pixel", y = "pixel") +
    theme_v
}

## Cells of this sample whose centroids fall in a window. `flat` already carries
## exact canvas px for every cell, so these fallbacks cost no extra I/O.
.win_cells <- function(rect) {
  q <- rect@coords[1L, ]
  flat[X_px >= q$xmin & X_px <= q$xmax & Y_px >= q$ymin & Y_px <= q$ymax,
       .(key, X_px, Y_px)]
}

## THE one definition of how a canvas-px panel is oriented. Panels 1-3 are in canvas
## px, where X_px runs OPPOSITE to x_slide_mm; panel 0 (the locator) is in mm. Without
## this the locator reads mirrored against the panels beside it.
##
## It has to be scale_x_reverse(), NOT reversed coord limits: ggplot2 SORTS the limits
## passed to coord_fixed/coord_cartesian, so coord_fixed(xlim = c(max, min)) silently
## does nothing (verified on ggplot2 4.0.3 -- panel x.range came back ascending). Only
## the scale transform actually flips the axis. It composes fine with geom_raster, which
## the image panel uses.
.coord_px <- function(rect) {
  q <- rect@coords[1L, ]
  list(scale_x_reverse(),
       coord_fixed(xlim = c(q$xmin, q$xmax), ylim = c(q$ymin, q$ymax)))
}

.coord_win <- function(rect) {
  list(.coord_px(rect), scale_y_reverse(), labs(x = "pixel", y = "pixel"), theme_v)
}

## Panel 2 for windows too wide for polygons: one point per cell, neighbourhood
## colours, uninvolved cells light grey underneath.
.nb_point_plot <- function(rect, focal_nb, min_dist_px, size) {
  w <- merge(.win_cells(rect), lc[, .(key, nb, scored)], by = "key", all.x = TRUE)
  inv <- w[scored %in% TRUE]                     # only scored nb get colours/hulls
  uns <- w[!is.na(nb) & !scored %in% TRUE]       # called but unscored
  nb_cols <- .nb_adjacent_colors(inv$nb, inv$X_px, inv$Y_px, min_dist_px)
  inv[, col := unname(nb_cols[as.character(nb)])]
  ## Hulls on centroids here -- no boundary polygons are loaded at this width.
  hulls <- .nb_hulls(inv$X_px, inv$Y_px, inv$nb, inv$col)

  p <- ggplot() +
    geom_point(data = w[is.na(nb)], aes(X_px, Y_px), colour = UNINV, size = size * 0.6) +
    geom_point(data = uns, aes(X_px, Y_px), colour = UNSCORED, size = size * 0.8) +
    geom_point(data = inv, aes(X_px, Y_px, colour = col), size = size) +
    scale_colour_identity() +
    .hull_layers(hulls, linewidth = 0.4)
  if (!is.na(focal_nb) && nrow(inv[nb == focal_nb]))
    p <- p + geom_point(data = inv[nb == focal_nb], aes(X_px, Y_px),
                        shape = 21, colour = "black", fill = NA,
                        size = size * 2.2, stroke = 0.4)
  list(p = p + .coord_win(rect) +
         labs(title = sprintf("2. %d cells by neighbourhood (%d nb; focal circled) -- centroids",
                              nrow(w), length(nb_cols))),
       hulls = hulls)
}

## Panel 3 for windows too wide to read transcripts: the same kappa/lambda
## question at cell level instead of transcript level.
.call_point_plot <- function(rect, size, hulls = NULL) {
  q <- rect@coords[1L, ]
  w <- lc[X_px >= q$xmin & X_px <= q$xmax & Y_px >= q$ymin & Y_px <= q$ymax]
  p <- ggplot(w, aes(X_px, Y_px, colour = call)) +
    geom_point(size = size) +
    scale_colour_manual(values = CALL_COLS, name = "light chain") +
    .hull_layers(hulls, linewidth = 0.4)
  p + .coord_win(rect) +
    labs(title = sprintf("3. %s light-chain-high cells by call (transcripts too wide to read)",
                         format(nrow(w), big.mark = ",")))
}

## The three views for one window. Panel titles stay short -- the full
## neighbourhood label goes on the figure once, via plot_annotation().
.views <- function(rect, inset = NULL, gp = NULL, focal_nb = NA_integer_, detail = TRUE) {
  img   <- cosmx_image(SAMPLE_ID, rect) %>% add_um_scale_bar()
  ## replaces the coord TiffPlotR set, to flip x; y stays as TiffPlotR had it
  p_img <- img@plots$colorized + .coord_px(rect) + theme_v + labs(title = "1. image")
  adj_px <- ADJ_MULT * (RADIUS / 1000) / COSMX_MM_PER_PX

  if (!detail) {
    ## Too wide for polygons or transcripts -- see MAX_DETAIL_FOVS.
    cp      <- .nb_point_plot(rect, focal_nb, adj_px, size = 0.35)
    p_cells <- cp$p
    ## the hulls are the most legible thing on the panel, so carry them onto the
    ## image and the call panel too
    p_img   <- p_img + .hull_layers(cp$hulls, linewidth = 0.4)
    p_genes <- .call_point_plot(rect, size = 0.25, hulls = cp$hulls)
    if (!is.null(inset)) {
      p_img   <- p_img   %>% shape_annotate(inset, color = "white", fill = NA, linewidth = 0.7)
      p_cells <- p_cells %>% shape_annotate(inset, color = "black", fill = NA, linewidth = 0.7)
      p_genes <- p_genes %>% shape_annotate(inset, color = "black", fill = NA, linewidth = 0.7)
    }
    return(list(loc = .locator_plot(rect), img = p_img,
                cells = p_cells, genes = p_genes))
  }

  poly <- cosmx_cells_in_rect(SAMPLE_ID, rect)
  if (is.null(poly)) {
    p_cells <- p_img + labs(title = "2. no cells in window")
    p_genes <- .gene_plot(.gp_subset(gp, rect), NULL, rect, size = 1.6)
  } else {
    ## Restrict the palette to the neighbourhoods actually matched in this
    ## window: cosmx_color_cells() stops if a matched value has no palette entry.
    hit <- cosmx_match_ids(lc$cell_id, poly@coords$name)
    win <- lc[match(hit, lc$cell_id)][scored %in% TRUE]
    nb_cols <- .nb_adjacent_colors(win$nb, win$X_px, win$Y_px, min_dist_px = adj_px)

    ## Cells with no neighbourhood (not light-chain-high, or dropped by the QC
    ## that step 04 inherited) are unmatched, so na_color paints them light grey.
    ## The "not found in metadata" message below is expected, not a join failure.
    poly_nb <- cosmx_color_cells(poly, as.data.frame(lc), "nb_chr",
                                 id_col = "cell_id",
                                 palette = c(nb_cols, unscored = UNSCORED,
                                             plasma = PLASMA_COL, halo = HALO_COL),
                                 na_color = UNINV, border = TRUE)
    focal_names <- if (is.na(focal_nb)) character(0) else
      poly_nb@coords$name[cosmx_poly_key(poly_nb@coords$name) %in%
                          ## %in%, not ==: nb is NA for the excluded plasma/halo cells,
                          ## and a logical index containing NA injects NA elements into a
                          ## base-R subset, which then breaks the id parsers
                          cosmx_seu_key(lc$cell_id[lc$nb %in% focal_nb])]
    cp <- .nb_fill_plot(poly_nb, rect, focal_names)
    p_cells <- cp$p +
      labs(title = sprintf("2. %d cells by neighbourhood (%d nb; focal nb outlined)",
                           nrow(poly@coords), length(unique(win$nb))))

    ## Same hulls over the image and the transcripts, so all three panels share one
    ## frame of reference. Transcript dots are half their former size and use the
    ## paler GENE_COLS, so the hull outlines stay the dominant edge.
    p_img   <- p_img + .hull_layers(cp$hulls, linewidth = 0.7)
    p_genes <- .gene_plot(.gp_subset(gp, rect), poly, rect,
                          size = if (is.null(inset)) 0.8 else 0.55,
                          hulls = cp$hulls)
  }
  gp_i <- .gp_subset(gp, rect)
  p_genes <- p_genes +
    labs(title = sprintf("3. %s IGKC / IGLC1,2 transcripts",
                         if (is.null(gp_i)) "no" else format(nrow(gp_i), big.mark = ",")))

  ## on the zoomed-out view, outline where the tight window sat
  if (!is.null(inset)) {
    p_img   <- p_img   %>% shape_annotate(inset, color = "white", fill = NA, linewidth = 0.7)
    p_cells <- p_cells %>% shape_annotate(inset, color = "black", fill = NA, linewidth = 0.7)
    p_genes <- p_genes %>% shape_annotate(inset, color = "black", fill = NA, linewidth = 0.7)
  }
  list(loc = .locator_plot(rect), img = p_img, cells = p_cells, genes = p_genes)
}

## ------------------ 3. render ----
for (r in seq_len(nrow(sel))) {
  s  <- sel[r]
  mem <- lc[nb %in% s$nb & !is.na(X_px)]
  if (!nrow(mem)) { warning("no located members for nb ", s$nb); next }

  cx <- mean(range(mem$X_px)); cy <- mean(range(mem$Y_px))
  dx <- diff(range(mem$X_px)); dy <- diff(range(mem$Y_px))
  side_tight <- max(dx, dy) * PAD
  side_tight <- max(side_tight, MIN_SIDE)

  ## Cross-check the window against step 04, which reached these cells through
  ## the mm route instead. Its span_mm is the max RADIAL distance from the
  ## centroid (lightchain_neighborhoods.R:119), so reproduce that, not the bbox.
  span_px_mm <- max(sqrt((mem$X_px - mean(mem$X_px))^2 +
                         (mem$Y_px - mean(mem$Y_px))^2)) * COSMX_MM_PER_PX
  cat(sprintf("\n--- %s nb %d (%s): %d members, span %.3f mm (table %.3f mm)\n",
              s$group, s$nb, s$class, nrow(mem), span_px_mm, s$span_mm))
  if (abs(span_px_mm - s$span_mm) > 0.005)
    warning(sprintf("nb %d span disagrees with the table: %.3f vs %.3f mm",
                    s$nb, span_px_mm, s$span_mm))

  lab_base <- sprintf("r%dum %s nb%d  k=%.2f (%dk/%dl/%dt)  q=%.1e  %s%s",
                      RADIUS, s$group, s$nb, s$kappa_frac,
                      s$n_kappa, s$n_lambda, s$n_tie, s$q, s$class,
                      if ("n_plasma_near" %in% names(s))
                        sprintf("  | %d plasma cells nearby", s$n_plasma_near) else "")

  sides <- side_tight * ZOOMS
  rects <- lapply(names(ZOOMS), function(z)
    cosmx_px_rect(max(0, cx - sides[z] / 2), min(CANVAS[["x"]], cx + sides[z] / 2),
                  max(0, cy - sides[z] / 2), min(CANVAS[["y"]], cy + sides[z] / 2),
                  name = sprintf("nb%d_%s", s$nb, z)))
  names(rects) <- names(ZOOMS)

  ## How much detail each zoom can afford, decided from the FOV count alone -- that
  ## is known before any csv is opened.
  n_fov <- vapply(rects, function(r) {
    fv <- cosmx_fovs_in_rect(SAMPLE_ID, r); if (is.null(fv)) 0L else nrow(fv)
  }, 1L)
  detail <- n_fov <= MAX_DETAIL_FOVS
  cat(sprintf("    zooms: %s\n", paste(sprintf("%s=%dpx/%dFOV/%s", names(ZOOMS),
              round(sides), n_fov, ifelse(detail, "detail", "points")), collapse = "  ")))

  ## Read transcripts ONCE, on the widest window that can afford them, and subset
  ## for the narrower ones -- each intersecting FOV's transcript csv is ~120 MB.
  gp <- NULL; gp_rect <- NULL
  wid <- names(ZOOMS)[detail]
  if (length(wid)) {
    gp_rect <- rects[[wid[length(wid)]]]
    gp <- cosmx_gene_points(SAMPLE_ID, gp_rect, names(GENE_COLS), palette = GENE_COLS)
  }

  for (z in names(ZOOMS)) {
    tight <- z == names(ZOOMS)[1]
    rect  <- rects[[z]]
    q     <- rect@coords[1L, ]
    side  <- max(q$xmax - q$xmin, q$ymax - q$ymin)   # after clamping to the canvas
    v <- .views(rect, inset = if (tight) NULL else rects[[1]], gp = gp,
                focal_nb = s$nb, detail = detail[[z]])

    pg <- v$loc + v$img + v$cells + v$genes +
      plot_layout(nrow = 1, widths = c(0.55, 1, 1, 1)) +
      plot_annotation(
        title = lab_base,
        subtitle = sprintf(
          "%s -- %.0f px (%.3f mm) window, %d FOVs, %s%s | neighbourhood colours arbitrary but adjacent differ; IGKC %s, IGLC1/2 %s",
          z, side, side * COSMX_MM_PER_PX, n_fov[[z]],
          if (detail[[z]]) "boundaries + transcripts"
          else sprintf("centroids (over %d FOVs)", MAX_DETAIL_FOVS),
          if (tight) "" else "; outline marks the tight window",
          GENE_HI, GENE_LO))
    out_f <- file.path(OUT_DIR, sprintf("%s_nb%d_%s.png", s$group, s$nb, z))
    ggsave(out_f, plot = pg, width = 22, height = 7.5, dpi = 140, limitsize = FALSE)
    cat("wrote ", out_f, "\n", sep = "")
  }
}

writeLines(format(Sys.time()), file.path(OUT_DIR, "views.done"))
cat("\ndone r", RADIUS, "um\n", sep = "")
