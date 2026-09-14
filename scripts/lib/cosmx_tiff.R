## ============================================================================
## cosmx_tiff.R -- map CosMx/Seurat coordinates into assembled OME-TIFF pixel
## space, and plot image + cell boundaries + transcript points for a window.
##
## Load with (from the project root):
##   source("scripts/lib/cosmx_tiff.R")
##
## ---------------------------------------------------------------------------
## HISTORY
## Consolidated from three scripts in the v1 project (P3_CosMx), which are left
## untouched there as the provenance record for the v1 figures:
##   assemble_cosmx/16_cosmx_tiff_mapping.R  -> constants, ids, placement, coords
##   20_cell_boundaries.R                    -> resources, boundaries, colouring
##   21_gene_points.R                        -> transcript points
## In v1 the FOV-placement computation existed in three copies, the FOV-overlap
## test in two, and the metadata key-match in two; all are single here. Names
## gained a cosmx_ prefix, sample_id became the first argument everywhere, and
## cells_in_rect()'s three-way polymorphic return split into three functions.
##
## ---------------------------------------------------------------------------
## COORDINATE SPACES -- three of them, easy to confuse:
##   mm         slide millimetres. Seurat img$centroids@coords; fov_positions.csv
##              X_mm/Y_mm; the x_slide_mm/y_slide_mm columns in output CSVs.
##   FOV-local  pixels within one FOV. Seurat meta.data x_FOV_px/y_FOV_px;
##              boundary CSV x_local/y_local; transcript CSV x/y; flat-file
##              CenterX_local_px/CenterY_local_px.
##   canvas     pixels in the assembled OME-TIFF -- named X_px/Y_px throughout.
##              This is what TiffRect and all plotting here use.
##
## The canvas transform (matches assemble_cosmx/13_build_ometiff.py exactly):
##   x0   = round((X_mm[fov] - min(X_mm)) / COSMX_MM_PER_PX)   # FOV placement
##   y0   = round((Y_mm[fov] - min(Y_mm)) / COSMX_MM_PER_PX)
##   X_px = x0 + (COSMX_FOV_PX - 1) - x_local     # tiles were mirrored in x
##   Y_px = y0 +                      y_local     # y is not mirrored
##
## Two ways to reach canvas px:
##   1. EXACT, per cell, from FOV-local px  -> cosmx_local_to_px/cosmx_cells_to_px
##   2. FITTED global mm->px affine          -> cosmx_fit_mm_affine/cosmx_mm_to_px
## Use (1) whenever FOV-local px are available. (2) exists to place ARBITRARY mm
## coordinates (e.g. a region bounding box recorded in mm) and reports residuals
## so the fit can be trusted.
##
## ---------------------------------------------------------------------------
## DO NOT "SIMPLIFY" THESE -- each was measured, and each has a silent failure
## mode. Full rationale in the plan; the short version:
##   * COSMX_MM_PER_PX / COSMX_FOV_PX are instrument constants.
##   * round(), and min(X_mm) taken AFTER the per-sample filter.
##   * (COSMX_FOV_PX - 1L) in the transform, but x0 + COSMX_FOV_PX for the FOV
##     bounding rects -- NOT an off-by-one, see cosmx_fovs_in_rect().
##   * cell ids are parsed from the END of the underscore split.
##   * polygon names are UNPADDED ("F448_C12"); padding breaks every join.
##   * cosmx_channel_colors() must return a list; channel name order is
##     positional and reordering silently recolours every image.
##   * cosmx_fit_mm_affine() fits two separable 1-D lm's PER SAMPLE with a free
##     intercept. RB1's mm frame is offset from IB1's by ~1 mm on both axes, so
##     a shared or analytic intercept displaces every RB1 figure by ~2 FOVs
##     while still reporting sub-pixel residuals.
##   * cosmx_cells_to_px() and cosmx_fit_mm_affine() filter rows DIFFERENTLY on
##     purpose; unifying them changes which cells enter the lm.
##   * cosmx_mm_to_px() reads its input columns positionally.
##   * .cosmx_gquants_to_precalc()'s min_q branch -- see the note there.
##
## Validated by v1's assemble_cosmx/12_calib_orientation.py,
## 14_calib_cellmap.py and 17_validate_mapping.R, and now additionally by:
##   scripts/sbatch/verify_cosmx_tiff.R        -- V3 asserts
##     X_px + CenterX_global_px == max(x0) + COSMX_FOV_PX - 1 exactly for every
##     one of the 4.4M cells in the flat-file metadata of both samples.
##   scripts/sbatch/verify_cosmx_tiff_seurat.R -- the Seurat-facing half; every
##     cell of IB1 (400k) and RB1 (1.4M) maps identically to the flat-file route,
##     and both match v1's functions exactly.
##
## Measured affines (from that second script), which is why the intercept is
## fitted per sample rather than derived from the canvas:
##   IB1  ax=-8313.869 bx=124851.0  ay=-8313.879 by=166305.0  resid=(0.00, 0.02)
##   RB1  ax=-8313.909 bx=124851.6  ay=-8313.876 by=166305.0  resid=(0.12, 0.12)
## IB1's intercepts equal its own canvas extent (124851, 166305) -- its mm origin
## sits on the canvas edge. RB1's do NOT equal its extent (116537, 157517): its mm
## frame is offset by ~8314 px in x (1.0002 mm) and ~8788 px in y (1.057 mm). So an
## intercept derived analytically from each sample's canvas would put every RB1
## figure about two FOVs off, with the residuals still reporting sub-pixel.
## ============================================================================

suppressPackageStartupMessages({
  library(methods)
  library(TiffPlotR)
  library(dplyr)
  library(tidyr)
})
## NOTE: v1 guarded these with `if (!"pkg" %in% loadedNamespaces()) library(pkg)`,
## which is not an attach guard -- a loaded-but-not-attached namespace passes it
## while leaving %>%, filter() and pivot_longer() unavailable. v1 also attached
## the whole tidyverse; consuming scripts must now attach ggplot2 / readr /
## stringr themselves if they need them.

## ---------------------------------------------------------------------------
## CONFIG
## ---------------------------------------------------------------------------
COSMX_MM_PER_PX    <- 120.280945e-6   # mm per full-res pixel (ImPixelSize_nm header)
COSMX_FOV_PX       <- 4256L           # FOV side, in px
COSMX_POS_CSV      <- "input/fov_positions.csv"
COSMX_TIFF_DIR     <- "/gpfs2/scratch/pi/avolaric/P3_CosMx_assemble_FOVs"
COSMX_GQUANT_CACHE <- "input/tiff_gquants_v2.Rds"
## NOTE: scripts/02_spatial_QC.R independently defines FOV_PX <- 4256L for its
## FOV-local edge distances. Keep the two literals in sync; verify V1 checks it.

## Memo for the two things that would otherwise be re-parsed on every call in a
## loop over windows: the ~200 KB placement CSV and the 227 KB quantile cache.
## Keyed on (path, mtime), so editing either file invalidates it. Invisible to
## callers; clear with cosmx_forget().
.cosmx_memo <- new.env(parent = emptyenv())

cosmx_forget <- function() {
  rm(list = ls(.cosmx_memo, all.names = TRUE), envir = .cosmx_memo)
  invisible(NULL)
}

.cosmx_memo_key <- function(prefix, path, ...) {
  paste(prefix, path, file.mtime(path), ..., sep = "|")
}

## ---------------------------------------------------------------------------
## RESOURCES
## ---------------------------------------------------------------------------

## Assembled OME-TIFFs, as a named list (names are sample ids: IB1, RB1).
cosmx_tiff_paths <- function(dir = COSMX_TIFF_DIR) {
  tiffs <- dir(dir, pattern = "tiff$", full.names = TRUE)
  names(tiffs) <- tiffs %>% basename %>% sub(".ome.tiff$", "", .)
  as.list(tiffs)
}

## Channel identities, in acquisition order. This order is POSITIONAL -- it maps
## onto channel indices 1..5 in the TIFF, so reordering silently recolours every
## image.
##   Channel-B PanCK | Channel-G G | Channel-Y Membrane | Channel-R CD45 | Channel-U DNA
cosmx_channel_names <- function() {
  c("PanCK", "G", "Membrane", "CD45", "DNA")
}

## Display colour per channel, same order as cosmx_channel_names().
## MUST stay a list: fetchTiffData.colorized() stops on a non-list.
cosmx_channel_colors <- function() {
  c("blue", "green", "yellow", "red", "violet") %>% as.list
}

## Reduce collected channel quantiles to the per-resolution min/max normalisation
## table that fetchTiffData() expects.
##
## Kept deliberately separate from the cache handling so it can be tested on the
## login node without touching a TIFF.
##
## The two subtleties, both load-bearing:
##  * filter(signal_value > 0) runs BEFORE the grouping. On this data the 30th
##    percentile is 0 in all 100 (sample x resolution x channel) groups, so it is
##    filtered out and `min_q %in% quantile_value` is always FALSE -- meaning
##    min_value is "the smallest non-zero quantile", not "the min_q quantile".
##    Collapsing the branch to signal_value[quantile_value == min_q] sets every
##    min_value to 0 and blows out the low end of every channel in every figure.
##  * `resolution` stays a factor and `channel` becomes numeric, because that is
##    what TiffPlotR's .handle_precalc_max() consumes: it selects rows by fetch
##    resolution, and stops on character channel names that aren't in
##    channel_names.
.cosmx_gquants_to_precalc <- function(tiff_gquants, min_q = 0.30, max_q = 0.99) {
  convert_q_matrix <- function(x) {
    pvals <- rownames(x) %>% sub("%", "", .) %>% as.numeric
    x <- as.data.frame(x)
    rownames(x) <- NULL
    x$quantile_value <- pvals / 100
    pivot_longer(x, -quantile_value, names_to = "channel", values_to = "signal_value")
  }

  gq_df <- lapply(tiff_gquants, function(samp_gquant) {
    sig_q <- samp_gquant$signal_quantiles
    names(sig_q) <- seq_along(sig_q)          # list index IS the resolution level
    lapply(sig_q, convert_q_matrix) %>% bind_rows(.id = "resolution")
  }) %>% bind_rows(.id = "sample_id")
  gq_df$resolution <- gq_df$resolution %>% as.numeric() %>% factor

  precalc <- gq_df %>%
    filter(signal_value > 0) %>%
    group_by(sample_id, resolution, channel) %>%
    summarise(
      max_value = signal_value[quantile_value == max_q],
      min_value = {
        if (min_q %in% quantile_value) signal_value[quantile_value == min_q] else min(signal_value)
      },
      .groups = "drop_last"
    )
  precalc$channel <- as.numeric(precalc$channel)

  split(precalc, precalc$sample_id)
}

## Per-sample precalc_max tables, as a named list of
## data.frame(sample_id, resolution, channel, max_value, min_value).
##
## Reads the collected quantiles from `cache` (227 KB). The 87 GB / 48 GB TIFFs
## do not need to be touched. refresh = TRUE recollects them, which is a batch
## job -- see the commented block below.
##
## min_q/max_q define the contrast of every figure. Do not change them casually:
## TiffPlotR's own quantiles_to_precalc_max() defaults (0.90 / 0.998) move IB1's
## DNA min_value from 140 to 2663, clipping essentially all nuclear signal to
## black. It also aggregates away the `resolution` column, which is why the
## reduction above is ours rather than a call into TiffPlotR.
cosmx_precalc_max <- function(cache = COSMX_GQUANT_CACHE,
                              min_q = 0.30, max_q = 0.99,
                              refresh = FALSE,
                              tiffs = NULL) {
  if (!refresh && !file.exists(cache))
    stop("no quantile cache at '", cache, "'. Stage it from the v1 project, or ",
         "pass refresh = TRUE to recollect (a batch job over the OME-TIFFs).",
         call. = FALSE)

  if (refresh) {
    if (is.null(tiffs)) tiffs <- cosmx_tiff_paths()
    ## Recollection is the expensive, well-tested half of TiffPlotR, so use it --
    ## but feed its output to OUR reduction above, to keep the contrast fixed.
    ## Not exercised since the v1 assembly; expect --mem=150G --time=12:00:00.
    tiff_gquants <- lapply(tiffs, TiffPlotR::gather_channel_quantiles)
    saveRDS(tiff_gquants, cache)
    cosmx_forget()
  }

  key <- .cosmx_memo_key("precalc", cache, min_q, max_q)
  if (!is.null(.cosmx_memo[[key]])) return(.cosmx_memo[[key]])
  out <- .cosmx_gquants_to_precalc(readRDS(cache), min_q = min_q, max_q = max_q)
  .cosmx_memo[[key]] <- out
  out
}

## ---------------------------------------------------------------------------
## CELL IDS AND JOIN KEYS
## ---------------------------------------------------------------------------

## The canonical CosMx cell id is 'c_<slide>_<fov>_<cell>', but Seurat cells are
## often sample-prefixed ('IB1_c_<slide>_<fov>_<cell>') after a merge with
## add.cell.ids. Read from the END (fov = 2nd-to-last token, cell = last) so any
## leading prefix is tolerated.
.cosmx_id_tokens <- function(cell_id) strsplit(as.character(cell_id), "_", fixed = TRUE)

cosmx_fov_from_id <- function(cell_id)
  as.integer(vapply(.cosmx_id_tokens(cell_id), function(z) z[length(z) - 1L], character(1)))

cosmx_cell_from_id <- function(cell_id)
  as.integer(vapply(.cosmx_id_tokens(cell_id), function(z) z[length(z)], character(1)))

## Join key "fov_cell" from a CosMx/Seurat cell id (prefix-robust).
##   "c_2_368_802" and "IB1_c_2_368_802" both -> "368_802"
cosmx_seu_key <- function(cell_id)
  paste(cosmx_fov_from_id(cell_id), cosmx_cell_from_id(cell_id), sep = "_")

## Join key "fov_cell" from a TiffPolygon name "F<fov>_C<cell>".
##   "F368_C802" -> "368_802".  Unmatched names give "NA_NA", not an error.
## as.integer() makes this robust to zero-padded names even though
## cosmx_fov_boundaries() deliberately does not pad -- see its note.
cosmx_poly_key <- function(poly_names) {
  pm <- regmatches(poly_names, regexec("F([0-9]+)_C([0-9]+)", poly_names))
  paste(as.integer(vapply(pm, function(z) z[2], character(1))),
        as.integer(vapply(pm, function(z) z[3], character(1))), sep = "_")
}

## Index of each polygon into a vector of cell ids. NA where absent.
.cosmx_key_hits <- function(cell_ids, poly_names)
  match(cosmx_poly_key(poly_names), cosmx_seu_key(cell_ids))

## Cell ids aligned to polygon names -- i.e. order-safe by poly_names.
cosmx_match_ids <- function(cell_ids, poly_names)
  cell_ids[.cosmx_key_hits(cell_ids, poly_names)]

## ---------------------------------------------------------------------------
## FOV PLACEMENT AND COORDINATE TRANSFORMS
## ---------------------------------------------------------------------------

## Where each FOV of one sample sits on the assembled canvas.
##
## Returns data.frame(fov, x0, y0, file): the canvas-px origin of the FOV, plus
## the morphology TIF path it came from. Every per-FOV file (cell boundaries,
## transcript coordinates) is derived from `file`, which is why it is returned --
## in v1 it was dropped here and the whole computation was reimplemented twice.
cosmx_fov_placement <- function(sample_id, pos_csv = COSMX_POS_CSV) {
  key <- .cosmx_memo_key("placement", pos_csv, sample_id)
  if (!is.null(.cosmx_memo[[key]])) return(.cosmx_memo[[key]])

  d <- read.csv(pos_csv, stringsAsFactors = FALSE)
  d <- d[d$sample == sample_id, , drop = FALSE]
  if (!nrow(d)) stop("no FOVs for sample ", sample_id, " in ", pos_csv)
  d$fov <- as.integer(sub(".*_F0*([0-9]+)\\.TIF$", "\\1", d$file))
  xmin <- min(d$X_mm); ymin <- min(d$Y_mm)
  out <- data.frame(fov  = d$fov,
                    x0   = round((d$X_mm - xmin) / COSMX_MM_PER_PX),
                    y0   = round((d$Y_mm - ymin) / COSMX_MM_PER_PX),
                    file = d$file,
                    stringsAsFactors = FALSE)
  ## A duplicate fov would make the two consumers disagree: match() takes the
  ## first, the %in% subset in cosmx_fovs_in_rect() takes both.
  stopifnot(!anyDuplicated(out$fov))

  .cosmx_memo[[key]] <- out
  out
}

## EXACT mapping: FOV-local px -> canvas px. Returns data.frame(X_px, Y_px).
##   fov, x_fov_px, y_fov_px: numeric vectors of equal length
cosmx_local_to_px <- function(sample_id, fov, x_fov_px, y_fov_px,
                              pos_csv = COSMX_POS_CSV) {
  ## This is the single chokepoint for the validated transform after the v1
  ## dedup, so refuse to run it with clobbered constants.
  stopifnot(identical(COSMX_FOV_PX, 4256L),
            identical(COSMX_MM_PER_PX, 120.280945e-6))
  pl <- cosmx_fov_placement(sample_id, pos_csv)
  i  <- match(fov, pl$fov)
  if (anyNA(i)) warning(sum(is.na(i)), " cells reference a FOV not in the image set")
  data.frame(X_px = pl$x0[i] + (COSMX_FOV_PX - 1L) - x_fov_px,
             Y_px = pl$y0[i] +                       y_fov_px)
}

## Which metadata column names a cell's sample. Returns NA_character_ to mean
## "do not filter" (single-sample object), and stops rather than silently
## mapping every cell with one sample's placement -- v1 skipped the filter
## whenever Run_Tissue_name was absent, which is the case for merged objects
## here (they carry sample_id / merge_sample instead).
.cosmx_sample_col <- function(sample_col, cols) {
  if (length(sample_col) == 1L && is.na(sample_col)) return(NA_character_)
  if (!is.null(sample_col)) {
    if (!sample_col %in% cols)
      stop("sample_col '", sample_col, "' not in metadata", call. = FALSE)
    return(sample_col)
  }
  cand <- c("Run_Tissue_name", "sample_id", "merge_sample")
  found <- cand[cand %in% cols]
  if (!length(found))
    stop("no sample column found (looked for ", paste(cand, collapse = ", "),
         "). Pass sample_col = \"<column>\", or sample_col = NA to map all cells ",
         "of a single-sample object.", call. = FALSE)
  found[1L]
}

## All cells of one sample in canvas px, from the FOV-local columns of a Seurat
## object's meta.data. Returns data.frame(cell_id, fov, X_px, Y_px).
cosmx_cells_to_px <- function(sample_id, seu,
                              xcol = "x_FOV_px", ycol = "y_FOV_px",
                              sample_col = NULL, pos_csv = COSMX_POS_CSV) {
  md   <- seu@meta.data
  cid  <- rownames(md)
  fov  <- cosmx_fov_from_id(cid)
  keep <- fov > 0 & !is.na(md[[xcol]]) & !is.na(md[[ycol]])

  sample_col <- .cosmx_sample_col(sample_col, colnames(md))
  if (!is.na(sample_col)) keep <- keep & md[[sample_col]] == sample_id

  md <- md[keep, , drop = FALSE]; cid <- cid[keep]; fov <- fov[keep]
  px <- cosmx_local_to_px(sample_id, fov, md[[xcol]], md[[ycol]], pos_csv)
  data.frame(cell_id = cid, fov = fov, X_px = px$X_px, Y_px = px$Y_px)
}

## Fit a global mm->px affine (X_px = ax*x_mm + bx ; Y_px = ay*y_mm + by) from
## this sample's cells, using the EXACT FOV-local px as ground truth. Returns
## list(sample_id, ax, bx, ay, by, resid_px = c(x=, y=), n).
##
## Two separable 1-D fits: no rotation or shear term, because the placement is
## axis-aligned. Both slopes are expected NEGATIVE (~ -8313.9 px/mm) and the
## residual SD sub-pixel; the guards below say so out loud, since a bad fit is
## otherwise only visible as a shifted image.
##
## The intercept is free and PER SAMPLE. That is not incidental -- RB1's mm frame
## is offset from IB1's by ~1.00 mm in x and ~1.06 mm in y, so any shared or
## analytically-derived intercept displaces every RB1 figure by ~2 FOVs while
## still reporting sub-pixel residuals.
##
## Note the row filter differs from cosmx_cells_to_px()'s on purpose: no fov > 0
## test, no sample filter (`image` already selects the sample's cells).
cosmx_fit_mm_affine <- function(sample_id, seu, image = sample_id,
                                xcol = "x_FOV_px", ycol = "y_FOV_px",
                                pos_csv = COSMX_POS_CSV) {
  if (!image %in% names(seu@images))
    stop("no image '", image, "' in object; have: ",
         paste(names(seu@images), collapse = ", "), call. = FALSE)

  co    <- seu@images[[image]]$centroids
  cells <- co@cells
  mm    <- co@coords                      # columns: x, y (mm)
  if (!all(cells %in% rownames(seu@meta.data)))
    stop("centroid names do not match meta.data rownames -- merged-object ",
         "prefix mismatch?", call. = FALSE)

  md  <- seu@meta.data[cells, , drop = FALSE]
  fov <- cosmx_fov_from_id(cells)
  ok  <- !is.na(md[[xcol]]) & !is.na(md[[ycol]]) & !is.na(fov)
  px  <- cosmx_local_to_px(sample_id, fov[ok], md[[xcol]][ok], md[[ycol]][ok], pos_csv)

  fx <- lm(px$X_px ~ mm[ok, 1]); fy <- lm(px$Y_px ~ mm[ok, 2])
  aff <- list(sample_id = sample_id,
              ax = unname(coef(fx)[2]), bx = unname(coef(fx)[1]),
              ay = unname(coef(fy)[2]), by = unname(coef(fy)[1]),
              resid_px = c(x = sd(resid(fx)), y = sd(resid(fy))),
              n = sum(ok))
  message(sprintf("[%s] n=%d  slope_x=%.2f slope_y=%.2f px/mm (|expect|=%.1f)",
                  sample_id, aff$n, aff$ax, aff$ay, 1 / COSMX_MM_PER_PX))
  message(sprintf("[%s] intercept_x=%.1f intercept_y=%.1f  resid SD: x=%.2f y=%.2f px",
                  sample_id, aff$bx, aff$by, aff$resid_px[1], aff$resid_px[2]))

  expect <- 1 / COSMX_MM_PER_PX
  if (aff$ax > 0 || aff$ay > 0)
    warning("expected both slopes negative, got ax=", signif(aff$ax, 6),
            " ay=", signif(aff$ay, 6))
  if (abs(abs(aff$ax) - expect) / expect > 0.01 ||
      abs(abs(aff$ay) - expect) / expect > 0.01)
    warning("slope is >1% off the expected ", signif(expect, 6), " px/mm")
  if (max(aff$resid_px) > 1)
    warning("residual SD exceeds 1 px (x=", signif(aff$resid_px[1], 3),
            " y=", signif(aff$resid_px[2], 3), ") -- do not trust this affine")
  aff
}

## Apply a fitted affine to arbitrary mm coordinates. `mm_xy` is a 2-column
## matrix/data.frame read POSITIONALLY (column 1 = x, column 2 = y) -- callers
## pass a bare centroids@coords matrix whose names are not usable.
cosmx_mm_to_px <- function(mm_xy, affine) {
  mm_xy <- as.matrix(mm_xy)
  stopifnot(ncol(mm_xy) == 2L)
  data.frame(X_px = affine$ax * mm_xy[, 1] + affine$bx,
             Y_px = affine$ay * mm_xy[, 2] + affine$by)
}

## A canvas-px TiffRect from any two opposite corners, bounds normalised.
##
## Use this instead of calling TiffRect() with mapped mm corners directly: both
## affine slopes are negative, so the mm minimum maps to the px MAXIMUM on both
## axes and a naive TiffRect(x_min_px, x_max_px, ...) is inverted. v1 handled
## this by swapping the arguments at the call site, which was correct but
## invisible.
cosmx_px_rect <- function(x1, x2, y1, y2, name = "rect") {
  TiffRect(xmin = pmin(x1, x2), xmax = pmax(x1, x2),
           ymin = pmin(y1, y2), ymax = pmax(y1, y2), name = name)
}

## ---------------------------------------------------------------------------
## WINDOWS: WHICH FOVS, WHICH CELLS
## ---------------------------------------------------------------------------

## FOVs (with placement) whose canvas box intersects a query rect.
## Returns data.frame(fov, x0, y0, file), or NULL if none -- silently, so
## callers can phrase their own message.
##
## On xmax = x0 + COSMX_FOV_PX rather than + COSMX_FOV_PX - 1: this is NOT an
## off-by-one against the (COSMX_FOV_PX - 1L) in the transform. Pixel indices run
## [x0, x0+4255]; as the continuous area they cover that is [x0, x0+4256), and
## rect_test_overlap() compares with strict `>`, i.e. treats rects as half-open.
## Narrowing it to +COSMX_FOV_PX-1 makes a rect starting inside the last pixel
## column of a FOV miss that FOV, silently dropping boundary cells and
## transcripts along every FOV seam. The current form errs toward inclusion: at
## worst one extra FOV is read, and its cells and transcripts are then removed by
## the filters downstream. Verify V5 pins this behaviour down.
cosmx_fovs_in_rect <- function(sample_id, rect, pos_csv = COSMX_POS_CSV) {
  stopifnot(is(rect, "TiffRect"), nrow(rect@coords) == 1L)
  pl <- cosmx_fov_placement(sample_id, pos_csv)
  fov_rects <- TiffRect(xmin = pl$x0, xmax = pl$x0 + COSMX_FOV_PX,
                        ymin = pl$y0, ymax = pl$y0 + COSMX_FOV_PX,
                        name = as.character(pl$fov))
  hits <- rect_test_overlap(fov_rects, rect, subset = TRUE)
  if (is.null(hits)) return(NULL)
  pl[pl$fov %in% as.integer(hits@coords$name), c("fov", "x0", "y0", "file"), drop = FALSE]
}

## CellBoundaries csv path for a FOV, derived from its morphology tif path:
##   .../CellStatsDir/Morphology2D/..._F00742.TIF
##   .../CellStatsDir/FOV00742/CellBoundaries_F00742.csv
.cosmx_boundary_path <- function(tif, fov) {
  csd <- dirname(dirname(tif))                       # CellStatsDir
  file.path(csd, sprintf("FOV%05d", fov), sprintf("CellBoundaries_F%05d.csv", fov))
}

## Cell boundaries for a set of FOVs, as one uncoloured TiffPolygon in canvas px.
## Returns NULL if nothing could be loaded. Cells with fewer than 3 vertices are
## dropped. Boundary vertices use the same FOV-local convention as centroids.
##
## Polygon names are "F<fov>_C<cell>", UNPADDED. Do not switch to %05d: the join
## keys are built by pasting the raw captured digits, so padded names would give
## "00448_00802" against cosmx_seu_key()'s "448_802" -- a total join failure that
## renders as "every cell is na_color" rather than as an error.
##
## Example
##   poly <- cosmx_fov_boundaries("IB1", c(448, 449))
##   range(unlist(poly$x)); range(unlist(poly$y))
cosmx_fov_boundaries <- function(sample_id, fovs, pos_csv = COSMX_POS_CSV) {
  pl <- cosmx_fov_placement(sample_id, pos_csv)

  polys <- list()
  for (f in fovs) {
    i <- match(f, pl$fov)
    if (is.na(i)) { warning("no placement for FOV ", f); next }
    cb <- .cosmx_boundary_path(pl$file[i], f)
    if (!file.exists(cb)) { warning("missing CellBoundaries: ", cb); next }
    b <- read.csv(cb)
    b$X <- pl$x0[i] + (COSMX_FOV_PX - 1L) - b$x_local
    b$Y <- pl$y0[i] +                       b$y_local
    b$nm <- sprintf("F%d_C%d", f, b$cellID)
    # drop degenerate cells (< 3 vertices)
    keep <- b$nm %in% names(which(table(b$nm) >= 3L))
    b <- b[keep, , drop = FALSE]
    if (!nrow(b)) next
    polys[[as.character(f)]] <- TiffPolygon(x = b$X, y = b$Y, name = b$nm)
  }
  if (!length(polys)) return(NULL)
  if (length(polys) == 1L) return(polys[[1]])
  do.call(c, unname(polys))   # unname so the TiffPolygon c() method binds `x`
}

## Cells whose boundaries overlap a canvas-px rect, as an UNCOLOURED TiffPolygon
## (or NULL). Colour it afterwards with either:
##   cosmx_color_cells(poly, meta, column)         # by a metadata value
##   TiffPlotR::assign_adjacent_colors(poly)       # just make neighbours distinct
##
## v1's cells_in_rect() did the fetching, the FOV listing and both colouring
## modes in one function, returning NULL, an integer vector or a TiffPolygon
## depending on flags -- and did not forward fill/border, so callers wanting a
## border had to bypass it anyway. Ask cosmx_fovs_in_rect() for the FOV list.
cosmx_cells_in_rect <- function(sample_id, rect, pos_csv = COSMX_POS_CSV) {
  stopifnot(is(rect, "TiffRect"), nrow(rect@coords) == 1L)
  fv <- cosmx_fovs_in_rect(sample_id, rect, pos_csv)
  if (is.null(fv)) { message("no FOVs intersect rect"); return(NULL) }
  fovs <- fv$fov
  message(length(fovs), " FOVs intersect rect: ", paste(sort(fovs), collapse = ", "))

  poly <- cosmx_fov_boundaries(sample_id, fovs, pos_csv)
  if (is.null(poly)) return(NULL)
  poly <- shape_overlaps_rect(poly, rect, subset = TRUE)
  if (is.null(poly)) { message("no cells overlap rect"); return(NULL) }
  message(nrow(poly@coords), " cells overlap rect")
  poly
}

## Colour cell polygons by a metadata column of a Seurat object or data.frame.
##
## Polygons are named "F<fov>_C<cell>"; metadata cells are matched on (fov, cell)
## parsed from the CosMx cell id "c_<slide>_<fov>_<cell>" (rownames by default,
## or a column named via id_col). The value->colour mapping is delegated to
## TiffPlotR::assign_value_colors(), which handles both categorical and
## continuous columns.
##
##   poly      TiffPolygon from cosmx_cells_in_rect() / cosmx_fov_boundaries()
##   meta      Seurat object or data.frame of metadata
##   column    metadata column to colour by
##   id_col    where to read cell ids from: rownames (default) or a column name
##   palette, continuous, na_color, fill, border  passed to assign_value_colors()
##
## Example
##   rect <- cosmx_px_rect(54250, 56750, 30250, 32750, name = "win")
##   poly <- cosmx_cells_in_rect("IB1", rect)
##   qc   <- read.csv("output/02_filter_QC/qc_spatial_cells.csv.gz")
##   poly <- cosmx_color_cells(poly, qc, "qc_spatial_staining", id_col = "cell")
##   attr(shape_meta(poly), "color_map")            # value -> colour legend
##   ggplot2::ggplot() %>% shape_annotate(poly)
##
##   # continuous column, drawn as an outline too
##   poly <- cosmx_color_cells(poly, seu, "Mean.G", border = TRUE)
##
##   # over the image
##   cosmx_image("IB1", rect) %>% shape_annotate(poly)
cosmx_color_cells <- function(poly, meta, column, id_col = NULL,
                              palette = NULL, continuous = NULL,
                              na_color = "grey80", fill = TRUE, border = FALSE) {
  md <- if (isS4(meta) && methods::is(meta, "Seurat")) meta@meta.data else as.data.frame(meta)
  if (!column %in% colnames(md)) stop("column '", column, "' not in metadata", call. = FALSE)

  ids  <- if (is.null(id_col)) rownames(md) else as.character(md[[id_col]])
  pkey <- poly@coords$name
  hit  <- .cosmx_key_hits(ids, pkey)
  vals <- md[[column]][hit]
  n_na <- sum(is.na(hit))
  if (n_na) message(n_na, "/", length(pkey), " cells not found in metadata (", na_color, ")")
  assign_value_colors(poly, vals, palette = palette, continuous = continuous,
                      na_color = na_color,
                      fill = fill, border = border)
}

## ---------------------------------------------------------------------------
## TRANSCRIPT POINTS
## ---------------------------------------------------------------------------

## Transcript coordinate csv for a FOV, resolved from its morphology tif path.
## Unlike .cosmx_boundary_path() this globs, so it can return 0 or several paths;
## callers must check length() and take [1].
.cosmx_transcript_path <- function(tif, fov) {
  base <- dirname(dirname(dirname(tif)))     # .../<run e.g. 20251218_234107_S1>
  Sys.glob(file.path(base, "AnalysisResults", "*", sprintf("FOV%05d", fov),
                     sprintf("*FOV%05d__complete_code_cell_target_call_coord.csv", fov)))
}

## Build a named gene->colour palette covering all present genes.
## Requested `genes` order is kept first; a named `palette` overrides those genes
## by name; an unnamed `palette` is applied positionally in that order; anything
## left over gets an automatic palette.
.cosmx_gene_palette <- function(present, genes, palette) {
  lev <- c(intersect(genes, present), setdiff(present, genes))
  n <- length(lev)
  auto <- if (n <= 8L)
    c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#a65628","#f781bf","#999999")[seq_len(max(n,1L))]
  else grDevices::hcl.colors(n, "Dark 3")
  pal <- stats::setNames(auto, lev)
  if (!is.null(palette)) {
    if (!is.null(names(palette))) {
      common <- intersect(names(palette), lev)
      pal[common] <- palette[common]
    } else {
      m <- min(length(palette), n)
      pal[seq_len(m)] <- palette[seq_len(m)]
    }
  }
  pal
}

## Transcript positions for genes of interest within a canvas-px rect.
##
##   sample_id        "IB1" / "RB1"
##   rect             single-row TiffRect in canvas pixel space
##   genes            target names, e.g. c("MS4A1","CD3D","IGHG1/2")
##   palette          named gene->colour, or unnamed (applied in gene order),
##                    or NULL for an automatic palette
##   endogenous_only  drop non-Endogenous codes (negative / system probes)
##
## Returns data.frame(X_px, Y_px, target, fov, CellId, CellComp, color) carrying
## a named `palette` attribute, or NULL. Ready for cosmx_add_genes().
##
## The rect filter here is an inclusive point-in-box test, deliberately: these
## are points, not shapes, and a transcript sitting exactly on the boundary
## should be drawn.
##
## Example
##   rect <- cosmx_px_rect(54250, 56750, 30250, 32750, name = "win")
##   gp <- cosmx_gene_points("IB1", rect, c("MS4A1", "CD3D", "IGHG1/2"),
##                           palette = c(MS4A1 = "cyan", CD3D = "yellow"))
##   attr(gp, "palette")
##   cosmx_add_genes(ggplot2::ggplot(), gp, size = 0.6)
cosmx_gene_points <- function(sample_id, rect, genes, palette = NULL,
                              endogenous_only = TRUE, pos_csv = COSMX_POS_CSV) {
  fv <- cosmx_fovs_in_rect(sample_id, rect, pos_csv)
  if (is.null(fv)) { message("no FOVs intersect rect"); return(NULL) }
  message(nrow(fv), " FOVs intersect rect: ", paste(sort(fv$fov), collapse = ", "))
  q <- rect@coords[1L, ]
  parts <- list()
  for (i in seq_len(nrow(fv))) {
    f <- fv$fov[i]
    p <- .cosmx_transcript_path(fv$file[i], f)
    if (!length(p) || !file.exists(p[1])) { warning("no transcript file for FOV ", f); next }
    tx <- data.table::fread(p[1], select = c("target", "x", "y", "CellId", "CellComp", "codeclass"))
    tx <- tx[tx$target %in% genes]
    if (endogenous_only) tx <- tx[tx$codeclass == "Endogenous"]
    if (!nrow(tx)) next
    X <- fv$x0[i] + (COSMX_FOV_PX - 1L) - tx$x
    Y <- fv$y0[i] +                       tx$y
    keep <- X >= q$xmin & X <= q$xmax & Y >= q$ymin & Y <= q$ymax
    if (!any(keep)) next
    parts[[as.character(f)]] <- data.frame(
      X_px = X[keep], Y_px = Y[keep], target = tx$target[keep], fov = f,
      CellId = tx$CellId[keep], CellComp = tx$CellComp[keep], stringsAsFactors = FALSE)
  }
  if (!length(parts)) { message("no transcripts for those genes in rect"); return(NULL) }
  gp <- do.call(rbind, parts)
  pal <- .cosmx_gene_palette(unique(gp$target), genes, palette)
  gp$color <- unname(pal[gp$target])
  attr(gp, "palette") <- pal
  message(nrow(gp), " transcripts across ", length(unique(gp$target)), " genes")
  gp
}

## ---------------------------------------------------------------------------
## PLOTTING
## ---------------------------------------------------------------------------

## The colorized multichannel image for a canvas-px rect, with this project's
## channel names, channel colours and per-sample normalisation filled in.
## Returns a TiffPlotData; the plot itself is at @plots$colorized.
##
## `precalc` accepts either the full named list from cosmx_precalc_max() or a
## single sample's data.frame. No channel coercion is needed here -- TiffPlotR's
## .handle_precalc_max() already selects the fetch resolution's rows.
##
## Example
##   img <- cosmx_image("IB1", rect) %>% add_um_scale_bar()
##   img@plots$colorized
cosmx_image <- function(sample_id, rect,
                        precalc = NULL, tiffs = NULL,
                        channel_names = cosmx_channel_names(),
                        channel_colors = cosmx_channel_colors(),
                        max_pixels = 800, ...) {
  if (is.null(tiffs))   tiffs   <- cosmx_tiff_paths()
  if (is.null(precalc)) precalc <- cosmx_precalc_max()

  tf <- tiffs[[sample_id]]
  if (is.null(tf))
    stop("no tiff for sample '", sample_id, "'; have: ",
         paste(names(tiffs), collapse = ", "), call. = FALSE)

  pm <- if (is.data.frame(precalc)) precalc else precalc[[sample_id]]
  if (is.null(pm))
    stop("no precalc_max for sample '", sample_id, "'; have: ",
         paste(names(precalc), collapse = ", "), call. = FALSE)

  fetchTiffData.colorized(tf, rect = rect,
                          channel_names = channel_names,
                          channel_colors = channel_colors,
                          precalc_max = pm,
                          max_pixels = max_pixels, ...)
}

## Add transcript points to a ggplot, coloured by gene, with a legend.
##
## Uses a `colour` scale on purpose, so it composes with boundary polygons that
## are filled via a `fill` scale. Point coords are canvas px, the same space as
## the boundaries and the image.
##
## Example
##   p <- cosmx_add_genes(img@plots$colorized, gp, size = 1.6)
cosmx_add_genes <- function(p, gp, size = 0.5, alpha = 1, shape = 16, legend = TRUE) {
  if (is.null(gp)) { message("no gene points to add"); return(p) }
  pal <- attr(gp, "palette")
  gp$target <- factor(gp$target, levels = names(pal))
  p +
    ggplot2::geom_point(data = gp,
                        ggplot2::aes(x = .data$X_px, y = .data$Y_px, colour = .data$target),
                        size = size, alpha = alpha, shape = shape) +
    ggplot2::scale_colour_manual(values = pal, drop = FALSE,
                                 guide = if (legend) "legend" else "none",
                                 name = "gene")
}
