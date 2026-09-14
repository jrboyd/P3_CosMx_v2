## ============================================================================
## 15_verify_imagej_export.R -- verify the ImageJ export end to end.
##
## The ROI writer (scripts/lib/imagej_roi.py) is hand-rolled, so it needs an
## independent reader to check it against. TiffPlotR::readImageJRois() is exactly
## that: a different codebase, written before this one, parsing the int16 path
## that ImageJ itself falls back on. If the corners survive that round trip, the
## bytes are right.
##
## Checks:
##   V1  every overview ROI recovers its corners from the zip to < 1 px
##   V2  every crop ROI likewise, and lands inside its crop's bounds
##   V3  the crop manifest covers every region exactly once
##   V4  the written TIFFs report the pixel size Fiji will use for its scale bar
##   V5  scaled overview coords are consistent with canvas coords / 2^level
##
##   sbatch --mem=16G --time=1:00:00 scripts/sbatch/run_rscript.sh \
##     scripts/15_verify_imagej_export.R IB1 5
## ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(TiffPlotR)
})

args   <- commandArgs(trailingOnly = TRUE)
SAMPLE <- if (length(args) >= 1) args[1] else "IB1"
LEVEL  <- if (length(args) >= 2) as.integer(args[2]) else 5L

OUT_DIR   <- "output/15_imagej"
CROPS_DIR <- file.path(OUT_DIR, "crops")
TOL       <- 1.0   # px; the ROI header stores integer corners, so <1 px is exact

fail <- 0L
ok   <- function(fmt, ...) cat(sprintf(paste0("  ok    ", fmt, "\n"), ...))
bad  <- function(fmt, ...) { fail <<- fail + 1L
                             cat(sprintf(paste0("  FAIL  ", fmt, "\n"), ...)) }

cat("============================================================\n")
cat("step 15 verify -- ", SAMPLE, " L", LEVEL, "\n", sep = "")
cat("============================================================\n")

REG <- fread(file.path(OUT_DIR, sprintf("%s_regions_canvas_px.csv", SAMPLE)))

## Recover (x, y, w, h) from whatever shape class readImageJRois() returned.
.roi_boxes <- function(r) {
  out <- list()
  if (!is.null(r$rects) && nrow(r$rects@coords)) {
    d <- as.data.table(r$rects@coords)
    out[[length(out) + 1L]] <- d[, .(name, rx = xmin, ry = ymin,
                                     rw = xmax - xmin, rh = ymax - ymin)]
  }
  if (!is.null(r$ellipses) && nrow(r$ellipses@coords)) {
    d <- as.data.table(r$ellipses@coords)
    out[[length(out) + 1L]] <- d[, .(name, rx = x0 - radius_x, ry = y0 - radius_y,
                                     rw = 2 * radius_x, rh = 2 * radius_y)]
  }
  if (!length(out)) return(data.table(name = character(), rx = numeric(),
                                      ry = numeric(), rw = numeric(), rh = numeric()))
  rbindlist(out)
}

## Names carry a metric suffix, so match on the id prefix instead.
.match_ids <- function(names, ids) {
  vapply(names, function(n) {
    hit <- ids[startsWith(n, ids)]
    if (!length(hit)) NA_character_ else hit[which.max(nchar(hit))]
  }, character(1))
}

## ---- V1: overview ROIs round-trip ---------------------------------------------
cat("\nV1 overview ROI round-trip\n")
ov_zip <- file.path(OUT_DIR, sprintf("%s_overview_L%d_RoiSet.zip", SAMPLE, LEVEL))
ov_csv <- file.path(OUT_DIR, sprintf("%s_regions_L%d.csv", SAMPLE, LEVEL))
if (!file.exists(ov_zip)) {
  cat("  skip  no ", ov_zip, " (overview not generated?)\n", sep = "")
} else {
  E <- fread(ov_csv)
  B <- .roi_boxes(readImageJRois(ov_zip))
  if (nrow(B) != nrow(E)) bad("zip has %d ROIs, csv has %d", nrow(B), nrow(E))
  B[, id := .match_ids(name, E$id)]
  if (anyNA(B$id)) bad("%d ROI name(s) do not map to a region id", sum(is.na(B$id)))
  M <- merge(B, E, by = "id")
  if (nrow(M) != nrow(E)) bad("merged %d of %d regions", nrow(M), nrow(E))
  d <- M[, pmax(abs(rx - x), abs(ry - y),
                abs(rw - width), abs(rh - height))]
  if (max(d) > TOL) {
    bad("max corner error %.2f px > %.1f", max(d), TOL)
    print(M[which.max(d), .(id, rx, x, ry, y, rw, width, rh, height)])
  } else {
    ok("%d ROIs round-trip, max corner error %.2f px", nrow(M), max(d))
  }
  ## V5 folded in here: the scaled coords must equal canvas coords / 2^level
  S <- merge(E[, .(id, x, y)], REG[, .(id, x1, y1)], by = "id")
  ds <- S[, max(pmax(abs(x - x1 / 2^LEVEL), abs(y - y1 / 2^LEVEL)))]
  if (ds > 0.02) bad("V5 scaled coords differ from canvas/2^%d by %.3f px", LEVEL, ds)
  else ok("V5 scaled coords match canvas px / 2^%d (max %.4f px)", LEVEL, ds)
}

## ---- V2 + V3: crops -----------------------------------------------------------
cat("\nV2 crop ROI round-trip, V3 manifest coverage\n")
man_f <- file.path(CROPS_DIR, sprintf("%s_crop_manifest.csv", SAMPLE))
if (!file.exists(man_f)) {
  cat("  skip  no ", man_f, " (crops not generated?)\n", sep = "")
} else {
  MAN <- fread(man_f)
  worst <- 0; n_roi <- 0L; oob <- 0L
  for (i in seq_len(nrow(MAN))) {
    m  <- MAN[i]
    zf <- file.path(CROPS_DIR, sub("\\.tif$", "_RoiSet.zip", m$crop_file))
    if (!file.exists(zf)) { bad("missing %s", zf); next }
    B   <- .roi_boxes(readImageJRois(zf))
    ids <- strsplit(m$region_ids, " ")[[1]]
    if (nrow(B) != length(ids))
      bad("%s: %d ROIs vs %d manifest ids", basename(zf), nrow(B), length(ids))
    B[, id := .match_ids(name, ids)]
    E <- REG[id %in% ids]
    M <- merge(B, E, by = "id")
    ## crop-local expectation: canvas corner minus the window origin
    d <- M[, pmax(abs(rx - (x1 - m$x0)), abs(ry - (y1 - m$y0)),
                  abs(rw - (x2 - x1)),   abs(rh - (y2 - y1)))]
    worst <- max(worst, if (length(d)) max(d) else 0)
    n_roi <- n_roi + nrow(M)
    oob   <- oob + M[, sum(rx < 0 | ry < 0 |
                           rx + rw > m$width | ry + rh > m$height)]
  }
  if (worst > TOL) bad("crop ROI max corner error %.2f px > %.1f", worst, TOL)
  else ok("%d crop ROIs across %d crops round-trip, max error %.2f px",
          n_roi, nrow(MAN), worst)
  if (oob > 0) bad("%d crop ROI(s) fall outside their crop bounds", oob)
  else ok("every crop ROI lies inside its crop")

  all_ids <- unlist(strsplit(MAN$region_ids, " "))
  dup <- all_ids[duplicated(all_ids)]
  miss <- setdiff(REG$id, all_ids)
  extra <- setdiff(all_ids, REG$id)
  ## Only the exported groups are expected in the manifest.
  miss <- intersect(miss, REG[group %in% unique(REG[id %in% all_ids, group]), id])
  if (length(dup))   bad("V3 %d region(s) appear in more than one crop: %s",
                         length(dup), paste(head(dup, 5), collapse = " "))
  if (length(miss))  bad("V3 %d exported region(s) missing from every crop: %s",
                         length(miss), paste(head(miss, 5), collapse = " "))
  if (length(extra)) bad("V3 %d manifest id(s) are not regions", length(extra))
  if (!length(dup) && !length(miss) && !length(extra))
    ok("V3 manifest covers each of %d exported regions exactly once", length(all_ids))
}

## ---- V4: what Fiji will read out of the written TIFFs ------------------------
## Not via RBioFormats: its coreMetadata has no pixel-size field for a plain
## ImageJ TIFF. tiff::readTIFF(payload = FALSE) reads the tags only (no pixel
## data) and hands back x.resolution in pixels-per-unit, plus the ImageJ
## description string that carries the unit and channel count.
cat("\nV4 pixel size and ImageJ metadata in the written TIFFs\n")
UMPX0 <- 120.280945e-6 * 1000

.check_tif <- function(f, exp_um, exp_ch) {
  a <- suppressWarnings(tiff::readTIFF(f, payload = FALSE))
  got_um <- 1 / a$x.resolution[1]
  if (abs(got_um - exp_um) / exp_um > 0.01)
    bad("%s: %.5f um/px, expected %.5f", basename(f), got_um, exp_um)
  else
    ok("%s: %.4f um/px (expected %.4f)", basename(f), got_um, exp_um)

  desc <- a$description[1]
  if (!grepl("unit=um", desc, fixed = TRUE))
    bad("%s: ImageJ header does not say unit=um", basename(f))
  else if (!grepl("mode=composite", desc, fixed = TRUE))
    bad("%s: ImageJ header is not composite mode", basename(f))
  else if (!grepl(sprintf("channels=%d", exp_ch), desc))
    bad("%s: ImageJ header does not report %d channels", basename(f), exp_ch)
  else
    ok("%s: unit=um, composite, channels=%d", basename(f), exp_ch)
}

ov_tif <- file.path(OUT_DIR, sprintf("%s_overview_L%d.tif", SAMPLE, LEVEL))
if (file.exists(ov_tif)) .check_tif(ov_tif, UMPX0 * 2^LEVEL, 5L)
if (file.exists(man_f)) {
  cf <- file.path(CROPS_DIR, fread(man_f)$crop_file[1])
  if (file.exists(cf)) .check_tif(cf, UMPX0, 5L)
}

cat("\n============================================================\n")
if (fail == 0L) cat("ALL CHECKS PASSED\n") else cat(fail, " CHECK(S) FAILED\n", sep = "")
cat("============================================================\n")
if (fail > 0L) quit(status = 1L)
