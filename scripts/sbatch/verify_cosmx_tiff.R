## ============================================================================
## Assertion suite for scripts/lib/cosmx_tiff.R
##
## Runs on the login node: <= 1 CPU, well under 4 GB, no readRDS of any project
## Seurat object, no TIFF access. V3 is the expensive one at ~90 s.
##
##   cd /gpfs1/pi/avolaric/files_jrboyd/P3_CosMx_v2
##   R_LIBS_USER=/gpfs1/home/j/r/jrboyd/R/x86_64-pc-linux-gnu-library/4.5 \
##     /gpfs1/sw/rh9/pkgs/stacks/gcc/13.3.0/R/4.5.1/bin/Rscript --vanilla \
##     scripts/sbatch/verify_cosmx_tiff.R
##
## V1 constants and sacred text        V6 per-FOV path derivation
## V2 old-vs-new regression (the gate) V7 gene palette truth table
## V3 whole-dataset exactness (gate)   V8 precalc_max fixture
## V4 cell ids and join keys           V9 invariants a tidy-up would break
## V5 rect normalisation, FOV edge
##
## Pass `only=V3` (or a comma list) as an argument to run a subset.
## ============================================================================

if (!file.exists("scripts/lib/cosmx_tiff.R"))
  stop("run from the project root: cd /gpfs1/pi/avolaric/files_jrboyd/P3_CosMx_v2",
       call. = FALSE)

LIB    <- "scripts/lib/cosmx_tiff.R"
V1_DIR <- "/gpfs1/pi/avolaric/files_jrboyd/P3_CosMx"
V1_POS <- file.path(V1_DIR, "fov_positions.csv")
FLAT   <- "/netfiles/volaric_research/DLBCL/CosMx/export1/flatFiles/%s/%s_metadata_file.csv.gz"
SAMPLES <- c("IB1", "RB1")

args <- commandArgs(trailingOnly = TRUE)
only <- NULL
if (length(args)) {
  hit <- grep("^only=", args, value = TRUE)
  if (length(hit)) only <- strsplit(sub("^only=", "", hit[1]), ",")[[1]]
}

source(LIB)

n_pass <- 0L
check <- function(label, expr) {
  ok <- isTRUE(tryCatch(expr, error = function(e) structure(FALSE, msg = conditionMessage(e))))
  if (!ok) {
    res <- tryCatch(expr, error = function(e) conditionMessage(e))
    stop("FAIL ", label, if (is.character(res)) paste0(" -- ", res) else "", call. = FALSE)
  }
  n_pass <<- n_pass + 1L
  cat("  ok  ", label, "\n", sep = "")
}
run <- function(id, fun) {
  if (!is.null(only) && !id %in% only) return(invisible(NULL))
  cat("\n== ", id, " ==\n", sep = "")
  fun()
}

## ---------------------------------------------------------------------------
## V1 -- constants, and the sacred expressions present verbatim
## ---------------------------------------------------------------------------
run("V1", function() {
  check("COSMX_FOV_PX is 4256L",       identical(COSMX_FOV_PX, 4256L))
  check("COSMX_MM_PER_PX unchanged",   identical(COSMX_MM_PER_PX, 120.280945e-6))
  check("pos csv default is in v2",    identical(COSMX_POS_CSV, "input/fov_positions.csv"))
  check("staged pos csv exists",       file.exists(COSMX_POS_CSV))
  check("staged quantile cache exists", file.exists(COSMX_GQUANT_CACHE))

  src <- readLines(LIB, warn = FALSE)
  ## the x mirror, in both the centroid and the boundary transform
  check("x mirror present twice",
        sum(grepl("(COSMX_FOV_PX - 1L) - ", src, fixed = TRUE)) >= 3L)
  ## the half-open FOV extent must NOT have been "fixed" to +FOV_PX-1
  check("FOV rects stay half-open",
        any(grepl("xmax = pl$x0 + COSMX_FOV_PX", src, fixed = TRUE)) &&
        !any(grepl("xmax = pl$x0 + COSMX_FOV_PX - 1", src, fixed = TRUE)))
  check("FOV regex unchanged",
        any(grepl('".*_F0*([0-9]+)\\\\.TIF$"', src, fixed = TRUE)))
  check("round() on both placement axes",
        sum(grepl("round((d$X_mm - xmin) / COSMX_MM_PER_PX)", src, fixed = TRUE),
            grepl("round((d$Y_mm - ymin) / COSMX_MM_PER_PX)", src, fixed = TRUE)) == 2L)

  ## the instrument constant is duplicated in 02_spatial_QC.R on purpose; keep in sync
  qc <- readLines("scripts/02_spatial_QC.R", warn = FALSE)
  check("02_spatial_QC.R agrees on 4256", any(grepl("FOV_PX\\s*<-\\s*4256L", qc)))
})

## ---------------------------------------------------------------------------
## V2 -- old vs new, identical() not all.equal(). The gate on rename + dedup.
##
## The v1 files are sys.source'd into their own env. Note 20_/21_ each source()
## 16_ at top level with local = FALSE, so v1's NM_MM / FOV / POS_CSV_DEFAULT land
## in globalenv; harmless, since no v1 name collides with a cosmx_ name.
## ---------------------------------------------------------------------------
run("V2", function() {
  v1 <- new.env(parent = globalenv())
  suppressMessages({
    sys.source(file.path(V1_DIR, "assemble_cosmx/16_cosmx_tiff_mapping.R"), v1)
    sys.source(file.path(V1_DIR, "20_cell_boundaries.R"), v1)
    sys.source(file.path(V1_DIR, "21_gene_points.R"), v1)
  })

  set.seed(1)
  for (s in SAMPLES) {
    o <- v1$load_fov_placement(s, V1_POS)
    n <- cosmx_fov_placement(s, V1_POS)
    check(paste0(s, ": placement fov/x0/y0 identical"),
          identical(o$fov, n$fov) && identical(o$x0, n$x0) && identical(o$y0, n$y0))
    check(paste0(s, ": placement gained the file column"),
          is.character(n$file) && nrow(n) == nrow(o))

    ## exact transform on synthetic FOV-local input
    f  <- sample(o$fov, 300, replace = TRUE)
    xl <- runif(300, 0, COSMX_FOV_PX - 1); yl <- runif(300, 0, COSMX_FOV_PX - 1)
    po <- v1$fov_local_to_tiff_px(f, xl, yl, s, V1_POS)
    pn <- cosmx_local_to_px(s, f, xl, yl, V1_POS)
    check(paste0(s, ": local->px identical on 300 synthetic points"),
          identical(po$X_px, pn$X_px) && identical(po$Y_px, pn$Y_px))

    ## FOV overlap: v1's fovs_in_rect vs the now-single implementation
    pl <- n
    for (k in c(1L, ceiling(nrow(pl) / 2), nrow(pl))) {
      r <- cosmx_px_rect(pl$x0[k] - 500, pl$x0[k] + 5000,
                         pl$y0[k] - 500, pl$y0[k] + 5000, name = "w")
      fo <- v1$fovs_in_rect(r, s, V1_POS)
      fn <- cosmx_fovs_in_rect(s, r, V1_POS)
      check(paste0(s, ": fovs_in_rect identical around row ", k),
            identical(fo$fov, fn$fov) && identical(fo$x0, fn$x0) &&
            identical(fo$y0, fn$y0) && identical(fo$file, fn$file))
    }
  }

  ## boundary loader: the inline placement copy is gone, geometry must not move
  fovs <- c(448, 449)
  suppressMessages({
    bo <- v1$load_fov_boundaries("IB1", fovs, V1_POS)
    bn <- cosmx_fov_boundaries("IB1", fovs, V1_POS)
  })
  check("boundaries: same cell count",  identical(nrow(bo@coords), nrow(bn@coords)))
  check("boundaries: same names",       identical(bo@coords$name, bn@coords$name))
  check("boundaries: identical x",      identical(unlist(bo$x), unlist(bn$x)))
  check("boundaries: identical y",      identical(unlist(bo$y), unlist(bn$y)))

  ## key builders and the shared match
  ids <- c("c_2_368_802", "IB1_c_2_368_802", "RB1_c_1_1_1")
  pn  <- c("F368_C802", "F1_C1")
  check("seu key identical to v1",  identical(v1$seu_id2key(ids), cosmx_seu_key(ids)))
  check("poly key identical to v1", identical(v1$cosmx_id2key(pn), cosmx_poly_key(pn)))
  check("match identical to v1",
        identical(v1$match_cell_id_from_fov(ids, pn), cosmx_match_ids(ids, pn)))

  ## The precalc reduction, against v1's own function. v1 hardcoded a
  ## cwd-relative cache path, so run it from the directory holding the cache.
  pv1 <- local({
    wd <- getwd(); on.exit(setwd(wd), add = TRUE)
    setwd(dirname(COSMX_GQUANT_CACHE))
    tryCatch(suppressMessages(suppressWarnings(v1$getCosMxPrecalcMax())),
             error = function(e) e)
  })
  check("v1 precalc ran", !inherits(pv1, "error") && identical(names(pv1), SAMPLES))
  pn2 <- cosmx_precalc_max()
  for (s in SAMPLES)
    check(paste0(s, ": precalc identical to v1"),
          identical(as.data.frame(pv1[[s]]), as.data.frame(pn2[[s]])))
})

## ---------------------------------------------------------------------------
## V3 -- whole-dataset exactness against NanoString's own global pixel frame.
##
## The assembly mirrored x and stacked FOVs, so canvas px and the vendor's
## global px run in opposite directions on x; their SUM must therefore be a
## single constant across every cell. That one identity validates round(), the
## (FOV_PX - 1) mirror, the whole X_mm/Y_mm table and the per-sample origin at
## once. No TIFF, no Seurat.
## ---------------------------------------------------------------------------
run("V3", function() {
  suppressPackageStartupMessages(library(data.table))
  setDTthreads(1)
  expect <- list(IB1 = c(124851, 166305), RB1 = c(116537, 157517))

  for (s in SAMPLES) {
    md <- fread(sprintf(FLAT, s, s),
                select = c("fov", "CenterX_local_px", "CenterY_local_px",
                           "CenterX_global_px", "CenterY_global_px"))
    pl <- cosmx_fov_placement(s, COSMX_POS_CSV)
    px <- cosmx_local_to_px(s, md$fov, md$CenterX_local_px, md$CenterY_local_px,
                            COSMX_POS_CSV)
    cx <- unique(px$X_px + md$CenterX_global_px)
    cy <- unique(px$Y_px + md$CenterY_global_px)
    check(sprintf("%s: x sum is one constant over %s cells", s, format(nrow(md), big.mark = ",")),
          length(cx) == 1L)
    check(paste0(s, ": y sum is one constant"), length(cy) == 1L)
    check(paste0(s, ": x constant == max(x0) + COSMX_FOV_PX - 1"),
          cx == max(pl$x0) + COSMX_FOV_PX - 1L)
    check(paste0(s, ": y constant == max(y0)"), cy == max(pl$y0))
    check(paste0(s, ": canvas extent matches the recorded value"),
          cx == expect[[s]][1] && cy == expect[[s]][2])
  }
  rm(md); invisible(gc(FALSE))
})

## ---------------------------------------------------------------------------
## V4 -- cell ids and join keys
## ---------------------------------------------------------------------------
run("V4", function() {
  ids <- c("c_2_368_802", "IB1_c_2_368_802", "RB1_c_1_1_1")
  check("fov parsed from the end",  identical(cosmx_fov_from_id(ids),  c(368L, 368L, 1L)))
  check("cell parsed from the end", identical(cosmx_cell_from_id(ids), c(802L, 802L, 1L)))
  check("seu key is prefix-robust",
        identical(cosmx_seu_key(ids), c("368_802", "368_802", "1_1")))
  check("poly key from F<fov>_C<cell>",
        identical(cosmx_poly_key(c("F368_C802", "F1_C1")), c("368_802", "1_1")))
  check("poly key tolerates padding",
        identical(cosmx_poly_key("F00448_C00802"), "448_802"))
  check("unmatched poly name gives NA_NA",
        identical(cosmx_poly_key("nonsense"), "NA_NA"))
  check("match is order-safe by poly names",
        identical(cosmx_match_ids(ids, c("F368_C802", "F1_C1")),
                  c("c_2_368_802", "RB1_c_1_1_1")))
  check("keys join across the two id styles",
        identical(cosmx_seu_key("IB1_c_2_448_12"), cosmx_poly_key("F448_C12")))
  check("boundary names are never zero-padded",
        !any(grepl("F%05d_C%05d", readLines(LIB, warn = FALSE), fixed = TRUE)))
})

## ---------------------------------------------------------------------------
## V5 -- rect normalisation, and the half-open FOV edge
## ---------------------------------------------------------------------------
run("V5", function() {
  r <- cosmx_px_rect(1000, 200, 3000, 500, name = "z")
  check("px rect normalises inverted bounds",
        r$xmin == 200 && r$xmax == 1000 && r$ymin == 500 && r$ymax == 3000)
  check("px rect leaves ordered bounds alone", {
    r2 <- cosmx_px_rect(200, 1000, 500, 3000)
    r2$xmin == 200 && r2$xmax == 1000 && r2$ymin == 500 && r2$ymax == 3000
  })

  pl <- cosmx_fov_placement("IB1")
  i  <- which(pl$fov == 400L)
  ## a rect starting inside the FOV's LAST pixel column must still hit it
  inside <- cosmx_px_rect(pl$x0[i] + COSMX_FOV_PX - 0.5, pl$x0[i] + COSMX_FOV_PX + 10,
                          pl$y0[i] + 10, pl$y0[i] + 20)
  check("last pixel column of a FOV is inside it",
        400L %in% cosmx_fovs_in_rect("IB1", inside)$fov)
  ## a rect starting exactly at the half-open boundary must not
  outside <- cosmx_px_rect(pl$x0[i] + COSMX_FOV_PX, pl$x0[i] + COSMX_FOV_PX + 10,
                           pl$y0[i] + 10, pl$y0[i] + 20)
  hit <- cosmx_fovs_in_rect("IB1", outside)
  check("touching the half-open edge is outside",
        is.null(hit) || !(400L %in% hit$fov))
  ## a rect well away from the sample returns NULL, silently
  check("no overlap returns NULL",
        is.null(cosmx_fovs_in_rect("IB1", cosmx_px_rect(-1e6, -9e5, -1e6, -9e5))))
})

## ---------------------------------------------------------------------------
## V6 -- per-FOV path derivation (stat calls and globs only)
## ---------------------------------------------------------------------------
run("V6", function() {
  for (s in SAMPLES) {
    pl <- cosmx_fov_placement(s)
    k  <- pl[unique(round(seq(1, nrow(pl), length.out = 20L))), , drop = FALSE]
    check(paste0(s, ": morphology TIFs exist"), all(file.exists(k$file)))
    check(paste0(s, ": derived CellBoundaries csvs exist"),
          all(file.exists(mapply(.cosmx_boundary_path, k$file, k$fov))))
    check(paste0(s, ": exactly one transcript csv per FOV"),
          all(mapply(function(f, v) length(.cosmx_transcript_path(f, v)),
                     k$file, k$fov) == 1L))
  }
})

## ---------------------------------------------------------------------------
## V7 -- gene palette truth table
## ---------------------------------------------------------------------------
run("V7", function() {
  p <- .cosmx_gene_palette(present = c("CD3D", "MS4A1"),
                           genes   = c("MS4A1", "CD3D"),
                           palette = c(MS4A1 = "cyan"))
  check("requested gene order is kept",  identical(names(p), c("MS4A1", "CD3D")))
  check("named palette overrides by name", p[["MS4A1"]] == "cyan")
  check("unnamed genes get the auto palette", p[["CD3D"]] == "#377eb8")

  p2 <- .cosmx_gene_palette(c("A", "B"), c("A", "B"), palette = c("red", "blue"))
  check("unnamed palette is positional", p2[["A"]] == "red" && p2[["B"]] == "blue")

  p3 <- .cosmx_gene_palette(LETTERS[1:12], LETTERS[1:12], NULL)
  check("more than 8 genes falls through to hcl.colors",
        length(p3) == 12L && identical(unname(p3), grDevices::hcl.colors(12, "Dark 3")))

  check("requested genes sort before unrequested",
        identical(names(.cosmx_gene_palette(c("B", "A"), "A", NULL)), c("A", "B")))
})

## ---------------------------------------------------------------------------
## V8 -- precalc_max fixture. Pins the contrast of every figure.
## ---------------------------------------------------------------------------
run("V8", function() {
  pm <- cosmx_precalc_max()
  check("one table per sample", identical(names(pm), c("IB1", "RB1")))
  check("10 resolutions x 5 channels each", all(vapply(pm, nrow, 1L) == 50L))
  check("resolution column retained",
        all(vapply(pm, function(d) "resolution" %in% names(d), TRUE)))
  check("channel is numeric",
        all(vapply(pm, function(d) is.numeric(d$channel), TRUE)))
  check("column order unchanged from v1",
        all(vapply(pm, function(d)
          identical(names(d), c("sample_id", "resolution", "channel",
                                "max_value", "min_value")), TRUE)))

  at <- function(d, res) {
    z <- as.data.frame(d[d$resolution == res, , drop = FALSE])
    z[order(z$channel), ]
  }
  chk <- function(d, res, mx, mn) {
    z <- at(d, res)
    all(z$max_value == mx) && all(z$min_value == mn)
  }
  check("IB1 resolution 1", chk(pm$IB1, 1, c(1800, 360, 7184, 14464, 10268),
                                          c( 156, 104,  304,   284,   140)))
  check("IB1 resolution 3", chk(pm$IB1, 3, c(1509, 318, 7715,  9101,  6362),
                                          c( 152, 109,  160,   237,   131)))
  check("RB1 resolution 1", chk(pm$RB1, 1, c(1824, 440, 4924,  5940,  9348),
                                          c( 152, 104,  240,   252,   132)))
  check("RB1 resolution 3", chk(pm$RB1, 3, c( 978, 283, 3304,  3165,  8715),
                                          c( 148, 107,  148,   165,   141)))

  ## the min_q branch: q30 is 0 everywhere, so it is filtered out before the
  ## grouping and min_value is "smallest non-zero quantile". Collapsing the
  ## branch would set every min_value to 0.
  gq <- readRDS(COSMX_GQUANT_CACHE)
  pv <- as.numeric(sub("%", "", rownames(gq$IB1$signal_quantiles[[1]]))) / 100
  check("both 0.30 and 0.99 are present in the cache", any(pv == 0.30) && any(pv == 0.99))
  check("q30 is zero, so the else branch is the live one",
        all(gq$IB1$signal_quantiles[[1]][pv == 0.30, ] == 0))
  check("min_value is never 0", all(vapply(pm, function(d) all(d$min_value > 0), TRUE)))

  ## memo must not change the answer
  check("memoised call is identical", identical(pm, cosmx_precalc_max()))
})

## ---------------------------------------------------------------------------
## V9 -- invariants a plausible tidy-up would break
## ---------------------------------------------------------------------------
run("V9", function() {
  cc <- cosmx_channel_colors()
  check("channel colours are a list", is.list(cc))
  check("five channel colours",       length(cc) == 5L)
  check("channel names and order",
        identical(cosmx_channel_names(), c("PanCK", "G", "Membrane", "CD45", "DNA")))
  check("names and colours line up",  length(cc) == length(cosmx_channel_names()))

  suppressPackageStartupMessages(library(ggplot2))
  p <- suppressMessages(cosmx_add_genes(ggplot(), NULL))
  check("NULL gene points is a no-op", length(p$layers) == 0L)

  check("sample column opt-out is honoured",
        is.na(.cosmx_sample_col(NA, c("a", "b"))))
  check("sample column is auto-detected",
        identical(.cosmx_sample_col(NULL, c("a", "sample_id")), "sample_id"))
  check("Run_Tissue_name wins when present",
        identical(.cosmx_sample_col(NULL, c("sample_id", "Run_Tissue_name")),
                  "Run_Tissue_name"))
  check("missing sample column stops",
        inherits(tryCatch(.cosmx_sample_col(NULL, c("a", "b")), error = function(e) e),
                 "error"))

  check("mm_to_px rejects a non-2-column input",
        inherits(tryCatch(cosmx_mm_to_px(matrix(1:9, 3), list(ax = 1, bx = 0, ay = 1, by = 0)),
                          error = function(e) e), "error"))
  check("mm_to_px reads columns positionally", {
    aff <- list(ax = -2, bx = 10, ay = 3, by = -5)
    r <- cosmx_mm_to_px(data.frame(whatever = 1, other = 2), aff)
    r$X_px == 8 && r$Y_px == 1
  })
})

cat("\n", n_pass, " checks passed\n", sep = "")
