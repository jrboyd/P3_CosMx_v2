## Pathology review bundle, stage 1 of 3: choose the cells.
##
##   sbatch --mem=16G scripts/sbatch/run_rscript.sh pathology_review_app/build/01_select_cells.R test
##   (or `full`; pathology_review_app/build/submit_build.sh runs all three stages)
##
## Writes pathology_review_app/output/<size>/selection.csv: every cell the renderer may
## use, in priority order within its stratum, with `need` = how many of that stratum the
## bundle wants. Strata carry BACKUPS (the renderer drops windows that are off-canvas,
## blank, or overlap a window already taken), so `need` is a target, not a row count.
##
## Design (CLAUDE.md, steps 28-31):
##   IB1 candidates  the step-29 ROI cells, the first k slots of each of the six groups
##                   (high/medium/low x random/top-score), so the review ties back to
##                   output/29_bright_cells/roi_png/.
##   RB1 candidates  RB1 has 9 cells >= 50 and none >= 100, so: 50-99 (all available),
##                   10-49 top HRS-leaning score, 10-49 random.
##   controls        same sample, chemo < 10, defect-free matched tier (the population
##                   the sources came from), each matched to a randomly chosen candidate
##                   of its own sample on area (+/-10%) and depth0 (+/-15%). Without them
##                   a "yes" rate has no reference (step 31's design).
##   range cells     16 random matched-tier cells per slide; the renderer sets each
##                   slide's display range from their windows, so the range does not
##                   depend on which cells were reviewed or on the bundle size.
##
## Repeats and the blinded order are assigned in stage 3, after image QC.
suppressMessages(library(data.table))
source("scripts/lib/cosmx_tiff.R")

SIZE <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "test"
DESIGN <- list(
  test = list(ib1_per_group = 4,  ib1_ctrl = 10, rb1_5099 = 3, rb1_top = 2, rb1_rand = 1, rb1_ctrl = 5),
  full = list(ib1_per_group = 20, ib1_ctrl = 40, rb1_5099 = 9, rb1_top = 6, rb1_rand = 5, rb1_ctrl = 20)
)[[SIZE]]
stopifnot(!is.null(DESIGN))
OUT <- file.path("pathology_review_app/output", SIZE)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
set.seed(38)
BACKUP <- 4                                    # candidates drawn per needed slot
SLIDE <- c(IB1 = 2L, RB1 = 1L)                 # c_2_ = IB1, c_1_ = RB1 (CLAUDE.md)
FLAT <- "input/CosMx/export1/flatFiles"

src <- fread("output/29_bright_cells/source_cells.csv.gz")
cols <- c("cell_id", "samp", "fov_id", "X_px", "Y_px", "chemo", "hrs_score", "k3",
          "Area.um2", "depth0")

## ---- IB1 candidates: the step-29 ROI cells, in slot order ------------------------
idx <- fread("output/29_bright_cells/roi_png/roi_png_index.csv")
ib1 <- merge(idx[, .(cell_id, src_group = group, priority = slot)], src[, ..cols], by = "cell_id")
stopifnot(nrow(ib1) == nrow(idx), all(ib1$samp == "IB1"))
ib1[, `:=`(kind = "candidate", stratum = paste0("IB1_", src_group), need = DESIGN$ib1_per_group)]

## ---- RB1 candidates ---------------------------------------------------------------
rb <- src[samp == "RB1"]
pick <- function(d, n, ord) {
  d <- d[ord][seq_len(min(.N, n * BACKUP))]
  d[, priority := seq_len(.N)]
}
r_hi  <- pick(rb[chemo >= 50], DESIGN$rb1_5099, sample(nrow(rb[chemo >= 50])))
lo    <- rb[chemo < 50]
r_top <- pick(lo, DESIGN$rb1_top, order(-lo$hrs_score))
lo2   <- lo[!cell_id %in% r_top$cell_id]
r_rnd <- pick(lo2, DESIGN$rb1_rand, sample(nrow(lo2)))
r_hi[,  `:=`(stratum = "RB1_50-99",      src_group = "RB1_50-99",      need = DESIGN$rb1_5099)]
r_top[, `:=`(stratum = "RB1_10-49_hrs",  src_group = "RB1_10-49_hrs",  need = DESIGN$rb1_top)]
r_rnd[, `:=`(stratum = "RB1_10-49_rand", src_group = "RB1_10-49_rand", need = DESIGN$rb1_rand)]
rb1 <- rbind(r_hi, r_top, r_rnd)[, c(cols, "src_group", "priority", "stratum", "need"), with = FALSE]
rb1[, kind := "candidate"]
cand <- rbind(ib1[, names(rb1), with = FALSE], rb1)

## ---- controls and range cells -----------------------------------------------------
reg <- fread(cmd = "zcat output/16_regions/matched/cell_regions.csv.gz",
             select = c("cell_id", "samp", "fov_id", "is_defect", "tier_matched",
                        "Area.um2", "nCount_RNA", "chemo_ct"))
reg <- reg[tier_matched == TRUE & is_defect == FALSE]
reg[, depth0 := nCount_RNA - chemo_ct]

match_controls <- function(s, n_need) {
  ## reference = the candidates the bundle will actually show first (priority <= need)
  ref <- cand[samp == s & priority <= need]
  pool <- reg[samp == s & chemo_ct < 10 & depth0 > 0]
  got <- character(0)
  out <- vector("list", n_need * BACKUP)
  for (i in seq_along(out)) {
    r <- ref[sample(.N, 1)]
    ok <- pool[abs(log(Area.um2 / r$Area.um2)) < 0.10 & abs(log(depth0 / r$depth0)) < 0.15 &
               !cell_id %in% got]
    if (!nrow(ok)) stop("no control match for ", r$cell_id)
    c1 <- ok[sample(.N, 1)]
    got <- c(got, c1$cell_id)
    out[[i]] <- c1[, .(cell_id, samp, fov_id, chemo = chemo_ct, Area.um2, depth0,
                       matched_to = r$cell_id)]
  }
  d <- rbindlist(out)
  d[, `:=`(kind = "control", stratum = paste0(s, "_control"), src_group = paste0(s, "_control"),
           need = n_need, priority = seq_len(.N), hrs_score = NA_real_, k3 = NA_integer_)]
}
ctrl <- rbind(match_controls("IB1", DESIGN$ib1_ctrl), match_controls("RB1", DESIGN$rb1_ctrl))

rng <- reg[, .SD[sample(.N, 16)], by = samp][, .(cell_id, samp, fov_id, chemo = chemo_ct,
                                                 Area.um2, depth0)]
rng[, `:=`(kind = "range", stratum = paste0(samp, "_range"), src_group = NA_character_,
           need = 16L, priority = seq_len(.N), hrs_score = NA_real_, k3 = NA_integer_,
           matched_to = NA_character_), by = samp]

## canvas px from FOV-local centroids (the analytic transform; CLAUDE.md "Images")
nc <- rbind(ctrl, rng, fill = TRUE)
xy <- rbindlist(lapply(names(SLIDE), function(s) {
  m <- fread(file.path(FLAT, s, paste0(s, "_metadata_file.csv.gz")),
             select = c("fov", "cell_ID", "CenterX_local_px", "CenterY_local_px"))
  m[, cell_id := sprintf("c_%d_%d_%d", SLIDE[[s]], fov, cell_ID)]
  m <- m[cell_id %in% nc[samp == s]$cell_id]
  p <- cosmx_local_to_px(s, m$fov, m$CenterX_local_px, m$CenterY_local_px)
  m[, .(cell_id, X_px = p$X_px, Y_px = p$Y_px)]
}))
nc <- merge(nc, xy, by = "cell_id")
stopifnot(nrow(nc) == nrow(ctrl) + nrow(rng), !anyNA(nc$X_px))

cand[, matched_to := NA_character_]
sel <- rbind(cand, nc, use.names = TRUE)
setorder(sel, samp, kind, stratum, priority)
fwrite(sel, file.path(OUT, "selection.csv"))

cat(sprintf("%s: wrote %d rows to %s\n", SIZE, nrow(sel), OUT))
print(sel[, .(rows = .N, need = need[1], med_area = median(Area.um2), med_depth0 = median(depth0),
              med_chemo = median(chemo)), by = .(samp, kind, stratum)])
