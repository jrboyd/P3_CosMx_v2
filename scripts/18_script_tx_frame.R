## ============================================================================
## 18_script_tx_frame.R -- transcripts in the DONOR-ALIGNED frame, plus the per-cell statistic.
##
## THE FRAME. Every focal cell is re-expressed with its own centroid at the origin and its
## A-marker donor straight UP: the donor centroid lands at (0, donorA_dist) for every cell,
## so pooling thousands of cells into one 2D density is meaningful. Leakage shows up as
## A-marker density displaced toward +y relative to the cell's other transcripts.
##
## WHY THE ROTATION IS DONE IN A NEGATED-Y FRAME. Canvas px runs opposite to slide mm on
## BOTH axes (CLAUDE.md, THE AXIS RULE), so rotating raw px would produce a MIRROR of the
## tissue, not a rotation of it. Negating y first makes the frame right-handed and physical,
## which matters the moment anyone reads left-vs-right in a panel. The rotation is asserted
## below: the donor must land on the +y axis to floating-point precision, or the frame is
## wrong and nothing downstream is interpretable.
##
## THE STATISTIC IS WITHIN-CELL AND NEEDS NO PERMUTATION LOOP. For a cell holding N
## transcripts of which k are family-A markers, ask: is the mean projection of the k
## A-markers onto the donor direction larger than the mean projection of all N? Under the
## null "the k A-markers are a random subset of this cell's N transcripts", the sampling
## distribution of that mean is exactly the sampling-without-replacement one --
##   E = mubar,  Var = (sigma^2 / k) * (N - k) / (N - 1)
## -- so the z-score is analytic and exact in mean and variance. No permutation is drawn,
## nothing is estimated from other cells, and the null cannot be broken by cell size, cell
## shape, elongation, local density, sequencing depth or k. That is the whole reason for
## working within-cell: step 17 built three BETWEEN-cell nulls and all three failed their
## own positive control.
##
## THREE DIRECTIONS, SAME CELLS, SAME TRANSCRIPTS. z_A is the test. z_C projects onto a
## donor of a third family (not the argmax call, not the cluster call) and z_rand onto a
## random direction. A cell-shape artefact moves all three; leakage moves only z_A.
##
## Usage: sbatch --mem=96G --time=3:00:00 scripts/sbatch/run_rscript.sh \
##          scripts/18_script_tx_frame.R IB1
## ============================================================================
suppressPackageStartupMessages({ library(data.table); setDTthreads(2) })
if (!file.exists("scripts/lib/lineage.R")) stop("run from the project root", call. = FALSE)
source("scripts/lib/lineage.R")

args   <- commandArgs(trailingOnly = TRUE)
SAMPLE <- if (length(args) >= 1) args[1] else "IB1"
OUT    <- "output/18_leakage"
TX_F   <- file.path(OUT, sprintf("%s_tx_raw.csv.gz",     SAMPLE))
FOC_F  <- file.path(OUT, sprintf("%s_focal_cells.csv.gz", SAMPLE))
CENT_F <- file.path(OUT, sprintf("%s_centroids.csv.gz",   SAMPLE))
for (f in c(TX_F, FOC_F, CENT_F)) if (!file.exists(f)) stop("missing ", f, call. = FALSE)
sc <- function(n, d) { fwrite(d, file.path(OUT, n)); cat("wrote ", file.path(OUT, n),
                                                         " (", nrow(d), " rows)\n", sep = "") }
UMPX      <- 120.280945e-6 * 1000
D_EXTRACT <- 25
FAMS      <- names(LINEAGE_SETS_COARSE)
GENE2FAM  <- setNames(rep(names(LINEAGE_SETS_COARSE), lengths(LINEAGE_SETS_COARSE)),
                      unlist(LINEAGE_SETS_COARSE, use.names = FALSE))

cat("============================================================\n")
cat("step 18 -- donor-aligned transcript frame: ", SAMPLE, "\n", sep = "")
cat("============================================================\n")

foc <- fread(FOC_F)
tx  <- fread(TX_F)
cat(sprintf("focal cells %d | transcripts read %d | distinct cells in tx %d\n",
            nrow(foc), nrow(tx), uniqueN(tx$cell_id)))
tx[, fam := GENE2FAM[target]]
tx[is.na(fam), fam := ""]
cat("\n=== transcripts by compartment ===\n"); print(tx[, .N, by = CellComp][order(-N)])
cat("\n=== marker-family transcripts ===\n");  print(tx[fam != "", .N, by = fam][order(-N)])

## ---- geometry helpers -------------------------------------------------------
## dx/dy: transcript minus a reference centroid, in physical um (y negated: see header).
## rot: rotate by phi = pi/2 - theta so that the direction `theta` becomes +y.
rot_xy <- function(dx, dy, theta) {
  phi <- pi/2 - theta
  list(x = dx*cos(phi) - dy*sin(phi), y = dx*sin(phi) + dy*cos(phi))
}

## ---- 1. focal cells' own transcripts ---------------------------------------
FK <- c("cell_id","A","cx_px","cy_px","theta_A","theta_C","theta_rand")
f1 <- tx[foc[, ..FK], on = .(cell_id), nomatch = NULL]
cat(sprintf("\nfocal transcripts: %d over %d cells\n", nrow(f1), uniqueN(f1$cell_id)))
f1[, `:=`(dx = (x_px - cx_px) * UMPX, dy = -(y_px - cy_px) * UMPX)]
f1[, c("xr","yr") := rot_xy(dx, dy, theta_A)]
f1[, is_A := as.integer(fam == A)]

## ASSERTION: the donor must land on the +y axis. If it does not, the frame is a mirror or
## the angle convention is inverted, and every density map built from it is wrong.
##
## The donor centroid is re-read from the centroid table rather than carried over from the
## selection script, so this checks the whole chain -- donor id -> its centroid -> theta_A
## -> the rotation -- instead of only the rotation formula. A mirror would show up here as
## an error of order the cell spacing (um), not of order floating point.
ce <- fread(CENT_F, select = c("cell_id","cx_px","cy_px"))
ce[, `:=`(dX_um = cx_px * UMPX, dY_um = -cy_px * UMPX)][, c("cx_px","cy_px") := NULL]
chk <- foc[is.finite(donorA_dist), .(cell_id, X_um, Y_um, donorA_id, donorA_dist, theta_A)]
chk <- ce[, .(donorA_id = cell_id, dX_um, dY_um)][chk, on = .(donorA_id), nomatch = NULL]
dchk <- rot_xy(chk$dX_um - chk$X_um, chk$dY_um - chk$Y_um, chk$theta_A)
err_x <- max(abs(dchk$x)); err_y <- max(abs(dchk$y - chk$donorA_dist))
cat(sprintf("\nframe check -- donor lands at (0, d) for %d cells: max |x| = %.3g um, max |y - d| = %.3g um\n",
            nrow(chk), err_x, err_y))
if (err_x > 1e-4 || err_y > 1e-4) stop("rotation is wrong: donor does not land on +y", call. = FALSE)
rm(ce, chk, dchk); invisible(gc())

## ---- 2. the per-cell statistic, three directions ---------------------------
## proj(theta) is the projection of each transcript onto the unit vector at angle theta.
## The z-score is the sampling-without-replacement null described in the header.
zstat <- function(p, isA, N) {
  k <- sum(isA); if (k == 0L || k >= N) return(c(NA_real_, NA_real_, NA_real_))
  mu <- mean(p); s2 <- sum((p - mu)^2) / (N - 1)
  if (s2 <= 0) return(c(NA_real_, NA_real_, NA_real_))
  mA <- sum(p * isA) / k
  v  <- (s2 / k) * ((N - k) / (N - 1))
  c(mA - mu, (mA - mu) / sqrt(v), mu)
}
f1[, c("x_px","y_px","cx_px","cy_px") := NULL]
f1[, `:=`(pA    = dx*cos(theta_A)    + dy*sin(theta_A),
          pC    = dx*cos(theta_C)    + dy*sin(theta_C),
          pRand = dx*cos(theta_rand) + dy*sin(theta_rand))]
f1[, c("dx","dy","theta_A","theta_C","theta_rand") := NULL]
invisible(gc())
cat("\ncomputing per-cell polarisation ...\n")
pol <- f1[, {
  N <- .N; a <- is_A
  rA <- zstat(pA, a, N); rC <- zstat(pC, a, N); rR <- zstat(pRand, a, N)
  .(n_tx = N, k_A = sum(a),
    delta_A = rA[1], z_A = rA[2], delta_C = rC[1], z_C = rC[2],
    delta_rand = rR[1], z_rand = rR[2],
    ## compartment split, A-markers vs everything else in the SAME cell. A gene-matched
    ## comparison is done downstream -- compartment bias is gene-dependent, so these
    ## counts are raw material, not a result.
    nucA = sum(a == 1L & CellComp == "Nuclear"), cytA = sum(a == 1L & CellComp == "Cytoplasm"),
    memA = sum(a == 1L & CellComp == "Membrane"),
    nucO = sum(a == 0L & CellComp == "Nuclear"), cytO = sum(a == 0L & CellComp == "Cytoplasm"),
    memO = sum(a == 0L & CellComp == "Membrane"))
}, by = cell_id]
f1[, c("pA","pC","pRand") := NULL]; invisible(gc())
pol <- foc[pol, on = .(cell_id)]
sc(sprintf("%s_polarisation.csv.gz", SAMPLE), pol)

cat("\n=== median z by argmax call x agreement (scoreable, adjacent, k_A >= 1) ===\n")
print(pol[scoreable == TRUE & adjacent == TRUE & k_A >= 1,
          .(cells = .N, z_A = round(median(z_A, na.rm = TRUE), 4),
            z_C = round(median(z_C, na.rm = TRUE), 4),
            z_rand = round(median(z_rand, na.rm = TRUE), 4)),
          by = .(A, agree)][order(A, agree)])

## ---- 3. donor transcripts, in the FOCAL cell's frame -----------------------
## Only the donor's MARKER transcripts: they are the context that makes a density map
## readable (the donor should appear as a lobe of family-A signal at +y), and carrying all
## ~200 transcripts per donor for every focal cell it serves would triple the table for
## background that the focal cell already supplies.
dn <- foc[is.finite(donorA_dist) & donorA_dist <= D_EXTRACT,
          .(cell_id, A, agree, cx_px, cy_px, theta_A, donorA_id)]
d1 <- tx[fam != "", .(donorA_id = cell_id, target, fam, CellComp, z, x_px, y_px)][
        dn, on = .(donorA_id), nomatch = NULL, allow.cartesian = TRUE]
d1[, `:=`(dx = (x_px - cx_px) * UMPX, dy = -(y_px - cy_px) * UMPX)]
d1[, c("xr","yr") := rot_xy(dx, dy, theta_A)]
d1[, is_A := as.integer(fam == A)]
cat(sprintf("\ndonor marker transcripts placed in focal frames: %d\n", nrow(d1)))

## ---- 4. the frame table ----------------------------------------------------
## TWO FILES, not one with a `role` column. The focal block alone is 77M rows; rbind()-ing
## the donor block onto it would double peak memory, and fwrite(append=TRUE) onto a .gz is
## not a supported combination. Downstream reads the focal file and only joins the donor
## file when a figure wants the neighbour drawn in.
KEEP <- c("cell_id","target","fam","is_A","CellComp","z","xr","yr")
f1[, `:=`(xr = round(xr, 2), yr = round(yr, 2))]
d1[, `:=`(xr = round(xr, 2), yr = round(yr, 2))]
cat("\n=== frame rows ===\n")
print(rbind(f1[, .(block = "focal", .N), by = is_A], d1[, .(block = "donor", .N), by = is_A])[order(block, -N)])
for (b in list(list("focal", f1), list("donor", d1))) {
  fn <- file.path(OUT, sprintf("%s_tx_frame_%s.csv.gz", SAMPLE, b[[1]]))
  fwrite(b[[2]][, ..KEEP], fn)
  cat("wrote ", fn, " (", nrow(b[[2]]), " rows, ",
      round(file.size(fn)/2^30, 2), " GiB)\n", sep = "")
}
cat("\ndone\n")
