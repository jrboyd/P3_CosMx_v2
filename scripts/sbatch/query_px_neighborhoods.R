## ============================================================================
## Ask whether given canvas-pixel coordinates land in a scored light-chain
## neighbourhood, at any radius, in either table set.
##
## Works entirely in canvas px via the exact per-cell route, so no mm affine is
## involved: every light-chain-high cell is placed with cosmx_local_to_px() and the
## query point is compared against those cells directly.
##
## Reports, per query point x radius x table set:
##   in_nb        the scored neighbourhood whose members enclose the point
##                (point within span_px of the member centroid), or none
##   nearest_nb   nearest scored neighbourhood regardless, with its distance
##   nearest_cell distance to the closest light-chain-high cell of any kind
##
##   sbatch --mem=32G --time=0:30:00 scripts/sbatch/run_rscript.sh \
##     scripts/sbatch/query_px_neighborhoods.R IB1 112000,167700 99500,58500 63500,2000
##
## Usage: Rscript scripts/sbatch/query_px_neighborhoods.R <sample_id> <x,y> [<x,y> ...]
## ============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L)
  stop("Usage: Rscript scripts/sbatch/query_px_neighborhoods.R <sample_id> <x,y> [<x,y> ...]",
       call. = FALSE)
SAMPLE_ID <- args[1]
Q <- do.call(rbind, lapply(args[-1], function(a) {
  v <- as.numeric(strsplit(a, ",")[[1]])
  if (length(v) != 2L) stop("bad coordinate '", a, "', want x,y", call. = FALSE)
  data.frame(X_px = v[1], Y_px = v[2])
}))
Q$label <- sprintf("(%.0f, %.0f)", Q$X_px, Q$Y_px)

if (!file.exists("scripts/lib/cosmx_tiff.R"))
  stop("run from the project root", call. = FALSE)
source("scripts/lib/cosmx_tiff.R")
suppressPackageStartupMessages(library(data.table))
setDTthreads(1)

FLAT   <- "/netfiles/volaric_research/DLBCL/CosMx/export1/flatFiles/%s/%s_metadata_file.csv.gz"
RADII  <- c(20, 30, 40, 60)
SETS   <- c("output/04_merge", "output/07_lightchain")
MM     <- COSMX_MM_PER_PX

## canvas extent, so an out-of-range query is called out rather than silently
## matched to the nearest edge cell
pl <- cosmx_fov_placement(SAMPLE_ID)
ext <- c(x = max(pl$x0) + COSMX_FOV_PX - 1L, y = max(pl$y0) + COSMX_FOV_PX - 1L)
cat(sprintf("%s canvas extent: 0..%d x 0..%d px\n", SAMPLE_ID, ext["x"], ext["y"]))
for (r in seq_len(nrow(Q))) {
  oob <- Q$X_px[r] < 0 || Q$X_px[r] > ext["x"] || Q$Y_px[r] < 0 || Q$Y_px[r] > ext["y"]
  cat(sprintf("  %s  %s\n", Q$label[r],
              if (oob) "*** OUTSIDE the canvas ***" else "inside the canvas"))
}

## exact canvas px for every cell of this sample
flat <- fread(sprintf(FLAT, SAMPLE_ID, SAMPLE_ID),
              select = c("cell", "fov", "CenterX_local_px", "CenterY_local_px"))
px <- cosmx_local_to_px(SAMPLE_ID, flat$fov, flat$CenterX_local_px, flat$CenterY_local_px)
flat[, `:=`(key = cosmx_seu_key(cell), X_px = px$X_px, Y_px = px$Y_px)]
setkey(flat, key)

out <- list()
for (set in SETS) {
  for (R in RADII) {
    nb_f <- file.path(set, sprintf("lightchain_neighborhoods_r%dum.csv", R))
    lc_f <- file.path(set, sprintf("lightchain_cells_r%dum.csv", R))
    if (!file.exists(nb_f) || !file.exists(lc_f)) next
    nbs <- fread(nb_f)
    lc  <- fread(lc_f)
    lc[, key := cosmx_seu_key(cell_id)]
    lc <- merge(lc, flat[, .(key, X_px, Y_px)], by = "key", all.x = TRUE)
    lc <- lc[!is.na(X_px)]

    ## per-neighbourhood centroid and radius, in px, from the members themselves
    cen <- lc[nb %in% nbs$nb, .(cx = mean(X_px), cy = mean(Y_px),
                                span_px = max(sqrt((X_px - mean(X_px))^2 +
                                                   (Y_px - mean(Y_px))^2)),
                                n = .N), by = nb]
    cen <- merge(cen, nbs[, .(nb, kappa_frac, class, artifact_suspect, n_cells)], by = "nb")

    for (r in seq_len(nrow(Q))) {
      d_nb   <- sqrt((cen$cx - Q$X_px[r])^2 + (cen$cy - Q$Y_px[r])^2)
      inside <- which(d_nb <= cen$span_px)
      i_near <- which.min(d_nb)
      d_cell <- min(sqrt((lc$X_px - Q$X_px[r])^2 + (lc$Y_px - Q$Y_px[r])^2))
      pick   <- if (length(inside)) inside[which.min(d_nb[inside])] else NA_integer_
      out[[length(out) + 1L]] <- data.table(
        query        = Q$label[r],
        set          = basename(set),
        radius_um    = R,
        in_nb        = if (is.na(pick)) NA_integer_ else cen$nb[pick],
        in_class     = if (is.na(pick)) NA_character_ else cen$class[pick],
        in_kappa     = if (is.na(pick)) NA_real_ else cen$kappa_frac[pick],
        nearest_nb   = cen$nb[i_near],
        nearest_mm   = round(d_nb[i_near] * MM, 3),
        nearest_kappa = cen$kappa_frac[i_near],
        nearest_class = cen$class[i_near],
        nearest_flag  = cen$artifact_suspect[i_near],
        cell_mm      = round(d_cell * MM, 3))
    }
  }
}
res <- rbindlist(out)

cat("\n================ per radius and table set ================\n")
print(res, nrows = 200)

cat("\n================ verdict per query ================\n")
for (q in unique(res$query)) {
  z <- res[query == q]
  hit <- z[!is.na(in_nb)]
  cat("\n", q, "\n", sep = "")
  cat(sprintf("  closest light-chain-high cell: %.3f mm away (min over all sets)\n",
              min(z$cell_mm)))
  if (nrow(hit)) {
    cat("  INSIDE a scored neighbourhood:\n")
    print(hit[, .(set, radius_um, in_nb, in_kappa, in_class)])
  } else {
    cat("  not inside any scored neighbourhood. Nearest:\n")
    print(z[order(nearest_mm)][1:min(4, .N),
            .(set, radius_um, nearest_nb, nearest_mm, nearest_kappa, nearest_class, nearest_flag)])
  }
}
fwrite(res, "output/07_lightchain/query_px_neighborhoods.csv")
cat("\nwrote output/07_lightchain/query_px_neighborhoods.csv\n")
