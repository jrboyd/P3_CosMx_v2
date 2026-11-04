## Step 34, polygon stage: the vendor cell outlines for the step-31 cells, in canvas px.
## Reads input/CosMx/export1/flatFiles/IB1/IB1-polygons.csv.gz (one row per vertex, FOV-local
## px) through an awk prefilter on (fov, cellID), maps to canvas px with the analytic
## transform, and checks each polygon's centroid against the cell's recorded centroid.
##   Rscript scripts/34_extract_polygons.R [cells.csv] [out_dir]
source("scripts/lib/cosmx_tiff.R")
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(trailingOnly = TRUE)
CELLS <- if (length(a) >= 1) a[1] else "output/31_radial_cd45/cells.csv"
OUT   <- if (length(a) >= 2) a[2] else "output/34_nuclear_size"
POLY  <- "input/CosMx/export1/flatFiles/IB1/IB1-polygons.csv.gz"
cells <- fread(CELLS)
cells[, c("fov", "cell_ID") := .(as.integer(tstrsplit(cell_id, "_")[[3]]), as.integer(tstrsplit(cell_id, "_")[[4]]))]
keys <- file.path(OUT, "poly_keys.tmp")
fwrite(cells[, .(fov, cell_ID)], keys, col.names = FALSE)
pg <- fread(cmd = sprintf("zcat %s | awk -F, 'NR==FNR{k[$1\",\"$2]; next} FNR==1 || (($1\",\"$2) in k)' %s -", POLY, keys))
file.remove(keys)
setnames(pg, c("fov", "cellID"), c("fov", "cell_ID"))
pg[, cell_id := sprintf("c_2_%d_%d", fov, cell_ID)]       # c_2_ = IB1 (CLAUDE.md)
pg[, order := seq_len(.N), by = cell_id]
xy <- cosmx_local_to_px("IB1", pg$fov, pg$x_local_px, pg$y_local_px)
pg[, `:=`(X_px = xy$X_px, Y_px = xy$Y_px)]
chk <- merge(pg[, .(cx = mean(X_px), cy = mean(Y_px), n_vertex = .N), by = cell_id], cells[, .(cell_id, X_px, Y_px)], by = "cell_id")
chk[, off_um := sqrt((cx - X_px)^2 + (cy - Y_px)^2) * COSMX_MM_PER_PX * 1000]
cat(sprintf("polygons: %d of %d cells; vertices median %d; vertex-mean vs recorded centroid: median %.2f um, 95th pct %.2f um\n",
            uniqueN(pg$cell_id), nrow(cells), as.integer(median(chk$n_vertex)), median(chk$off_um), quantile(chk$off_um, .95)))
stopifnot(uniqueN(pg$cell_id) >= 0.98 * nrow(cells), median(chk$off_um) < 2)
fwrite(pg[order(cell_id, order), .(cell_id, order, X_px = round(X_px, 1), Y_px = round(Y_px, 1))], file.path(OUT, "polygons_canvas.csv"))
fwrite(chk[, .(cell_id, n_vertex, off_um = round(off_um, 2))], file.path(OUT, "polygon_check.csv"))
