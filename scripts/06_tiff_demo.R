## ============================================================================
## 06_tiff_demo.R -- end-to-end smoke test and worked example for
## scripts/lib/cosmx_tiff.R.
##
## Renders one window of an assembled OME-TIFF four ways:
##   A  image only, with a um scale bar
##   B  image + cell boundaries, neighbours coloured distinctly
##   C  image + transcript points for a couple of genes
##   D  cell boundaries coloured by a per-cell metadata column
##
## Everything here is driven from fov_positions.csv and a flat CSV -- no Seurat
## object is loaded, so the only reason this needs a compute node is the TIFF
## fetch itself.
##
##   sbatch --mem=64G --time=1:00:00 --job-name=cosmx-demo \
##     scripts/sbatch/run_rscript.sh scripts/06_tiff_demo.R IB1 400 2000
##
## Usage: Rscript scripts/06_tiff_demo.R <sample_id> [fov] [px_size]
##   sample_id  IB1 | RB1
##   fov        FOV number to centre the window on   (default 400)
##   px_size    window side, in canvas px            (default 4256, one whole FOV)
## ============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (!length(args))
  stop("Usage: Rscript scripts/06_tiff_demo.R <sample_id> [fov] [px_size]", call. = FALSE)
SAMPLE_ID <- args[1]
FOV       <- if (length(args) >= 2) as.integer(args[2]) else 400L
PX_SIZE   <- if (length(args) >= 3) as.numeric(args[3]) else 4256

if (!file.exists("scripts/lib/cosmx_tiff.R"))
  stop("run from the project root: cd /gpfs1/pi/avolaric/files_jrboyd/P3_CosMx_v2",
       call. = FALSE)
source("scripts/lib/cosmx_tiff.R")
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

## ------------------ CONFIG ----
OUT_DIR   <- "output/06_tiff_demo"
FIG_DIR   <- file.path(OUT_DIR, "figures")
QC_CSV    <- "output/02_filter_QC/qc_spatial_cells.csv.gz"
## A per-cell column, so panel D shows structure rather than one flat colour.
## This one also double-checks the mapping: it flags cells near a FOV edge in
## FOV-LOCAL px, and because x is mirrored on the canvas, the "low_x" band must
## come out on the RIGHT of the image and "high_x" on the left.
QC_COLUMN <- "qc_spatial_fovpx"
GENES     <- c(IGKC = "red", "IGLC1/2" = "cyan")
TAG       <- sprintf("%s_fov%d_%dpx", SAMPLE_ID, FOV, round(PX_SIZE))

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
done_f <- file.path(OUT_DIR, paste0(TAG, ".done"))
out_f  <- file.path(FIG_DIR, paste0(TAG, ".png"))

theme_demo <- theme(panel.background = element_blank(),
                    panel.grid       = element_blank(),
                    legend.position  = "bottom")

## ------------------ 1. the window ----
## Centre a square window on the chosen FOV. cosmx_px_rect() normalises the
## bounds, so it does not matter which corner is which.
pl <- cosmx_fov_placement(SAMPLE_ID)
i  <- match(FOV, pl$fov)
if (is.na(i))
  stop("FOV ", FOV, " is not in ", SAMPLE_ID, " (have ", min(pl$fov), "-",
       max(pl$fov), ")", call. = FALSE)

cx <- pl$x0[i] + (COSMX_FOV_PX - 1L) / 2
cy <- pl$y0[i] + (COSMX_FOV_PX - 1L) / 2
rect <- cosmx_px_rect(cx - PX_SIZE / 2, cx + PX_SIZE / 2,
                      cy - PX_SIZE / 2, cy + PX_SIZE / 2, name = TAG)
message(sprintf("[%s] FOV %d centre = (%.0f, %.0f) canvas px; window %.0f px",
                SAMPLE_ID, FOV, cx, cy, PX_SIZE))

## ------------------ 2. panel A: the image ----
## cosmx_image() fills in the tiff path, channel names/colours and this sample's
## precalc_max normalisation, which is what fixes the contrast across figures.
img   <- cosmx_image(SAMPLE_ID, rect) %>% add_um_scale_bar()
p_img <- img@plots$colorized + theme_demo + labs(title = "image")

## ------------------ 3. panel B: cell boundaries ----
## cosmx_cells_in_rect() returns UNCOLOURED polygons. Colour them however you
## like afterwards -- here just so neighbours are distinguishable. Annotating
## with fill = "#00000000" draws outlines only, leaving the image visible.
poly <- cosmx_cells_in_rect(SAMPLE_ID, rect)
if (is.null(poly)) stop("no cells in this window -- try another FOV", call. = FALSE)

poly_adj <- assign_adjacent_colors(poly, fill = TRUE, border = TRUE)
p_cells  <- (p_img %>% shape_annotate(poly_adj, fill = "#00000000")) +
  guides(color = "none", fill = "none") +
  labs(title = sprintf("%d cell boundaries", nrow(poly@coords)))

## ------------------ 4. panel C: transcript points ----
## Points are canvas px, the same space as the image and the boundaries.
gp      <- cosmx_gene_points(SAMPLE_ID, rect, names(GENES), palette = GENES)
p_genes <- cosmx_add_genes(p_img, gp, size = 1.2) +
  labs(title = if (is.null(gp)) "transcripts (none found)"
               else sprintf("%d transcripts", nrow(gp)))

## ------------------ 5. panel D: boundaries coloured by metadata ----
## The join key is (fov, cell) WITHIN a sample -- ids like "c_2_1_1" occur in both
## samples, so the metadata table must be filtered to this sample first.
qc <- data.table::fread(QC_CSV, select = c("cell", "sample_id", "fov", QC_COLUMN))
qc <- as.data.frame(qc[qc$sample_id == SAMPLE_ID])
message(nrow(qc), " metadata rows for ", SAMPLE_ID)

poly_qc <- cosmx_color_cells(poly, qc, QC_COLUMN, id_col = "cell", border = TRUE)
p_meta  <- (p_img %>% shape_annotate(poly_qc, alpha = 0.55)) +
  labs(title = QC_COLUMN)

## The value -> colour mapping is on the coloured polygons, but shape_annotate()
## does not surface it as a legend (the panel shows the channel legend instead),
## so log it. See pkg_tiffplotR_features.txt in the v1 project.
cmap <- attr(shape_meta(poly_qc), "color_map")   # data.frame(value, fill)
if (is.data.frame(cmap))
  message(QC_COLUMN, " colours: ",
          paste(cmap$value, cmap$fill, sep = "=", collapse = "  "))

## ------------------ 6. assemble ----
pg <- (p_img + p_cells) / (p_genes + p_meta)
ggsave(out_f, plot = pg, width = 14, height = 12, dpi = 150)
message("wrote ", out_f)

writeLines(format(Sys.time()), done_f)
