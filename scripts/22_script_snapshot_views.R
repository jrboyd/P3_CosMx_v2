#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------------------
## step 22e -- the snapshots: tight crops of real cell pairs, on real morphology.
##
## WHAT EACH PANEL SHOWS, and why each layer is there.
##   image            the 5-channel composite, so a "cell" can be checked against morphology
##   cell outlines    every segmented boundary in the window -- this is what makes the
##                    leakage question answerable by eye: two positive cells 6 um apart
##                    either have two nuclei and two boundaries or they do not
##   transcripts      the LIGAND (cyan) and RECEPTOR (magenta) transcript positions, which
##                    is the actual evidence; CellComp is carried so an UNASSIGNED transcript
##                    sitting between the two cells is visible rather than inferred
##   focal markers    the two cells the pair is made of, ringed and labelled with their
##                    counts, cell type and the centre-to-centre distance
##
## COLOUR. Ligand cyan, receptor magenta -- neither is a channel colour (PanCK blue, G green,
## membrane yellow, CD45 red, DNA violet) and neither is a reserved disease hue (FL orange
## #eb6834, Hodgkin blue #0b3d7a), so nothing here collides with a convention.
##
## COORDINATES. Canvas px throughout, obtained ANALYTICALLY via cosmx_local_to_px() from the
## flat metadata's FOV-local centres -- not from the fitted mm->px affine, which is only
## needed when the source is slide mm. No axis reversal is applied anywhere: TiffPlotR's
## colorized plot establishes the canvas-px coordinate system already and every overlay
## (shape_annotate, cosmx_add_genes) is drawn in that same space. This follows
## scripts/06_tiff_demo.R verbatim, which is the worked example precisely so that THE AXIS
## RULE does not have to be re-derived per script.
##
## TWO VARIANTS PER PAIR, and the labelling is not decorative:
##   cherry  the clearest examples, ranked on min(ligand, receptor) counts. NOT typical.
##   repr    a uniform draw from the same candidate pool, <=2 per FOV. Judge on this one.
## Each figure's subtitle states which it is, the candidate pool size it was drawn from, and
## step 22's enrichment and q for that pair, so a panel cannot be read without its context.
##
## Usage: sbatch --mem=64G --time=6:00:00 scripts/sbatch/run_rscript.sh \
##          scripts/22_script_snapshot_views.R [tier] [pair-substring]
## ---------------------------------------------------------------------------------------
if (!file.exists("scripts/lib/cosmx_tiff.R")) stop("run from the project root", call. = FALSE)
source("scripts/lib/cosmx_tiff.R")
suppressPackageStartupMessages({library(ggplot2); library(patchwork); library(data.table)})
setDTthreads(1); options(width = 200)

AG    <- commandArgs(TRUE)
TIER  <- if (length(AG) >= 1) AG[1] else "matched"
ONLY  <- if (length(AG) >= 2) AG[2] else NA_character_
SNAP  <- file.path("output/22_lr", TIER, "snapshots")
FIG   <- file.path(SNAP, "figures"); dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
FLAT  <- "/netfiles/volaric_research/DLBCL/CosMx/export1/flatFiles/%s/%s_metadata_file.csv.gz"

LIG_COL <- "#00e5ff"; REC_COL <- "#ff2fd0"
PAD_UM  <- 22          # context beyond the two cells
MIN_UM  <- 55          # never crop tighter than this
INK <- "#111111"

P <- fread(file.path(SNAP, "snapshot_picks.csv"))
TG <- fread(file.path(SNAP, "snapshot_targets.csv"))
CP <- fread(file.path(SNAP, "snapshot_candidate_pool.csv"))
if (!is.na(ONLY)) P <- P[grepl(ONLY, pair, fixed = TRUE)]
stopifnot(nrow(P) > 0)

## ---- canvas px for the focal cells, analytically ---------------------------------------
px <- rbindlist(lapply(unique(P$samp), function(SM) {
  fl <- fread(sprintf(FLAT, SM, SM), select = c("cell","fov","CenterX_local_px","CenterY_local_px"))
  ## cosmx_seu_key() must be applied to BOTH sides -- it maps "c_2_100_1007" to the
  ## within-sample join key "100_1007". Keying only the flat side and comparing against raw
  ## cell_ids matched nothing and silently dropped every panel.
  need <- unique(cosmx_seu_key(c(P[samp == SM, L_cell], P[samp == SM, R_cell])))
  fl[, cid := cosmx_seu_key(cell)]
  fl <- fl[cid %in% need]
  q <- cosmx_local_to_px(SM, fl$fov, fl$CenterX_local_px, fl$CenterY_local_px)
  data.table(samp = SM, cid = fl$cid, fov = fl$fov, X_px = q$X_px, Y_px = q$Y_px)
}))
setkey(px, samp, cid)
gx <- function(sm, id, col) px[.(sm, id)][[col]]
P[, `:=`(L_key = cosmx_seu_key(L_cell), R_key = cosmx_seu_key(R_cell))]
P[, `:=`(L_X = gx(samp, L_key, "X_px"), L_Y = gx(samp, L_key, "Y_px"),
         R_X = gx(samp, R_key, "X_px"), R_Y = gx(samp, R_key, "Y_px")), by = seq_len(nrow(P))]
stopifnot(P[, sum(!is.na(L_X))] > 0)
miss <- P[is.na(L_X) | is.na(R_X)]
if (nrow(miss)) { cat("dropping", nrow(miss), "panels with no canvas coordinate\n")
                  P <- P[!is.na(L_X) & !is.na(R_X)] }
UMPX <- COSMX_MM_PER_PX * 1000        # um per canvas px

theme_snap <- theme_void(base_size = 8) +
  theme(plot.title = element_text(size = 6.0, colour = INK, hjust = 0, margin = margin(b = 1)),
        legend.position = "none",
        plot.margin = margin(1, 1, 1, 1))

one_panel <- function(r) {
  half_um <- max(MIN_UM, r$dist_um + 2*PAD_UM) / 2
  half    <- half_um / UMPX
  cx <- (r$L_X + r$R_X)/2; cy <- (r$L_Y + r$R_Y)/2
  rect <- cosmx_px_rect(cx - half, cx + half, cy - half, cy + half,
                        name = sprintf("%s_%s_%02d", r$pair, r$variant, r$panel))
  img <- tryCatch(cosmx_image(r$samp, rect) %>% add_um_scale_bar(), error = function(e) NULL)
  if (is.null(img)) return(NULL)
  p <- img@plots$colorized + theme_snap
  ## every boundary in the window: outlines only, so the image stays visible
  poly <- tryCatch(cosmx_cells_in_rect(r$samp, rect), error = function(e) NULL)
  if (!is.null(poly)) p <- (p %>% shape_annotate(poly, fill = "#00000000",
                                                 color = "#ffffff55")) +
                             guides(color = "none", fill = "none")
  ## transcripts
  pal <- c(LIG_COL, REC_COL); names(pal) <- c(r$ligand, r$receptor)
  gp <- tryCatch(cosmx_gene_points(r$samp, rect, c(r$ligand, r$receptor), palette = pal),
                 error = function(e) NULL)
  if (!is.null(gp) && nrow(gp)) {
    ## Per-FOV transcript files mark an unassigned transcript with CellId 0 (and an empty
    ## CellComp). The legend claims the white rings mean "assigned to no cell", so the claim
    ## is checked in the log rather than assumed.
    gp$unassigned <- is.na(gp$CellId) | gp$CellId == 0 | gp$CellComp == "" |
                     gp$CellComp == "None"
    .TALLY$tx <<- .TALLY$tx + nrow(gp)
    .TALLY$un <<- .TALLY$un + sum(gp$unassigned)
    .TALLY$comp <<- unique(c(.TALLY$comp, gp$CellComp))
    p <- cosmx_add_genes(p, gp, size = 1.5)
    ## ring the unassigned ones: they are the leakage-relevant subset
    if (any(gp$unassigned)) p <- p + geom_point(
      data = gp[gp$unassigned, , drop = FALSE], aes(X_px, Y_px),
      shape = 21, colour = "#ffffff", fill = NA, stroke = .35, size = 2.6,
      inherit.aes = FALSE)
    n_l <- sum(gp$target == r$ligand); n_r <- sum(gp$target == r$receptor)
    n_u <- sum(gp$unassigned)
  } else { n_l <- 0L; n_r <- 0L; n_u <- 0L }
  ## the two focal cells
  foc <- data.frame(X = c(r$L_X, r$R_X), Y = c(r$L_Y, r$R_Y),
                    col = c(LIG_COL, REC_COL))
  p <- p +
    geom_point(data = foc, aes(X, Y), shape = 3, size = 2.2, stroke = .5,
               colour = foc$col, inherit.aes = FALSE) +
    geom_point(data = foc, aes(X, Y), shape = 21, size = 7, stroke = .5,
               colour = foc$col, fill = NA, inherit.aes = FALSE) +
    annotate("segment", x = r$L_X, xend = r$R_X, y = r$L_Y, yend = r$R_Y,
             colour = "#ffffffcc", linewidth = .25, linetype = "22")
  ## Titles are kept short on purpose: the first version spelled out both gene names and
  ## 10-character cell types, which overran the panel width and collided with the neighbouring
  ## panel's title. Gene names are in the figure header, so the panel needs only counts.
  ab <- function(x) { x <- sub("Stromal/FDC","Strom",x); x <- sub(" cell","",x)
                      x <- sub("Endothelial","Endo",x); substr(x, 1, 6) }
  ## Clip to the window. cosmx_cells_in_rect() returns every cell OVERLAPPING the rect, so
  ## boundaries run past the image edge and render as grey flaps outside it. Limits are given
  ## SORTED, so this clips only -- it does not flip an axis (see THE AXIS RULE; orientation is
  ## already set by TiffPlotR's scales, never by coords).
  q <- rect@coords[1L, ]
  p <- p + coord_cartesian(xlim = sort(c(q$xmin, q$xmax)),
                           ylim = sort(c(q$ymin, q$ymax)), expand = FALSE)
  p + labs(title = sprintf("%d. %.0fum  L%d %s > R%d %s  tx%d/%d%s",
        r$panel, r$dist_um, r$L_count, ab(r$L_type), r$R_count, ab(r$R_type),
        n_l, n_r, if (n_u) sprintf(" %du", n_u) else ""))
}

for (pr in unique(P$pair)) for (vr in unique(P[pair == pr, variant])) {
  .TALLY <- list(tx = 0L, un = 0L, comp = character())
  D <- P[pair == pr & variant == vr][order(panel)]
  tg <- TG[pair == pr][1]; cp <- CP[pair == pr][1]
  cat(sprintf("\n=== %s / %s: %d panels\n", pr, vr, nrow(D)))
  ps <- lapply(seq_len(nrow(D)), function(i) {
    o <- tryCatch(one_panel(D[i]), error = function(e) { cat("  panel", i, "failed:",
                    conditionMessage(e), "\n"); NULL })
    if (is.null(o)) plot_spacer() else o
  })
  ## Wrapped explicitly: the unwrapped version ran off the right edge of the page.
  sub <- paste(strwrap(sprintf(paste0(
    "%s in %s %s, band %s um. Step 22: enrichment %.3f, q = %.3g. %s of %d candidate ",
    "pairs in %d FOVs. Cyan = %s transcripts, ringed cell is the ligand cell; magenta = ",
    "%s. Dashed line joins the pair; thin white outlines are segmentation boundaries. ",
    "%d transcripts drawn, %d assigned to no cell (white ring). NOTE: leakage assigns a ",
    "transcript to the WRONG cell, which no annotation can reveal -- judge that from ",
    "whether the two ringed cells have separate nuclei and boundaries. Panel titles: ",
    "distance, ligand count and cell type > receptor count and cell type, transcripts drawn."),
    pr, tg$samp, tg$bin_region, tg$band, tg$enrich_fov, tg$q,
    if (vr == "cherry") "CHERRY-PICKED, the clearest available" else
                        "REPRESENTATIVE, a uniform draw",
    cp$candidates, cp$fovs, tg$ligand, tg$receptor, .TALLY$tx, .TALLY$un),
    width = 168), collapse = "\n")
  cat(sprintf("  transcripts drawn %d, unassigned %d, CellComp seen: %s\n",
              .TALLY$tx, .TALLY$un, paste(sort(.TALLY$comp), collapse = "|")))
  g <- wrap_plots(ps, nrow = 3, ncol = 3) +
    plot_annotation(
      title = sprintf("%s  -  %s  (%s)", pr, toupper(vr), tg$status),
      subtitle = sub,
      theme = theme(plot.title = element_text(face = "bold", size = 12),
                    plot.subtitle = element_text(size = 7.0, colour = "#4a4844")))
  f <- file.path(FIG, sprintf("fig_22_snap_%s_%s.png", gsub("[^A-Za-z0-9]", "", pr), vr))
  ggsave(f, g, width = 9.6, height = 10.4, dpi = 190, bg = "white", limitsize = FALSE)
  cat("wrote", f, "\n")
}
writeLines(format(Sys.time()), file.path(SNAP, "snapshot_views.done"))
cat("\nOK\n")
