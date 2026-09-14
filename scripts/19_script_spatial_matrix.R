## ============================================================================
## 19_script_spatial_matrix.R -- cell type VS cell type, spatially, as a matrix.
##
## The small-multiple map (celltype_spatial.png) answers "where is each type concentrated"
## one type at a time, and each panel is scaled to its own ceiling, so colour cannot be
## compared between panels. That makes the one question you actually want to ask of a tissue
## -- do these two types occupy the SAME places or different ones -- unanswerable from it.
##
## This is the pairwise version. Each panel holds one ordered pair and each bin carries
##      balance = (k_row - k_col) / (k_row + k_col)
## so +1 means the bin is entirely the ROW type, -1 entirely the COLUMN type, and 0 an even
## mix. Two types that co-localise give a pale panel; two that segregate give a panel split
## into red and blue territories. The statistic is a ratio, so it is immune to the overall
## cellularity that dominates a raw density map, and being bounded it needs no clamping.
##
## COLOURS follow the project's diverging rule -- blue to red, red the positive end -- so
## RED IS THE ROW TYPE and BLUE IS THE COLUMN TYPE. The matrix is antisymmetric by
## construction: the lower triangle is the upper triangle with the colours swapped. It is
## drawn in full anyway, because being able to read one row across is worth the redundancy.
## The diagonal is a type against itself and is left blank.
##
## Bins with fewer than MIN_PAIR cells of the two types combined are background: a balance
## computed over two cells is noise, and at +/-1 it would be the loudest thing on the figure.
##
## Facet strips carry each type's own ceiling (its 98th percentile bin fraction), the same
## number shown in celltype_spatial.png, so the two figures can be read together.
##
## THE with_qc VARIANT adds the over-fixation defect as an eleventh group. Those 1,647,768
## cells carry no cell-type label -- they are `unresolved`, and correctly so -- but they are
## not absent from the tissue, and every cell-type panel above is drawn on a section whose
## periphery is largely them. Putting the defect in as a row and a column shows, per cell
## type, whether its territory overlaps the damaged rim or avoids it. A type that sits
## against the defect is a type whose spatial claims are weakest, and that is not visible
## from a figure the defect has been removed from.
##
## Usage: sbatch --mem=96G --cpus-per-task=2 --time=2:00:00 \
##          scripts/sbatch/run_rscript.sh scripts/19_script_spatial_matrix.R [celltypes|with_qc]
## ============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); setDTthreads(2)
})
OUT <- "output/19_celltypes"; FIG <- file.path(OUT,"figures")
sf <- function(n,p,w,h){ggsave(file.path(FIG,n),p,width=w,height=h,dpi=150,bg="#fcfcfb",limitsize=FALSE);cat("wrote",file.path(FIG,n),"\n")}
## TWO ENCODINGS, AND `bivariate` IS THE BETTER ONE.
##   balance    (k_row - k_col)/(k_row + k_col) on a blue-yellow-red diverging scale. Purely
##              RELATIVE, and that is its flaw: it discards magnitude. A bin with 50 row-type
##              cells and no column-type reads +1, and so does a bin with 8 and none -- the
##              same saturated red. It answers "which of the two dominates here" and nothing
##              about how much of either is present.
##   bivariate  each type drives its own colour channel, mixed subtractively from white.
##              row-high/col-zero, row-middling/col-zero and both-high are three different
##              colours, so magnitude and balance are both legible. This is what a two-colour
##              row/column scheme should be.
args <- commandArgs(trailingOnly = TRUE)
MODE <- if (length(args) >= 1) args[1] else "celltypes"
ENC  <- if (length(args) >= 2) args[2] else "bivariate"
if (!ENC %in% c("balance","bivariate")) stop("encoding must be balance or bivariate", call.=FALSE)
if (!MODE %in% c("celltypes","with_qc")) stop("mode must be celltypes or with_qc", call.=FALSE)
BIN_MM <- 0.10; MIN_PAIR <- 8L
DEFECT <- "over-fixation defect"
SURF<-"#fcfcfb"; INK<-"#0b0b0b"; INK2<-"#52514e"; MUTED<-"#898781"
## THE MIDPOINT IS YELLOW, NOT PALE GREY, AND THE BACKGROUND IS DARKER THAN BEFORE.
## A diverging scale normally runs through near-white, but here near-white is also what an
## EMPTY bin looks like, and those two mean opposite things -- "these two types are evenly
## mixed here" against "there are no cells of either type here". With a #f0efec midpoint on
## an #e8e7e1 background they were indistinguishable, which is worst exactly where the
## interesting answer lives: a pale panel is the signature of co-localisation. The ends stay
## blue-to-red per the project's diverging rule; only the midpoint hue and the background
## tone change.
BG      <- "#cbcac3"   # no data -- neutral grey, deliberately darker than the midpoint
MID_COL <- "#f7efa8"   # even mix -- pale yellow, unmistakably a VALUE and not background
ROW_COL <- "#d03b3b"   # red  = the ROW type
COL_COL <- "#184f95"   # blue = the COLUMN type
ORD <- c("B cell","Plasmablast","Plasma cell","T cell","NK cell","Myeloid","pDC",
         "Endothelial","Stromal/FDC","Hodgkin chemokine niche")

## The bins are recomputed here rather than reused, because the defect cells are outside the
## high-confidence object and so absent from spatial_celltype_bins.csv.gz by construction.
a  <- fread(file.path(OUT,"all_cells_celltypes.csv.gz"), na.strings=c("NA",""),
            select=c("cell_id","samp","confidence","celltype_cluster","celltype_projected"))
md <- fread("output/16_regions/cells_meta_all.csv.gz",
            select=c("cell_id","x_slide_mm","y_slide_mm","is_defect"))
d <- md[a, on=.(cell_id)]
d[, cls := fifelse(!is.na(celltype_cluster), celltype_cluster, celltype_projected)]
if (MODE == "with_qc") {
  d[is.na(cls) & is_defect == TRUE, cls := DEFECT]
  ORD <- c(ORD, DEFECT)
}
d <- d[!is.na(cls) & cls %in% ORD & !is.na(x_slide_mm)]
d[, `:=`(bx = floor(x_slide_mm/BIN_MM), by_ = floor(y_slide_mm/BIN_MM))]
tot <- d[, .(n_bin = .N), by=.(samp, bx, by_)]
g <- d[, .(k = .N), by=.(samp, bx, by_, cls)]
g <- tot[g, on=.(samp,bx,by_)][, `:=`(frac = k/n_bin, x = bx*BIN_MM, y = by_*BIN_MM)]
cat(sprintf("mode %s | cells %d | bin-by-class rows %d | classes %d\n",
            MODE, nrow(d), nrow(g), uniqueN(g$cls)))
print(d[, .(cells = .N), by=.(cls, samp)][order(cls, samp)])

## each type's ceiling, matching celltype_spatial.png
ceil <- g[frac > 0, .(hi = quantile(frac, .98)), by=.(samp, cls)]
ceil[, strip := sprintf("%s\n<=%.0f%%", cls, 100*hi)]

for (S in c("IB1","RB1")) {
  z <- g[samp == S]
  cl <- intersect(ORD, unique(z$cls))
  ## every ordered pair; k is the count of that type in the bin, so the balance is exact
  w <- dcast(z, bx + by_ + x + y ~ cls, value.var="k", fill=0L)
  ## each type's own ceiling on the FRACTION, so minor types stay visible -- the same
  ## 98th-percentile convention the other spatial figures use. Defined BEFORE the pair loop
  ## because the bivariate scaling needs it there.
  cs <- ceil[samp == S]
  hi <- setNames(cs$hi, cs$cls)
  nb <- unique(z[, .(bx, by_, n_bin)])
  w  <- nb[w, on = .(bx, by_)]
  P <- rbindlist(lapply(cl, function(A) rbindlist(lapply(cl, function(B) {
    if (A == B) return(NULL)
    d <- w[, .(x, y, n_bin, kA = get(A), kB = get(B))]
    d <- d[kA + kB >= MIN_PAIR]
    if (!nrow(d)) return(NULL)
    d[, `:=`(a = pmin((kA/n_bin)/hi[[A]], 1), b = pmin((kB/n_bin)/hi[[B]], 1))]
    d[, .(row_t = A, col_t = B, x, y, bal = (kA - kB)/(kA + kB), a, b)]
  }))))
  ## subtractive mix from white: (0,0) white, (1,0) row colour, (0,1) column colour,
  ## (1,1) a dark blend. No-data bins stay the panel background, a grey that no mix produces.
  crgb <- grDevices::col2rgb(ROW_COL)/255; ccgb <- grDevices::col2rgb(COL_COL)/255
  P[, hex := grDevices::rgb(
      pmax(1 - a*(1-crgb[1]) - b*(1-ccgb[1]), 0),
      pmax(1 - a*(1-crgb[2]) - b*(1-ccgb[2]), 0),
      pmax(1 - a*(1-crgb[3]) - b*(1-ccgb[3]), 0))]
  P <- cs[, .(cls, strip)][P, on = c(cls = "row_t")]; setnames(P, c("cls","strip"), c("row_t","row_s"))
  P <- cs[, .(cls, strip)][P, on = c(cls = "col_t")]; setnames(P, c("cls","strip"), c("col_t","col_s"))
  P[, row_s := factor(row_s, levels = cs[match(cl, cls), strip])]
  P[, col_s := factor(col_s, levels = cs[match(cl, cls), strip])]
  cat(sprintf("\n%s: %d pair-bins over %d pairs\n", S, nrow(P), uniqueN(paste(P$row_t,P$col_t))))
  ## The encoding-dependent parts are PARENTHESISED single terms, so the shared chain that
  ## follows (facets, coord, labs, theme) applies to both. Written as a bare if/else the
  ## chain binds only to the else branch, and the bivariate figure silently came out as one
  ## unfaceted default-themed panel.
  lyr <- if (ENC == "bivariate") geom_tile(aes(fill = hex), width=BIN_MM, height=BIN_MM)
         else                    geom_tile(aes(fill = bal), width=BIN_MM, height=BIN_MM)
  scl <- if (ENC == "bivariate") scale_fill_identity()
         else scale_fill_gradient2(low = COL_COL, mid = MID_COL, high = ROW_COL, midpoint = 0,
                limits = c(-1,1), name = NULL, breaks = c(-1,0,1),
                labels = c("all column type","even mix","all row type"))
  p <- ggplot(P, aes(x, y)) + lyr + scl +
    facet_grid(row_s ~ col_s, switch = "y") + coord_fixed() +
    labs(title = sprintf("%s -- do two cell types occupy the same places?%s", S,
                         if (MODE == "with_qc") "  (over-fixation defect included as a group)" else ""),
         subtitle = if (ENC == "bivariate") paste0(
             "Each type drives its OWN colour channel: RED = the row type, BLUE = the column type, each scaled to its own ceiling (in the strips).\n",
             "Saturation carries HOW MUCH -- pale red is a little of the row type, deep red a lot. DARK PURPLE is both present together, i.e. co-localisation.\n",
             "WHITE is a bin holding cells but little of either. This fixes the diverging version, where 8-vs-0 and 50-vs-0 were the same red.")
           else paste0("Each bin: (row - column) / (row + column) over the two types' cell counts. RED = the row type, BLUE = the column type,\n",
                           "yellow = an even mix. NOTE: purely relative -- it discards magnitude, so 8-vs-0 and 50-vs-0 both read as full red.\n",
                           "Strips carry each type's own ceiling (98th-percentile bin fraction), as in celltype_spatial.png."),
         x = "slide x (mm)", y = "slide y (mm)",
         caption = sprintf("%g um bins; a bin needs %d cells of the pair combined. GREY = neither type present there -- a colour no mix produces, so it cannot be confused with %s.\nThe matrix is antisymmetric: the lower triangle is the upper with the two colours swapped. Diagonal blank. Plasmablast is absent -- at 5,056 cells it never reaches %d in a bin.",
                           BIN_MM*1000, MIN_PAIR,
                           if (ENC == "bivariate") "white (present, little of either)" else "yellow (present, evenly mixed)",
                           MIN_PAIR)) +
    theme_minimal(base_size = 9) +
    theme(plot.background=element_rect(fill=SURF,colour=NA),
          panel.background=element_rect(fill=BG,colour=NA), panel.grid=element_blank(),
          panel.spacing=unit(1,"pt"),
          axis.text=element_blank(), axis.ticks=element_blank(),
          axis.title=element_text(colour=INK2,size=8),
          strip.text.x=element_text(colour=INK2,size=5.5,face="bold",lineheight=1.05),
          strip.text.y.left=element_text(colour=INK2,size=5.5,face="bold",angle=0,lineheight=1.05),
          strip.placement="outside",
          plot.title=element_text(colour=INK,face="bold",size=12),
          plot.subtitle=element_text(colour=INK2,size=7.5),
          plot.caption=element_text(colour=MUTED,size=6,hjust=0),
          legend.position="bottom", legend.key.width=unit(40,"pt"), legend.key.height=unit(7,"pt"),
          legend.text=element_text(colour=INK2,size=7))
  side <- length(cl)
  ## 2D colour key: the whole (row, column) space as a square, so the reader can decode any
  ## panel colour without guessing.
  if (ENC == "bivariate") {
    kg <- CJ(a = seq(0,1,length.out=25), b = seq(0,1,length.out=25))
    kg[, hex := grDevices::rgb(
        pmax(1 - a*(1-crgb[1]) - b*(1-ccgb[1]), 0),
        pmax(1 - a*(1-crgb[2]) - b*(1-ccgb[2]), 0),
        pmax(1 - a*(1-crgb[3]) - b*(1-ccgb[3]), 0))]
    key <- ggplot(kg, aes(b, a, fill = hex)) + geom_raster() + scale_fill_identity() +
      coord_fixed() +
      labs(x = "more of the COLUMN type (blue) ->", y = "more of the ROW type (red) ->",
           title = "colour key: both channels at once") +
      theme_minimal(8) +
      theme(plot.background=element_rect(fill=SURF,colour=NA), panel.grid=element_blank(),
            axis.text=element_blank(), axis.ticks=element_blank(),
            axis.title=element_text(colour=INK2,size=7),
            plot.title=element_text(colour=INK2,size=7,face="bold"))
    ## Composed as its OWN panel beneath the matrix rather than as an inset. An inset is
    ## placed in figure-relative coordinates, and at bottom-left it landed on top of the
    ## last row's strip label -- and any hand-tuned position would move again the moment
    ## with_qc adds an eleventh row. coord_fixed on the key keeps it square inside a short
    ## strip, so it cannot be distorted either.
    p <- p / key + patchwork::plot_layout(heights = c(1, 0.085))
  }
  sf(sprintf("celltype_spatial_matrix%s%s_%s.png",
             if (MODE == "with_qc") "_with_qc" else "",
             if (ENC == "bivariate") "_bivar" else "", S),
     p, 1.55*side + 1.5, 1.55*side + 2.0)
}
cat("\ndone\n")
