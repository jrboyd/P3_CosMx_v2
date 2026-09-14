## ============================================================================
## 19_script_spatial_density.R -- where each cell type and each QC class sits in the tissue.
##
## TWO FIGURES:
##   celltype_spatial.png  one panel per cell type per sample -- where is each concentrated
##   qc_spatial.png        one panel per QC / confidence class -- where the flags fall and
##                         where the rescue actually recovered cells
##
## WHAT IS PLOTTED IS A LOCAL FRACTION, NOT A RAW DENSITY. A raw density map of one cell type
## cannot separate "no cells of this type here" from "no tissue here", and it also just
## redraws overall cellularity: the densest part of the section looks like the densest part
## for every type. Each bin instead carries the SHARE of its cells that carry the label, so
## the map answers "concentrated" rather than "numerous". Bins under MIN_CELLS are drawn as
## background, because a fraction over four cells is noise.
##
## COORDINATES ARE SLIDE mm AND THE AXES ARE PLAIN. CLAUDE.md's axis rule in one line:
## canvas px reverses BOTH axes, slide mm reverses NEITHER. Working in mm therefore needs no
## scale_*_reverse, and adding one would mirror the tissue without erroring -- which is
## exactly how step 16 shipped two mirrored panels. The samples have different extents, so
## each is drawn on its own coord_fixed rather than sharing a scale.
##
## THE CELL TYPE USED IS THE COMBINED ONE: the cluster label where it exists (high
## confidence), the projected label where the rescue passed its stratum's gate. That is the
## point of the rescue -- it is what makes the maps dense enough to read -- and the QC figure
## shows exactly which parts of the tissue that added.
##
## Usage: sbatch --mem=96G --cpus-per-task=2 --time=2:00:00 \
##          scripts/sbatch/run_rscript.sh scripts/19_script_spatial_density.R
## ============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); setDTthreads(2)
})
OUT <- "output/19_celltypes"; FIG <- file.path(OUT,"figures")
sf <- function(n,p,w,h){ggsave(file.path(FIG,n),p,width=w,height=h,dpi=170,bg="#fcfcfb",limitsize=FALSE);cat("wrote",file.path(FIG,n),"\n")}
sc <- function(n,d){fwrite(d,file.path(OUT,n));cat("wrote",file.path(OUT,n),"\n")}
BIN_MM    <- 0.10   # 100 um bins
MIN_CELLS <- 20L    # cells in a bin before its fraction is drawn
CLAMP_Q   <- c(0.02, 0.98)

SURF<-"#fcfcfb"; INK<-"#0b0b0b"; INK2<-"#52514e"; MUTED<-"#898781"; BG<-"#e8e7e1"
ORD <- c("B cell","Plasmablast","Plasma cell","T cell","NK cell","Myeloid","pDC",
         "Endothelial","Stromal/FDC","Hodgkin chemokine niche")

a  <- fread(file.path(OUT,"all_cells_celltypes.csv.gz"), na.strings=c("NA",""))
md <- fread("output/16_regions/cells_meta_all.csv.gz",
            select=c("cell_id","x_slide_mm","y_slide_mm","is_defect","tier_matched","nCount_RNA"))
d <- md[a, on=.(cell_id)]
d[, celltype_final := fifelse(!is.na(celltype_cluster), celltype_cluster, celltype_projected)]
cat(sprintf("cells %d | with a final label %d (%.1f%%)\n", nrow(d),
            d[!is.na(celltype_final), .N], 100*d[, mean(!is.na(celltype_final))]))
## the low-information floor: non-defect AND tier-matched, yet outside the high-confidence
## object -- the only thing that can exclude such a cell is the floor added in 19a
d[, qc_class := fifelse(is_defect, "over-fixation defect",
                 fifelse(!tier_matched, "failed spatial QC tier",
                  fifelse(confidence == "high", "high confidence", "low-information floor")))]
cat("\n=== QC classes ===\n"); print(d[, .N, by=.(qc_class, confidence)][order(-N)])

d[, `:=`(bx = floor(x_slide_mm/BIN_MM), by_ = floor(y_slide_mm/BIN_MM))]
tot <- d[, .(n_bin = .N), by=.(samp, bx, by_)]
cat(sprintf("\nbins with >= %d cells: %d of %d\n", MIN_CELLS,
            tot[n_bin >= MIN_CELLS, .N], nrow(tot)))

frac_map <- function(col, lev) {
  x <- d[!is.na(get(col)), .(k = .N), by=.(samp, bx, by_, cls = get(col))]
  x <- tot[x, on=.(samp,bx,by_)][n_bin >= MIN_CELLS]
  x[, frac := k/n_bin]
  ## a bin with tissue but none of this class must read as zero, not as missing
  full <- CJ(idx = seq_len(nrow(tot[n_bin >= MIN_CELLS])), cls = lev, unique = TRUE)
  base <- tot[n_bin >= MIN_CELLS][, idx := .I]
  full <- base[full, on=.(idx)]
  x <- x[full, on=.(samp,bx,by_,cls)]
  x[is.na(frac), frac := 0]
  x[, cls := factor(cls, levels = lev)]
  x[, `:=`(x = bx*BIN_MM, y = by_*BIN_MM)]
  x[]
}

## EACH PANEL IS SCALED TO ITS OWN TYPE, with the absolute ceiling written into the strip so
## nothing is lost by it. One scale across all types is dominated by whichever is most
## abundant -- B cells reach ~80% of a bin in RB1, which flattened every other RB1 panel to
## background and made the figure unreadable for exactly the minor populations it exists to
## locate. "Where is this type concentrated" is a within-type question.
##
## geom_TILE, not geom_raster: the grid has holes wherever a bin fell below MIN_CELLS, and
## geom_raster assumes a complete lattice -- given gaps it SHIFTS pixels to close them, which
## silently distorts a spatial map.
panel <- function(g, S, ttl, sub_, pal_hi) {
  z <- copy(g[samp == S])
  z[, hi := quantile(frac[frac > 0], CLAMP_Q[2], na.rm = TRUE), by = cls]
  z[!is.finite(hi) | hi <= 0, hi := 1]
  z[, rel := pmin(frac / hi, 1)]
  key <- z[, .(hi = hi[1]), by = cls]
  key[, strip := sprintf("%s
up to %.0f%%", as.character(cls), 100 * hi)]
  setorder(key, cls)
  z <- key[, .(cls, strip)][z, on = .(cls)]
  z[, strip := factor(strip, levels = key$strip)]
  ggplot(z, aes(x, y, fill = rel)) + geom_tile(width = BIN_MM, height = BIN_MM) +
    scale_fill_gradientn(colours = c(BG, "#cfd8e6", pal_hi), limits = c(0, 1),
                         oob = scales::squish, labels = scales::percent,
                         name = "share of the bin's cells, as a fraction of this type's own ceiling") +
    facet_wrap(~ strip, nrow = 2) + coord_fixed() +
    labs(title = ttl, subtitle = sub_, x = "slide x (mm)", y = "slide y (mm)") +
    theme_minimal(base_size = 10) +
    theme(plot.background=element_rect(fill=SURF,colour=NA),
          panel.background=element_rect(fill=BG,colour=NA), panel.grid=element_blank(),
          axis.text=element_text(colour=MUTED,size=7), axis.title=element_text(colour=INK2,size=8),
          strip.text=element_text(colour=INK2,face="bold",size=7,lineheight=1.1),
          plot.title=element_text(colour=INK,face="bold"),
          plot.subtitle=element_text(colour=INK2,size=8), legend.position="bottom",
          legend.key.width=unit(30,"pt"), legend.key.height=unit(7,"pt"))
}

## ---- 1. cell types -----------------------------------------------------------
LEV <- intersect(ORD, unique(d$celltype_final))
G <- frac_map("celltype_final", LEV); sc("spatial_celltype_bins.csv.gz", G)
p <- panel(G,"IB1","IB1 -- where each cell type is concentrated",
           sprintf("Share of each %g um bin's cells carrying the label, each panel scaled to its OWN type (absolute ceiling in the strip).\nBins under %d cells are background. Slide mm, plain axes -- mm needs no reversal, canvas px would need both.",
                   BIN_MM*1000, MIN_CELLS), "#7a1b1b") /
     panel(G,"RB1","RB1 -- where each cell type is concentrated",
           "same scale construction, computed within RB1", "#7a1b1b")
sf("celltype_spatial.png", p, 17, 15)

## ---- 2. QC and rescue --------------------------------------------------------
QLEV <- c("high confidence","failed spatial QC tier","over-fixation defect","low-information floor")
Q <- frac_map("qc_class", QLEV); sc("spatial_qc_bins.csv.gz", Q)
CLEV <- c("high","rescued","unresolved")
C <- frac_map("confidence", CLEV)
q <- panel(Q,"IB1","IB1 -- where the QC flags fall","share of each bin's cells in each QC class","#1b4f7a") |
     panel(Q,"RB1","RB1 -- where the QC flags fall","","#1b4f7a")
r <- panel(C,"IB1","IB1 -- where the rescue recovered cells",
           "high = clustered directly; rescued = projected label passed its stratum gate; unresolved = no label","#0b6b4f") |
     panel(C,"RB1","RB1 -- where the rescue recovered cells","","#0b6b4f")
sf("qc_spatial.png", (q / r) + plot_layout(heights=c(1,1)), 20, 14)

cat("\n=== QC class by sample (% of the sample's cells) ===\n")
print(dcast(d[, .N, by=.(samp,qc_class)][, pct := round(100*N/sum(N),1), by=samp],
            qc_class ~ samp, value.var="pct"))
cat("\n=== confidence by sample ===\n")
print(dcast(d[, .N, by=.(samp,confidence)][, pct := round(100*N/sum(N),1), by=samp],
            confidence ~ samp, value.var="pct"))
cat("\ndone\n")
